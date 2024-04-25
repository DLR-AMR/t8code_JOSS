---
title: 't8code - modular adaptive mesh refinement in the exascale era'
tags:
  - C
  - C++
  - adaptive mesh refinement
  - exascale
  - hypbrid meshes
  - modularity
authors:
  - name: Johannes Holke
    # orcid: 0000-0000-0000-0000
    corresponding: true # (This is how to denote the corresponding author)
    equal-contrib: true
    affiliation: 1
  - name: Johannes Markert
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 1
  - name: David Knapp
    equal-contrib: true
    affiliation: 1
  - given-names: Lukas Dreyer
    affiliation: 1
  - given-names: Sandro Elsweijer
    affiliation: 1
  - given-names: Niklas Böing
    affiliation: 1
  - given-names: Chiara Hergl
    affiliation: 1
  - given-names: Prasanna Ponnusamy
    affiliation: 1
  - given-names: Achim Basermann
    affiliation: 1
affiliations:
 - name: German Aerospace Center (DLR), Institute for Software Technology, Cologne, Germany
   index: 1
date: 06 March 2024
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

We present our software library t8code for scalable dynamic adaptive mesh
refinement (AMR) officially released in 2022 [@Holke_t8code_2022]. t8code is
written in C/C++, open source, and readily available at
[www.dlr-amr.github.io/t8code](www.dlr-amr.github.io/t8code). The library provides fast and memory
efficient parallel algorithms for dynamic AMR to handle tasks such as mesh
adaptation, load-balancing, ghost computation, feature search and more. t8code
can manage meshes with over one trillion mesh elements [@holke_optimized_2021]
and scales up to one million parallel processes [@holke_scalable_2018]. It is
intended to be used as mesh management backend in scientific and engineering
simulation codes paving the way towards high-performance applications of the
upcoming exascale era.

# Statement of need

Adaptive mesh refinement has been established as a successful approach
for scientific and engineering simulations over the past decades
[@TEUNISSEN2019106866; @10.1145/1268776.1268779; @doi:10.1137/0733054;
@doi:10.1137/0715049]. By modifying the mesh resolution locally according to
problem specific indicators, the computational power is efficiently
concentrated where needed and the overall memory usage is reduced by orders of
magnitude. However, managing adaptive meshes and associated data is a very
challenging task, especially for parallel codes. Implementing fast and scalable
AMR routines generally leads to a large development overhead motivating the
need for external mesh management libraries like t8code.

Currently, t8code's AMR routines support a wide range of element types:
vertices, lines, quadrilaterals, triangles, hexahedra, tetrahedra, prisms, and
pyramids. The latter having a $1:10$ refinement rule with tetrahedra emerging
as child elements [@Knapp20].  Additionally, implementation of other refinement
patterns and element shapes is possible according to the specific requirements
of the application. t8code aims to provide a comprehensive mesh management
framework for a wide range of use cases in science and engineering
applications.

# Fundamental Concepts

t8code is based on the concept of tree-based adaptive mesh refinement.
Starting point is an unstructured input mesh, which we call coarse mesh that
describes the geometry of the computational domain. The coarse mesh elements
are refined recursively in a structured pattern, resulting in refinement trees
of which we store only minimal information of the finest elements (the leafs of
the tree). We call this resulting fine mesh the forest.

By enumerating the children in the refinement pattern we obtain a space-filling
curve (SFC) logic. Via these SFCs, all elements in a refinement tree are assigned an
index and are stored in linear order of these indices. Information such as
coordinates or element neighbors do not need to be stored explicitly, but can
be recovered from the index and the appropriate information of the coarse
elements. The forest mesh can be distributed, that is, at any time, each
parallel process only stores a unique portion of the forest mesh, the
boundaries of which are calculated from the SFC indices; see
\autoref{fig:SpaceFillingCurves}.

While in the past being successfully applied to quadrilateral
and hexahedral meshes [@burstedde_p4est_2011; @weinzierl_peano_2019],
t8code extends these SFC techniques in a modular fashion, such that arbitrary
element shapes are supported. We achieve this modularity through a novel
decoupling approach that separates high-level (mesh global) algorithms from
low-level (element local) implementations. All high-level algorithms can then
be applied to different implementations of element shapes and refinement
patterns. A mix of different element shapes in the same mesh is also
supported.

![Left: Quad-tree of an exemplary forest mesh consisting of two trees
($\text{k}_{\text{0}}$, $\text{k}_{\text{1}}$) distributed over three parallel
processes p0 to p2. The SFC is represented by a black curve tracing only the
finest elements (leaf nodes) of each tree. Right: Sketch of the associated
mixed shape mesh refined up to level three.\label{fig:SpaceFillingCurves}](pics/t8code_sfc_hybrid.png)

# Performance

t8code supports distributed coarse meshes of arbitrary size and complexity,
which we tested for up to 370 million input elements [@burstedde_coarse_2017].
Moreover, we present some of our benchmark results from various
performance studies conducted on the JUQUEEN [@juqueen_fz_juelich] and the
JUWELS [@juwels_fz_juelich] supercomputers at the Jülich Supercomputing
Center. t8code's ghost and partition routines are exceptionally fast with
proper scaling of up to 1.1 trillion mesh elements; see
\autoref{tab:t8code_runtimes}, [@holke_optimized_2021]. 
Furthermore, in a prototype code [@Dreyer2021] implementing a high-order
discontinuous Galerkin method (DG) for advection-diffusion equations on
dynamically adaptive hexahedral meshes we obverve a 12 times speed-up compared
to non-AMR meshes with only an overall 15\% runtime contribution of
t8code; see \autoref{fig:t8code_runtimes}. 

+----------------+-------------------+--------------------+--------+-----------+
| \# process     | \# elements       | \# elem. / process | Ghost  | Partition |
+:==============:+:=================:+:==================:+:======:+:=========:+
|     49,152     | 1,099,511,627,776 |     22,369,621     | 2.08 s |   0.73 s  |
+----------------+-------------------+--------------------+--------+-----------+
|     98,304     | 1,099,511,627,776 |     11,184,811     | 1.43 s |   0.33 s  |
+================+===================+====================+========+===========+
| Table 1: Runtimes on JUQUEEN for the ghost layer and partitioning operations |
| for a distributed mesh consisting of 1.1 trillion elements.                  |
| \label{tab:t8code_runtimes}                                                  |
+================+===================+====================+========+===========+

![Runtimes on JUQUEEN of the different components of our DG prototype code
coupled with t8code. Note that all features associated with dynamical mesh
adaptation utilize only around 15\% of the total runtime largely independent of
the number of processes.\label{fig:t8code_runtimes}
](pics/t8code_runtimes_2.png)

# Conclusion

In this note, we introduce our open source AMR library t8code. We give a brief
overview of the fundamental design principles and high-level operations. Due to
the high modularity, t8code can be easily extended for a wide range of use
cases. Performance results confirm that t8code is a solid choice for mesh
management in high-performance applications in the upcoming exascale era.

For further information beyond this short note and also for code examples, we
refer to our
[Documentation](https://dlr-amr.github.io/t8code/pages/documentation.html) and
[Wiki](https://github.com/DLR-AMR/t8code/wiki) reachable via our homepage
[www.dlr-amr.github.io/t8code](www.dlr-amr.github.io/t8code) and our technical
publications on t8code [@holke_scalable_2018; @burstedde_coarse_2017;
@holke_optimized_2021; @burstedde_tetrahedral_2016; @Knapp20;
@Becker_hanging_faces; @elsweijer_curved_2021; @Dreyer2021;
@Lilikakis_removing; @Holke_t8code_2022].

# References
