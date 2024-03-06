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
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

In this note we present our software library t8code for
scalable dynamic adaptive mesh refinement (AMR) officially released in
2022 [@Holke_t8code_2022]. t8code is written in C/C++, open source, and
readily available at www.dlr-amr.github.io/t8code. The
library provides fast and memory efficient parallel algorithms for dynamic AMR
to handle tasks such as mesh adaptation, load-balancing, ghost computation,
feature search and more. t8code can manage meshes with over one
trillion mesh elements [@holke_optimized_2021] and scales up to one
million parallel processes [@holke_scalable_2018]. It is intended to
be used as mesh management back end in scientific and engineering simulation
codes paving the way towards high-performance applications of the upcoming
exascale era.

# Introduction

AMR has been established as a successful approach for scientific and
engineering simulations over the past decades [@TEUNISSEN2019106866;
@10.1145/1268776.1268779; @doi:10.1137/0733054; @doi:10.1137/0715049]. By
modifying the mesh resolution locally according to problem specific indicators,
the computational power is efficiently concentrated where needed and the
overall memory usage is reduced by orders of magnitude. However, managing
adaptive meshes and associated data is a very challenging task, especially for
parallel codes. Implementing fast and scalable AMR routines generally leads to
a large development overhead motivating the need for external mesh management
libraries like t8code.

t8code is written in C/C++, open source, and the latest release can be obtained at
https://dlr-amr.github.io/t8code [@Holke_t8code_2022].
It uses efficient space-filling curves (SFC) to manage the data in structured
refinement trees. While in the past being successfully applied to quadrilateral
and hexahedral meshes [@burstedde_p4est_2011; @weinzierl_peano_2019],
t8code extends these SFC techniques in a modular fashion, such that arbitrary
element shapes are supported. We achieve this modularity through a novel
decoupling approach that separates high-level (mesh global) algorithms from
low-level (element local) implementations. All high-level algorithms can then
be applied to different implementations of element shapes and refinement
patterns. A mix of different element shapes in the same mesh is also
supported.

Currently, t8code provides implementations of Morton type SFCs with
$1:2^d$ refinement for vertices ($d=0$), lines ($d=1$), quadrilaterals,
triangles ($d=2$), hexahedra, tetrahedra, prisms, and pyramids ($d=3$). The
latter having a $1:10$ refinement rule with tetrahedra emerging as child
elements [@Knapp20]. Additionally, implementation of other refinement
patterns and SFCs is possible according to the specific requirements of the
application.

The purpose of this note is to provide a brief overview and a first point of
entrance for software developers working on codes storing data on (distributed)
meshes. 

<!---
# The structure is as follows: \Secref{sec:foundations} gives a brief
# outline of the fundamental algorithms, \Secref{sec:interface} presents the
# interface, \Secref{sec:modularity} emphasizes the modularity of \tetcode while
# \Secref{sec:results} shows some performance results. Finally, in
# \Secref{sec:conclusion} we draw a conclusion and give a brief outlook.
-->

For further information beyond this short note and also for code examples, we
refer to our Documentation and Wiki [@Holke_t8code_2022] and our other
technical papers on t8code
[@holke_scalable_2018; @burstedde_coarse_2017; @holke_optimized_2021; @burstedde_tetrahedral_2016;
@Knapp20; @Becker_hanging_faces; @elsweijer_curved_2021; @Dreyer2021; @Lilikakis_removing].


# Fundamental Concepts

t8code is based on the concept of tree-based adaptive mesh refinement.
Starting point is an unstructured input mesh, which we call coarse mesh that
describes the geometry of the computational domain. The coarse mesh elements
are refined recursively in a structured pattern, resulting in refinement trees
of which we store only minimal information of the finest elements (the leafs of
the tree). We call this resulting fine mesh the forest.

By enumerating the children in the refinement pattern we obtain a space-filling
curve logic. Via these SFCs, all elements in a refinement tree are assigned an
index and are stored in linear order of these indices. Information such as
coordinates or element neighbors do not need to be stored explicitly, but can
be recovered from the index and the appropriate information of the coarse
elements. The less elements the input mesh has, the more memory and runtime are
saved through the SFC logic. t8code supports distributed coarse meshes of
arbitrary size and complexity, which we tested for up to 370 million input
elements~ [@burstedde_coarse_2017].

The forest mesh is distributed, that is, at any time, each parallel process
only stores a unique portion of the forest mesh, the boundaries of which are
calculated from the SFC indices; see Fig. \autoref{fig:SpaceFillingCurves}.

![Left: Quad-tree of an exemplary \forest mesh consisting of two trees
($\text{k}_{\text{0}}$, $\text{k}_{\text{1}}$) distributed over three parallel
processes P0 to P2. The SFC is represented by a black curve tracing only the
finest elements (leaf nodes) of each tree. Right: Sketch of the associated
triangular mesh refined up to level three.\label{fig:SpaceFillingCurves}](pics/forestmesh.pdf)

# References
