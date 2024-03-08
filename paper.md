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
elements [@burstedde_coarse_2017].

The forest mesh is distributed, that is, at any time, each parallel process
only stores a unique portion of the forest mesh, the boundaries of which are
calculated from the SFC indices; see \autoref{fig:SpaceFillingCurves}.

![Left: Quad-tree of an exemplary forest mesh consisting of two trees
($\text{k}_{\text{0}}$, $\text{k}_{\text{1}}$) distributed over three parallel
processes P0 to P2. The SFC is represented by a black curve tracing only the
finest elements (leaf nodes) of each tree. Right: Sketch of the associated
triangular mesh refined up to level three.\label{fig:SpaceFillingCurves}](pics/forestmesh.pdf)

# Interfacing with t8code

In this section we discuss the main interface of t8code and how an
application would use it. While t8code offers various ways to interact with
meshes and data, we restrict ourselves to the most important functionality
here.

Every application is different and comes with their own requirements, data, and
adaptation criteria. In order to support a wide variety of use cases, our core
philosophy for t8code is to impose as few assumptions and to offer as much
freedom as possible. We cater for this by applying the Hollywood principle:
"Don't call us, we'll call you!". Whenever an application needs to interact
with the mesh, e.g., adapting the mesh, interpolating data, etc., we offer
suitable callback handlers.

The application developer implements custom callback functions and registers
them via the t8code application programming interface (API). Any mesh
specific details on how to access individual elements in the forest is opaque
to the application and internally handled by t8code in an efficient manner.
Of course, any typical application using hierarchical meshes needs to store
data on the elements of a forest. This data might correspond to some simulated
state variables, e.g., fluid velocity and temperature in a CFD simulation. In
accordance to our core philosophy, the data is only loosely coupled with
t8code's data structures. In order to properly access the application data in
the callbacks, the data simply needs to be provided as a consecutive array with
one entry per element enumerated in SFC order. For parallel applications,
access to neighboring elements across parallel zones (ghost layer) is provided
in a similar fashion.

## An example application

In the following section, we want to discuss the most important high-level
operations implemented in t8code. For this, consider a 3D numerical solver
application that traces a flow bubble moving around a rotating cylinder. The
application runs in parallel and the mesh is dynamically adapted in (almost)
every time step resolving the moving bubble with higher resolution than the
surrounding domain. These perpetual mesh changes constantly require the flow
state data to be interpolated from one adaption step to the next. A
visualization of such a setup might look like \autoref{fig:curved_advection}.

![Meshed region of fluid flow around a rotating cylinder. The green blob
corresponds to a bubble that is transported within the moving fluid. The mesh
is particularly refined along the boundary of the bubble. Colors encode the
element's distance from the bubble.
\label{fig:curved_advection}
](pics/curved_hybrid.png)

The standard way to implement such an application is to use the following
high-level t8code operations: New, Adapt, Balance, Interpolate,
Partition, Ghost, Iterate. This is also illustrated in the flowchart in
\autoref{fig:flowchart}. Next, we give more details about the different
operations:

New
: Construct a new, uniformly refined mesh from a coarse geometry mesh. This
: mesh is already distributed across the parallel processes. This step is usually
: only carried out once during the preprocessing phase.
   
Adapt
: Decide for each element whether to refine, coarsen, or pass according to the
: results of a criterion provided by a custom adaption callback.

Balance
: Establishesa 2:1 balance condition, meaning that afterwards the refinement
: levels of neighboring elements are either the same or differ by at most $\pm
: 1$. Note, this operation only refines elements, never coarsens them.
: Applications are free to decide whether they require the balance condition or
: not.

Interpolate
: Interpolate data from one forest mesh to another. For each element that was
: refined, coarsened or remained the same, an application provided callback is
: executed deciding how to map the data onto the new mesh.

Partition
: Re-partition the mesh across all parallel processes, such that each process
: has the same computational load (e.g. element count). Due to the SFC logic,
: this operation is very efficient and may be carried out in each time step.

Partition Data
: Redistribute any user defined data from the original mesh to the
: re-partitioned one. Input is an array with one entry for each element of the
: original forest containing the application data, output is an array with one
: entry for each element of the re-partitioned forest, containing the same data
: (that may previously have been on a different process).

Ghost
: Compute a list of all ghost elements of the current process. Ghosts are
: elements that are neighbors to elements of the process, but do not belong to
: the process itself.

Ghost exchange
: Transfer application specific data across all ghost elements. Input is an
: array of application data with one filled entry for each local element and one
: unfilled entry for each ghost. On output the entries at the ghost elements will
: be filled with the corresponding values from the neighbor processes.

Iterate
: Iterates through the mesh, providing face neighbor information for each
: element passed as an argument to the callback. In our example application, it
: is used to carry out the advection step of the bubble.

Search
: May be used additionally for extra tasks, such as searching for particles, or
: identifying flow features. It hierarchically iterates through the mesh and
: executes a callback function on all elements that match a given criterion.
: Leveraging the SFC tree logic, \texttt{Search} omits large chunks of the mesh
: if they do not match the criterion. Hence, it does not necessarily inspect each
: individual element and therefore performs much faster than a linear
: search [@holke_optimized_2021; @BursteddeSearch20].

![Flowchart of a typical simulation code which interacts with t8code.
Information about the different operations can be found in the text.
\label{fig:flowchart}
](pics/curved_hybrid.png)

# Modularity \& Extensibility

A distinct feature of t8code compared to similar AMR libraries is its high
modularity achieved by decoupling high-level from low-level algorithms and
coming along with it the support for arbitrary element shapes and refinement
patterns. It also allows to combine different element shapes within the same
mesh (hybrid meshes).

All high-level operations use the low-level algorithms only as a black box. For
example, mesh adaption routines iterate through the mesh and when necessary call
low-level algorithms for retrieving the children or the parent to
refine or coarsen an element. In order to implement the logic of
the adaption, however, no knowledge of the implementation details of these
low-level functions is required.

Thus, for each individual tree we can simply replace the underlying
implementation of the low-level algorithms (e.g. from tetrahedra to hexahedra)
without affecting the high-level functionality. We achieve this by
encapsulating all shape-specific element operations such as parent/child
computation, face-neighbor computation, SFC index computation and more in an
abstract C++ base class. The different element shapes and refinement
patterns are then specializations of this base class. Hence, t8code can be
easily extended - also by application developers - to support other refinement
patterns and SFCs.

Moreover, this very high degree of modularity allows us to support an even
wider range of non-standard additions. For example, the insertion of
sub-elements to resolve hanging nodes [@Becker_hanging_faces]
in quadrilateral meshes. Each quad element that has a hanging node is
subdivided into a set of several triangles eliminating the hanging node.

Furthermore, we added support for holes in the mesh by selectively
deleting elements [@Lilikakis_removing]. This feature can be used to
incorporate additional geometry information into the mesh. Similar to marking
elements as getting refined or coarsened, we can additionally mark elements as
getting removed. These elements will be eliminated completely from the SFC
reducing the overall memory footprint.

Addionally, we support curved hexahedra with geometry-informed
AMR [@elsweijer_curved_2021]. Thus, information such as element volumes,
face areas, or positions of interpolation/quadrature points in high order
meshes can be calculated exactly with respect to the actual geometry. Another
use case is to start with a very coarse input mesh and geometrically refine the
mesh maxing out the performance benefits of tree-based AMR.

![Strong scaling on JUWELS with tetrahedral elements. We plot the runtimes of
Ghost and Partition routines with a refinement band from levels 8 to
10 after four time steps. Hence, the forest mesh consists of approximately
$1.91$ billion tetrahedra. As observed in the plot, we achieve perfect scaling
for the Ghost algorithm in the number $G/P$ of ghosts per process. The runtime of
Partition is below 0.1 seconds even for the largest run. More details can
be found in [@holke_optimized_2021].
\label{fig:t8code_strong_scaling}](pics/t8code_strong_scaling.png)

# Performance

In this section we present some of our benchmark results from various
performance studies conducted on the JUQUEEN [@juqueen_fz_juelich] and the
JUWELS [@juwels_fz_juelich] supercomputers at the Jülich Supercomputing
Center. t8code's Ghost and Partition routines are exceptionally fast with
proper scaling of up to 1.1 trillion mesh elements; see
\autoref{tab:t8code_runtimes}, [@holke_optimized_2021].  In
\autoref{fig:t8code_strong_scaling} we show a strong scaling result for a
tetrahedral mesh achieving ideal strong scaling efficiency for the Ghost algorithm.
Furthermore, in a prototype code [@Dreyer2021] implementing a high-order
discontinuous Galerkin method (DG) for advection-diffusion equations on
dynamically adaptive hexahedral meshes we obverve a 12 times speed-up compared
to non-AMR meshes with only an overall 10 to 15\% runtime contribution of
t8code; see autoref{fig:t8code_runtimes}. 

+----------------+-------------------+--------------------+--------+-----------+
| \# process     | \# elements       | \# elem. / process | Ghost  | Partition |
+:==============:+:=================:+:==================:+:======:+:=========:+
|     49,152     | 1,099,511,627,776 |      22,369,621    | 2.08 s |   0.73 s  |
+----------------+-------------------+--------------------+--------+-----------+
|     98,304     | 1,099,511,627,776 |     11,184,811     | 1.43 s |   0.33 s  |
+================+===================+====================+========+===========+
| Runtimes on JUQUEEN for the ghost layer and partitioning operations for a    |
| distributed mesh consisting of 1.1 trillion elements.                        |
| \label{tab:t8code_runtimes}                                                  |
+================+===================+====================+========+===========+

![Runtimes on JUQUEEN of the different components of our DG prototype code
coupled with t8code. Note that all features associated with dynamical mesh
adaptation utilize only around 15\% of the total runtime largely independent of
the number of processes.\label{fig:t8code_runtimes}
](pics/t8code_runtimes_2.png)

# Conclusion

In this note we introduce our open source AMR library t8code. We give a brief
overview of the fundamental algorithms and data structures, namely our modular
SFC, and outline a general usage pattern when an application interacts with the
library. Due to its high modularity, t8code can be easily extended for a wide
range of use cases. Performance results confirm that t8code is a solid choice
for mesh management in high-performance applications in the upcoming exascale
era.

Future efforts will include an integration of our techniques into simulation
use cases with in-depth performance and accuracy evaluations. Additionally, we
strive to extend all presented features to all element shapes and space
dimensions. Other possible extensions that we plan to research in the near
future are mesh adaption of prism boundary layers and the support for
an-isotropic refinement.

# References
