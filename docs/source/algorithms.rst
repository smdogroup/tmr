Algorithms
**********

The following sections describe the algorithms that are implemented with
:class:`TMR` at a high-level. 

Serial meshing capabilities
===========================
The 2D and mapped 2D meshing capabilities in :class:`TMR` are built around the
blossom-quad algorithm :cite:`Remacle:2012:blossom-quad`. This algorithm
generally produces good quality quadrilateral meshes with some restrictions on
the input geometry and the mesh size. The maximum mesh size is approximately
:math:`\frac{1}{2}` the minimum geometric feature size or edge length, or
:math:`\frac{1}{4}` the minimum hole size. This requirement arises because
there is a minimum of 2 elements per edge or 4 elements around a hole. Sharp
corners with small acute angles limit local mesh quality and lead to the
introduction of local kite elements.

Blossom-quad first builds a triangular surface mesh which is then converted to a
quadrilateral mesh through an optimal weighted matching using the blossom
algorithm. :class:`TMR` uses a generalization of Rebay's algorithm for efficient
unstructured triangular mesh generation
:cite:`Rebay:1993:unstructured-triangle`. This algorithm generates a triangular
surface mesh using a frontal method that updates a Delaunay triangularization
based on the Bowyer--Watson algorithm :cite:`Shewchuk:2012:Delaunay-notes`. Once
the triangular mesh is generated, it is smoothed using a Laplacian
technique. The recombination of the triangular mesh into a quad mesh is then
performed using the weighted matching algorithm.  The weights between adjacent
triangles are computed based on a combined quad quality, which is a function of
the maximum interior angle in the quadrilateral. The weights are a nonlinear
function of this quality metric, where angles close to :math:`180^{\circ}` are
most heavily penalized. Adjacent triangular elements along the boundary are also
linked with imaginary quads with a large weight. 

The matching algorithm selects which quads to combine in order to create the
mesh with the best overall quality based on the weights. The quad mesh is then
post-processed to remove poor quality quadrilateral elements with specific
topologies. Finally, the resulting quad mesh is smoothed using a quad-specific
algorithm :cite:`Giuliani:1982:quad-smoothing`. The smoothing algorithm
minimizes a combination of the squeeze and distortion in the mesh. 

Mapped Quadtree and octree meshing
==================================
Quadtree and octree meshing techniques provide an intermediate step between
fully structured meshes and unstructured meshes.  There are several ways that an
quadtree or octree can be stored.  The most natural way is to store the data in
a tree structure where the leaves of the tree represent elements in the mesh.
However, this type of data structure method requires extra overhead and memory
needed to traverse through the depth of the tree.

In :class:`TMR`, the leaves of the quadtree or octree are stored in an array data
structure directly, without referencing the nodes within the tree hierarchy.
This approach is based on the work of :cite:`BursteddeWilcoxGhattas11` and
:cite:`IsaacBursteddeWilcoxEtAl15`. Within this data structure, the leaves
represent either elements or nodes within the mesh. To mesh realistic
geometries, the method employs a semi-structured approach in which a global
unstructured quadrilateral or hexahedral super-mesh defines the connections
between quadtree or octrees. Each quadtree or octree is then parametrically
mapped to the geometry using an abstract geometry layer described below.

Each quadrant or octant within a local quadtree is uniquely identified based on
its local :math:`x`-:math:`y`-:math:`z` coordinates and its refinement
level. The coordinates correspond to the lower left hand corner of the quadrant
or octant and the level is used to determine the length of the side of an
element.This length can be used to determine the local coordinates of
neighboring elements. The tree can be traversed by creating parent or child
quadrants.

To enable local changes in the mesh refinement, adjacent quadrants may have
different levels, producing elements with different side-lengths. Dependent or
hanging nodes must be added when the mesh refinement level changes between
adjacent quadrants. To ensure finite-element compatibility between adjacent
elements, dependent nodes must be constrained along an edge. To ensure
compatibility, :class:`TMR` imposes that the relative difference in levels
between two adjacent quadrants cannot exceed one. This is referred to as 2:1
balancing.  After a refinement step, the global mesh must be balanced so that
both intra- and inter-quadtree quadrants are 2:1 balanced. Inter-quadtree
operations are performed by mapping the quadrants from their locally-aligned
coordinate system to the coordinate system of an adjacent quadtree along common
shared edges or corners.

In addition to the element mesh, the semi-structured method must also create and
uniquely order the nodes. A local Morton encoding on each tree to facilitate
these operations. The order of the quadrilaterals within the global quadtree
mesh is based on the unstructured global mesh. The ownership of shared corners,
edges, and faces is determined in advance to avoid duplicating node numbers.
During the node ordering process, hanging nodes are labeled and their
corresponding independent neighbors are determined so that they can be
eliminated using compatibility constraints during finite-element analysis.

Interface to CAD
================
The interface to CAD geometry is provided through an intermediate interface in
:class:`TMR`. An implementation of this interface is provided for OpenCascade
which enables models to be loaded directly from OpenCascade or intermediate
tools. For instance, OpenCascade can load geometry that has been exported in a
STEP format.

:class:`TMR` has internal definitions for both geometry and topological
entities.  The geometric entities include nodes, curves, and surfaces, while the
topological entities include vertex, edge, edge loop, face, and volume objects.
The topological entities describe the logical relationships between their
underlying geometric representations. Each vertex is associated with a node,
each edge with a curve, and each face with surface. Each edge contains a first
and second vertex denoting the start and end point of the underlying curve. An
edge loop contains a series of edges that are connected together to form a
closed loop on a surface. Each face is bounded by a series of edge loops that
define the extent of the face and any holes or cutouts that form the surface.
The volume object contains a series of faces that bound the volume to create a
watertight surface.

The orientation of geometry objects within the model are of critical importance.
Incorrect orientation information can produce a model with an incompatible
topology that does not reflect the underlying geometry. Each edge has a natural
orientation defined by its first and second vertices. However, within an edge
loop, the orientation of a particular edge may be reversed relative to its
natural orientation. Therefore, the edge loop object stores both the list of
edges that form the loop and their orientations relative to the natural edge
orientation.

Faces in :class:`TMR` are defined in their natural orientation. This means that
the parametric areas computed are positive such that the parametric system is
right-handed. The orientation of an edge loop on a face is always defined such
that the material lies to the left of an edge within an edge loop when walking
the loop in its positive orientation (taking into account the relative
orientation in the edge loop object). The natural orientation may vary from the
orientation used within the CAD package and so :class:`TMR` stores a flag to
indicate whether the orientation of the surface is flipped. The surfaces that
create a volume have an orientation, defined with an outward normal
direction. The volume object contains a list of surfaces and their orientations
relative to the natural surface orientation stored in :class:`TMR` that
indicates whether the surface is outward facing.

The meshes within :class:`TMR` are generated and stored on a component-level
basis. Meshing proceeds through the hierarchy of topological entities from
vertices, edges, faces to volumes. Each component stores its own portion of the
mesh using its natural orientation. However, when the mesh is extracted, the
orientation is converted to the CAD-based orientation by flipping the
orientation when it does not align.

.. bibliography:: ../refs.bib
