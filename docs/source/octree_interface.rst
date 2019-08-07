Quad/octree-level classes
=========================

The quadtree and octree level classes are used to create, refine and manipulate
the semi-structured mesh in parallel. The quadrants or octants are stored in a
distributed manner across processors. The following section describes the
functions for the forest of octree operations. Analogous functionality is
defined for the forest of quadtree object as well. The first step in creating a
:class:`~tmr.TMR.OctForest` object is typically creating a :class:`~tmr.TMR.Topology`
which defines a model with mapped 2D and 3D geometry that consists entirely of
non-degenerate quadrilateral surfaces with four non-degenerate edges and
hexahedral volumes with 6 enclosing surfaces. This type of object can be
created using the constructor as follows:

* .. autoclass:: tmr.TMR.Topology

The :class:`~tmr.TMR.Topology` object defines additional functionality that is
generally required only for lower-level operations. Once the topology object has
been created, the :class:`~tmr.TMR.OctForest` can be initialized. It has the
following functions:

* .. autoclass:: tmr.TMR.OctForest
    :members:

Similarly, :class:`~tmr.TMR.QuadForest` can be initialized and has the similar functions:

* .. autoclass:: tmr.TMR.QuadForest
    :members:

Typical Usage
-------------
The typical usage for a :class:`~tmr.TMR.OctForest` would consist of the following:

#. Create the object and call :func:`~tmr.TMR.OctForest.setTopology` to set the
   super mesh
#. Create an initial element mesh by calling :func:`~tmr.TMR.OctForest.createTrees`

#. Create a refined mesh by duplicating and refining the octree forest by
   calling :func:`~tmr.TMR.OctForest.duplicate` followed by
   :func:`~tmr.TMR.OctForest.refine`. Note that the length of the integer array
   passed to refine must be equal to the number of elements in the original mesh.

#. Balance the mesh and create nodes by calling :func:`~tmr.TMR.OctForest.balance`
   then :func:`~tmr.TMR.OctForest.createNodes`

#. Create a sequence of coarser lower-order meshes by repeated calling
   :func:`~tmr.TMR.OctForest.duplicate` and :func:`~tmr.TMR.OctForest.coarsen`

#. Create interpolants between meshes for multigrid methods by calling
   :func:`~tmr.TMR.OctForest.createInterpolation`  

See :doc:`example` for an example of this usage.
