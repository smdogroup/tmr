Mesh generation level classes
=============================

:class:`~tmr.TMR.Mesh` provides the primary interface to the meshing
functionality in :class:`TMR` and is used to create quadrilateral and hexahedral
meshes. The actual meshing operations are performed on the root processor with
rank 0, but all connectivity and location information is broadcast to all
processors.

:class:`~tmr.TMR.Mesh` attempts to mesh all :class:`tmr.TMR.Face` and :class:`tmr.TMR.Volume`
objects contained in :class:`~tmr.TMR.Model`. If only a surface mesh is desired, create
a new :class:`~tmr.TMR.Model` class without the volume objects.

:class:`TMR` only supports swept volume mesh generation. This requires the
specification of additional information to indicate what surfaces should be
linked and what direction should be swept. This information is provided by
indicating source and target geometry entities. When a source/target pair is
set, the mesh connectivity is copied from the source to the target. A swept
volume can only be created if:

#. :class:`~tmr.TMR.Volume` contains a source/target :class:`~tmr.TMR.Face` object
   pair that share the same number of edge loops and each source/target edge
   pair has the same number of nodes.

#. All other :class:`~tmr.TMR.Face` objects are structured with the same number
   of nodes along the swept direction.

To ensure these conditions apply, it is often necessary to set source/target
pairs for faces and edges.

The primary mesh-level functions are as follows:

* .. autoclass:: tmr.TMR.Mesh
      :members:

:class:`~tmr.TMR.MeshOptions` class defines a number of options that modify
the meshing algorithm. These include the following:

* .. autoclass:: tmr.TMR.MeshOptions
      :members:

:class:`~tmr.TMR.ElementFeatureSize` class is a base class that can be used to define the feature
size of the elements in the mesh.

* .. autoclass:: tmr.TMR.ElementFeatureSize

The mesh itself is created by meshing the vertices, edges, faces, and volumes
within the mesh.  These meshes are stored in the following classes:

The :class:`~tmr.TMR.EdgeMesh` object is created by first computing the number of
nodes along its length, and then globally ordering the nodes. The number of
nodes is determined based on the following criteria:

#. If the first and second vertices are different, then at least 3 nodes are
   created along the edge.
       
#. If the first and second vertices are the same and the edge is not
   degenerate, then it forms a loop back on itself and at least 5 nodes are
   created (double counting the first and last node number).

#. If the edge is a target edge, the number of nodes is taken from the source
   edge.

#. Otherwise, the number of nodes is selected as an odd number that most
   closely matches the spacing requested along the edge.

* .. autoclass:: tmr.TMR.EdgeMesh
          
* .. autoclass:: tmr.TMR.FaceMesh
   
* .. autoclass:: tmr.TMR.VolumeMesh
