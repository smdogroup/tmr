Geometry classes
================

TMREntity
---------
Almost all classes in :class:`TMR` inherit from :class:`~tmr.TMR.Entity`. Those that do
not are low-level data structures that are used within high-level
operations. :class:`~tmr.TMR.Entity` defines a reference counting scheme for
all objects allocated on the heap. This is handled automatically through the
Python interface.

Loading a model
---------------

The primary geometry-level classes used within :class:`TMR` consist of entity-level classes 
(:class:`~tmr.TMR.Vertex`, :class:`~tmr.TMR.Edge`, and :class:`~tmr.TMR.Face`) and container-level classes 
(:class:`~tmr.TMR.EdgeLoop`, :class:`~tmr.TMR.Volume`, and :class:`~tmr.TMR.Model`).
The member functions for these classes are described below.

The geometry of a model is contained within the class :class:`~tmr.TMR.Model` and
is often created by a call to the function :func:`~tmr.TMR.LoadModel`. This function
loads the model data from a STEP, IGES or EGADS file and creates the TMR versions of
the geometric entities and containers that define the model. This call attempts
to remove any extraneous geometry defined in the file but not contained
within the model. This function creates and initializes the internal topology of
the model, which can then be used in subsequent meshing operations:

* .. autofunction:: tmr.TMR.LoadModel

Alternatively, a :class:`~tmr.TMR.Model` instance can be created by direct in-memory
wrapping of EGADS objects through egads4py:

* .. autofunction:: tmr.TMR.ConvertEGADSModel

A direct connection between the geometry and octree or quadtree levels is
provided through :class:`~tmr.TMR.Topology` that defines the topology for
mapped-quadrilateral and hexahedral geometries. This class can be created
directly, but it is most common to create this using the mesh-level classes.

Geometry interface
------------------

The interface classes to the underlying geometry consist of:

* .. autoclass:: tmr.TMR.Vertex
     :members:

* .. autoclass:: tmr.TMR.Edge
      :members:

* .. autoclass:: tmr.TMR.Face
      :members:

The other objects define collections of geometric entities:

* .. autoclass:: tmr.TMR.EdgeLoop
      :members:
          
* .. autoclass:: tmr.TMR.Volume
      :members:

* .. autoclass:: tmr.TMR.Model
      :members:
