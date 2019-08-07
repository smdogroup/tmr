STEP import and meshing
=======================

The following example loads in a CAD model from a STEP file, and creates a volume mesh.
The code consists of the following steps:

* Import the CAD model from a STEP geometry file into a TMR-compatible objects using ``TMR.LoadModel`` 
* Extract lists of the volumes, faces and edges defined from the file
* Set source and target edges and faces. In the source-target relationship, target edges will have the same number of nodes as their source. Target faces will have the same surface mesh topology as their source faces. Source and target faces must be in the same volume object.
* The call ``mesh = TMR.Mesh(comm, geo)`` generates a new mesh. If the model input contains vertices, edges and faces, then a surface mesh is created, if it also contains volumes, then a volume mesh is created.

Listed below is an example of using TMR to generate a quadtree/octree mesh

.. code-block:: python
    
    from mpi4py import MPI
    from tmr import TMR

    # Set the communicator
    comm = MPI.COMM_WORLD

    # Load the model from the STEP file
    geo = TMR.LoadModel('geometry.stp')

    # Get the volumes
    vols = geo.getVolumes()

    # Get the edges/faces from the geometry
    faces = geo.getFaces()
    edges = geo.getEdges()

    # Set the source/target relationships
    faces[4].setSource(vols[0], faces[5])
    edges[8].setSource(edges[5])

    # Create the geometry
    mesh = TMR.Mesh(comm, geo)

    # Mesh the part
    opts = TMR.MeshOptions()
    opts.num_smoothing_steps = 0

    # Mesh the geometry with the given target size
    htarget = 4.0
    mesh.mesh(htarget, opts=opts)

It is sometimes necessary to identify the numbering of vertex, edges, faces, and volumes within a file.
A convenient way to visualize this information is to create a tecplot file with labels.
An example of generating this information is shown below:

.. code-block:: python
    
    from mpi4py import MPI
    from tmr import TMR

    # Set the communicator
    comm = MPI.COMM_WORLD

    # Load the model from the STEP file
    geo = TMR.LoadModel('geometry.stp')
    
    # Output the geometry model with the numbering convention
    geo.writeModelToTecplot('tecplot_surfaces.dat')
