Example
*******

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
