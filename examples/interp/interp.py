import os
import numpy as np
from mpi4py import MPI
from tmr import TMR
from tacs import TACS

# Set the communicator
comm = MPI.COMM_WORLD

# The fine octree forest
fine = None

stepfile = "beam.stp"
if os.path.isfile(stepfile):
    # Load the geometry model
    geo = TMR.LoadModel(stepfile)

    # Mark the boundary condition faces
    faces = geo.getFaces()
    volumes = geo.getVolumes()
    faces[4].setSource(volumes[0], faces[5])

    # Create the mesh
    mesh = TMR.Mesh(comm, geo)

    # Set the meshing options
    opts = TMR.MeshOptions()
    opts.frontal_quality_factor = 1.25
    opts.num_smoothing_steps = 10
    opts.write_mesh_quality_histogram = 0

    # Create the surface mesh
    htarget = 4.0
    mesh.mesh(htarget, opts)

    # Create a model from the mesh
    model = mesh.createModelFromMesh()

    # Create the corresponding mesh topology from the mesh-model
    topo = TMR.Topology(comm, model)
    fine = TMR.OctForest(comm)
    fine.setTopology(topo)
else:
    conn = np.array(
        [[0, 1, 3, 4, 6, 7, 9, 10], [8, 11, 2, 5, 7, 10, 1, 4]], dtype=np.intc
    )
    fine = TMR.OctForest(comm)
    fine.setConnectivity(conn)

# Create the fine mesh
fine.createTrees(0)
refine = np.zeros(len(fine.getOctants()), dtype=np.intc)
if comm.rank == 0:
    refine[0] = 4
fine.refine(refine)
fine.balance(1)
fine.repartition()
fine.setMeshOrder(4, TMR.GAUSS_LOBATTO_POINTS)
fine.createNodes()

# Create a refinement array
octants = fine.getOctants()
refine = np.zeros(len(octants), dtype=np.intc)
for i in range(len(octants)):
    if i % 7 == 0:
        refine[i] = 2
    elif i % 5 == 0:
        refine[i] = -1

# Make the coarse tree finer than the fine tree for testing purposes
coarse = fine.duplicate()
coarse.refine(refine)
coarse.balance(0)
coarse.repartition()
coarse.setMeshOrder(2, TMR.GAUSS_LOBATTO_POINTS)
coarse.createNodes()

coarse_range = coarse.getNodeRange()
nc = coarse_range[comm.rank + 1] - coarse_range[comm.rank]
coarse_map = TACS.NodeMap(comm, nc)

fine_range = fine.getNodeRange()
nf = fine_range[comm.rank + 1] - fine_range[comm.rank]
fine_map = TACS.NodeMap(comm, nf)

# Create the two interpolations fine -> coarse and coarse -> fine
interp = TACS.VecInterp(coarse_map, fine_map, 1)
fine.createInterpolation(coarse, interp)
interp.initialize()

interp2 = TACS.VecInterp(fine_map, coarse_map, 1)
coarse.createInterpolation(fine, interp2)
interp2.initialize()
