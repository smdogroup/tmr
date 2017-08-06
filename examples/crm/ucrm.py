from mpi4py import MPI
from tmr import TMR

# Create the communicator
comm = MPI.COMM_WORLD

# Load the model from the STEP file
geo = TMR.LoadModel('first-section.stp')

# Get the volumes
vols = geo.getVolumes()

# Get the edges/faces from the geometry
faces = geo.getFaces()
edges = geo.getEdges()

# Set the master/target relationships
faces[4].setMaster(vols[0], faces[5])
edges[8].setMaster(edges[5])

# Create the geometry
mesh = TMR.Mesh(comm, geo)

# Mesh the part
opts = TMR.MeshOptions()
opts.num_smoothing_steps = 0

# Mesh the geometry with the given target size
htarget = 4.0
mesh.mesh(htarget, opts=opts)

# Write the mesh to a bdf file
mesh.writeToVTK("volume-mesh.vtk")
