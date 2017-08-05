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

# Set the master relations
for i in xrange(len(faces)):
    faces[i].writeToVTK('face%d.vtk'%(i))

faces[4].setMaster(faces[5])
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
mesh.writeToBDF("volume-mesh.bdf", 'quads')
mesh.writeToVTK("volume-mesh.vtk", 'quads')
