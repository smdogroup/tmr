from mpi4py import MPI
from tmr import TMR

# The communicator
comm = MPI.COMM_WORLD

# Load the geometry model
geo = TMR.LoadModel('beam.stp')

# Mark the boundary condition faces
verts = geo.getVertices()
faces = geo.getFaces()
volumes = geo.getVolumes()
faces[3].setAttribute('fixed')
faces[4].setSource(volumes[0], faces[5])
verts[4].setAttribute('pt1')
verts[3].setAttribute('pt2')

# Set the boundary conditions for the problem
bcs = TMR.BoundaryConditions()
bcs.addBoundaryCondition('fixed')

# Create the mesh
mesh = TMR.Mesh(comm, geo)

# Set the meshing options
opts = TMR.MeshOptions()

# Create the surface mesh
mesh.mesh(4.0, opts)

# Create a model from the mesh
model = mesh.createModelFromMesh()

# Create the corresponding mesh topology from the mesh-model 
topo = TMR.Topology(comm, model)

# Create the quad forest and set the topology of the forest
forest = TMR.OctForest(comm)
forest.setTopology(topo)

# Create the trees
forest.createTrees(0)

# Get the octants
octants = forest.getOctants()

# TMR.MAX_LEVEL
lev = 8
hmax = 1 << TMR.MAX_LEVEL
h = 1 << (TMR.MAX_LEVEL - lev)

for i in range(len(octants)):
    oc = octants[i]
    oc.level = lev
    oc.x = hmax - h
    oc.y = hmax - h
    oc.z = hmax - h
    octants[i] = oc

forest.writeForestToVTK('pre_forest.vtk')
forest.balance(1)
forest.writeForestToVTK('forest.vtk')
