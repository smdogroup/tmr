from mpi4py import MPI
from tmr import TMR
from tacs import TACS, elements, constitutive
import numpy as np
import argparse
import os

def createTACS(comm, forest, order=2):
    '''Create the TACSAssembler object from the forest'''

    # Create the nodes
    forest.createNodes(order)
    forest.createDepNodeConn()

    # Create the connectivity
    quads = forest.createMeshConn()

    # Set the connectivity
    num_elems = quads.shape[0]

    # Get the dependent node connectivity
    ptr, conn, weights = forest.getDepNodeConn()
    ndep = len(ptr)-1

    # Get the owned node range
    r = forest.getNodeRange()
    owned_nodes = r[comm.rank+1] - r[comm.rank]
    
    # Create the TACS assembler object   
    per_node = 6
    assembler = TACS.Assembler.create(comm, per_node,
                                      owned_nodes, num_elems, ndep)

    # Set the connectivity
    quad_ptr = np.arange(0, order*order*num_elems+1, order*order,
                         dtype=np.intc)
    assembler.setElementConnectivity(quads.flatten(), quad_ptr)

    # Set the dependent node connectivity
    assembler.setDependentNodes(ptr, conn, weights)

    # Create the elements
    elems = []
    for i in xrange(num_elems):
        # Set constitutive properties
        rho = 2500.0 # density, kg/m^3
        E = 70e9 # elastic modulus, Pa
        nu = 0.3 # poisson's ratio
        kcorr = 5.0 / 6.0 # shear correction factor
        ys = 350e6 # yield stress, Pa
        min_thickness = 0.002
        max_thickness = 0.20
        thickness = 0.02
        
        stiff = constitutive.isoFSDT(rho, E, nu, kcorr, ys, thickness, i,
                                     min_thickness, max_thickness)
        if order == 2:
            element = elements.MITCShell(2, stiff)
        elif order == 3:
            element = elements.MITCShell(3, stiff)
        elems.append(element)

    # Set the elements
    assembler.setElements(elems)

    # Initialize the assembler object
    assembler.initialize()

    return assembler

# Set the communicator
comm = MPI.COMM_WORLD

# Create an argument parser to read in arguments from the commnad line
p = argparse.ArgumentParser()
p.add_argument('--htarget', type=float, default=10.0)
p.add_argument('--filename', type=str, default=None, help='STEP file name')
p.add_argument('--output', type=str, default='surface-mesh.vtk',
               help='output file name')
p.add_argument('--forest_output', type=str, default='forest-mesh.vtk',
               help='forest output file name')
args = p.parse_args()

# Get the value of the filename
filename = args.filename
if not os.path.isfile(filename):
    raise ValueError('File %s does not exist'%(filename))

# Set the value of the target length scale in the mesh
htarget = args.htarget

# Load the geometry model
geo = TMR.LoadModel(filename)

# Create a model by discarding the volumes
verts = geo.getVertices()
edges = geo.getEdges()
faces = geo.getFaces()
geo = TMR.Model(verts, edges, faces)

# Create the new mesh
mesh = TMR.Mesh(comm, geo)

# Set the meshing options
opts = TMR.MeshOptions()
opts.frontal_quality_factor = 1.25
opts.num_smoothing_steps = 10
opts.write_mesh_quality_histogram = 1
opts.triangularize_print_iter = 50000

# Create the surface mesh
mesh.mesh(htarget, opts)

# Write the surface mesh to a file
mesh.writeToVTK(args.output)

# Create a model from the mesh
model = mesh.createModelFromMesh()

# Create the corresponding mesh topology from the mesh-model 
topo = TMR.Topology(comm, model)

# Create the quad forest and set the topology of the forest
forest = TMR.QuadForest(comm)
forest.setTopology(topo)

# Create random trees and balance the mesh. Print the output file
nlevels = 2
# forest.createTrees(nlevels)
forest.createRandomTrees(nrand=10, min_lev=0, max_lev=7)
t0 = MPI.Wtime()
forest.balance(1)
if comm.rank == 0:
    print 'Balance time ', MPI.Wtime() - t0
filename = os.path.splitext(args.forest_output)[0] + '%d.vtk'%(comm.rank)
forest.writeForestToVTK(filename)

# Create the forests
forests = [ forest ]
t0 = MPI.Wtime()
assemblers = [ createTACS(comm, forest) ]
if comm.rank == 0:
    print 'TACS creation time ', MPI.Wtime() - t0

for i in xrange(nlevels):
    forests.append(forests[-1].coarsen())
    t0 = MPI.Wtime()
    forests[-1].balance(1)
    if comm.rank == 0:
        print 'Coarse balance time ', MPI.Wtime() - t0
    assemblers.append(createTACS(comm, forests[-1]))

# Create the multigrid object
mg = TMR.createMg(assemblers, forests)
mat = mg.getMat()




