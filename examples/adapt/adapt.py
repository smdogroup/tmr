from mpi4py import MPI
from tmr import TMR
from tacs import TACS, elements, constitutive
import numpy as np
import argparse
import os

class CreateMe(TMR.QuadCreator):
    def __init__(self, bcs):
        TMR.QuadCreator.__init__(bcs)

    def createElement(self, order, quad):
        # Set constitutive properties
        rho = 0.102 # density, lb/in^3
        E = 10e6 # elastic modulus, psi
        nu = 0.3 # poisson's ratio
        kcorr = 5.0/6.0 # shear correction factor
        ys = 73e3 # yield stress, psi
        min_thickness = 0.2
        max_thickness = 1.0
        thickness = 0.5 # in
        
        stiff = constitutive.isoFSDT(rho, E, nu, kcorr, ys, 
                                     thickness, quad.face,
                                     min_thickness, max_thickness)
        
        theta = 30.0*np.pi/180.0
        stiff.setRefAxis(np.array([np.sin(theta), np.cos(theta), 0.0]))
        element = None
        if order == 2:
            element = elements.MITCShell(2, stiff)
        elif order == 3:
            element = elements.MITCShell(3, stiff)

        return element

# Set the communicator
comm = MPI.COMM_WORLD

# Create an argument parser to read in arguments from the commnad line
p = argparse.ArgumentParser()
p.add_argument('--htarget', type=float, default=10.0)
p.add_argument('--filename', type=str, default='wingbox_solid1.stp', 
               help='STEP file name')
p.add_argument('--output', type=str, default='surface-mesh.vtk',
               help='output file name')
p.add_argument('--forest_output', type=str, default='forest-mesh.vtk',
               help='forest output file name')
p.add_argument('--order', type=int, default=3, help='mesh order = 2 or 3')
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

geo.writeModelToTecplot('model.dat')

# Set the names of the geometric objects that will be fixed
faces[28].setAttribute('fixed')
for index in [67, 72, 74, 75]:
    edges[index].setAttribute('fixed')
for index in [44, 46, 47, 45]:
    verts[index].setAttribute('fixed')

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

# Set the element order
order = args.order

# Create random trees and balance the mesh. Print the output file
if order == 3:
    nlevels = 2
elif order == 2:
    nlevels = 3

# Create the forest
forest.createTrees(nlevels+1)
forest.balance(1)

# Make the creator class
bcs = TMR.BoundaryConditions()
bcs.addBoundaryCondition('fixed')
creator = CreateMe(bcs)

# Create the forests
forests = [ forest ]
assemblers = [ creator.createTACS(order, forest) ]

if order == 3:
    order = 2
    forests.append(forests[-1].duplicate())
    forests[-1].balance(1)
    assemblers.append(creator.createTACS(order, forests[-1]))

for i in xrange(nlevels-1):
    forests.append(forests[-1].coarsen())
    forests[-1].balance(1)
    assemblers.append(creator.createTACS(order, forests[-1]))

# Create the multigrid object
mg = TMR.createMg(assemblers, forests)

# Create a solution vector
ans = assemblers[0].createVec()
res = assemblers[0].createVec()
alpha = 1.0
beta = 0.0
gamma = 0.0
mg.assembleJacobian(alpha, beta, gamma, res)

# Factor the preconditioner
mg.factor()

# Set up a arbitrary load
v = res.getArray()
v[2::6] += 1.0
assemblers[0].applyBCs(res)

subspace = 100
gmres = TACS.KSM(mg.getMat(), mg, subspace, isFlexible=1)
gmres.setMonitor(comm, 'GMRES', 1)
gmres.solve(res, ans)
assemblers[0].setVariables(ans)

# Output for visualization 
flag = (TACS.ToFH5.NODES |
        TACS.ToFH5.DISPLACEMENTS |
        TACS.ToFH5.STRAINS)
f5 = TACS.ToFH5(assemblers[0], TACS.PY_SHELL, flag)
f5.writeToFile('visualization.f5')




