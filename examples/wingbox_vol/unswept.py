from mpi4py import MPI
from tmr import TMR
from paropt import ParOpt
from tacs import TACS, elements, constitutive, functions
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

class CreateMe(TMR.OctTopoCreator):
    def __init__(self, bcs, filt):
        TMR.OctTopoCreator.__init__(bcs, filt)

    def createElement(self, order, octant, index, weights):
        '''Create the element'''
        rho = 2600.0
        E = 70e9
        nu = 0.3
        stiff = TMR.OctStiffness(rho, E, nu, index, weights, q=5.0)
        elem = elements.Solid(2, stiff)
        return elem

def addFaceTraction(order, forest, attr, assembler, tr):
    trac = []
    for findex in range(6):
        trac.append(elements.Traction3D(order, findex, tr[0], tr[1], tr[2]))

    # Retrieve octants from the forest
    octants = forest.getOctants()
    face_octs = forest.getOctsWithAttribute(attr)
    aux = TACS.AuxElements()

    for i in range(len(face_octs)):
        index = octants.findIndex(face_octs[i])
        if index is not None:
            aux.addElement(index, trac[face_octs[i].tag])

    return aux

def createTopoProblem(forest, order=2, nlevels=2):
    # Create the forest
    forests = []
    filters = []
    assemblers = []
    varmaps = []
    vecindices = []

    # Create the trees, rebalance the elements and repartition
    forest.balance(1)
    forest.repartition()
    forests.append(forest)

    # Create the filter
    filtr = forest.coarsen()
    filtr.balance(1)
    filters.append(filtr)
    
    # Make the creator class
    creator = CreateMe(bcs, filters[-1])
    assemblers.append(creator.createTACS(order, forest))
    varmaps.append(creator.getMap())
    vecindices.append(creator.getIndices())

    for i in xrange(nlevels-1):
        forest = forests[-1].coarsen()
        forest.balance(1)
        forest.repartition()
        forests.append(forest)

        # Create the filter
        filtr = forest.coarsen()
        filtr.balance(1)
        filters.append(filtr)

        # Make the creator class
        creator = CreateMe(bcs, filters[-1])
        assemblers.append(creator.createTACS(order, forest))
        varmaps.append(creator.getMap())
        vecindices.append(creator.getIndices())

    # Create the multigrid object
    mg = TMR.createMg(assemblers, forests)

    # Create the topology optimization problem
    problem = TMR.TopoProblem(assemblers, filters, varmaps, vecindices, mg)

    return assemblers[0], problem

def initGeo():
    filename = 'geom/unswept.stp'
    geo = TMR.LoadModel(filename)

    # Create a model by discarding the volumes
    verts = geo.getVertices()
    edges = geo.getEdges()
    faces = geo.getFaces()
    vols = geo.getVolumes()

    # set top edges to bottom edges
    edges[3].setSource(edges[1])
    edges[11].setSource(edges[10])
    edges[9].setSource(edges[7])
    edges[6].setSource(edges[4])
    # set side edges left to right
    edges[2].setSource(edges[0])
    edges[8].setSource(edges[0])
    edges[5].setSource(edges[2])
    # set top face to bottom face
    faces[4].setSource(vols[0], faces[5])

    geo = TMR.Model(verts, edges, faces, vols)

    return geo

# Set the communicator
comm = MPI.COMM_WORLD

# Create the geometry
geo = initGeo()

# Mark the boundary condition faces
faces = geo.getFaces()
volumes = geo.getVolumes()
faces[0].setAttribute('fixed') # fix root
faces[5].setAttribute('surface') # top face

# Set the boundary conditions for the problem
bcs = TMR.BoundaryConditions()
bcs.addBoundaryCondition('fixed')

# Create the geometry    
mesh = TMR.Mesh(comm, geo)

# # Mesh the wing box
opts = TMR.MeshOptions()
opts.write_mesh_quality_histogram = 1
#opts.num_smoothing_steps = 10

htarget = 0.3
mesh.mesh(htarget, opts=opts)
mesh.writeToVTK('unswept-mesh.vtk', 'hex')
# Create a model from the mesh
model = mesh.createModelFromMesh()

# Create the corresponding mesh topology from the mesh-model 
topo = TMR.Topology(comm, model)

# Create the quad forest and set the topology of the forest
forest = TMR.OctForest(comm)
forest.setTopology(topo)

# Create the trees, rebalance the elements and repartition
nlevels = 2
order = 2
forest.createTrees(nlevels-1)
assembler, problem = createTopoProblem(forest,
                                       nlevels=nlevels)
aux = addFaceTraction(order, forest, 'surface', assembler,
                      [0.0, 0.0, 1000.0])

assembler.zeroVariables()
force = assembler.createVec()
assembler.setAuxElements(aux)
assembler.assembleRes(force)
force.scale(-1.0)
forces = [force]

problem.setLoadCases(forces)
funcs = [functions.StructuralMass(assembler)]
initial_mass = assembler.evalFunctions(funcs)
m_fixed = 0.15*initial_mass/0.95
print m_fixed
problem.addConstraints(0, funcs, [-m_fixed], [-1.0/m_fixed])
problem.setObjective([1.0])
problem.initialize()
problem.setPrefix('./results_unswept/')

max_bfgs = 20
opt = ParOpt.pyParOpt(problem, max_bfgs, ParOpt.BFGS)
opt.setOutputFrequency(1)
opt.setOutputFile('results_unswept/opt.dat')
opt.checkGradients(1e-6)
#opt.optimize()

# Output for visualization
flag = (TACS.ToFH5.NODES |
        TACS.ToFH5.DISPLACEMENTS |
        TACS.ToFH5.STRAINS |
        TACS.ToFH5.STRESSES |
        TACS.ToFH5.EXTRAS)
f5 = TACS.ToFH5(assembler, TACS.PY_SOLID, flag)
f5.writeToFile('unswept.f5')
