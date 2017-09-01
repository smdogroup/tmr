from mpi4py import MPI
from tmr import TMR
from paropt import ParOpt
from tacs import TACS, elements, constitutive, functions
import numpy as np
import argparse
import os

class CreateMe(TMR.OctTopoCreator):
    def __init__(self, bcs, filt):
        TMR.OctTopoCreator.__init__(bcs, filt)

    def createElement(self, order, octant, index, weights):
        '''Create the element'''
        rho = 2600.0
        E = 70e9
        nu = 0.3
        stiff = TMR.OctStiffness(rho, E, nu, index, weights, q=8.0)
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

    for i in xrange(nlevels):
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

comm = MPI.COMM_WORLD

# Load the geometry model
geo = TMR.LoadModel('bracket_solid.stp')

# Mark the boundary condition faces
faces = geo.getFaces()
volumes = geo.getVolumes()
faces[15].setAttribute('fixed')
faces[4].setSource(volumes[0], faces[1])
faces[4].setAttribute('surface')

# Set the boundary conditions for the problem
bcs = TMR.BoundaryConditions()
bcs.addBoundaryCondition('fixed')

# Create the mesh
mesh = TMR.Mesh(comm, geo)

# Set the meshing options
opts = TMR.MeshOptions()
opts.frontal_quality_factor = 1.25
opts.num_smoothing_steps = 10
opts.write_mesh_quality_histogram = 1

# Create the surface mesh
htarget = 4.0
mesh.mesh(htarget, opts)

# Write the surface mesh to a file
mesh.writeToVTK('bracket_mesh.vtk')

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
forest.createTrees(nlevels)
assembler, problem = createTopoProblem(forest,
                                       nlevels=nlevels)
aux1 = addFaceTraction(order, forest, 'surface', assembler,
                       [1.0, 1.0, 1.0])
aux2 = addFaceTraction(order, forest, 'surface', assembler,
                       [0.0, 0.0, 1.0])

assembler.zeroVariables()
force1 = assembler.createVec()
assembler.setAuxElements(aux1)
assembler.assembleRes(force1)
force1.scale(-1.0)

assembler.zeroVariables()
force2 = assembler.createVec()
assembler.setAuxElements(aux2)
assembler.assembleRes(force2)
force2.scale(-1.0)

# Set the load cases
forces = [force1, force2]
problem.setLoadCases(forces)

# Compute the volume of the bracket
r = 7.5
a = 50.0
t = 25.0
vol = (3*a*a*t - 3*np.pi*r*r*t)*2600
vol_frac = 0.2

# Set the constraints
funcs = [functions.StructuralMass(assembler)]
initial_mass = assembler.evalFunctions(funcs)
m_fixed =  vol_frac*vol

problem.addConstraints(0, funcs, [-m_fixed], [-1.0/m_fixed])
problem.addConstraints(1, [], [], [])
problem.setObjective([1.0, 1.0])

# Initialize the problem and set the prefix
problem.initialize()
problem.setPrefix('./results/')

max_bfgs = 20
opt = ParOpt.pyParOpt(problem, max_bfgs, ParOpt.BFGS)
opt.setOutputFrequency(1)
opt.setOutputFile("paropt_output.out")
opt.optimize()

print assembler.evalFunctions(funcs)/initial_mass

# Output for visualization
flag = (TACS.ToFH5.NODES |
        TACS.ToFH5.DISPLACEMENTS |
        TACS.ToFH5.STRAINS |
        TACS.ToFH5.STRESSES |
        TACS.ToFH5.EXTRAS)
f5 = TACS.ToFH5(assembler, TACS.PY_SOLID, flag)
f5.writeToFile('bracket.f5')
