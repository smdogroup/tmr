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
        stiff = TMR.OctStiffness(rho, E, nu, index, weights, q=5.0)
        elem = elements.Solid(2, stiff)
        return elem

comm = MPI.COMM_WORLD

# Load the geometry model
geo = TMR.LoadModel('bracket_solid.stp')

# Mark the boundary condition faces
faces = geo.getFaces()
volumes = geo.getVolumes()
faces[15].setAttribute('fixed')
faces[4].setSource(volumes[0], faces[1])

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

# Create the forest
forests = []
filters = []
assemblers = []
varmaps = []
vecindices = []

# Create the trees, rebalance the elements and repartition
nlevels = 1
order = 2
forest.createTrees(nlevels)
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

force = assemblers[0].createVec()
force.getArray()[::] = 1.0
assemblers[0].applyBCs(force)
forces = [force]

problem.setLoadCases(forces)
funcs = [functions.StructuralMass(assemblers[0])]
m_fixed = 10.0
problem.addConstraints(0, funcs, [-m_fixed], [1.0])
problem.setObjective([1.0])
problem.initialize()
problem.setPrefix('')

max_bfgs = 20
opt = ParOpt.pyParOpt(problem, max_bfgs, ParOpt.BFGS)
opt.optimize()
