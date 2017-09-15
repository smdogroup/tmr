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

def addVertexLoad(comm, order, forest, attr, assembler, F):
    # Retrieve octants from the forest
    octants = forest.getOctants()
    node_octs = forest.getNodesWithAttribute(attr)
    force = assembler.createVec()
    f_array = force.getArray()
    node_range = forest.getNodeRange()
    mpi_rank = comm.Get_rank()
    for i in range(len(node_octs)):
        if (node_octs[i].tag >= node_range[mpi_rank]) and \
               (node_octs[i].tag < node_range[mpi_rank+1]): 
            index = node_octs[i].tag-node_range[mpi_rank]
            
            f_array[3*index] -= F[0]
            f_array[3*index+1] -= F[1]
            f_array[3*index+2] -= F[2]
    return force

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

    return assemblers[0], problem, filters[0], varmaps[0]

comm = MPI.COMM_WORLD

# Load the geometry model
geo = TMR.LoadModel('beam.stp')
geo.writeModelToTecplot('model.dat')

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
opts.frontal_quality_factor = 1.25
opts.num_smoothing_steps = 10
opts.write_mesh_quality_histogram = 1

# Create the surface mesh
htarget = 4.0
mesh.mesh(htarget, opts)

# Write the surface mesh to a file
mesh.writeToVTK('beam_mesh.vtk')

# Create a model from the mesh
model = mesh.createModelFromMesh()

# Create the corresponding mesh topology from the mesh-model 
topo = TMR.Topology(comm, model)

# Create the quad forest and set the topology of the forest
forest = TMR.OctForest(comm)
forest.setTopology(topo)

# Compute the volume of the bracket
r = 10.0
a = 50.0
vol = r*r*a
vol_frac = 0.15

initial_mass = vol*2600
m_fixed =  vol_frac*initial_mass

# Set the max nmber of iterations
max_iterations = 4

# Set parameters for later usage
min_level = 1
max_level = 9
nlevels = 3
order = 2
forest.createTrees(nlevels)

# The old filter/map classes
old_filtr = None
old_varmap = None

# Values from the previuos iteration
old_dvs = None
old_zl = None
old_zu = None
old_z = 0.0

# Set the values of the objective array
obj_array = [1.0e3]

for ite in xrange(max_iterations):
    # Create the TACSAssembler and TMRTopoProblem instance
    if ite == 1:
        nlevels += 1
    assembler, problem, filtr, varmap = createTopoProblem(forest,
                                                          nlevels=nlevels-1)
    
    # Set the constraint type
    funcs = [functions.StructuralMass(assembler)]
    # Add the point loads to the vertices
    force1 = addVertexLoad(comm,order, forest, 'pt1', assembler,
                           [0.0, -1000., 0.0])
    force2 = addVertexLoad(comm,order, forest, 'pt2', assembler,
                           [0.0, 0.0, -1000.])
    
    force1.axpy(1.0, force2)
    force1.scale(-1.0)
    
    # Set the load cases
    forces = [force1]
    problem.setLoadCases(forces)
    
    # Do the interpolation
    if old_filtr:
        interp = TACS.VecInterp(old_varmap, varmap)
        filtr.createInterpolation(old_filtr, interp)
        interp.initialize()

        # PVec objects
        new_dvs = problem.createDesignVec()
        new_zl = problem.createDesignVec()
        new_zu = problem.createDesignVec()
        
        # Do the interpolation
        new_dvs_vec = problem.convertPVecToVec(new_dvs)
        new_zl_vec = problem.convertPVecToVec(new_zl)
        new_zu_vec = problem.convertPVecToVec(new_zu)

        old_dvs_vec = problem.convertPVecToVec(old_dvs)
        old_zl_vec = problem.convertPVecToVec(old_zl)
        old_zu_vec = problem.convertPVecToVec(old_zu)

        # Interpolate the design variables, lagrange multipliers
        interp.mult(old_dvs_vec, new_dvs_vec)
        interp.mult(old_zl_vec, new_zl_vec)
        interp.mult(old_zu_vec, new_zu_vec)

        # Set the new design variables
        problem.setInitDesignVars(new_dvs)
    
    problem.addConstraints(0, funcs, [-m_fixed], [-1.0/m_fixed])
    problem.setObjective(obj_array)
    
    # Initialize the problem and set the prefix
    problem.initialize()
    problem.setPrefix('./results_mem/')
    
    max_bfgs = 20
    opt = ParOpt.pyParOpt(problem, max_bfgs, ParOpt.BFGS)
    max_opt_iters = 250
    opt.setMaxMajorIterations(max_opt_iters)
    problem.setIterationCounter(max_opt_iters*ite)
    opt.setAbsOptimalityTol(1e-6)
    opt.setOutputFrequency(1)
    opt.setOutputFile('./results_mem/paropt_output%d.out'%(ite))
        
    # Get new barrier parameters
    if ite > 0:
        # Do the interpolation
        x_old, z, zw, zl, zu = opt.getOptimizedPoint()
        
        # Set the values of the new multipliers
        opt.setInitStartingPoint(0)
        z[0] = old_z
        zl.copyValues(new_zl)
        zu.copyValues(new_zu)
        
        # Reset the bounds and input the new barrier parameter
        new_barrier = opt.getComplementarity()
        opt.setInitBarrierParameter(new_barrier)

    # Optimize the new point
    opt.optimize()
    
    # Get the final values of the design variables
    x_old, z, zw, zl, zu = opt.getOptimizedPoint()
    
    # Set the old design variable vector to what was once the new
    # design variable values. Copy the values from getOptimizedPoint
    old_z = z[0]
    old_dvs = x_old
    old_zl = zl
    old_zu = zu
    
    # Set the old filter/variable map for the next time through the
    # loop so that we can interpolate design variable values
    old_varmap = varmap
    old_filtr = filtr

    # Output for visualization
    flag = (TACS.ToFH5.NODES |
            TACS.ToFH5.EXTRAS)
    f5 = TACS.ToFH5(assembler, TACS.PY_SOLID, flag)
    f5.writeToFile('./results_mem/beam%d.f5'%(ite))

    # Create refinement array
    num_elems = assembler.getNumElements()
    refine = np.ones(num_elems, dtype=np.int32)

    # Refine based solely on the value of the density variable
    elems = assembler.getElements()
    
    for i in xrange(num_elems):        
        c = elems[i].getConstitutive()
        if c is not None:
            rho = c.getDVOutputValue(0, np.zeros(3,dtype=float))
            # Refine things differently depending on whether
            # the density is above or below a threshold
            if rho > 0.5:
                refine[i] = int(1)
            elif rho < 0.05:
                refine[i] = int(-1)
                
    # Refine the forest
    forest.refine(refine)

