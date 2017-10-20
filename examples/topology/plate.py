from mpi4py import MPI
from tmr import TMR
from paropt import ParOpt
from tacs import TACS, elements, constitutive, functions
import numpy as np
import argparse
import os

class CreateMe(TMR.QuadTopoCreator):
    def __init__(self, bcs, filt):
        TMR.QuadTopoCreator.__init__(bcs, filt)

    def createElement(self, order, quadrant, index, weights):
        '''Create the element'''
        rho = 2600.0
        E = 70e9
        nu = 0.3
        stiff = TMR.QuadStiffness(rho, E, nu, index, weights, q=5.0)
        elem = elements.PlaneQuad(2, stiff)
        return elem

def addVertexLoad(comm, order, forest, attr, assembler, F):
    # Retrieve octants from the forest
    quadrants = forest.getQuadrants()
    node_octs = forest.getNodesWithAttribute(attr)
    force = assembler.createVec()
    f_array = force.getArray()
    node_range = forest.getNodeRange()
    mpi_rank = comm.Get_rank()
    for i in range(len(node_octs)):
        if (node_octs[i].tag >= node_range[mpi_rank]) and \
               (node_octs[i].tag < node_range[mpi_rank+1]): 
            index = node_octs[i].tag-node_range[mpi_rank]
            
            f_array[2*index] -= F[0]
            f_array[2*index+1] -= F[1]
            
    return force

def addEdgeTraction(order, forest, attr, assembler, tr):
    trac = []
    
    for findex in range(4):
        trac.append(elements.PSQuadTraction(findex, tr[0]*np.ones(2), \
                                            tr[1]*np.ones(2)))

    # Retrieve octants from the forest
    quadrants = forest.getQuadrants()
    edge_quads = forest.getQuadsWithAttribute(attr)
    aux = TACS.AuxElements()

    for i in range(len(edge_quads)):
        index = quadrants.findIndex(edge_quads[i])
        if index is not None:
            aux.addElement(index, trac[edge_quads[i].tag])

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

# Create an argument parser to read in arguments from the commnad line
p = argparse.ArgumentParser()
p.add_argument('--prefix', type=str, default='./results')
p.add_argument('--vol_frac', type=float, default=0.15)
p.add_argument('--htarget', type=float, default=1.5)
p.add_argument('--max_opt_iters', type=int, default=250)
p.add_argument('--opt_abs_tol', type=float, default=1e-6)
p.add_argument('--opt_barrier_frac', type=float, default=0.25)
p.add_argument('--opt_barrier_power', type=float, default=1.0)
p.add_argument('--output_freq', type=int, default=1)
p.add_argument('--init_depth', type=int, default=2)
p.add_argument('--mg_levels', type=int, nargs='+', default=[2, 2, 3, 3])
p.add_argument('--max_lbfgs', type=int, default=10)
p.add_argument('--hessian_reset', type=int, default=10)
args = p.parse_args()

# Set the parameter to use paropt or MMA
use_paropt = False

# The communicator
comm = MPI.COMM_WORLD

# Load the geometry model
geo = TMR.LoadModel('beam2d.stp')

# Mark the boundary condition faces
verts = geo.getVertices()
edges = geo.getEdges()
volumes = geo.getVolumes()
edges[1].setAttribute('fixed')
verts[0].setAttribute('pt1')

# Set the boundary conditions for the problem
bcs = TMR.BoundaryConditions()
bcs.addBoundaryCondition('fixed')

# Create the mesh
mesh = TMR.Mesh(comm, geo)

# Set the meshing options
opts = TMR.MeshOptions()
opts.write_mesh_quality_histogram=1

# Create the surface mesh
mesh.mesh(args.htarget, opts)

# Create a model from the mesh
model = mesh.createModelFromMesh()

# Create the corresponding mesh topology from the mesh-model 
topo = TMR.Topology(comm, model)

# Create the quad forest and set the topology of the forest
forest = TMR.QuadForest(comm)
forest.setTopology(topo)

# Compute the volume of the plate
r = 10.0
a = 30.0
vol = r*a
vol_frac = args.vol_frac

# Set the fixed mass
density = 2600.0
initial_mass = vol*density
m_fixed = vol_frac*initial_mass

# Set the max nmber of iterations
mg_levels = args.mg_levels
max_iterations = 2#len(mg_levels)

# Set parameters for later usage
order = 2
forest.createTrees(args.init_depth)

# The old filter/map classes
old_filtr = None
old_varmap = None

# Values from the previous iteration
old_dvs = None
old_z = 0.0
olz_zl = None
old_zu = None

# The volumes of the elements in the filter mesh
filtr_volumes = None

# Set the values of the objective array
obj_array = [1.0e2]

for ite in xrange(max_iterations):
    # Create the TACSAssembler and TMRTopoProblem instance
    nlevs = mg_levels[ite]
    assembler, problem, filtr, varmap = createTopoProblem(forest, 
                                                          nlevels=nlevs)
   
    # Write out just the mesh - for visualization
    flag = (TACS.ToFH5.NODES)

    f5 = TACS.ToFH5(assembler, TACS.PY_PLANE_STRESS, flag)
    f5.writeToFile(os.path.join(args.prefix, 'plate_mesh%d.f5'%(ite)))
    
    # Set the constraint type
    funcs = [functions.StructuralMass(assembler)]
    # Add the point loads to the vertices
    force1 = addVertexLoad(comm, order, forest, 'pt1', assembler,
                           [0.0, -1000.])

    force1.scale(-1.0)    
    # Set the load cases
    forces = [force1]
    problem.setLoadCases(forces)

    # Set the mass constraint
    problem.addConstraints(0, funcs, [-m_fixed], [-1.0/m_fixed])
    problem.setObjective(obj_array)
    
    # Initialize the problem and set the prefix
    problem.initialize()
    problem.setPrefix(args.prefix)

    if use_paropt:
        # Create the ParOpt problem
        opt = ParOpt.pyParOpt(problem, args.max_lbfgs, ParOpt.BFGS)

        # Set the norm type to use
        opt.setNormType(ParOpt.L1_NORM)

        # Set parameters
        opt.setMaxMajorIterations(args.max_opt_iters)
        opt.setHessianResetFreq(args.hessian_reset)
        problem.setIterationCounter(args.max_opt_iters*ite)
        opt.setAbsOptimalityTol(args.opt_abs_tol)
        opt.setBarrierFraction(args.opt_barrier_frac)
        opt.setBarrierPower(args.opt_barrier_power)
        opt.setOutputFrequency(args.output_freq)
        opt.setOutputFile(os.path.join(args.prefix, 
                                       'paropt_output%d.out'%(ite)))
        
        # If the old filter exists, we're on the second iteration
        if old_filtr:
            # Create the interpolation
            interp = TACS.VecInterp(old_varmap, varmap)
            filtr.createInterpolation(old_filtr, interp)
            interp.initialize()
        
            # Get the optimization variables for the new optimizer
            x, z, zw, zl, zu = opt.getOptimizedPoint()

            # Do not try to estimate the new multipliers
            opt.setInitStartingPoint(0)
        
            # Set the values of the new mass constraint multipliers
            z[0] = old_z

            # Do the interpolation
            old_vec = problem.convertPVecToVec(old_dvs)
            x_vec = problem.convertPVecToVec(x)
            interp.mult(old_vec, x_vec)

            # Set the new design variables
            problem.setInitDesignVars(x)

            # Compute the new filter areas
            filtr_areas = problem.createAreaVec()
            areas = filtr_areas.getArray()
            # Do the interpolation of the multipliers
            old_vec = problem.convertPVecToVec(old_zl)
            zl_vec = problem.convertPVecToVec(zl)
            interp.mult(old_vec, zl_vec)
            zl[:] *= areas

            # Do the interpolation
            old_vec = problem.convertPVecToVec(old_zu)
            zu_vec = problem.convertPVecToVec(zu)
            interp.mult(old_vec, zu_vec)
            zu[:] *= areas
        
            # Reset the complementarity 
            new_barrier = opt.getComplementarity()
            opt.setInitBarrierParameter(new_barrier)
        else:
            filtr_areas = problem.createAreaVec()

        # Optimize the new point
        opt.optimize()
    
        # Get the final values of the design variables
        x, z, zw, zl, zu = opt.getOptimizedPoint()
    
        # Set the old design variable vector to what was once the new
        # design variable values. Copy the values from
        # getOptimizedPoint
        old_z = z[0]
        old_dvs = x
        old_zl = zl
        old_zu = zu
        
        # Divide the bound multipliers by their associated volumes
        areas = filtr_areas.getArray()
        zl[:] = zl[:]/areas
        zu[:] = zu[:]/areas
    else:
        # Here we use MMA, and only worry about interpolating the
        # variables
        max_mma_iters = 5

        # Set the ParOpt problem into MMA
        mma = ParOpt.pyMMA(problem)
        mma.setPrintLevel(2)
        mma.setOutputFile(os.path.join(args.prefix, 
                                       'mma_output%d.out'%(ite)))

        # Create the ParOpt problem
        opt = ParOpt.pyParOpt(mma, args.max_lbfgs,
                              ParOpt.NO_HESSIAN_APPROX)

        # Set parameters
        opt.setMaxMajorIterations(args.max_opt_iters)
        opt.setHessianResetFreq(args.hessian_reset)
        opt.setAbsOptimalityTol(args.opt_abs_tol)
        opt.setBarrierFraction(args.opt_barrier_frac)
        opt.setBarrierPower(args.opt_barrier_power)
        opt.setOutputFrequency(args.output_freq)
        opt.setUseDiagHessian(1)
        opt.setMajorIterStepCheck(1)

        # Set the output file name
        filename = os.path.join(args.prefix, 'paropt%04d.out'%(ite))
        opt.setOutputFile(filename)

        if ite == 0:
            # Set the starting point using the mass fraction
            x = mma.getOptimizedPoint()
            x[:] = 0.5
            problem.setInitDesignVars(x)

        # If the old filter exists, we're on the second iteration
        if old_filtr:
            # Create the interpolation
            interp = TACS.VecInterp(old_varmap, varmap)
            filtr.createInterpolation(old_filtr, interp)
            interp.initialize()
        
            # Get the optimization variables for the new optimizer
            x = mma.getOptimizedPoint()

            # Do the interpolation
            old_vec = problem.convertPVecToVec(old_dvs)
            x_vec = problem.convertPVecToVec(x)
            interp.mult(old_vec, x_vec)

            # Set the initial design variable values
            problem.setInitDesignVars(x)
            
        # Initialize the subproblem
        mma.initializeSubProblem()
        opt.resetDesignAndBounds()
        
        # Enter the optimization loop
        for i in range (max_mma_iters):
            opt.setInitBarrierParameter(10.0)
            opt.optimize()

            # Write the solution out to a file
            if i % args.output_freq == 0:
                itr = max_mma_iters*ite + i
                filename = os.path.join(args.prefix, 
                                        'beam_topo%04d.f5'%(itr))
                flag = (TACS.ToFH5.NODES |
                        TACS.ToFH5.EXTRAS)
                f5 = TACS.ToFH5(assembler, TACS.PY_PLANE_STRESS, flag)
                f5.writeToFile(filename)

            opt.checkGradients(1e-7)
            
            # Get the optimized point
            x, z, zw, zl, zu = opt.getOptimizedPoint()           
            mma.setMultipliers(z, zw)
            mma.initializeSubProblem(x)
            opt.resetDesignAndBounds()
            
            # Compute the KKT error
            l1_norm, linfty_norm, infeas = mma.computeKKTError()
            if comm.rank == 0:
                print 'z = ', z
                print 'l1_norm = ', l1_norm
                print 'infeas = ', infeas

            if l1_norm < 1e-3 and infeas < 1e-6:
                break

        # Set the old values of the variables
        old_dvs = mma.getOptimizedPoint()
            
    # Set the old filter/variable map for the next time through the
    # loop so that we can interpolate design variable values
    old_varmap = varmap
    old_filtr = filtr

    # Output for visualization
    flag = (TACS.ToFH5.NODES |
            TACS.ToFH5.EXTRAS)
    f5 = TACS.ToFH5(assembler, TACS.PY_PLANE_STRESS, flag)
    f5.writeToFile(os.path.join(args.prefix, 'tacs_output%d.f5'%(ite)))
    # Create refinement array
    num_elems = assembler.getNumElements()
    refine = np.ones(num_elems, dtype=np.int32)
    # Refine based solely on the value of the density variable
    elems = assembler.getElements()
    for i in xrange(num_elems):        
        c = elems[i].getConstitutive()
        if c is not None:
            rho = c.getDVOutputValue(0, np.zeros(2, dtype=float))
            # Refine things differently depending on whether the
            # density is above or below a threshold
            if rho > 0.25:
                refine[i] = int(1)
            elif rho < 0.05:
                refine[i] = int(-1)

    # Refine the forest
    forest.refine(refine)
   
