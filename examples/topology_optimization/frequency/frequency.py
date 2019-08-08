from __future__ import print_function
from mpi4py import MPI
from tmr import TMR
from paropt import ParOpt
from tacs import TACS, elements, constitutive, functions
import numpy as np
import argparse
import os

class OctCreator(TMR.OctTopoCreator):
    def __init__(self, bcs, filt, props):
        TMR.OctTopoCreator.__init__(bcs, filt)
        self.props = props

    def createElement(self, order, octant, index, weights):
        '''Create the element'''
        stiff = TMR.OctStiffness(self.props, index, weights)
        elem = elements.Solid(2, stiff)
        return elem

class CreatorCallback:
    def __init__(self, bcs, props):
        self.bcs = bcs
        self.props = props

    def creator_callback(self, forest):
        filtr = forest.duplicate()
        filtr.coarsen()
        creator = OctCreator(self.bcs, filtr, self.props)
        return creator, filtr

def create_forest(comm, depth, htarget):
    # Load the geometry model
    geo = TMR.LoadModel('beam.stp')

    # Mark the boundary condition faces
    verts = geo.getVertices()
    edges = geo.getEdges()
    faces = geo.getFaces()

    faces[3].setName('fixed')
    faces[4].setSource(volumes[0], faces[5])
    verts[4].setName('pt1')
    verts[3].setName('pt2')

    # Set the boundary conditions for the problem
    bcs = TMR.BoundaryConditions()
    bcs.addBoundaryCondition('fixed')
    
    # Create the mesh
    mesh = TMR.Mesh(comm, geo)
    
    # Set the meshing options
    opts = TMR.MeshOptions()
    
    # Create the surface mesh
    mesh.mesh(args.htarget, opts)
    
    # Create a model from the mesh
    model = mesh.createModelFromMesh()
    
    # Create the corresponding mesh topology from the mesh-model 
    topo = TMR.Topology(comm, model)
    
    # Create the quad forest and set the topology of the forest
    forest = TMR.OctForest(comm)
    forest.setTopology(topo)

    # Create the trees
    forest.createTrees(depth)
    
    return forest
    
# Create an argument parser to read in arguments from the commnad line
p = argparse.ArgumentParser()
# Output options
p.add_argument('--prefix', type=str, default='./results')
p.add_argument('--output_freq', type=int, default=1)

# Mesh parameters
p.add_argument('--htarget', type=float, default=4.0)
p.add_argument('--init_depth', type=int, default=1)
p.add_argument('--order', type=int, default=2)
p.add_argument('--refine_flag', type=bool, nargs='+', default=[True])

# Solver parameters
p.add_argument('--mg_levels', type=int, nargs='+', default=[2])
p.add_argument('--use_decrease_order', action='store_true', default=True)

# Frequency constraint parameters
p.add_argument('--num_eigs', type=int, default=9)
p.add_argument('--num_recycle', type=int, default=9)
p.add_argument('--use_jd', action='store_true', default=False)

# Optimizer settings
p.add_argument('--max_opt_iters', type=int, nargs='+',
               default=[250])
p.add_argument('--opt_abs_tol', type=float, default=1e-6)
p.add_argument('--opt_barrier_frac', type=float, default=0.25)
p.add_argument('--opt_barrier_power', type=float, default=1.0)
p.add_argument('--qn_subspace', type=int, default=2)
p.add_argument('--tr_penalty', type=float, default=15.0)

# Optimization parameters
p.add_argument('--vol_frac', type=float, default=0.15)
p.add_argument('--q_penalty', type=float, nargs='+', default=[8.0])
p.add_argument('--ks_weight', type=float, nargs='+', default=[25.0])
p.add_argument('--eps', type=float, nargs='+', default=[0.1])

args = p.parse_args()

# Set the communicator
comm = MPI.COMM_WORLD

# Print out all of the arguments to the command line
if comm.rank == 0:
    for arg in vars(args):
        print('%-20s'%(arg), getattr(args, arg))

# Set the optimization parameters
optimization_options = {
    # Parameters for the trust region method
    'tr_init_size': 0.01,
    'tr_max_size': 0.05,
    'tr_min_size': 1e-6,
    'tr_eta': 0.25,
    'tr_penalty_gamma': args.tr_penalty,

    # Parameters for the interior point method (used to solve the
    # trust region subproblem)
    'max_qn_subspace': args.qn_subspace,
    'tol': 1e-8,
    'maxiter': 500,
    'norm_type': 'L1',
    'barrier_strategy': 'Complementarity fraction',
    'start_strategy': 'Affine step'}

prefix = args.prefix
optimization_options['output_file'] = os.path.join(prefix, 'output_file.dat')
optimization_options['tr_output_file'] = os.path.join(prefix, 'tr_output_file.dat')

# Set the max nmber of outer iterations
# (for continuation strategies and/or mesh refinement)
mg_levels = args.mg_levels
max_iterations = len(mg_levels)

# Set the material properties
rho = [2600.0]
E = [70e9]
nu = [0.3]

# Set the number of variables per node
vars_per_node = 1
if (len(rho) > 1):
    vars_per_node = 1+len(rho)

# Create the stiffness properties object
props = TMR.StiffnessProperties(rho, E, nu, k0=1e-3,
                                eps=0.2, q=5.0)

# Compute the volume of the bracket
r = 10.0
a = 50.0
vol = r*r*a
initial_mass = vol*np.average(rho)
m_fixed = args.vol_frac*initial_mass

# Set the variables that we want to output for visualization
flag = (TACS.ToFH5.NODES |
        TACS.ToFH5.EXTRAS)

# Set the values of the objective array
obj_array = [ 1.0e2 ]

# Set the original filter to NULL
orig_filter = None
xopt = None

# Start the optimization
for step in range(max_iterations):
    nlevs = mg_levels[step]

    # Allocate the creator callback function
    obj = CreatorCallback(bcs, props)

    # Define the filter type
    filter_type = 'conform'

    # Create the problem and filter objects
    problem = TopOptUtils.createTopoProblem(forest,
                                            obj.creator_callback, filter_type, nlevels=nlevels, lowest_order=2)
    
    # Extract the filter to interpolate design variables
    filtr = problem.getFilter()

    # Get the assembler object we just created
    assembler = problem.getAssembler()

    # Create the load vectors, combine them, and set them
    T = 1e3
    force1 = computeVertexLoad('pt1', forest, assembler, [0.0, -T, 0.0])
    force2 = computeVertexLoad('pt2', forest, assembler, [0.0, 0.0, -T])
    force1.axpy(1.0, force2)
    assembler.reorderVec(force1)
    problem.setLoadCases([force1])
    
    # Add the fixed mass constraint
    # (m_fixed - m(x))/m_fixed >= 0.0
    funcs = [functions.StructuralMass(assembler)]
    problem.addConstraints(0, funcs, [-m_fixed], [-1.0/m_fixed])
    
    # Add the natural frequency constraint
    omega_min = (0.5/np.pi)*(2e4)**0.5
    freq_opts = {'use_jd':args.use_jd, 'num_eigs':args.num_eigs,
                 'num_recycle':args.num_recycle, 'track_eigen_iters':nlevs}
    TopOptUtils.addNaturalFrequencyConstraint(problem, omega_min, freq_opts)

    # Set the compliance objective
    problem.setObjective(obj_array)

    # Initialize the problem and set the prefix
    problem.initialize()
    problem.setPrefix(args.prefix)

    problem.setIterationCounter(sum(args.max_opt_iters[:step]))

    if orig_filter is not None:
        # Create one of the new design vectors
        x = problem.createDesignVec()
        TopOptUtils.interpolateDesignVec(orig_filter, xopt, filtr, x)
        problem.setInitDesignVars(x)

    # Set the new original filter
    orig_filter = filtr

    # Set the max number of iterations
    optimization_options['maxiter'] = args.max_opt_iters[step]
    
    # Optimize the problem
    opt = TopOptUtils.TopologyOptimizer(problem, optimization_options)
    xopt = opt.optimize()
    
    # Create refinement array
    num_elems = assembler.getNumElements()
    refine = np.zeros(num_elems, dtype=np.int32)

    # Refine based solely on the value of the density variable
    elems = assembler.getElements()
    
    for i in range(num_elems):        
        c = elems[i].getConstitutive()
        if c is not None:
            if vars_per_node == 1:
                density = c.getDVOutputValue(0, np.zeros(3, dtype=float))
            else:
                density = 1.0 - c.getDVOutputValue(2, np.zeros(3, dtype=float))

            # Refine things differently depending on whether the
            # density is above or below a threshold
            if density >= 0.5:
                refine[i] = int(1)
            elif density < 0.05:
                refine[i] = int(-1)
    
    # Refine the forest
    forest.refine(refine, min_lev=0)

    # Output for visualization
    f5 = TACS.ToFH5(assembler, TACS.PY_SOLID, flag)
    f5.writeToFile('beam{0}.f5'.format(step))
