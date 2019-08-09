"""
Compliance minimization with mass and frequency constraints

This example demonstrates:

1) Creating meshes using TMR.Creator methods
2) Setting up a topology optimization problem with frequency constraints
3) Design-feature based AMR
4) Use of a Lagrange filter

For a full list of arguments type:

python frequency.py --help
"""

from __future__ import print_function
from mpi4py import MPI
from tmr import TMR, TopOptUtils
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
    """
    Create an initial forest for analysis and optimization.

    This code loads in the model, sets names, meshes the geometry and creates
    an OctForest from the mesh. The forest is populated with octrees with
    the specified depth.

    Args:
        comm (MPI_Comm): MPI communicator
        depth (int): Depth of the initial trees
        htarget (float): Target global element mesh size

    Returns:
        OctForest: Initial forest for topology optimization
    """
    # Load the geometry model
    geo = TMR.LoadModel('../cantilever/cantilever.stp')

    # Mark the boundary condition faces
    verts = geo.getVertices()
    edges = geo.getEdges()
    faces = geo.getFaces()
    volumes = geo.getVolumes()

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

def create_problem(forest, bcs, props, nlevels):
    """
    Create a TopoProblem instance for mass and frequency constrained compliance
    minimization.

    This problem takes in a forest at the current refinement level, boundary
    condition information storing the names of the geometric entities to apply
    Dirichlet boundary conditions, the material properties and the number of
    multigrid levels.

    Args:
        forest (OctForest): Forest object
        bcs (BoundaryConditions): Boundary condition object
        props (StiffnessProperties): Material properties object
        nlevels (int): number of multigrid levels

    Returns:
        TopoProblem: Topology optimization problem instance
    """

    # Allocate the creator callback function
    obj = CreatorCallback(bcs, props)

    # Define the filter type
    filter_type = 'conform'

    # Create the problem and filter objects
    problem = TopOptUtils.createTopoProblem(forest, obj.creator_callback,
        filter_type, nlevels=nlevels, lowest_order=2)

    # Get the assembler object we just created
    assembler = problem.getAssembler()

    # Set the load
    P = 1.0e3
    force = TopOptUtils.computeVertexLoad('pt1', forest, assembler, [0, P, 0])
    temp = TopOptUtils.computeVertexLoad('pt2', forest, assembler, [0, 0, P])
    force.axpy(1.0, temp)

    problem.setLoadCases([force])

    # Compute the fixed mass target
    lx = 50.0 # mm
    ly = 10.0 # mm
    lz = 10.0 # mm
    vol = lx*ly*lz
    vol_frac = 0.25
    density = 2600.0
    m_fixed = vol_frac*(vol*density)

    # Add the fixed mass constraint
    # (m_fixed - m(x))/m_fixed >= 0.0
    funcs = [functions.StructuralMass(assembler)]
    problem.addConstraints(0, funcs, [-m_fixed], [-1.0/m_fixed])

    # Add the natural frequency constraint
    omega_min = (0.5/np.pi)*(2e4)**0.5
    freq_opts = {'use_jd':args.use_jd, 'num_eigs':args.num_eigs,
                 'num_recycle':args.num_recycle, 'track_eigen_iters':nlevels}
    TopOptUtils.addNaturalFrequencyConstraint(problem, omega_min, **freq_opts)

    # Set the compliance objective
    problem.setObjective(obj_array)

    # Initialize the problem and set the prefix
    problem.initialize()

    return problem

# Create an argument parser to read in arguments from the commnad line
p = argparse.ArgumentParser()
# Output options
p.add_argument('--prefix', type=str, default='./results')
p.add_argument('--output_freq', type=int, default=1)

# Mesh parameters
p.add_argument('--htarget', type=float, default=5.0)
p.add_argument('--init_depth', type=int, default=4)
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
p.add_argument('--max_opt_iters', type=int, nargs='+', default=[250])
p.add_argument('--opt_abs_tol', type=float, default=1e-6)
p.add_argument('--opt_barrier_frac', type=float, default=0.25)
p.add_argument('--opt_barrier_power', type=float, default=1.0)
p.add_argument('--qn_subspace', type=int, default=2)
p.add_argument('--tr_penalty', type=float, default=15.0)

# Optimization parameters
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
    'maxiter': 50,
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

# Create the stiffness properties object
props = TMR.StiffnessProperties(rho, E, nu, k0=1e-3, eps=0.2, q=5.0)

# Set the values of the objective array
obj_array = [ 1.0e2 ]

# Set the boundary conditions for the problem
bcs = TMR.BoundaryConditions()
bcs.addBoundaryCondition('fixed', [0, 1, 2], [0.0, 0.0, 0.0])

# Create the initial forest
forest = create_forest(comm, args.init_depth, args.htarget)

# Set the original filter to NULL
orig_filter = None
xopt = None

# Start the optimization
for step in range(max_iterations):
    # Create the TopoProblem instance
    nlevels = mg_levels[step]
    problem = create_problem(forest, bcs, props, nlevels)
    problem.setPrefix(args.prefix)

    # Extract the filter to interpolate design variables
    filtr = problem.getFilter()

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

    # Refine based solely on the value of the density variable
    assembler = problem.getAssembler()
    lower = 0.05
    upper = 0.5
    TopOptUtils.densityBasedRefine(forest, assembler, lower=lower, upper=upper)

    # Repartition the mesh
    forest.repartition()

    # Output for visualization
    flag = (TACS.ToFH5.NODES | TACS.ToFH5.EXTRAS)
    f5 = TACS.ToFH5(assembler, TACS.PY_SOLID, flag)
    f5.writeToFile('beam{0}.f5'.format(step))
