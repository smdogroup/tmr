from __future__ import print_function
from mpi4py import MPI
from tmr import TMR
from paropt import ParOpt
from tacs import TACS, elements, constitutive, functions
import numpy as np
import argparse
import os

class QuadConformCreator(TMR.QuadConformTopoCreator):
    def __init__(self, bcs, forest, order, interp, props=None):
        # Set the interpolation for the new filter
        TMR.QuadConformTopoCreator.__init__(bcs, forest, order, interp)
        
        # Create the array of properties
        self.props = props
        
    def createElement(self, order, quadrant, index, filtr):
        '''Create the element'''
        stiff = TMR.ThermoQuadStiffness(self.props,
                                        index, None, filtr)
        elem = elements.PSThermoelasticQuad(order, stiff)
        return elem

class CreatorCallback:
    def __init__(self, bcs, props):
        self.bcs = bcs
        self.props = props

    def creator_callback(self, forest):

        order = forest.getMeshOrder()-1
        interp = TMR.BERNSTEIN_POINTS

        filtr = forest.duplicate()
        filtr.coarsen()
        creator = OctCreator(self.bcs, filter, self.props)
        return creator, filtr

def create_forest(comm, depth, htarget):
    # Load the geometry model
    geo = TMR.LoadModel('biclamped_traction.stp')

    # Mark the boundary condition faces
    verts = geo.getVertices()
    edges = geo.getEdges()
    faces = geo.getFaces()

    edges[1].setName('fixed')
    edges[9].setName('fixed')
    edges[4].setName('traction')

    # Set the boundary conditions for the problem
    bcs = TMR.BoundaryConditions()
    bcs.addBoundaryCondition('fixed', [0,1,2], [0.0,0.0,0.])

    # Create the mesh
    mesh = TMR.Mesh(comm, geo)

    # Set the meshing options
    opts = TMR.MeshOptions()

    # Create the surface mesh
    mesh.mesh(htarget, opts)

    # Create a model from the mesh
    model = mesh.createModelFromMesh()

    # Create the corresponding mesh topology from the mesh-model 
    topo = TMR.Topology(comm, model)

    # Create the quad forest and set the topology of the forest
    forest = TMR.QuadForest(comm)
    forest.setTopology(topo)

    # Set parameters for later usage
    forest.createTrees(depth)

    return forest

# Set the optimization parameters
optimization_options = {
    # Parameters for the trust region method
    'tr_init_size': 0.01,
    'tr_max_size': 0.05,
    'tr_min_size': 1e-6,
    'tr_eta': 0.25,
    'tr_penalty_gamma': 20.0,

    # Parameters for the interior point method (used to solve the
    # trust region subproblem)
    'max_qn_subspace': 2,
    'tol': 1e-8,
    'maxiter': 500,
    'norm_type': 'L1',
    'barrier_strategy': 'Complementarity fraction',
    'start_strategy': 'Affine step'}

prefix = 'results'
optimization_options['output_file'] = os.path.join(prefix, 'output_file.dat')
optimization_options['tr_output_file'] = os.path.join(prefix, 'tr_output_file.dat')

# Create an argument parser to read in arguments from the command line
p = argparse.ArgumentParser()
p.add_argument('--prefix', type=str, default='./results')
p.add_argument('--vol_frac', type=float, default=0.25)
p.add_argument('--htarget', type=float, default=2.5e-3)
p.add_argument('--max_opt_iters', type=int, nargs='+',
               default=[250])
p.add_argument('--opt_abs_tol', type=float, default=1e-8)
p.add_argument('--opt_barrier_frac', type=float, default=0.25)
p.add_argument('--opt_barrier_power', type=float, default=1.0)
p.add_argument('--output_freq', type=int, default=1)
p.add_argument('--init_depth', type=int, default=1)
p.add_argument('--mg_levels', type=int, nargs='+', default=[2])
p.add_argument('--tr_penalty', type=float, default=15.0)
p.add_argument('--qn_subspace', type=int, default=2)
p.add_argument('--use_L1', action='store_true', default=True)
p.add_argument('--use_Linf', action='store_true', default=False)
p.add_argument('--use_decrease_order', action='store_true', default=True)
p.add_argument('--order', type=int, default=4)
p.add_argument('--q_penalty', type=float, default=8.0)
p.add_argument('--omega', type=float, default=1.0)
args = p.parse_args()

# Print out all of the arguments to the command line
if comm.rank == 0:
    for arg in vars(args):
        print('%-20s'%(arg), getattr(args, arg))

# Set the max nmber of iterations
mg_levels = args.mg_levels
max_iterations = len(mg_levels)

# Set the communicator
comm = MPI.COMM_WORLD

# Set the qn subspace size
optimization_options['max_qn_subspace'] = args.qn_subspace

# Compute the volume of the bracket
r = 0.06
a = 0.04
vol = r*a
vol_frac = args.vol_frac

# Set the values of the objective array
obj_array = [ 1e0 ]
thickness = 1.0e-2
rho = [2600*thickness, 1300.*thickness]
E = [70e9*thickness, 35e9*thickness]
nu = [0.3, 0.3]
aT = [23.5e-6, 23.5e-6*0.5]
kcond = [130.0*thickness, 65.0*thickness]
ys = [450e6*thickness, 275e6*thickness]

# Set the fixed mass
max_density = sum(rho)/len(rho)
initial_mass = vol*max_density
m_fixed = vol_frac*initial_mass

# Set the number of variables per node
vars_per_node = 1
if (len(rho) > 1):
    vars_per_node = 1+len(rho)

# Create the stiffness properties object
props = TMR.QuadStiffnessProperties(rho, E, nu, _aT=aT, _kcond=kcond,
                                    k0=1e-3, eps=0.2,
                                    q=args.q_penalty, qtemp=0.0,
                                    qcond=0.0)

# Create the tractions
T = 2.5e6
thermal_tractions = []
tx = np.ones(order)*tr[0]
ty = np.ones(order)*tr[1]
for findex in range(4):
    thermal_tractions.append(elements.PSThermoQuadTraction(findex, tx, ty))


# Create the problem and filter object

problem, filter_obj = TopOptUtils.createTopoProblem(forest,
    obj.creator_callback, filter_type, nlevels=nlevels, lowest_order=3)


time_array = np.zeros(sum(args.max_opt_iters[:]))
t0 = MPI.Wtime()

for step in range(max_iterations):
    # Create the TACSAssembler and TMRTopoProblem instance
    nlevels = mg_levels[step]

    # Allocate the creator callback function
    obj = CreatorCallback(bcs, props)

    filter_type = 'conform'
    args.order
    
    problem, filter_obj = TopOptUtils.createTopoProblem(forest,
        ob.creator_callback, filter_type, nlevels=nlevels, lowest_order=3)

    # Set the constraint type
    funcs = [functions.StructuralMass(assembler)]

    # Get the assembler object we just created
    assembler = filter_obj.getAssembler()

    # Allocate a thermal traction boundary condition
    force1 = TopOptUtils.computeTractionLoad('traction', forest, assembler,
        thermal_tractions)
    
    # Set the load cases
    problem.setLoadCases([force1])

    # Set the mass constraint
    # (m_fixed - m(x))/m_fixed >= 0.0    
    problem.addConstraints(0, funcs, [-m_fixed], [-1.0/m_fixed])
    problem.setObjective(obj_array, [functions.Compliance(assembler)])
    
    # Initialize the problem and set the prefix
    problem.initialize()
    problem.setPrefix(args.prefix)

    optimization_options[maxiter] = rgs.max_opt_iters[step]

    problem.setIterationCounter(sum(args.max_opt_iters[:step]))

    if orig_filter:
        # Create one of the new design vectors
        x = problem.createDesignVec()
        TopOptUtils.interpolateDesignVec(orig_filter, xopt, filtr, x)
        problem.setInitDesignVars(x)

    # Optimize the problem
    opt = TopOptUtils.TopologyOptimizer(problem, optimization_options)
    xopt = opt.optimize()
 
    # Set the old filter/variable map for the next time through the
    # loop so that we can interpolate design variable values
    old_varmap = varmap
    old_filtr = filtr

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
            if density >= 0.1:
                refine[i] = int(1)
            elif density < 0.05:
                refine[i] = int(-1)
    
    # Refine the forest
    forest.refine(refine, min_lev=0)

# Do an averaging of all the values and write to text file
new_array=np.zeros(len(time_array))
comm.Allreduce(time_array, new_array,op=MPI.SUM)
new_array /= comm.size
if comm.rank == 0:
    np.savetxt(os.path.join(args.prefix,'time.out'), new_array, fmt='%1.6e')
