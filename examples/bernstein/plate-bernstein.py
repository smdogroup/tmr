from __future__ import print_function
from mpi4py import MPI
from tmr import TMR
from paropt import ParOpt
from tacs import TACS, elements, constitutive, functions
import numpy as np
import argparse
import os

class CreateMe(TMR.QuadConformTopoCreator):
    def __init__(self, bcs, forest, order, interp, props=None):
        # Set the interpolation for the new filter
        order = forest.getMeshOrder()-1
        interp = TMR.BERNSTEIN_POINTS
        TMR.QuadConformTopoCreator.__init__(bcs, forest, order, interp)
        
        # Create the array of properties
        self.props = props
        
    def createElement(self, order, quadrant, index, filtr):
        '''Create the element'''
        stiff = TMR.ThermoQuadStiffness(self.props,
                                        index, None, filtr)
        elem = elements.PSThermoelasticQuad(order, stiff)
        return elem
    
def addVertexLoad(comm, forest, attr, assembler, F):
    # Retrieve quadrant from the forest
    quadrants = forest.getQuadrants()
    node_quads = forest.getNodesWithName(attr)
    force = assembler.createVec()
    f_array = force.getArray()
    node_range = forest.getNodeRange()
    mpi_rank = comm.Get_rank()
    for i in range(len(node_quads)):
        if ((node_quads[i] >= node_range[mpi_rank]) and
            (node_quads[i] < node_range[mpi_rank+1])): 
            index = node_quads[i] - node_range[mpi_rank]
            f_array[3*index] -= F[0]
            f_array[3*index+1] -= F[1]
            
    return force

def addEdgeTraction(order, forest, attr, assembler, tr):
    trac = []
    tx = np.ones(order)*tr[0]
    ty = np.ones(order)*tr[1]
    for findex in range(4):
        trac.append(elements.PSThermoQuadTraction(findex, tx, ty))

    # Retrieve octants from the forest
    quadrants = forest.getQuadrants()
    edge_quads = forest.getQuadsWithName(attr)
    aux = TACS.AuxElements()
    
    for i in range(len(edge_quads)):
        aux.addElement(edge_quads[i].tag, trac[edge_quads[i].info])

    return aux

def createTopoProblem(props, forest, order=3, nlevels=2,
                      ordering=TACS.PY_MULTICOLOR_ORDER,
                      ite=1):
    # Create the forest
    forests = []
    filters = []
    assemblers = []
    varmaps = []
    vecindices = []

    # Create the trees, rebalance the elements and repartition
    forest.balance(1)
    forest.setMeshOrder(order)
    forest.repartition()
    forests.append(forest)
    use_bernstein=1

    # Make the creator class
    creator = CreateMe(bcs, forests[-1], props)
    assemblers.append(creator.createTACS(forest, ordering=ordering))
    filters.append(creator.getFilter())

    for i in range(nlevels-1):
        if args.use_decrease_order:
            forest = forests[-1].duplicate()
            forest.balance(1)        
            if order-1-i >= 3:
                forest.setMeshOrder(order-1-i)
            else:
                forest.setMeshOrder(3)
        else:
            forest = forests[-1].coarsen()
            forest.balance(1)
            forest.setMeshOrder(order)
        
        forests.append(forest)

        # Make the creator class
        creator = CreateMe(bcs, forests[-1], props)
        assemblers.append(creator.createTACS(forest, ordering=ordering))
        filters.append(creator.getFilter())

    # Create the multigrid object
    mg = TMR.createMg(assemblers, forests, omega=args.omega)
   
    # Create the topology optimization problem
    vars_per_node = 1
    nmats = props.getNumMaterials()
    if nmats > 1:
        vars_per_node = nmats+1

    # Set the filter type
    fltr = TMR.ConformFilter(assemblers, filters, vars_per_node=vars_per_node)
    problem = TMR.TopoProblem(fltr, mg)

    return assemblers[0], problem, filters[0], varmaps[0]

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

# The communicator
comm = MPI.COMM_WORLD

# Print out all of the arguments to the command line
if comm.rank == 0:
    for arg in vars(args):
        print('%-20s'%(arg), getattr(args, arg))

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
mesh.mesh(args.htarget, opts)

# Create a model from the mesh
model = mesh.createModelFromMesh()

# Create the corresponding mesh topology from the mesh-model 
topo = TMR.Topology(comm, model)

# Create the quad forest and set the topology of the forest
forest = TMR.QuadForest(comm)
forest.setTopology(topo)

# Compute the volume of the bracket
r = 0.06
a = 0.04
vol = r*a
vol_frac = args.vol_frac

# Set the max nmber of iterations
mg_levels = args.mg_levels
max_iterations = len(mg_levels)

# Set parameters for later usage
forest.createTrees(args.init_depth)

# The old filter/map classes
old_filtr = None
old_varmap = None

# Values from the previous iteration
old_dvs = None
old_z = 0.0
olz_zl = None
old_zu = None

# Set the values of the objective array
obj_array = [ 1e0 ]
thickness = 1.0e-2
rho = [2600*thickness, 1300.*thickness]
E = [70e9*thickness, 35e9*thickness]
nu = [0.3, 0.3]
aT = [23.5e-6, 23.5e-6*0.5]
kcond = [130.0*thickness, 65.0*thickness]
ys = [450e6*thickness, 275e6*thickness]

# Create the stiffness properties object
props = TMR.QuadStiffnessProperties(rho, E, nu, _aT=aT, _kcond=kcond,
                                    k0=1e-3, eps=0.2,
                                    q=args.q_penalty, qtemp=0.0,
                                    qcond=0.0)

# Set the fixed mass
max_density = sum(rho)/len(rho)
initial_mass = vol*max_density
m_fixed = vol_frac*initial_mass

# Set the number of variables per node
vars_per_node = 1
if (len(rho) > 1):
    vars_per_node = 1+len(rho)
    
time_array = np.zeros(sum(args.max_opt_iters[:]))
t0 = MPI.Wtime()

for step in range(max_iterations):
    # Create the TACSAssembler and TMRTopoProblem instance
    nlevs = mg_levels[step]
    assembler, problem, filtr, varmap = createTopoProblem(props, forest, 
                                                          nlevels=nlevs,
                                                          order=args.order)    
    # Set the constraint type
    funcs = [functions.StructuralMass(assembler)]
    T = 2.5e6
    aux = addEdgeTraction(forest.getMeshOrder(), forest, 'traction', assembler,
                          [0.0, -T])
    assembler.setAuxElements(aux)
    force1 = assembler.createVec()
    assembler.assembleRes(force1)
    force1.scale(-1.0)

    # Set the auxilary elements back to None
    assembler.setAuxElements(None)
    
    # Set the load cases
    forces = [force1]
    problem.setLoadCases(forces)

    # Set the mass constraint
    # (m_fixed - m(x))/m_fixed >= 0.0    
    problem.addConstraints(0, funcs, [-m_fixed], [-1.0/m_fixed])
    problem.setObjective(obj_array, [functions.Compliance(assembler)])
    
    # Initialize the problem and set the prefix
    problem.initialize()
    problem.setPrefix(args.prefix)

    # Create the quasi-Newton Hessian approximation
    qn_subspace_size = args.qn_subspace
    qn = ParOpt.LBFGS(problem, subspace=qn_subspace_size)
    problem.setIterationCounter(sum(args.max_opt_iters[:step]))
    # Trust region problem parameters
    tr_size = 0.01
    tr_max_size = 0.02
    tr_min_size = 1e-6
    eta = 0.25
    tr_penalty = args.tr_penalty
    
    # Create the trust region problem
    tr = ParOpt.pyTrustRegion(problem, qn, tr_size,
                              tr_min_size, tr_max_size, eta, tr_penalty)
    # Create the ParOpt problem
    opt = ParOpt.pyParOpt(tr, qn_subspace_size, ParOpt.NO_HESSIAN_APPROX)
    
    # Set the penalty parameter internally in the code. These must be
    # consistent between the trust region object and ParOpt.
    opt.setPenaltyGamma(tr_penalty)
    
    # Set parameters for ParOpt in the subproblem
    opt.setMaxMajorIterations(500)
    opt.setAbsOptimalityTol(args.opt_abs_tol)
    
    # Don't update the quasi-Newton method
    opt.setQuasiNewton(qn)
    opt.setUseQuasiNewtonUpdates(0)
    
    if args.use_L1:
        opt.setNormType(ParOpt.L1_NORM)
    elif args.use_Linf:
        opt.setNormType(ParOpt.INFTY_NORM)
    else:
        opt.setNormType(ParOpt.L2_NORM)

    opt.setStartingPointStrategy(ParOpt.AFFINE_STEP)
    opt.setBarrierStrategy(ParOpt.COMPLEMENTARITY_FRACTION)
    
    # Set the design variable bounds
    filename = os.path.join(args.prefix, 'paropt%d.out'%(step))
    opt.setOutputFile(filename)
    
    # If the old filter exists, we're on the second iteration
    if old_filtr:            
        # Create the interpolation
        interp = TACS.VecInterp(old_varmap, varmap, vars_per_node)
        filtr.createInterpolation(old_filtr, interp)
        interp.initialize()
        
        # Get the optimization variables for the new optimizer
        x, z, zw, zl, zu = opt.getOptimizedPoint()
        
        # Do the interpolation
        old_vec = problem.convertPVecToVec(old_x)
        x_vec = problem.convertPVecToVec(x)
        interp.mult(old_vec, x_vec)
        
        # Set the new design variables
        problem.setInitDesignVars(x)
        
    # Initialize the problem
    tr.initialize()
       
    # Iterate
    for i in range(args.max_opt_iters[step]):
        itr = sum(args.max_opt_iters[:step]) + i
        if comm.rank == 0:
            print("Iteration[%d]"%(itr))
        opt.setInitBarrierParameter(10.0)
        opt.resetDesignAndBounds()            
        opt.optimize()
        
        # Get the optimized point
        x, z, zw, zl, zu = opt.getOptimizedPoint()
        
        # Write the solution out to a file
        if i % args.output_freq == 0:
            # Compute the iteration counter
            itr = sum(args.max_opt_iters[:step]) + i
            
            # Write f5 object
            # Output for visualization
            flag = (TACS.ToFH5.NODES |
                    TACS.ToFH5.DISPLACEMENTS |
                    TACS.ToFH5.EXTRAS)
            f5 = TACS.ToFH5(assembler, TACS.PY_PLANE_STRESS, flag)
            f5.writeToFile(os.path.join(args.prefix, 'plate%04d.f5'%(itr))) 
    
            # Update the trust region method
            infeas, l1, linfty = tr.update(x, z, zw)
           
            t1_ind = MPI.Wtime()
            time_array[itr] = t1_ind-t0
            t0 = 1.0*t1_ind
            
        # Check for convergence using relatively strict tolerances
        if infeas < 1e-6 and l1 < 1e-4 and step > 1:
            break
            
    # Set the old values of the variables
    old_x, z, zw, zl, zu = opt.getOptimizedPoint()
 
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
