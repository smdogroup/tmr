from mpi4py import MPI
from tmr import TMR
from paropt import ParOpt
from tacs import TACS, elements, constitutive, functions
import numpy as np
import argparse
import os

class CreateMe(TMR.OctBernsteinTopoCreator):
    def __init__(self, bcs, forest, use_bernstein, props):
        TMR.OctBernsteinTopoCreator.__init__(bcs, forest, use_bernstein)        
        # Create the array of properties
        self.props = props
        
    def createElement(self, order, octant, index, filtr):
        '''Create the element'''
        stiff = TMR.ThermoOctStiffness(self.props,
                                       index, None, filtr)
        elem = elements.SolidThermo(order, stiff)
        
        return elem
     
def addVertexLoad(comm, forest, attr, assembler, F):
    # Retrieve octant from the forest
    octants = forest.getOctants()
    node_octs = forest.getNodesWithName(attr)
    force = assembler.createVec()
    f_array = force.getArray()
    node_range = forest.getNodeRange()
    mpi_rank = comm.Get_rank()
    for i in range(len(node_octs)):
        if ((node_octs[i] >= node_range[mpi_rank]) and
            (node_octs[i] < node_range[mpi_rank+1])): 
            index = node_octs[i] - node_range[mpi_rank]
            f_array[4*index] += F[0]
            f_array[4*index+1] += F[1]
            f_array[4*index+2] += F[2]
    return force

def addFaceTraction(order, forest, attr, assembler, tr):
    trac = []
    for findex in range(6):
        trac.append(elements.TACS3DThermoTraction(order, findex, tr[0], tr[1],
                                                  tr[2]))

    # Retrieve octants from the forest
    octants = forest.getOctants()
    face_octs = forest.getOctsWithName(attr)
    aux = TACS.AuxElements()

    for i in range(len(face_octs)):
        aux.addElement(face_octs[i].tag, trac[face_octs[i].info])

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
    forest.setMeshOrder(order)#, interp=TMR.GAUSS_LOBATTO_POINTS)
    forest.repartition()
    forests.append(forest)
    use_bernstein = args.use_bernstein
    # Make the creator class
    creator = CreateMe(bcs, forests[-1], use_bernstein, props)
    assemblers.append(creator.createTACS(forest, ordering=ordering))
    varmaps.append(creator.getMap())
    vecindices.append(creator.getIndices())
    filters.append(creator.getFilter())
    
    for i in xrange(nlevels-1):
        if args.use_decrease_order:
            forest = forests[-1].duplicate()
            forest.balance(1)
            if order-i-1 >= 3:
                forest.setMeshOrder(order-i-1)
            else: 
                forest.setMeshOrder(3)
        else:
            forest = forests[-1].coarsen()
            forest.balance(1)
            forest.setMeshOrder(order)#,interp=TMR.GAUSS_LOBATTO_POINTS)
        #forest.repartition()
        forests.append(forest)

        # Make the creator class
        creator = CreateMe(bcs, forests[-1], use_bernstein, props)
        assemblers.append(creator.createTACS(forest, ordering=ordering))
        varmaps.append(creator.getMap())
        vecindices.append(creator.getIndices())
        filters.append(creator.getFilter())

    # Create the multigrid object
    mg = TMR.createMg(assemblers, forests, omega=args.omega)
    
    # Create the topology optimization problem
    vars_per_node = 1
    nmats = props.getNumMaterials()
    if nmats > 1:
        vars_per_node = nmats+1
    problem = TMR.TopoProblem(assemblers, filters, mg,
                              varmaps, vecindices,  helmholtz_radius=0.25/100.,
                              vars_per_node=vars_per_node)

    return assemblers[0], problem, filters[0], varmaps[0]

# Create an argument parser to read in arguments from the commnad line
p = argparse.ArgumentParser()
p.add_argument('--prefix', type=str, default='./results')
p.add_argument('--vol_frac', type=float, default=0.3)
p.add_argument('--htarget', type=float, default=0.025)
p.add_argument('--max_opt_iters', type=int, nargs='+',
               default=[2000])
p.add_argument('--opt_abs_tol', type=float, default=1e-6)
p.add_argument('--opt_barrier_frac', type=float, default=0.25)
p.add_argument('--opt_barrier_power', type=float, default=1.0)
p.add_argument('--output_freq', type=int, default=1)
p.add_argument('--init_depth', type=int, default=1)
p.add_argument('--mg_levels', type=int, nargs='+', default=[2])
p.add_argument('--max_lbfgs', type=int, default=10)
p.add_argument('--hessian_reset', type=int, default=10)
p.add_argument('--tr_penalty', type=float, default=25.0)
p.add_argument('--qn_subspace', type=int, default=25)
p.add_argument('--use_paropt', action='store_true', default=True)
p.add_argument('--use_mma', action='store_true', default=False)
p.add_argument('--use_L1', action='store_true', default=False)
p.add_argument('--use_Linf', action='store_true', default=False)
p.add_argument('--use_decrease_order', action='store_true', default=False)
p.add_argument('--use_bernstein', action='store_true', default=False)
p.add_argument('--order', type=int, default=3)
p.add_argument('--q_penalty', type=float, default=8.0)
p.add_argument('--omega', type=float, default=1.0)
args = p.parse_args()

# Set the parameter to use paropt or MMA
use_paropt = False
if args.use_mma:
    use_paropt = False

# The communicator
comm = MPI.COMM_WORLD

# Print out all of the arguments to the command line
if comm.rank == 0:
    for arg in vars(args):
        print('%-20s'%(arg), getattr(args, arg))

# Load the geometry model
step_file = 'beam-hole.step'
bdf_file = 'beam-hole-coarse.bdf'
geo = TMR.LoadModel(step_file)

# Mark the boundary condition faces
verts = geo.getVertices()
edges = geo.getEdges()
faces = geo.getFaces()
volumes = geo.getVolumes()

faces[0].setName('fixed')
faces[6].setName('T1')
faces[7].setName('T1')
faces[3].setSource(volumes[0], faces[1])

# Load the face mesh from BDF file
external_mesh = TACS.MeshLoader(comm)
external_mesh.scanBDFFile(bdf_file)

# Broadcast the data
if comm.rank == 0:
   quad_ptr, quads, comps, Xpts = external_mesh.getConnectivity()
   quads = comm.bcast(quads, root=0)
   Xpts = comm.bcast(Xpts, root=0)
else:
   quads = comm.bcast(None, root=0)
   Xpts = comm.bcast(None, root=0)

conn = np.reshape(quads,(-1,4))
conn = conn.astype(np.intc)
Xpts = np.reshape(Xpts,(-1,3))
for i in range(conn.shape[0]):
   t = conn[i,2]
   conn[i,2] = conn[i,3]
   conn[i,3] = t
   
for i in range(conn.shape[0]):
    # Perform a check for the surface normal
    a = Xpts[conn[i,2]] - Xpts[conn[i,0]]
    b = Xpts[conn[i,3]] - Xpts[conn[i,1]]
    c = np.cross(a, b)

    if c[1] < 0.0:
        conn[i,:] = conn[i,::-1]

# For edge 0 and 2
A = Xpts[Xpts[:,0] == 0]
ext_e_mesh = TMR.EdgeMesh(comm, edges[0], A)
edges[0].setMesh(ext_e_mesh)
A[:,1] = 0.25
ext_e_mesh = TMR.EdgeMesh(comm, edges[2], A)
edges[2].setMesh(ext_e_mesh)

# For edge 5 and 10
A = Xpts[Xpts[:,2] == 0.25]
ext_e_mesh = TMR.EdgeMesh(comm, edges[5], A)
edges[5].setMesh(ext_e_mesh)
A[:,1] = 0.25
ext_e_mesh = TMR.EdgeMesh(comm, edges[10], A)
edges[10].setMesh(ext_e_mesh)

# # For edge 6 and 12
A = Xpts[Xpts[:,0] == 1.0]
ext_e_mesh = TMR.EdgeMesh(comm, edges[6], A)
edges[6].setMesh(ext_e_mesh)
A[:,1] = 0.25
ext_e_mesh = TMR.EdgeMesh(comm, edges[12], A)
edges[12].setMesh(ext_e_mesh)

# For edge 4 and 11
A = Xpts[Xpts[:,2] == 0.0]
ext_e_mesh = TMR.EdgeMesh(comm, edges[4], A)
edges[4].setMesh(ext_e_mesh)
A[:,1] = 0.25
ext_e_mesh = TMR.EdgeMesh(comm, edges[11], A)
edges[11].setMesh(ext_e_mesh)

hcircle = 0.3141592653589793/12
for e in [7, 8, 13, 14]:
    edge_mesh = TMR.EdgeMesh(comm, edges[e])
    edge_mesh.mesh(hcircle)
    edges[e].setMesh(edge_mesh)

# Create the surface mesh
ext_face_mesh = TMR.FaceMesh(comm, faces[1], Xpts, conn)
faces[1].setMesh(ext_face_mesh)

# Set the boundary conditions for the problem
bcs = TMR.BoundaryConditions()
bcs.addBoundaryCondition('fixed', [0,1,2,3],
                         [0.0,0.0,0.0,0.0])

# Create the mesh
mesh = TMR.Mesh(comm, geo)
opts = TMR.MeshOptions()
opts.frontal_quality_factor = 1.25
opts.num_smoothing_steps = 0
opts.triangularize_print_iter = 50000
opts.write_mesh_quality_histogram = 1
opts.reset_mesh_objects = 0

mesh.mesh(args.htarget, opts=opts)

# Create a model from the mesh
model = mesh.createModelFromMesh()

# Create the corresponding mesh topology from the mesh-model 
topo = TMR.Topology(comm, model)

# Create the oct forest and set the topology of the forest
forest = TMR.OctForest(comm)
forest.setTopology(topo)

# Compute the volume of the beam
vol = (1*.25-np.pi*.05**2)*.25
vol_frac = args.vol_frac

# Set the max nmber of iterations
mg_levels = args.mg_levels
max_iterations = len(mg_levels)
if len(args.max_opt_iters) < max_iterations:
    n = len(args.max_opt_iters) - max_iterations
    args.max_opt_iters.extend(n*[args.max_opt_iters[-1]])
    
# Set parameters for later usage
order = args.order
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
obj_array = [ 1e-2 ]
rho = [2600, 1300.]
E = [70e9, 35e9]
nu = [0.3, 0.3]
aT = [23.5e-6, 23.5e-6*0.5]
kcond = [130.0, 65.0]

# Create the stiffness properties object
props = TMR.StiffnessProperties(rho, E, nu, _aT=aT, _kcond=kcond, k0=1e-6,
                                q=args.q_penalty, qtemp=0.0,
                                qcond=args.q_penalty)

# Set the fixed mass
max_density = sum(rho)/len(rho)
initial_mass = vol*max_density
m_fixed = vol_frac*initial_mass

# Set the number of variables per node
vars_per_node = 1
if (len(rho) > 1):
    vars_per_node = 1+len(rho)

for step in xrange(max_iterations):
    # Create the TACSAssembler and TMRTopoProblem instance
    nlevs = mg_levels[step]
    assembler, problem, filtr, varmap = createTopoProblem(props, forest, 
                                                          nlevels=nlevs,
                                                          order=args.order)

    # Set the constraint type
    funcs = [functions.StructuralMass(assembler)]
    trac_mag = 58.0e5
    aux = addFaceTraction(forest.getMeshOrder(), forest, 'T1', assembler,
                          [0.0, 0.0, -trac_mag])
    force1 = assembler.createVec()
    assembler.setAuxElements(aux)
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
    use_tr = 0
    use_paropt = 1
    if use_tr:
        # Create the quasi-Newton Hessian approximation
        qn_subspace_size = args.qn_subspace
        qn = ParOpt.LBFGS(problem, subspace=qn_subspace_size)

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
            if comm.rank == 0:
                print("Iteration[%d]"%(sum(args.max_opt_iters[:step]) + i))
            opt.setInitBarrierParameter(10.0)
            opt.resetDesignAndBounds()
            opt.optimize()
            exit(0)
            # Get the optimized point
            x, z, zw, zl, zu = opt.getOptimizedPoint()

            # Write the solution out to a file
            if i % args.output_freq == 0:
                # Compute the iteration counter
                itr = sum(args.max_opt_iters[:step]) + i

                # Get the vector and convert it
                vec = problem.convertPVecToVec(x)

                for k in range(1, vars_per_node):
                    f = 'levelset05_var%d_binary%04d.bstl'%(k, itr)
                    filename = os.path.join(args.prefix, f)
                    TMR.writeSTLToBin(filename, filtr, vec, offset=k)
                
            # Update the trust region method
            infeas, l1, linfty = tr.update(x, z, zw)

            # Check for convergence using relatively strict tolerances
            if infeas < 1e-4 and l1 < 0.01:
                break

        # # Set the old values of the variables
        old_x, z, zw, zl, zu = opt.getOptimizedPoint()
    elif use_paropt:
        # Create the ParOpt problem
        opt = ParOpt.pyParOpt(problem, args.max_lbfgs, ParOpt.BFGS)

        # Set parameters
        opt.setMaxMajorIterations(args.max_opt_iters[step])
        opt.setHessianResetFreq(args.hessian_reset)
        problem.setIterationCounter(args.max_opt_iters[step]*step)
        opt.setAbsOptimalityTol(args.opt_abs_tol)
        opt.setBarrierFraction(args.opt_barrier_frac)
        opt.setBarrierPower(args.opt_barrier_power)
        opt.setOutputFrequency(args.output_freq)
        opt.setOutputFile(os.path.join(args.prefix, 
                                       'paropt_output%d.out'%(step)))

        opt.checkGradients(1e-6)
        # Output for visualization
        flag = (TACS.ToFH5.NODES |
                TACS.ToFH5.DISPLACEMENTS |
                TACS.ToFH5.EXTRAS)
        f5 = TACS.ToFH5(assembler, TACS.PY_SOLID, flag)
        f5.writeToFile(os.path.join(args.prefix, 'beam%d.f5'%(step)))
        
        exit(0)
    # Set the old filter/variable map for the next time through the
    # loop so that we can interpolate design variable values
    old_varmap = varmap
    old_filtr = filtr

    # Output for visualization
    flag = (TACS.ToFH5.NODES |
            TACS.ToFH5.DISPLACEMENTS |
            TACS.ToFH5.EXTRAS)
    f5 = TACS.ToFH5(assembler, TACS.PY_SOLID, flag)
    f5.writeToFile(os.path.join(args.prefix, 'beam%d.f5'%(step)))

    # Create refinement array
    num_elems = assembler.getNumElements()
    refine = np.zeros(num_elems, dtype=np.int32)

    # Refine based solely on the value of the density variable
    elems = assembler.getElements()
    
    for i in xrange(num_elems):        
        c = elems[i].getConstitutive()
        if c is not None:
            if vars_per_node == 1:
                density = c.getDVOutputValue(0, np.zeros(3, dtype=float))
            else:
                density = 1.0-c.getDVOutputValue(2, np.zeros(3, dtype=float))
            # Refine things differently depending on whether the
            # density is above or below a threshold
            if density >= 0.1:
                refine[i] = int(1)
            elif density < 0.05:
                refine[i] = int(-1)
    
    # Refine the forest
    forest.refine(refine)

