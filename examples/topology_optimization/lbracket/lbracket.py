
from mpi4py import MPI
from tmr import TMR, TopOptUtils
from paropt import ParOpt
from tacs import TACS, elements, constitutive, functions
import numpy as np
import argparse
import os
from scipy.optimize import minimize

class OptWeights:
    def __init__(self, diag, X, H):
        self.diag = diag
        self.X = X
        self.n = self.X.shape[0]

        # Compute the normalization
        if len(self.X.shape) == 1:
            self.delta = np.max(np.absolute(self.X - self.X[self.diag]))
        else:
            self.delta = np.sqrt(np.max(
                np.sum((self.X - self.X[self.diag,:])*(self.X - self.X[self.diag,:]), axis=1)))

        self.dim = 3
        if len(self.X.shape) == 1 or self.X.shape[1] == 1:
            self.dim = 1

            # Compute the constraint matrix
            A = np.zeros((2, self.n-1))

            # Populate the b vector
            b = np.zeros(2)
            b[1] = H[0,0]

            index = 0
            for i in range(self.n):
                if i != self.diag:
                    dx = (self.X[i] - self.X[self.diag])/self.delta

                    A[0,index] = dx
                    A[1,index] = 0.5*dx**2
                    index += 1

        elif self.X.shape[1] == 2:
            self.dim = 2

            # Compute the constraint matrix
            A = np.zeros((5, self.n-1))

            # Populate the b vector
            b = np.zeros(5)
            b[2] = H[0,0]
            b[3] = H[1,1]
            b[4] = 2.0*H[0,1]

            index = 0
            for i in range(self.n):
                if i != self.diag:
                    dx = (self.X[i,0] - self.X[self.diag,0])/self.delta
                    dy = (self.X[i,1] - self.X[self.diag,1])/self.delta

                    A[0,index] = dx
                    A[1,index] = dy
                    A[2,index] = 0.5*dx**2
                    A[3,index] = 0.5*dy**2
                    A[4,index] = dx*dy
                    index += 1
        else:
            # Compute the constraint matrix
            A = np.zeros((9, self.n-1))

            # Populate the b vector
            b = np.zeros(9)
            b[3] = H[0,0]
            b[4] = H[1,1]
            b[5] = H[2,2]
            b[6] = 2*H[1,2]
            b[7] = 2*H[0,2]
            b[8] = 2*H[0,1]

            index = 0
            for i in range(self.n):
                if i != self.diag:
                    dx = (self.X[i,0] - self.X[self.diag,0])/self.delta
                    dy = (self.X[i,1] - self.X[self.diag,1])/self.delta
                    dz = (self.X[i,2] - self.X[self.diag,2])/self.delta

                    A[0,index] = dx
                    A[1,index] = dy
                    A[2,index] = dz

                    A[3,index] = 0.5*dx**2
                    A[4,index] = 0.5*dy**2
                    A[5,index] = 0.5*dz**2

                    A[6,index] = dy*dz
                    A[7,index] = dx*dz
                    A[8,index] = dx*dy
                    index += 1

        self.b = b
        self.A = A

        return

    def obj_func(self, w):
        return 0.5*np.sum(w**2)

    def obj_func_der(self, w):
        return w

    def con_func(self, w):
        '''Compute the coefficients'''
        return np.dot(self.A, w) - self.b

    def con_func_der(self, w):
        '''Compute the derivative of the coefficients'''
        return self.A

    def set_alphas(self, w, alpha):
        alpha[:] = 0.0
        index = 0
        for i in range(self.n):
            if i != self.diag:
                alpha[i] = w[index]/self.delta**2
                alpha[self.diag] += w[index]/self.delta**2
                index += 1

        alpha[self.diag] += 1.0

        return

class Mfilter(TMR.HelmholtzPUFilter):
    def __init__(self, N, assemblers, filter, vars_per_node=1):
        return

    def getInteriorStencil(self, diag, X, alpha):
        '''
        Get an interior stencil point.
        '''
        self.r = 0.1
        H = self.r**2*np.eye(3)

        # Reshape the values in the matrix
        X = X.reshape((-1, 3))
        X = X[:,:2]
        n = X.shape[0]

        # Set up the optimization problem
        opt = OptWeights(diag, X, H)

        # Set the bounds and initial point
        w0 = np.ones(n-1)
        bounds = []
        for i in range(n-1):
            bounds.append((0, None))

        res = minimize(opt.obj_func, w0, jac=opt.obj_func_der, method='SLSQP', bounds=bounds,
                       constraints={'type': 'eq', 'fun': opt.con_func, 'jac': opt.con_func_der})

        # Set the optimized alpha values
        opt.set_alphas(res.x, alpha)

        return

    def getBoundaryStencil(self, diag, normal, X, alpha):
        '''
        Get a sentcil point on the domain boundary
        '''
        self.r = 0.1
        H = self.r**2*np.eye(2)

        # Reshape the values in the matrix
        X = X.reshape((-1, 3))
        X = X[:,:2]
        n = X.shape[0]

        # # Reduce the problem to a 2d problem on linearization of the
        # # the domain boundary. First, compute an arbitrary direction
        # # that is not aligned along the normal direction
        # index = np.argmin(np.absolute(normal))
        # t = np.zeros(3)
        # t[index] = 1.0

        # # Compute the in-plane directions (orthogonal to the normal direction)
        # t2 = np.cross(t, normal)
        # t1 = np.cross(normal, t2)

        # # Reduce the problem on the boundary
        # Xt = np.zeros(n, 2))
        # Xt[:,0] = np.dot(X - X[diag,:], t1)
        # Xt[:,1] = np.dot(X - X[diag,:], t2)

        t = np.array([normal[1], -normal[0]])
        Xt = np.dot(X - X[diag,:], t)

        # Set up the optimization problem
        opt = OptWeights(diag, Xt, H)

        # Set the bounds and initial point
        w0 = np.ones(n-1)
        bounds = []
        for i in range(n-1):
            bounds.append((0, None))

        res = minimize(opt.obj_func, w0, jac=opt.obj_func_der, method='SLSQP', bounds=bounds,
                       constraints={'type': 'eq', 'fun': opt.con_func, 'jac': opt.con_func_der})

        # Set the optimized alpha values
        opt.set_alphas(res.x, alpha)

        return

class QuadConformCreator(TMR.QuadConformTopoCreator):
    """
    This class is called to create a TACSAssembler class with an underlying
    filter mesh that conforms with Assembler mesh. The input to the class
    consists of the boundary conditions, the forest of quadtrees for the
    analysis mesh, the order of the conforming filter QuadForest mesh, and
    the type of interpolant to be used.
    """
    def __init__(self, bcs, forest, order=2, interp=TMR.BERNSTEIN_POINTS,
                 design_vars_per_node=1, props=None):
        # Store the properties
        self.props = props
        self.element = None

    def createTopoElement(self, order, filtr):
        """
        Create an element for the entire mesh
        """

        # Create the constitutive object - one for the entire mesh
        self.con = TMR.QuadConstitutive(props=self.props, forest=filtr)

        # Create the model
        self.model = elements.LinearElasticity2D(self.con)

        # Set the basis functions and create the element
        if order == 2:
            self.basis = elements.LinearQuadBasis()
        elif order == 3:
            self.basis = elements.QuadraticQuadBasis()
        elif order == 4:
            self.basis = elements.CubicQuadBasis()
        elif order == 5:
            self.basis = elements.QuarticQuadBasis()
        elif order == 6:
            self.basis = elements.QuinticQuadBasis()

        # Create the element type
        self.element = elements.Element2D(self.model, self.basis)

        return

    def createElement(self, order, quadrant, index, filtr):
        """
        Create the element for the specified quadrant

        Args:
            order (int): The order of the element
            quadrant (TMR.Quadrant): The quadrant to be build for this element
            index (list): The global numbers for the quadrant nodes
            filtr (QuadForest): The QuadForest for the filter

        Returns:
            TACS.Element: Element to place within the Assembler
        """
        if self.element is None:
            self.createTopoElement(order, filtr)
        return self.element

class CreatorCallback:
    def __init__(self, bcs, props):
        self.bcs = bcs
        self.props = props

    def creator_callback(self, forest):
        """
        Given the forest, instantiate a creator class that will populate a
        TACSAssembler object with the elements. This allocates the
        QuadConformCreator class above and also returns the QuadForest object
        associated with the filter. This is needed in the createTopoProblem
        call.
        """
        order = forest.getMeshOrder()
        interp = TMR.BERNSTEIN_POINTS
        dvs_per_node = self.props.getDesignVarsPerNode()
        creator = QuadConformCreator(self.bcs, forest, order=order,
                                     design_vars_per_node=dvs_per_node,
                                     interp=interp, props=self.props)
        filtr = creator.getFilter()
        return creator, filtr

def create_forest(comm, depth, htarget):
    """
    Create an initial forest for analysis and optimization

    This code loads in the model, sets names, meshes the geometry and creates
    a QuadForest from the mesh. The forest is populated with quadtrees with
    the specified depth.

    Args:
        comm (MPI_Comm): MPI communicator
        depth (int): Depth of the initial trees
        htarget (float): Target global element mesh size

    Returns:
        QuadForest: Initial forest for topology optimization
    """

    # Load in the L-bracket model
    geo = TMR.LoadModel('2d-bracket-fillet.stp')
    verts = geo.getVertices()
    edges = geo.getEdges()
    faces = geo.getFaces()

    geo = TMR.Model(verts, edges, faces)

    # Set the edges
    edges[5].setName('fixed')
    edges[1].setName('traction')

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

def filter_callback(assemblers, filters):
    N = 20
    r = 0.1
    return Mfilter(N, assemblers, filters)

def create_problem(forest, bcs, props, nlevels, iter_offset=0, m_fixed=0.0):
    """
    Create the TMRTopoProblem object and set up the topology optimization problem.

    This code is given the forest, boundary conditions, material properties and
    the number of multigrid levels. Based on this info, it creates the TMRTopoProblem
    and sets up the mass-constrained compliance minimization problem. Before
    the problem class is returned it is initialized so that it can be used for
    optimization.

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

    # Create a conforming filter
    filter_type = filter_callback

    # Characteristic length of the domain
    a = 0.1
    b = (2.0/5.0)*a
    area = a**2 - (a - b)**2
    r0 = 0.1*np.sqrt(area)

    # Create the problem and filter object
    problem = TopOptUtils.createTopoProblem(forest, obj.creator_callback, filter_type,
                                            nlevels=nlevels, lowest_order=2,
                                            r0=r0, N=50, use_galerkin=True,
                                            design_vars_per_node=1)

    # Get the assembler object we just created
    assembler = problem.getAssembler()

    # Get the basis object from one of the elements
    elems = assembler.getElements()
    basis = elems[0].getElementBasis()

    # Create the traction objects that will be used later..
    Ty = -2.5e6 # Traction force component in the y-direction
    vpn = elems[0].getVarsPerNode()
    trac = [0.0, Ty]
    tractions = []
    for findex in range(4):
        tractions.append(elements.Traction2D(vpn, findex, basis, trac))

    # Allocate a thermal traction boundary condition
    force1 = TopOptUtils.computeTractionLoad('traction', forest, assembler,
                                             tractions)

    # Set the load case
    problem.setLoadCases([force1])

    # Set the constraint functions
    funcs = [functions.StructuralMass(assembler)]

    # Set the mass constraint
    # (m_fixed - m(x))/m_fixed >= 0.0
    problem.addConstraints(0, funcs, [-m_fixed], [-1.0/m_fixed])

    # Set the values of the objective array
    obj_array = [ 1.0 ]
    ksfail = functions.KSFailure(assembler, 100.0)
    ksfail.setKSFailureType('discrete')
    problem.setObjective(obj_array, [ksfail])

    # Initialize the problem and set the prefix
    problem.initialize()

    # Set the callback for generating output
    cb = OutputCallback(assembler, iter_offset=iter_offset)
    problem.setOutputCallback(cb.write_output)

    return problem


class OutputCallback:
    def __init__(self, assembler, iter_offset=0):
        self.fig = None

        # Set the output file name
        flag = (TACS.OUTPUT_CONNECTIVITY |
                TACS.OUTPUT_NODES |
                TACS.OUTPUT_DISPLACEMENTS |
                TACS.OUTPUT_STRAINS |
                TACS.OUTPUT_EXTRAS)
        self.f5 = TACS.ToFH5(assembler, TACS.PLANE_STRESS_ELEMENT, flag)
        self.iter_offset = iter_offset

    def write_output(self, prefix, itr, oct_forest, quad_forest, x):
        self.f5.writeToFile('results/output%d.f5'%(itr + self.iter_offset))

# Set the optimization parameters
optimization_options = {
    # Parameters for the trust region method
    'tr_init_size': 0.01,
    'tr_max_size': 0.05,
    'tr_min_size': 1e-5,
    'tr_eta': 0.1,
    'tr_penalty_gamma': 5.0,
    'tr_write_output_freq': 1,
    'tr_infeas_tol': 1e-5,
    'tr_l1_tol': 1e-5,
    'tr_adaptive_gamma_update': False,
    'tr_linfty_tol': 0.0, # Don't use the l-infinity norm in the stopping criterion

    # Parameters for the interior point method (used to solve the
    # trust region subproblem)
    'max_qn_subspace': 5,
    'tol': 1e-8,
    'maxiter': 500,
    'norm_type': 'L1',
    'barrier_strategy': 'Complementarity fraction',
    'start_strategy': 'Affine step'}

# Create an argument parser to read in arguments from the command line
p = argparse.ArgumentParser()
p.add_argument('--prefix', type=str, default='./results')
p.add_argument('--vol_frac', type=float, default=0.3)
p.add_argument('--htarget', type=float, default=2.5e-3)
p.add_argument('--max_opt_iters', type=int, default=200)
p.add_argument('--init_depth', type=int, default=1)
p.add_argument('--mg_levels', type=int, default=3)
p.add_argument('--order', type=int, default=3)
p.add_argument('--q_penalty', type=float, default=8.0)
args = p.parse_args()

# Set the communicator
comm = MPI.COMM_WORLD

# Print out all of the arguments to the command line
if comm.rank == 0:
    for arg in vars(args):
        print('%-20s'%(arg), getattr(args, arg))

# Units are in m

# Create the first material properties object
rho = 2600.0
E = 70e9
nu = 0.3
ys = 100e6
mat = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)

# Set the fixed mass
a = 0.1
b = (2.0/5.0)*a
vol = a**2 - (a - b)**2
full_mass = vol*rho
m_fixed = args.vol_frac*full_mass

# Create the stiffness properties object
props = TMR.StiffnessProperties(mat, q=args.q_penalty,
                                qcond=args.q_penalty, eps=0.05, k0=1e-9)

# Set the boundary conditions for the problem
bcs = TMR.BoundaryConditions()
bcs.addBoundaryCondition('fixed', [0, 1], [0.0, 0.0])

# Create the initial forest
forest = create_forest(comm, args.init_depth, args.htarget)
forest.setMeshOrder(args.order, TMR.GAUSS_LOBATTO_POINTS)

# Create the TMRTopoProblem instance
problem = create_problem(forest, bcs, props, args.mg_levels, m_fixed=m_fixed)
problem.setPrefix(args.prefix)

# Check the gradient
problem.checkGradients(1e-6)

# Set parameters
optimization_options['maxiter'] = args.max_opt_iters
optimization_options['output_file'] = os.path.join(args.prefix,
                                                   'output_file.dat')
optimization_options['tr_output_file'] = os.path.join(args.prefix,
                                                      'tr_output_file.dat')

# Optimize the problem
opt = TopOptUtils.TopologyOptimizer(problem, optimization_options)
xopt = opt.optimize()
