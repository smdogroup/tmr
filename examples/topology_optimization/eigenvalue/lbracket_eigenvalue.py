"""
Solve the eigenvalue-constrained optimization problem

"""
from mpi4py import MPI
from tmr import TMR, TopOptUtils
from paropt import ParOpt
from tacs import TACS, elements, constitutive, functions
import numpy as np
import argparse
import os

# Import the cantilever example set up
import sys
sys.path.append('../lbracket/')
from lbracket import MFilterCreator, create_forest, CreatorCallback, OutputCallback

class FrequencyConstraint:
    def __init__(self, omega=1.0, eig_scale=1.0, num_eigenvalues=10,
                 max_jd_size=50, max_gmres_size=15, rho=100.0,
                 conscale=10.0, contype='semi-def'):

        self.omega = omega
        self.num_eigenvalues = num_eigenvalues
        self.max_jd_size = max_jd_size
        self.max_gmres_size = max_gmres_size
        self.rho = rho
        self.eig_scale = eig_scale
        self.conscale = conscale
        self.contype = contype

        # Set None objects for things that we'll allocate later
        self.fltr = None
        self.mg = None
        self.assembler = None
        self.oper = None
        self.jd = None
        self.W = None

        # Matrix for the semi-definite approach
        self.mat = None

        # Matrix for the normal method
        self.kmat = None
        self.mmat = None
        self.temp = None

        return

    def constraint(self, fltr, mg):
        """
        Evaluate the frequency constraint.
        """

        if self.fltr is None:
            self.mg = mg
            self.fltr = fltr
            self.assembler = self.fltr.getAssembler()
            self.comm = self.assembler.getMPIComm()

            if self.contype == 'semi-def':
                # Set the matrix for the operator
                self.mat = self.assembler.createMat()

                # Set up the operator with the given matrix and the
                # multigrid preconditioner
                self.oper = TACS.JDSimpleOperator(self.assembler, self.mat, self.mg)
            else:
                # Create the mass and stiffness matrices
                self.mmat = self.assembler.createMat()
                self.kmat = self.assembler.createMat()
                self.temp = self.assembler.createVec()

                # Creat ethe operator
                self.oper = TACS.JDFrequencyOperator(self.assembler,
                                                     self.kmat, self.mmat,
                                                     self.mg.getMat(), self.mg)

            # Create the eigenvalue solver and set the number of
            # recycling eigenvectors
            self.jd = TACS.JacobiDavidson(self.oper, self.num_eigenvalues,
                                          self.max_jd_size,
                                          self.max_gmres_size)
            self.jd.setTolerances(1e-6, 1e-8)
            self.jd.setRecycle(self.num_eigenvalues)

            # Create temporary vectors needed for the constraint
            # computation
            self.W = []
            for i in range(self.num_eigenvalues):
                self.W.append(self.assembler.createVec())

        if self.contype == 'semi-def':
            # Assemble the preconditioner at the point K - omega**2*M
            self.mg.assembleMatCombo(TACS.STIFFNESS_MATRIX, 1.0,
                                     TACS.MASS_MATRIX, -0.95*self.omega**2)
            self.mg.factor()

            # Assemble the exact matrix at K + lambda*G
            self.assembler.assembleMatCombo(TACS.STIFFNESS_MATRIX, 1.0,
                                            TACS.MASS_MATRIX, -self.omega**2,
                                            self.mat)
        else:
            # Assemble the matrices for the generalized eigenvalue problem
            self.assembler.assembleMatType(TACS.STIFFNESS_MATRIX, self.kmat)
            self.assembler.assembleMatType(TACS.MASS_MATRIX, self.mmat)
            self.mg.assembleMatType(TACS.STIFFNESS_MATRIX)
            self.mg.factor()

        # Solve the problem
        self.jd.solve(print_flag=True, print_level=0)

        # Check if the solve was successful, otherwise try again
        if self.jd.getNumConvergedEigenvalues() < self.num_eigenvalues:
            self.jd.solve(print_flag=True, print_level=1)

        # Extract the eigenvalues and eigenvectors
        self.eigs = np.zeros(self.num_eigenvalues)
        for i in range(self.num_eigenvalues):
            eig, error = self.jd.extractEigenvector(i, self.W[i])
            self.eigs[i] = eig

        # Scale the eigenvalues by the scaling factor
        self.eigs[:] /= self.eig_scale

        # Compute the value of the eigenvalue constraint
        eig_min = np.min(self.eigs)

        self.beta = 0.0
        self.eta = np.zeros(self.num_eigenvalues)
        for i, eig in enumerate(self.eigs):
            self.eta[i] = np.exp(-self.rho*(eig - eig_min))
            self.beta += self.eta[i]

        # Complete the computation of the eta values
        self.eta = self.eta/self.beta

        if self.contype == 'semi-def':
            # Compute the KS function of the minimum eigenvalues
            cval = self.eig_scale*(eig_min - np.log(self.beta)/self.rho)
        elif self.contype == 'ie-semi-def':
            # Compute the induced exponential function:
            # sum_{i} eig[i]*e^(-rho*eig[i])/ sum_{i} e^{-rho*eig[i]}
            cval = self.eig_scale*np.dot(self.eta, self.eigs)
        else:
            cval = self.eig_scale*((eig_min - np.log(self.beta)/self.rho) - self.omega**2)

        # Scale the constraint by the constraint scaling
        cval *= self.conscale

        if self.comm.rank == 0:
            print('KS frequency constraint: %25.10e'%(cval))

        return [cval]

    def constraint_gradient(self, fltr, mg, vecs):
        """
        Compute the constraint gradient
        """
        dcdx = vecs[0]
        dcdx.zeroEntries()

        for i in range(self.num_eigenvalues):
            # Compute the scalar factor to add to
            scale = self.eta[i]*self.conscale

            # Add the derivative of (K - omega**2*M) w.r.t. the
            # design variables to the constraint gradient
            if self.contype == 'semi-def':
                self.assembler.addMatDVSensInnerProduct(
                    scale, TACS.STIFFNESS_MATRIX,
                    self.W[i], self.W[i], dcdx)
                self.assembler.addMatDVSensInnerProduct(
                    -scale*self.omega**2, TACS.MASS_MATRIX,
                     self.W[i], self.W[i], dcdx)
            elif self.contype == 'ie-semi-def':
                cval = np.dot(self.eta, self.eigs)
                scale *= 1.0 + self.rho*self.eig_scale*(self.eigs[i] - cval)

                self.assembler.addMatDVSensInnerProduct(
                    scale, TACS.STIFFNESS_MATRIX,
                    self.W[i], self.W[i], dcdx)
                self.assembler.addMatDVSensInnerProduct(
                    -scale*self.omega**2, TACS.MASS_MATRIX,
                    self.W[i], self.W[i], dcdx)
            else:
                # Adjust the scaling factor
                self.mmat.mult(self.W[i], self.temp)
                scale /= self.temp.dot(self.W[i])

                self.assembler.addMatDVSensInnerProduct(
                    scale, TACS.STIFFNESS_MATRIX,
                    self.W[i], self.W[i], dcdx)
                self.assembler.addMatDVSensInnerProduct(
                    -scale*self.eigs[i], TACS.MASS_MATRIX,
                    self.W[i], self.W[i], dcdx)

        return

class BucklingConstraint:
    def __init__(self, load_factor=1.0, eig_scale=1.0, num_eigenvalues=10,
                 max_jd_size=50, max_gmres_size=15, rho=100.0):

        self.load_factor = load_factor
        self.num_eigenvalues = num_eigenvalues
        self.max_jd_size = max_jd_size
        self.max_gmres_size = max_gmres_size
        self.rho = rho
        self.eig_scale = eig_scale

        # Set None objects for things that we'll allocate later
        self.fltr = None
        self.mg = None
        self.assembler = None
        self.mat = None
        self.oper = None
        self.jd = None
        self.W = None
        self.dfdu = None
        self.adjoint = None

        return

    def constraint(self, fltr, mg):
        """
        Evaluate the buckling constraint.
        """

        if self.fltr is None:
            self.mg = mg
            self.fltr = fltr
            self.assembler = self.fltr.getAssembler()
            self.comm = self.assembler.getMPIComm()

            # Get the matrix associated with the pre conditioner
            self.mat = self.assembler.createMat()

            # Set up the operator with the given matrix and the multigrid preconditioner
            self.oper = TACS.JDSimpleOperator(self.assembler, self.mat, self.mg)

            # Create the eigenvalue solver and set the number of recycling eigenvectors
            self.jd = TACS.JacobiDavidson(self.oper, self.num_eigenvalues,
                                          self.max_jd_size,
                                          self.max_gmres_size)

            # Here, we're satisfied with
            eig_rtol = 1e-5
            eig_atol = 1e-5
            self.jd.setTolerances(eig_rtol, eig_atol)

            # Set the number of eigenvectors to recycle. Don't aggressively recycle,
            # b/c this can lead to a lack of robustness
            num_recycle = self.num_eigenvalues >> 1
            self.jd.setRecycle(self.num_eigenvalues)

            self.gmres = TACS.KSM(self.mg.getMat(), self.mg, 25, nrestart=20, isFlexible=1)
            self.gmres.setMonitor(self.comm)
            self.gmres.setTolerances(1e-10, 1e-30)

            # Create temporary vectors needed for the constraint
            # computation
            self.W = []
            for i in range(self.num_eigenvalues):
                self.W.append(self.assembler.createVec())

            self.vars = self.assembler.createVec()
            self.dfdu = self.assembler.createVec()
            self.adjoint = self.assembler.createVec()

        # We would normally solve for the solution path at this point
        # and set those variables into the assembler object. Here, we
        # know that the path variables have already been set during
        # the objective call, so we skip that step here. This is a bit
        # risky, but works because the constraint and constraint
        # gradient call always occur after the objective has been
        # called. This allows us to safely skip this computation here.
        self.assembler.getVariables(self.vars)

        # Assemble the preconditioner at the point K + lambda*G
        self.mg.assembleMatType(TACS.STIFFNESS_MATRIX)
        self.mg.factor()

        # Assemble the exact matrix at K + lambda*G
        self.assembler.assembleMatCombo(TACS.STIFFNESS_MATRIX, 1.0,
                                        TACS.GEOMETRIC_STIFFNESS_MATRIX, self.load_factor,
                                        self.mat)

        # Solve the problem
        self.jd.solve(print_flag=True, print_level=1)

        # Check if the solve was successful, otherwise try again
        nconv = self.jd.getNumConvergedEigenvalues()
        if nconv < self.num_eigenvalues:
            # Re-assemble and factor the preconditioner. Discard the part from
            # the geometric stiffness matrix
            self.mg.assembleMatType(TACS.STIFFNESS_MATRIX)
            self.mg.factor()

            # Temporarily set the number of recycling vectors to zero
            self.jd.solve(print_flag=True, print_level=1)

        # Extract the eigenvalues and eigenvectors
        self.eigs = np.zeros(self.num_eigenvalues)
        for i in range(self.num_eigenvalues):
            eig, error = self.jd.extractEigenvector(i, self.W[i])
            self.eigs[i] = eig

        # Scale the eigenvalues
        self.eigs /= self.eig_scale

        # Compute the value of the eigenvalue constraint
        eig_min = np.min(self.eigs)

        self.beta = 0.0
        self.eta = np.zeros(self.num_eigenvalues)
        for i, eig in enumerate(self.eigs):
            self.eta[i] = np.exp(-self.rho*(eig - eig_min))
            self.beta += self.eta[i]

        # Complete the computation of the eta values
        self.eta = self.eta/self.beta

        # Compute the KS function of the minimum eigenvalues
        cval = self.eig_scale*(eig_min - np.log(self.beta)/self.rho)

        if self.comm.rank == 0:
            print('KS stability constraint: %25.10e'%(cval))

        return [cval]

    def constraint_gradient(self, fltr, mg, vecs):
        """
        Compute the constraint gradient
        """
        dcdx = vecs[0]
        dcdx.zeroEntries()

        # Set variable values
        self.assembler.setVariables(self.vars)

        # Assemble and factor the stiffness matrix
        self.mg.assembleMatType(TACS.STIFFNESS_MATRIX)
        self.mg.factor()

        for i in range(self.num_eigenvalues):
            # Compute the scalar factor to add to
            scale = self.eta[i]

            # Add the derivative of (K + load_factor*G) w.r.t. the design
            # variables to the constraint gradient
            self.assembler.addMatDVSensInnerProduct(
                scale, TACS.STIFFNESS_MATRIX,
                self.W[i], self.W[i], dcdx)
            self.assembler.addMatDVSensInnerProduct(
                scale*self.load_factor, TACS.GEOMETRIC_STIFFNESS_MATRIX,
                self.W[i], self.W[i], dcdx)

            # Compute the derivative of d(W[i]^{T}*G*W[i])/d(path)
            self.assembler.evalMatSVSensInnerProduct(
                TACS.GEOMETRIC_STIFFNESS_MATRIX,
                self.W[i], self.W[i], self.dfdu)

            # Add the derivative of the path variables w.r.t. the
            # design variables
            self.gmres.solve(self.dfdu, self.adjoint)

            alpha = -scale*self.load_factor
            self.assembler.addAdjointResProducts([self.adjoint], [dcdx], alpha=alpha)

        return

def create_problem(forest, bcs, props, nlevels,
                   omega=10.0, load_factor=1.0, probtype='frequency',
                   num_eigenvalues=5, contype='semi-def',
                   r0_frac=0.05, N=20, iter_offset=0, m_fixed=0.0):
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
        omega (float): Natural frequency constraint value
        num_eigenvalues (int): Number of eigenvalues to include in the constraint
        contype (str): Type of constraint formulation
        r0_frac (float): Fraction of the characteristic domain length
        N (int): Number of iterations of the discrete filter
        iter_offset (int): Iteration offset counter
        m_fixed (float): Fixed mass fraction

    Returns:
        TopoProblem: Topology optimization problem instance
    """
    # Allocate the creator callback function
    obj = CreatorCallback(bcs, props)

    # Create a discrete M-filter
    mfilter = MFilterCreator(r0_frac, N)
    filter_type = mfilter.filter_callback

    # Create the problem and filter object
    problem = TopOptUtils.createTopoProblem(forest, obj.creator_callback, filter_type,
                                            nlevels=nlevels, lowest_order=2,
                                            use_galerkin=True,
                                            design_vars_per_node=1)

    # Get the assembler object we just created
    assembler = problem.getAssembler()

    # Get the basis object from one of the elements
    elems = assembler.getElements()
    basis = elems[0].getElementBasis()

    # Create the traction objects that will be used later..
    P = -5.0
    area = 0.005
    Ty = P/area # Traction force component in the y-direction
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

    # Set the objective
    problem.setObjective([0.1])

    # Set the constraint functions
    funcs = [functions.StructuralMass(assembler)]

    # Set the mass constraint
    # (m_fixed - m(x))/m_fixed >= 0.0
    mass_scale = 1.0
    problem.addConstraints(0, funcs, [-m_fixed], [-mass_scale/m_fixed])

    if probtype == 'frequency':
        # Add the frequency constraint
        eig_scale = 1.0
        if contype == 'semi-def':
            eig_scale = 1e-4

        freq = FrequencyConstraint(omega=omega, num_eigenvalues=num_eigenvalues,
                                   rho=10.0, eig_scale=eig_scale,
                                   conscale=1.0, contype=contype)
        problem.addConstraintCallback(1, freq.constraint, freq.constraint_gradient)
    else:
        # Add the buckling constraint callback
        eig_scale = 1e-4
        buckling = BucklingConstraint(load_factor=load_factor,
                                      num_eigenvalues=num_eigenvalues,
                                      rho=10.0, eig_scale=eig_scale)
        problem.addConstraintCallback(1, buckling.constraint,
                                      buckling.constraint_gradient)

    # Initialize the problem and set the prefix
    problem.initialize()

    return problem

def write_dvs_to_file(xopt, assembler, filename):
    # Output the original design variables before filtering
    rho_vec = assembler.createDesignVec()
    assembler.getDesignVars(rho_vec)
    x_vec = TMR.convertPVecToVec(xopt)
    assembler.setDesignVars(x_vec)

    # visualize
    flag = (TACS.OUTPUT_CONNECTIVITY |
            TACS.OUTPUT_NODES |
            TACS.OUTPUT_EXTRAS)
    f5 = TACS.ToFH5(assembler, TACS.PLANE_STRESS_ELEMENT, flag)
    f5.writeToFile(filename)

    # Set the tacs design vars back to the interpolated densities
    assembler.setDesignVars(rho_vec)

    return

# Set the optimization parameters
optimization_options = {
    # Parameters for the trust region method
    'tr_init_size': 0.01,
    'tr_max_size': 0.05,
    'tr_min_size': 1e-5,
    'tr_eta': 0.1,
    'tr_penalty_gamma': 100.0,
    'tr_penalty_gamma_list': [100.0, 1e6],
    'tr_adaptive_gamma_update': True, # Set whether to use an adaptive penalty or not
    'tr_penalty_gamma_max': 1e6, # Set the maximum penalty parameter
    'tr_write_output_freq': 1,
    'tr_infeas_tol': 1e-4,
    'tr_l1_tol': 1e-5,
    'tr_print_level': 2,
    'tr_linfty_tol': 0.0, # Don't use the l-infinity norm in the stopping criterion

    # Parameters for the interior point method (used to solve the
    # trust region subproblem)
    'qn_type': 'BFGS', # 'No Hessian approx'
    'max_qn_subspace': 5,
    'output_level': 2,
    'output_freq': 10,
    'tol': 1e-8,
    'maxiter': 500,
    'norm_type': 'L1',
    'barrier_strategy': 'Monotone',
    'start_strategy': 'Affine step'}

if __name__ == '__main__':
    # Create an argument parser to read in arguments from the command line
    p = argparse.ArgumentParser()
    p.add_argument('--prefix', type=str, default='./results')
    p.add_argument('--vol_frac', type=float, default=0.4)
    p.add_argument('--htarget', type=float, default=2.5e-3)
    p.add_argument('--max_opt_iters', type=int, default=200)
    p.add_argument('--init_depth', type=int, default=1)
    p.add_argument('--mg_levels', type=int, default=3)
    p.add_argument('--order', type=int, default=2)
    p.add_argument('--q_penalty', type=float, default=8.0)
    p.add_argument('--q_mass', type=float, default=3.0)
    p.add_argument('--N', type=int, default=10)
    p.add_argument('--r0_frac', type=float, default=0.05)
    p.add_argument('--use_project', action='store_true', default=False)
    p.add_argument('--use_simp', action='store_true', default=False)
    p.add_argument('--fs_type', type=str, default='none',
                   help='feature size refinement type: point, box, or None')
    p.add_argument('--probtype', type=str, default='frequency',
                   help='type of problem to use: frequency or buckling')
    p.add_argument('--contype', type=str, default='semi-def',
                   help='constraint type: semi-def or normal')
    p.add_argument('--num_eigenvalues', type=int, default=5,
                   help='number of eigenvalues use in the constraint')
    p.add_argument('--omega', type=float, default=10.0,
                   help='natural frequency constraint value')
    p.add_argument('--load_factor', type=float, default=1.0,
                   help='load factor value')
    args = p.parse_args()

    # Set the communicator
    comm = MPI.COMM_WORLD

    # Print out all of the arguments to the command line
    if comm.rank == 0:
        for arg in vars(args):
            print('%-20s'%(arg), getattr(args, arg))

    # Ensure that the prefix directory exists
    if comm.rank == 0 and not os.path.isdir(args.prefix):
        os.mkdir(args.prefix)

    # Set a barrier
    comm.Barrier()

    # Create the first material properties object
    rho = 2600.0
    E = 70e3
    nu = 0.3
    ys = 100.0
    mat = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)

    # Set the fixed mass
    a = 0.1
    b = (2.0/5.0)*a
    area = a**2 - (a - b)**2
    full_mass = area*rho
    m_fixed = args.vol_frac*full_mass

    # Set the characteristic domain length as the a-value
    domain_length = a

    # Create the stiffness properties object
    penalty_type = 'RAMP'
    if args.use_simp:
        penalty_type = 'SIMP'

    # Set the default to the lowest penalization q = 0.0 for RAMP, in
    # case the penalty type is SIMP, it gets reset to 1.0
    props = TMR.StiffnessProperties(mat, q=0.0, qmass=0.0, use_project=False,
                                    eps=0.05, k0=1e-4, penalty_type=penalty_type,
                                    beta=10.0)

    # Set the boundary conditions for the problem
    bcs = TMR.BoundaryConditions()
    bcs.addBoundaryCondition('fixed', [0, 1], [0.0, 0.0])

    # Create the initial forest
    forest = create_forest(comm, args.init_depth, args.htarget,
                           fs_type=args.fs_type,
                           filename='../lbracket/2d-bracket-fillet.stp')
    forest.writeToVTK(os.path.join(args.prefix, 'forest.vtk'))
    forest.setMeshOrder(args.order, TMR.GAUSS_LOBATTO_POINTS)

    # Set the original filter to NULL
    orig_filter = None
    xopt = None
    iter_offset = 0

    # Check what penalization strategy to use. Use a sequence of penalization
    # values which increase with iteration
    q_vals = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    if props.penalty_type == 'SIMP':
        q_vals = [1.0, 1.5, 2.0, 2.5, 3.0]

    max_iterations = len(q_vals)

    for step in range(max_iterations):
        # Set the penalization type
        props.stiffness_penalty_value = q_vals[step]

        # Create the TMRTopoProblem instance
        mg_levels = args.mg_levels
        if step > 0:
            mg_levels += 1
        problem = create_problem(forest, bcs, props, mg_levels,
                                 omega=args.omega, load_factor=args.load_factor,
                                 contype=args.contype, probtype=args.probtype,
                                 num_eigenvalues=args.num_eigenvalues,
                                 m_fixed=m_fixed, r0_frac=args.r0_frac, N=args.N,
                                 iter_offset=iter_offset)

        # Get the assembler object
        assembler = problem.getAssembler()

        # Set the callback for generating output
        cb = OutputCallback(args.prefix, assembler, iter_offset=iter_offset, freq=5)
        problem.setOutputCallback(cb.write_output)

        # Keep counting the total number of iterations
        iter_offset += args.max_opt_iters
        problem.setPrefix(args.prefix)

        # Check the gradient
        problem.checkGradients(1e-6)

        # Test the element implementation
        if comm.rank == 0:
            assembler.testElement(0, 2)

        # Extract the filter to interpolate design variables
        filtr = problem.getFilter()

        if orig_filter is not None:
            # Create one of the new design vectors
            x = problem.createDesignVec()
            TopOptUtils.interpolateDesignVec(orig_filter, xopt, filtr, x)
            problem.setInitDesignVars(x)

        # Set the new original filter
        orig_filter = filtr

        # Set parameters
        if step == max_iterations-1:
            optimization_options['maxiter'] = 4*args.max_opt_iters
        else:
            optimization_options['maxiter'] = args.max_opt_iters
        optimization_options['output_file'] = os.path.join(args.prefix,
                                                           'output_file%d.dat'%(step))
        optimization_options['tr_output_file'] = os.path.join(args.prefix,
                                                              'tr_output_file%d.dat'%(step))

        # Optimize the problem
        opt = TopOptUtils.TopologyOptimizer(problem, optimization_options)
        xopt = opt.optimize()

        # Refine based solely on the value of the density variable
        write_dvs_to_file(xopt, assembler, os.path.join(args.prefix, 'dv_output%d.f5'%(step)))

        if step == max_iterations-2:
            # Duplicate the forest before the refinement process - this is to avoid
            # refinement on the existing forest which will be needed again
            forest = forest.duplicate()

            # Perform refinement based on distance
            dist_file = os.path.join(args.prefix, 'distance_solution%d.f5'%(step))
            refine_distance = 0.025*domain_length
            TopOptUtils.targetRefine(forest, filtr, assembler, refine_distance,
                                     interface_lev=args.init_depth+1, interior_lev=args.init_depth,
                                     domain_length=domain_length, filename=dist_file)

            # Repartition the mesh
            forest.balance(1)
            forest.repartition()

        # If there is projection, use it after the first iteration only
        if args.use_project:
            props.use_project = True
