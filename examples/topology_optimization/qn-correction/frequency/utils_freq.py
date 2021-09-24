from tmr import TMR, TopOptUtils
from tacs import TACS, elements, functions
from egads4py import egads
import numpy as np
import openmdao.api as om
import os
import sys

sys.path.append('../eigenvalue')
from utils import OctCreator, CreatorCallback, MFilterCreator, OutputCallback

# Print colored text in terminal
try:
    from termcolor import colored
except:
    pass

class FrequencyConstr:
    """
    A class that evaluates the smallest eigenvalue, the objective is evaluated
    using an objective callback. We also add non-design mass to loaded nodes in
    order to form a well-posed mass minimization problem under frequency constriant

    this constraint takes the following form:
        c = ks >= 0
    """

    itr = 0

    def __init__(self, prefix, domain, forest, len0, AR, ratio,
                 iter_offset, lambda0, eig_scale=1.0, num_eigenvalues=10,
                 eig_method='jd', max_jd_size=100, max_gmres_size=30,
                 max_lanczos=60, ksrho=50, non_design_mass=5.0, jd_use_recycle=True,
                 lanczos_shift=-10):
        """
        Args:
            eig_scale: scale the eigenvalues internally in order to acquire better
                       KS approximation with smaller skrho
            num_eigenvalues: number of smallest eigenvalues to compute
            ksrho: KS parameter
        """

        # Check input
        if eig_method != 'jd' and eig_method != 'lanczos':
            raise ValueError("Invalid method for eigensolver.")
        self.eig_method = eig_method

        # Set objects
        self.forest = forest

        # Set up parameters
        self.prefix = prefix
        self.domain = domain
        self.iter_offset = iter_offset
        self.lx = len0*AR
        self.ly = len0
        self.lz = len0
        if domain == 'lbracket':
            self.ly = len0*ratio
        self.ratio = ratio
        self.lambda0 = lambda0
        self.eig_scale = eig_scale
        self.num_eigenvalues = num_eigenvalues
        self.max_jd_size = max_jd_size
        self.max_lanczos = max_lanczos
        self.max_gmres_size = max_gmres_size
        self.ksrho = ksrho
        self.non_design_mass = non_design_mass
        self.jd_use_recycle = jd_use_recycle

        self.lanczos_shift = lanczos_shift

        self.fltr = None
        self.mg = None
        self.assembler = None
        self.comm = None
        self.oper = None
        self.jd = None
        self.lanczos = None
        self.eig = None

        # TACS Vectors
        self.eigv = None
        self.rho = None
        self.rho_original = None
        self.update = None
        self.temp = None
        self.mvec = None

        # TACS Matrices
        self.mmat = None
        self.m0mat = None
        self.Amat = None

        # We keep track of failed qn correction
        self.curvs = []

        return

    def constraint(self, fltr, mg):
        """
        Evaluate the KS aggregation of the smallest eigenvalue for the generalized
        eigenvalue problem:

        ks = -1.0/rho * ln tr exp(-rho*A)

        """

        if self.fltr is None:
            self.mg = mg
            self.fltr = fltr
            self.assembler = self.fltr.getAssembler()
            self.svec = self.assembler.createDesignVec()
            self.comm = self.assembler.getMPIComm()
            self.rank = self.comm.rank

            # Initialize space for matrices and vectors
            self.mmat = self.assembler.createMat()
            self.m0mat = self.assembler.createMat()  # For non design mass
            self.Amat = self.assembler.createMat()
            self.eig = np.zeros(self.num_eigenvalues)
            self.eigv = []
            for i in range(self.num_eigenvalues):
                self.eigv.append(self.assembler.createVec())
            self.deig = []
            for i in range(self.num_eigenvalues):
                self.deig.append(self.assembler.createDesignVec())

            # Allocate vectors for qn correction
            self.rho = self.assembler.createDesignVec()
            self.rho_original = self.assembler.createDesignVec()
            self.update = self.assembler.createDesignVec()
            self.temp = self.assembler.createDesignVec()

            # Set up Jacobi-Davidson eigensolver
            if self.eig_method == 'jd':

                # Create the operator with given matrix and multigrid preconditioner
                self.oper = TACS.JDSimpleOperator(self.assembler, self.Amat, self.mg)

                # Create the eigenvalue solver and set the number of recycling eigenvectors
                self.jd = TACS.JacobiDavidson(self.oper, self.num_eigenvalues,
                                              self.max_jd_size, self.max_gmres_size)
                self.jd.setTolerances(eig_rtol=1e-6, eig_atol=1e-6, rtol=1e-6, atol=1e-12)
                self.jd.setThetaCutoff(0.01)

            # Set up Shift-and-invert Lanczos method
            else:
                assert(self.eig_method == 'lanczos')

                # Set up linear solver
                gmres_iters = 15
                nrestart = 0
                is_flexible = 0
                self.ksm_Amat = self.assembler.createMat()  # We need a separate matrix to apply shift
                ksm = TACS.KSM(self.ksm_Amat, self.mg, gmres_iters, nrestart, is_flexible)

                # Set up lanczos solver
                eig_tol = 1e-6
                ep_oper = TACS.EPShiftInvertOp(self.lanczos_shift, ksm)
                self.sep = TACS.SEPsolver(ep_oper, self.max_lanczos,
                                          TACS.SEP_FULL, self.assembler.getBcMap())
                self.sep.setTolerances(eig_tol, TACS.SEP_SMALLEST, self.num_eigenvalues)

            '''
            Create a non-design mass matrix
            '''
            self.mvec = self.assembler.createDesignVec()
            mvals = self.mvec.getArray()

            # Get nodal locations
            Xpts = self.forest.getPoints()

            # Note: the local nodes are organized as follows:
            # |--- dependent nodes -- | ext_pre | -- owned local -- | - ext_post -|

            # Get number of local nodes in the current processor
            n_local_nodes = Xpts.shape[0]

            # Get numbder of dependent nodes
            _ptr, _conn, _weights = self.forest.getDepNodeConn()

            # Get number of ext_pre nodes
            n_ext_pre = self.forest.getExtPreOffset()

            # Get numbder of own nodes:
            offset = n_ext_pre

            # # Loop over all owned nodes and set non-design mass values
            tol = 1e-6
            if self.domain == 'cantilever':
                xmin = self.lx - tol
                xmax = self.lx + tol
                ymin = 0.25*self.ly - tol
                ymax = 0.75*self.ly + tol
                zmin = 0.0*self.lz - tol
                zmax = 0.2*self.lz + tol

            elif self.domain == 'michell':
                xmin = self.lx - tol
                xmax = self.lx + tol
                ymin = 0.25*self.ly - tol
                ymax = 0.75*self.ly + tol
                zmin = 0.4*self.lz - tol
                zmax = 0.6*self.lz + tol

            elif self.domain == 'mbb':
                xmin = 0.0*self.lx - tol
                xmax = 0.2*self.lx + tol
                ymin = 0.25*self.ly - tol
                ymax = 0.75*self.ly + tol
                zmin = self.lz - tol
                zmax = self.lz + tol

            elif self.domain == 'lbracket':
                RATIO = self.ratio
                xmin = self.lx - tol
                xmax = self.lx + tol
                ymin = 0.25*self.ly - tol
                ymax = 0.75*self.ly + tol
                zmin = 0.5*RATIO*self.lz - tol
                zmax = 1.0*RATIO*self.lz + tol

            else:
                print("[Warning]Unsupported domain type for non-design mass!")

            for i in range(offset, n_local_nodes):
                x, y, z = Xpts[i]
                if xmin < x < xmax:
                    if ymin < y < ymax:
                        if zmin < z < zmax:
                            mvals[i-offset] = 1.0

            # Assemble a constant non-design mass matrix
            dv = self.assembler.createDesignVec()
            self.assembler.getDesignVars(dv)
            self.assembler.setDesignVars(self.mvec)
            self.assembler.assembleMatType(TACS.MASS_MATRIX, self.m0mat)
            self.m0mat.scale(self.non_design_mass)
            self.assembler.setDesignVars(dv)

        # Assemble the mass matrix
        self.assembler.assembleMatType(TACS.MASS_MATRIX, self.mmat)
        self.mmat.axpy(1.0, self.m0mat)
        self.assembler.applyMatBCs(self.mmat)

        # Reference the matrix associated with the multigrid preconditioner
        # We finally want:
        # mgmat = K - 0.95*lambda0*M
        # here we assemble the K part first
        self.assembler.assembleMatType(TACS.STIFFNESS_MATRIX, self.Amat)
        mgmat = self.mg.getMat()
        mgmat.copyValues(self.Amat)

        # Assemble A matrix for the simple eigenvalue problem and apply bcs
        # A = K - lambda0*M
        self.Amat.axpy(-self.lambda0, self.mmat)
        self.assembler.applyMatBCs(self.Amat)

        # Finish assembling the multigrid preconditioner
        mgmat.axpy(-0.95*self.lambda0, self.mmat)
        self.assembler.applyMatBCs(mgmat)

        # Factor the multigrid preconditioner
        self.mg.assembleGalerkinMat()
        self.mg.factor()

        """
        Solve the eigenvalue problem
        """
        if self.eig_method == 'jd':
            if self.jd_use_recycle:
                num_recycle = self.num_eigenvalues
            else:
                num_recycle = 0
            self.jd.setRecycle(num_recycle)
            if self.comm.rank == 0:
                print("[JD Recycle] JD solver is recycling {:d} eigenvectors...".format(num_recycle))
            self.jd.solve(print_flag=True, print_level=1)

            # Check if succeeded, otherwise try again
            nconvd = self.jd.getNumConvergedEigenvalues()
            if nconvd < self.num_eigenvalues:
                if self.comm.rank == 0:
                    print("[Warning] Jacobi-Davidson failed to converge"
                        " for the first run, starting rerun...")

                if self.jd_use_recycle:
                    self.jd.setRecycle(nconvd)

                # Modify the matrix associated with the preconditioner
                # so that the matrix is positive definite
                I = self.assembler.createMat()
                I.addDiag(1.0) # Create an identity matrix

                # Update mgmat so that it's positive definite
                eig0, err = self.jd.extractEigenvalue(0)
                if eig0 > 0:
                    if self.comm.rank == 0:
                        print("[mgmat] Smallest eigenvalue is already positive, don't update mgmat!")
                    else:
                        mgmat.axpy(-eig0, I)
                        self.assembler.applyMatBCs(mgmat)
                        self.mg.assembleGalerkinMat()
                        self.mg.factor()

                # Rerun the solver
                self.jd.solve(print_flag=True, print_level=1)
                nconvd = self.jd.getNumConvergedEigenvalues()

                # If it still fails, raise error, save fail f5 and exit
                if nconvd < self.num_eigenvalues:
                    msg = "No enough eigenvalues converged! ({:d}/{:d})".format(
                        nconvd, self.num_eigenvalues)

                    # set the unconverged eigenvector as state variable for visualization
                    for i in range(self.num_eigenvalues):
                        self.eig[i], error = self.jd.extractEigenvector(i, self.eigv[i])
                    self.assembler.setVariables(self.eigv[nconvd])

                    flag_fail = (TACS.OUTPUT_CONNECTIVITY |
                                TACS.OUTPUT_NODES |
                                TACS.OUTPUT_DISPLACEMENTS |
                                TACS.OUTPUT_EXTRAS)
                    f5_fail = TACS.ToFH5(self.assembler, TACS.SOLID_ELEMENT, flag_fail)
                    f5_fail.writeToFile(os.path.join(self.prefix, "fail.f5"))

                    raise ValueError(msg)

            # Extract eigenvalues and eigenvectors
            for i in range(self.num_eigenvalues):
                self.eig[i], error = self.jd.extractEigenvector(i, self.eigv[i])

        elif self.eig_method == 'lanczos':
            # For shift-and-invert lanczos, we need to apply shift to
            # both multigrid preconditioner and the underlying matrix
            # of the Krylov subspace solver
            self.ksm_Amat.copyValues(self.Amat)
            self.ksm_Amat.addDiag(-self.lanczos_shift)
            mgmat.addDiag(-self.lanczos_shift)
            self.sep.solve(self.comm, print_flag=True)
            for i in range(self.num_eigenvalues):
                self.eig[i], error = self.sep.extractEigenvector(i, self.eigv[i])

        else:
            raise ValueError("Invalid eig_method")

        # Debug: print out residuals
        debug_initialized = 0
        if debug_initialized == 0:
            debug_initialized = 1
            res = self.assembler.createVec()
            one = self.assembler.createVec()
            Av  = self.assembler.createVec()
            one_arr = one.getArray()
            one_arr[:] = 1.0
            debug_counter = 0

            residual  = np.zeros(self.num_eigenvalues)
            eigvec_l1 = np.zeros(self.num_eigenvalues)
            eigvec_l2 = np.zeros(self.num_eigenvalues)

        for i in range(self.num_eigenvalues):
            self.Amat.mult(self.eigv[i], Av)  # Compute Av

            res.copyValues(Av)
            res.axpy(-self.eig[i], self.eigv[i]) # Compute res = Av - lambda*v

            residual[i] = res.norm()
            eigvec_l1[i] = self.eigv[i].dot(one) # Compute l1 norm
            eigvec_l2[i] = self.eigv[i].dot(self.eigv[i])**0.5  # Compute l2 norm

        debug_counter += 1
        if self.assembler.getMPIComm().rank == 0:
            print("Optimization iteration:{:4d}".format(debug_counter))
            print("{:4s}{:15s}{:15s}{:15s}".format("No", "Eig Res", "Eigv l1 norm", "Eigv l2 norm"))
            for i in range(self.num_eigenvalues):
                print("{:4d}{:15.5e}{:15.5e}{:15.5e}".format(i, residual[i], eigvec_l1[i], eigvec_l2[i]))

        # Set first eigenvector as state variable for visualization
        self.assembler.setVariables(self.eigv[0])

        # Scale eigenvalues for a better KS approximation
        self.eig[:] *= self.eig_scale

        # Compute the minimal eigenvalue
        eig_min = np.min(self.eig)

        # Compute KS aggregation
        self.eta = np.exp(-self.ksrho*(self.eig - eig_min))
        self.beta = np.sum(self.eta)
        ks = (eig_min - np.log(self.beta)/self.ksrho)
        self.eta = self.eta/self.beta

        # Scale eigenvalue back
        self.eig[:] /= self.eig_scale

        # Print values
        if self.comm.rank == 0:
            print('{:30s}{:20.10e}'.format('[Constr] KS eigenvalue:', ks))
            print('{:30s}{:20.10e}'.format('[Constr] min eigenvalue:', eig_min))

        return [ks]


    def constraint_gradient(self, fltr, mg, vecs):
        """
        gradient of the spectral aggregate
        g = sum_k eta_k phi^T dAdx phi
        """

        # We only have one constraint
        dcdrho = vecs[0]

        # Zero out the gradient vector
        dcdrho.zeroEntries()

        for i in range(self.num_eigenvalues):

            # Compute the coefficient
            coeff = self.eta[i]*self.eig_scale

            # Compute gradient of eigenvalue
            self.deig[i].zeroEntries()
            self.assembler.addMatDVSensInnerProduct(
                coeff, TACS.STIFFNESS_MATRIX,
                self.eigv[i], self.eigv[i], self.deig[i])

            self.assembler.addMatDVSensInnerProduct(
                -coeff*self.lambda0, TACS.MASS_MATRIX,
                self.eigv[i], self.eigv[i], self.deig[i])

            # Make sure the vector is properly distributed over all processors
            self.deig[i].beginSetValues(op=TACS.ADD_VALUES)
            self.deig[i].endSetValues(op=TACS.ADD_VALUES)

            # Add the contribution
            dcdrho.axpy(1.0, self.deig[i])

        # Compute gradient norm
        norm = dcdrho.norm()
        if self.comm.rank == 0:
            print("{:30s}{:20.10e}".format('[Constr] gradient norm:', norm))

        self.dcdrho = dcdrho
        return

    def qn_correction(self, x, z, zw, s, y):
        """
        Update y:
        y <- y + z*F^T P Fs

        where:
        F: filter matrix
        P: Positive definite part of the constraint Hessian

        Note:
        x is raw design variable (unfiltered) and it's NOT equal to
        the design variable in the assembler:

        if:
        self.assembler.getDesignVars(dv)

        Then:
        dv == Fx

        Args:
            x (PVec): unfiltered design vector (not used because we get x directly from assembler)
            s (PVec): unfiltered update step
            y (PVec): y = Bs
            z, zw: multipliers for dense and sparse constraints, only z is used here
        """

        # Finite difference step length for computing second order
        # derivative of stiffness matrix
        h = 1e-8

        # Zero out the update vector for y
        self.update.zeroEntries()

        # Get current nodal density
        self.assembler.getDesignVars(self.rho)
        self.rho_original.copyValues(self.rho)

        # Apply filter to step vector and perturb rho
        self.svec.zeroEntries()
        self.fltr.applyFilter(TMR.convertPVecToVec(s), self.svec)
        self.rho.axpy(h, self.svec)

        # Zero out temp vector to store gradient information
        self.temp.zeroEntries()

        # set density = rho + h*s
        self.assembler.setDesignVars(self.rho)

        for i in range(self.num_eigenvalues):
            # Compute g(rho + h*s) for d2Kdx2
            coeff1 = self.eta[i]*self.eig_scale
            self.assembler.addMatDVSensInnerProduct(coeff1, TACS.STIFFNESS_MATRIX,
                self.eigv[i], self.eigv[i], self.temp)

            # Compute g(rho + h*s) for d2Mdx2
            coeff2 = -self.eta[i]*self.lambda0*self.eig_scale
            self.assembler.addMatDVSensInnerProduct(coeff2, TACS.MASS_MATRIX,
                self.eigv[i], self.eigv[i], self.temp)

        # set density = rho
        self.assembler.setDesignVars(self.rho_original)

        for i in range(self.num_eigenvalues):
            # Compute g(rho + h*s) - g(rho) for d2Kdx2
            coeff1 = self.eta[i]*self.eig_scale
            self.assembler.addMatDVSensInnerProduct(-coeff1, TACS.STIFFNESS_MATRIX,
                self.eigv[i], self.eigv[i], self.temp)

            # Compute g(rho + h*s) - g(rho) for d2Mdx2
            coeff2 = -self.eta[i]*self.lambda0*self.eig_scale
            self.assembler.addMatDVSensInnerProduct(-coeff2, TACS.MASS_MATRIX,
                self.eigv[i], self.eigv[i], self.temp)

        # Distribute the temp vector
        self.temp.beginSetValues(op=TACS.ADD_VALUES)
        self.temp.endSetValues(op=TACS.ADD_VALUES)

        # Add to the update
        self.update.axpy(1.0/h, self.temp)

        # Compute norms for output
        x_norm = x.norm()
        s_norm = s.norm()
        y_norm = y.norm()
        dy_norm = self.update.norm()
        curvature = self.svec.dot(self.update)
        if self.comm.rank == 0:
            if curvature < 0:
                try:
                    print(colored("curvature: {:20.10e}".format(curvature), "red"))
                except:
                    print("curvature: {:20.10e}".format(curvature))
            else:
                try:
                    print(colored("curvature: {:20.10e}".format(curvature), "green"))
                except:
                    print("curvature: {:20.10e}".format(curvature))

            print("norm(x):   {:20.10e}".format(x_norm))
            print("norm(s):   {:20.10e}".format(s_norm))
            print("norm(y):   {:20.10e}".format(y_norm))
            print("norm(dy):  {:20.10e}".format(dy_norm))

        # Update y
        y_wrap = TMR.convertPVecToVec(y)

        if curvature > 0:
            self.fltr.applyTranspose(self.update, self.update)
            y_wrap.axpy(z[0], self.update)

        self.curvs.append(curvature)

        return

    def getQnUpdateCurvs(self):
        return self.curvs

class MassObj:
    """
    Mass objective takes the following form:
        obj = m / m_fixed
    """

    def __init__(self, m_fixed, comm):

        self.m_fixed = m_fixed
        self.comm = comm
        self.rank = self.comm.Get_rank()

        self.assembler = None
        self.fltr = None
        self.mass_func = None

        return

    def objective(self, fltr, mg):
        if self.fltr is None:
            self.fltr = fltr
            self.assembler = self.fltr.getAssembler()
            self.mass_func = functions.StructuralMass(self.assembler)

        # Eval mass
        mass = self.assembler.evalFunctions([self.mass_func])[0]
        obj = mass / self.m_fixed
        if self.rank == 0:
            print("{:30s}{:20.10e}".format('[Obj] mass objective:',obj))

        return obj

    def objective_gradient(self, fltr, mg, dfdrho):
        # We only have one constraint
        dfdrho.zeroEntries()

        # Evaluate the mass gradient
        self.assembler.addDVSens([self.mass_func], [dfdrho], alpha=1.0/self.m_fixed)

        # Compute norm
        norm = dfdrho.norm()
        if self.rank == 0:
            print("{:30s}{:20.10e}".format('[Con] gradient norm:', norm))
        return

def create_problem(prefix, domain, forest, bcs, props, nlevels, lambda0, ksrho,
                   vol_frac=0.25, r0_frac=0.05, len0=1.0, AR=1.0, ratio=0.4,
                   density=2600.0, iter_offset=0, qn_correction=True,
                   non_design_mass=5.0, eig_scale=1.0, eq_constr=False,
                   num_eigenvalues=10, eig_method='jd', max_jd_size=100, jd_use_recycle=True,
                   max_gmres_size=30, max_lanczos=60, lanczos_shift=-10):
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
        density (float): Density to use for the mass computation
        iter_offset (int): iteration counter offset

    Returns:
        TopoProblem: Topology optimization problem instance
    """

    # Create the problem and filter object
    N = 20
    mfilter = MFilterCreator(r0_frac, N, a=len0)
    filter_type = mfilter.filter_callback
    obj = CreatorCallback(bcs, props)
    problem = TopOptUtils.createTopoProblem(forest, obj.creator_callback,
                                            filter_type, use_galerkin=True,
                                            nlevels=nlevels)

    # Get the assembler object we just created
    assembler = problem.getAssembler()

    # Compute the fixed mass target
    lx = len0*AR # mm
    ly = len0 # mm
    lz = len0 # mm
    if domain == 'lbracket':
        ly = len0*ratio
    vol = lx*ly*lz
    if domain == 'lbracket':
        S1 = lx*lz
        S2 = lx*lz*(1.0-ratio)**2
        vol = (S1-S2)*ly
    m_fixed = vol_frac*(vol*density)

    # Add objective callback
    obj_callback = MassObj(m_fixed, assembler.getMPIComm())
    problem.addObjectiveCallback(obj_callback.objective,
                                 obj_callback.objective_gradient)

    # Add constraint callback
    constr_callback = FrequencyConstr(prefix, domain, forest, len0, AR, ratio,
                                     iter_offset, lambda0,
                                     ksrho=ksrho,
                                     non_design_mass=non_design_mass,
                                     eig_scale=eig_scale,
                                     num_eigenvalues=num_eigenvalues,
                                     eig_method=eig_method,
                                     max_lanczos=max_lanczos,
                                     max_jd_size=max_jd_size,
                                     max_gmres_size=max_gmres_size,
                                     jd_use_recycle=jd_use_recycle,
                                     lanczos_shift=lanczos_shift)
    nineq = 1
    if eq_constr is True:
        nineq = 0
    problem.addConstraintCallback(1, nineq, constr_callback.constraint,
                                  constr_callback.constraint_gradient)

    # Use Quasi-Newton Update Correction if specified
    if qn_correction:
        problem.addQnCorrectionCallback(1, constr_callback.qn_correction)

    # Set output callback
    cb = OutputCallback(assembler, iter_offset=iter_offset)
    problem.setOutputCallback(cb.write_output)

    return problem, obj_callback, constr_callback
