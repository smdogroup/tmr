from tmr import TMR, TopOptUtils
from tacs import TACS, elements, functions
from paropt import ParOpt
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
                 max_lanczos=60, ksrho=50, jd_use_recycle=True,
                 lanczos_shift=-10.0, jd_use_Amat_shift=False):
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
        self.jd_use_recycle = jd_use_recycle

        self.lanczos_shift = lanczos_shift
        self.jd_use_Amat_shift = jd_use_Amat_shift
        self.old_min_eigval = 0.0

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
        # self.m0mat = None
        # self.k0mat = None
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
            # self.m0mat = self.assembler.createMat()  # For non design mass
            # self.k0mat = self.assembler.createMat()  # For non design stiffness
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
            self.update_vals = self.update.getArray()
            self.temp = self.assembler.createDesignVec()
            self.temp_vals = self.temp.getArray()

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

        # Assemble the mass matrix
        self.assembler.assembleMatType(TACS.MASS_MATRIX, self.mmat)
        # self.mmat.axpy(1.0, self.m0mat)
        self.assembler.applyMatBCs(self.mmat)

        # Reference the matrix associated with the multigrid preconditioner
        # We finally want:
        # mgmat = K - 0.95*lambda0*M
        # here we assemble the K part first
        self.assembler.assembleMatType(TACS.STIFFNESS_MATRIX, self.Amat)
        # self.Amat.axpy(1.0, self.k0mat)
        self.assembler.applyMatBCs(self.Amat)
        mgmat = self.mg.getMat()
        mgmat.copyValues(self.Amat)

        # Assemble A matrix for the simple eigenvalue problem and apply bcs
        # A = K - lambda0*M
        self.Amat.axpy(-self.lambda0, self.mmat)

        # We may shift the eigenvalues of A by adding value to diagonal entries:
        # A <- A - (old-1.0)*I
        # so that all eigenvalues for A are likely to be positive (and no smaller than 1.0)
        if self.jd_use_Amat_shift:
            self.Amat.addDiag(1.0 - self.old_min_eigval)

        # Apply bcs to A
        self.assembler.applyMatBCs(self.Amat)

        # Finish assembling the multigrid preconditioner
        mgmat.axpy(-0.95*self.lambda0, self.mmat)
        if self.jd_use_Amat_shift:
            mgmat.addDiag(1.0 - self.old_min_eigval)
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

                # Update mgmat so that it's positive definite
                eig0, err = self.jd.extractEigenvalue(0)
                if eig0 > 0:
                    if self.comm.rank == 0:
                        print("[mgmat] Smallest eigenvalue is already positive, don't update mgmat!")
                else:
                    mgmat.addDiag(-eig0)
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

            # Adjust eigenvalues and shift matrices back
            if self.jd_use_Amat_shift:
                for i in range(self.num_eigenvalues):
                    self.eig[i] += (self.old_min_eigval - 1.0)
                self.Amat.addDiag(self.old_min_eigval - 1.0)
                self.assembler.applyMatBCs(self.Amat)

                # Set the shift value for next optimization iteration
                self.old_min_eigval = self.eig[0]  # smallest eigenvalue

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

    def qn_correction(self, zero_idx, z, s, y):
        """
        Update y:
        y <- y + z*F^T P Fs

        where:
        F: filter matrix
        P: Positive definite part of the constraint Hessian

        Note:
        x, s correspond raw design variable (unfiltered) and it's NOT equal to
        the design variable in the assembler:

        if:
        self.assembler.getDesignVars(dv)

        Then:
        dv == Fx

        Inputs:
            zero_idx (int list): indices to-be-zeroed in order to compute
                                 the Hessian-vector product for reduced problem
            s (PVec): unfiltered update step
            z (array-like): multipliers for dense constraints


        Outputs:
            y (PVec): y = Bs
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
            # Zero entries in temp, if needed
            if zero_idx:
                self.temp_vals[zero_idx] = 0.0

            # Compute g(rho + h*s) for d2Kdx2
            coeff1 = self.eta[i]*self.eig_scale
            self.assembler.addMatDVSensInnerProduct(coeff1, TACS.STIFFNESS_MATRIX,
                self.eigv[i], self.eigv[i], self.temp)

            # Zero entries in temp, if needed
            if zero_idx:
                self.temp_vals[zero_idx] = 0.0

            # Compute g(rho + h*s) for d2Mdx2
            coeff2 = -self.eta[i]*self.lambda0*self.eig_scale
            self.assembler.addMatDVSensInnerProduct(coeff2, TACS.MASS_MATRIX,
                self.eigv[i], self.eigv[i], self.temp)

        # set density = rho
        self.assembler.setDesignVars(self.rho_original)

        for i in range(self.num_eigenvalues):
            # Zero entries in temp, if needed
            if zero_idx:
                self.temp_vals[zero_idx] = 0.0

            # Compute g(rho + h*s) - g(rho) for d2Kdx2
            coeff1 = self.eta[i]*self.eig_scale
            self.assembler.addMatDVSensInnerProduct(-coeff1, TACS.STIFFNESS_MATRIX,
                self.eigv[i], self.eigv[i], self.temp)

            # Zero entries in temp, if needed
            if zero_idx:
                self.temp_vals[zero_idx] = 0.0

            # Compute g(rho + h*s) - g(rho) for d2Mdx2
            coeff2 = -self.eta[i]*self.lambda0*self.eig_scale
            self.assembler.addMatDVSensInnerProduct(-coeff2, TACS.MASS_MATRIX,
                self.eigv[i], self.eigv[i], self.temp)

        # Distribute the temp vector
        self.temp.beginSetValues(op=TACS.ADD_VALUES)
        self.temp.endSetValues(op=TACS.ADD_VALUES)

        # Zero entries in temp, if needed
        if zero_idx:
            self.temp_vals[zero_idx] = 0.0

        # Add to the update
        self.update.axpy(1.0/h, self.temp)
        if zero_idx:
            self.update_vals[zero_idx] = 0.0

        # Compute norms for output
        rho_norm = self.rho.norm()
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

            print("norm(rho):   {:20.10e}".format(rho_norm))
            print("norm(s):   {:20.10e}".format(s_norm))
            print("norm(y):   {:20.10e}".format(y_norm))
            print("norm(dy):  {:20.10e}".format(dy_norm))

        # Update y
        y_wrap = TMR.convertPVecToVec(y)

        if curvature > 0:
            self.fltr.applyTranspose(self.update, self.update)
            if zero_idx:
                self.update_vals[zero_idx] = 0.0
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
                   density=2600.0, iter_offset=0,
                   eig_scale=1.0, eq_constr=False,
                   num_eigenvalues=10, eig_method='jd', max_jd_size=100, jd_use_recycle=True,
                   max_gmres_size=30, max_lanczos=60, lanczos_shift=-10, jd_use_Amat_shift=False):
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
                                     eig_scale=eig_scale,
                                     num_eigenvalues=num_eigenvalues,
                                     eig_method=eig_method,
                                     max_lanczos=max_lanczos,
                                     max_jd_size=max_jd_size,
                                     max_gmres_size=max_gmres_size,
                                     jd_use_recycle=jd_use_recycle,
                                     lanczos_shift=lanczos_shift,
                                     jd_use_Amat_shift=jd_use_Amat_shift)
    nineq = 1
    if eq_constr is True:
        nineq = 0
    problem.addConstraintCallback(1, nineq, constr_callback.constraint,
                                  constr_callback.constraint_gradient)

    # Set output callback
    cb = OutputCallback(assembler, iter_offset=iter_offset)
    problem.setOutputCallback(cb.write_output)

    return problem, obj_callback, constr_callback

class ReducedProblem(ParOpt.Problem):
    """
    A reduced problem by fixing some design variables in the original problem
    """

    def __init__(self, original_prob, fixed_dv_idx : list,
        fixed_dv_val=1.0, qn_correction_func=None):

        self.prob = original_prob
        self.assembler = self.prob.getAssembler()
        self.comm = self.assembler.getMPIComm()
        self.ncon = 1  # Hard-coded for now
        self.qn_correction_func = qn_correction_func

        # Allocate full-size vectors for the original problem
        self._x = self.prob.createDesignVec()
        self._g = self.prob.createDesignVec()
        self._A = []
        for i in range(self.ncon):
            self._A.append(self.prob.createDesignVec())

        # Allocate helper variables for the qn correction, if needed
        if self.qn_correction_func:
            self._s = self.prob.createDesignVec()
            self._y = self.prob.createDesignVec()

        # Get indices of fixed design variables, these indices
        # are with respect to the original full-sized problem
        self.fixed_dv_idx = fixed_dv_idx
        self.fixed_dv_val = fixed_dv_val

        # Compute the indices of fixed design variables, these indices
        # are with respect to the original full-sized problem
        self.free_dv_idx = [i for i in range(len(self._x)) if i not in self.fixed_dv_idx]
        self.nvars = len(self.free_dv_idx)

        super().__init__(self.comm, self.nvars, self.ncon)
        return

    def getNumCons(self):
        return self.ncon

    def getVarsAndBounds(self, x, lb, ub):

        # # Set bounds
        # x[:] = self._x0[self.free_dv_idx]
        # lb[:] = self._lb[self.free_dv_idx]
        # ub[:] = self._ub[self.free_dv_idx]
        x[:] = 0.95
        lb[:] = 1e-3
        ub[:] = 1.0
        return

    def evalObjCon(self, x):
        # Populate full-sized design variable
        self.reduDVtoDV(x, self._x)

        # Run analysis in full-sized problem
        fail, fobj, con = self.prob.evalObjCon(self.ncon, self._x)

        return fail, fobj, con

    def evalObjConGradient(self, x, g, A):
        # Populate full-sized design variable
        self.reduDVtoDV(x, self._x)

        # Run analysis in full-sized problem to get gradient
        fail = self.prob.evalObjConGradient(self._x, self._g, self._A)

        # Get reduced gradient and constraint jacobian
        self.DVtoreduDV(self._g, g)
        for i in range(self.ncon):
            self.DVtoreduDV(self._A[i], A[i])

        return fail

    def reduDVtoDV(self, reduDV, DV, fixed_val=None):
        '''
        Convert the reduced design vector to full-sized design vector
        '''
        if fixed_val is None:
            val = self.fixed_dv_val
        else:
            val = fixed_val

        if self.fixed_dv_idx:
            DV[self.fixed_dv_idx] = val
        DV[self.free_dv_idx] = reduDV[:]

        return

    def DVtoreduDV(self, DV, reduDV):
        '''
        Convert the full-sized design vector to reduced design vector
        '''
        reduDV[:] = DV[self.free_dv_idx]

        return

    def computeQuasiNewtonUpdateCorrection(self, x, z, zw, s, y):
        if self.qn_correction_func:
            # Populate full-sized vectors
            self.reduDVtoDV(s, self._s, fixed_val=0.0)
            self.reduDVtoDV(y, self._y, fixed_val=0.0)

            # Update y and copy back
            self.qn_correction_func(self.fixed_dv_idx, z, self._s, self._y)
            self.DVtoreduDV(self._y, y)
        return

def getFixedDVIndices(forest, domain, len0, AR, ratio):
    """
    Get indices for fixed design variables
    """

    fixed_dv_idx = []

    # Compute geometric parameters
    lx = len0*AR
    ly = len0
    lz = len0
    if domain == 'lbracket':
        ly = len0*ratio

    # Get nodal locations
    Xpts = forest.getPoints()

    # Note: the local nodes are organized as follows:
    # |--- dependent nodes -- | ext_pre | -- owned local -- | - ext_post -|

    # Get number of local nodes in the current processor
    n_local_nodes = Xpts.shape[0]

    # Get numbder of dependent nodes
    _ptr, _conn, _weights = forest.getDepNodeConn()

    # Get number of ext_pre nodes
    n_ext_pre = forest.getExtPreOffset()

    # Get numbder of own nodes:
    offset = n_ext_pre

    # # Loop over all owned nodes and set non-design mass values
    tol = 1e-6 # Make sure our ranges are inclusive
    depth = 0.1  # depth for non-design mass
    if domain == 'cantilever':
        xmin = (1-depth)*lx - tol
        xmax = lx + tol
        ymin = 0.25*ly - tol
        ymax = 0.75*ly + tol
        zmin = 0.0*lz - tol
        zmax = 0.2*lz + tol

    elif domain == 'michell':
        xmin = (1-depth)*lx - tol
        xmax = lx + tol
        ymin = 0.25*ly - tol
        ymax = 0.75*ly + tol
        zmin = 0.4*lz - tol
        zmax = 0.6*lz + tol

    elif domain == 'mbb':
        xmin = 0.0*lx - tol
        xmax = 0.2*lx + tol
        ymin = 0.25*ly - tol
        ymax = 0.75*ly + tol
        zmin = (1-depth)*lz - tol
        zmax = lz + tol

    elif domain == 'lbracket':
        xmin = (1-depth)*lx - tol
        xmax = lx + tol
        ymin = 0.25*ly - tol
        ymax = 0.75*ly + tol
        zmin = 0.5*ratio*lz - tol
        zmax = 1.0*ratio*lz + tol

    else:
        raise ValueError("[Error]Unsupported domain type for non-design mass!")

    for i in range(offset, n_local_nodes):
        x, y, z = Xpts[i]
        if xmin < x < xmax:
            if ymin < y < ymax:
                if zmin < z < zmax:
                    fixed_dv_idx.append(i-offset)

    return fixed_dv_idx

class ReduOmAnalysis(om.ExplicitComponent):
    '''
    This class wraps the analyses with openmdao interface such that
    the optimization can be run with different optimizers such as
    SNOPT and IPOPT.
    Note that the design/gradient vectors manipulated in this class
    are all global vectors. Local components can be queried by:

    local_vec = global_vec[start:end]

    where:
    start = self.offsets[rank]
    end = start + self.sizes[rank]
    '''

    def __init__(self, comm, redu_prob, x0):
        '''
        Args:
            comm (MPI communicator)
            problem (ParOpt.Problem)
            x0 (indexible array object)
        '''

        super().__init__()

        # Copy over parameters
        self.comm = comm
        self.problem = redu_prob

        # Compute sizes and offsets
        local_size = len(x0)
        sizes = [ 0 ]*comm.size
        offsets = [ 0 ]*comm.size
        sizes = comm.allgather(local_size)
        if comm.size > 1:
            for i in range(1,comm.size):
                offsets[i] = offsets[i-1] + sizes[i-1]
        self.sizes = sizes
        self.offsets = offsets

        # Get number of constraints
        self.ncon = self.problem.getNumCons()

        # Compute some indices and dimensions
        self.local_size = self.sizes[self.comm.rank]
        self.global_size = np.sum(self.sizes)
        self.start = self.offsets[self.comm.rank]
        self.end = self.start + self.local_size

        # Allocate paropt vectors
        self.x = self.problem.createDesignVec()
        self.g = self.problem.createDesignVec()
        self.A = []
        for i in range(self.ncon):
            self.A.append(self.problem.createDesignVec())

        return

    def setup(self):
        self.add_input('x', shape=(self.global_size,))
        self.add_output('obj', shape=1)
        self.add_output('con', shape=1)

        self.declare_partials(of='obj', wrt='x')
        self.declare_partials(of='con', wrt='x')

        return

    def compute(self, inputs, outputs):
        # Broadcase x from root to all processor
        # In this way we only use optimization result from
        # root and implicitly discard results from any other
        # optimizer to prevent potential inconsistency
        if self.comm.rank == 0:
            x = inputs['x']
        else:
            x = None
        x =self.comm.bcast(x, root=0)

        # Each processor only use its owned part to evaluate func and grad
        self.x[:] = x[self.start:self.end]
        fail, fobj, cons = self.problem.evalObjCon(self.x)

        if fail:
            raise RuntimeError("Failed to evaluate objective and constraints!")
        else:
            outputs['obj'] = fobj
            outputs['con'] = cons[0]

        # Barrier here because we don't do block communication
        self.comm.Barrier()

        return

    def compute_partials(self, inputs, partials):
        # Broadcase x from root to all processor
        # In this way we only use optimization result from
        # root and implicitly discard results from any other
        # optimizer to prevent potential inconsistency
        if self.comm.rank == 0:
            x = inputs['x']
        else:
            x = None
        x =self.comm.bcast(x, root=0)

        # Each processor only use its owned part to evaluate func and grad
        self.x[:] = x[self.start:self.end]
        fail = self.problem.evalObjConGradient(self.x, self.g, self.A)

        if fail:
            raise RuntimeError("Failed to evaluate objective and constraints!")
        else:
            global_g = self.comm.allgather(np.array(self.g))
            global_g = np.concatenate(global_g)
            global_A = []
            for i in range(self.ncon):
                global_A.append(self.comm.allgather(np.array(self.A[i])))
                global_A[i] = np.concatenate(global_A[i])

            partials['obj', 'x'] = global_g
            partials['con', 'x'] = global_A
        return

    def globalVecToLocalvec(self, global_vec, local_vec):
        '''
        Assign corresponding part of the global vector to local vector
        '''
        local_vec[:] = global_vec[self.start: self.end]
        return

class GeneralEigSolver:
    """
    This class checks the actual smallest generalized eigenvalue with the
    Jacobi-Davidson method
    """
    def __init__(self, problem, N=10, max_jd_size=100, max_gmres_size=30):
        """
        Args:
            problem (TMR.TopoProblem): An instance that contains all components we need.
            N (int): number of eigenpairs sought
        """

        self.N = N

        # Get references from problem instance
        self.assembler = problem.getAssembler()
        self.filter = problem.getTopoFilter()
        self.mg = problem.getMg()

        # Allocate space for K and M matrix
        self.kmat = self.assembler.createMat()
        self.mmat = self.assembler.createMat()

        # Allocate space for the nodal density, eigenvalues and eigenvectors
        self.rho = self.assembler.createDesignVec()
        self.evals = np.zeros(N)
        self.evecs = [self.assembler.createVec() for _ in range(N)]

        # Create operator for the generalized eigenvalue problem
        self.oper = TACS.JDFrequencyOperator(self.assembler,
                                             self.kmat, self.mmat,
                                             self.mg.getMat(), self.mg)

        # Set up the JD solver
        self.jd = TACS.JacobiDavidson(self.oper, N, max_jd_size, max_gmres_size)
        self.jd.setTolerances(eig_rtol=1e-6, eig_atol=1e-8, rtol=1e-12, atol=1e-15)
        self.jd.setRecycle(self.N)

        # Temp vectors to check residual
        self.MatVec = self.assembler.createVec()
        self.temp = self.assembler.createVec()
        self.res = np.zeros(N)

        return

    def compute(self, x):
        """
        Take in x, update the design variable in the assembler, assemble
        K, M and mg matrices, and solve the generalized eigenvalue problem.

        Args:
            x (PVec): the (unfiltered) design variable of the full (not reduced) problem
        """

        # Update the assembler with x
        self.filter.applyFilter(TMR.convertPVecToVec(x), self.rho)
        self.assembler.setDesignVars(self.rho)

        # Update matrices and factor the multigrid preconditioner
        self.assembler.assembleMatType(TACS.STIFFNESS_MATRIX, self.kmat)
        self.assembler.assembleMatType(TACS.MASS_MATRIX, self.mmat)
        self.mg.assembleMatType(TACS.STIFFNESS_MATRIX)
        self.mg.factor()

        # Solve and check success
        self.jd.solve(print_flag=True, print_level=0)
        assert( self.jd.getNumConvergedEigenvalues() >= self.N)

        # Extract eigenpairs
        for i in range(self.N): self.evals[i], err = self.jd.extractEigenvector(i, self.evecs[i])

        # Check residuals
        for i in range(self.N):
            self.kmat.mult(self.evecs[i], self.MatVec)  # Compute K*v
            self.temp.copyValues(self.MatVec)
            self.mmat.mult(self.evecs[i], self.MatVec)  # Compute M*v
            self.temp.axpy(-self.evals[i], self.MatVec)
            self.res[i] = self.temp.norm()

        return self.evals, self.evecs, self.res