from tmr import TMR, TopOptUtils
from tacs import TACS, elements, functions
from egads4py import egads
import numpy as np
import openmdao.api as om
import os
import sys

sys.path.append('../eigenvalue')
from utils import OctCreator, CreatorCallback, MFilterCreator, OutputCallback
from utils import MassConstr, FrequencyObj
from utils import OmAnalysis

sys.path.append('../frequency')
from utils_freq import FrequencyConstr

sys.path.append('../compliance')
from utils_comp import CompObj

# Print colored text in terminal
try:
    from termcolor import colored
except:
    print("[Unavailable module] termcolor is not installed!")

def create_problem(prefix, domain, forest, bcs, props, nlevels, lambda0, ksrho,
                   has_freq_constr=True,
                   vol_frac=0.25, r0_frac=0.05, len0=1.0, AR=1.0, ratio=0.4,
                   density=2600.0, iter_offset=0,
                   qn_correction_comp=True, qn_correction_freq=True,
                   comp_scale=1.0, eig_scale=1.0, eq_constr=False,
                   max_jd_size=100, max_gmres_size=30):
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
        vol_frac (float): Volume fraction for the mass constraint
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
        S2 = lx*lz*(1-ratio)**2
        vol = (S1-S2)*ly
    m_fixed = vol_frac*(vol*density)

    # Add objective callback
    obj_callback = CompObj(prefix, domain, forest, len0, AR, ratio,
                           iter_offset, comp_scale=comp_scale)
    problem.addObjectiveCallback(obj_callback.objective,
                                 obj_callback.objective_gradient)

    # Instantiate constraint callback classes
    massconstr = MassConstr(m_fixed, assembler.getMPIComm())
    freqconstr = FrequencyConstr(prefix, domain, forest, len0, AR, ratio,
                                 iter_offset, lambda0,
                                 ksrho=ksrho,
                                 non_design_mass=0.0,
                                 eig_scale=eig_scale,
                                 max_jd_size=max_jd_size,
                                 max_gmres_size=max_gmres_size)

    # Add constraint callback
    if has_freq_constr:
        nconstr = 2
        constr_callback_fun = lambda fltr, mg :[
            massconstr.constraint(fltr, mg)[0],
            freqconstr.constraint(fltr, mg)[0]]
        def constr_callback_grad(fltr, mg, vecs):
            massconstr.constraint_gradient(fltr, mg, [vecs[0]])
            freqconstr.constraint_gradient(fltr, mg, [vecs[1]])
            return
    else:
        nconstr = 1
        constr_callback_fun = massconstr.constraint
        constr_callback_grad = massconstr.constraint_gradient

    nineq = nconstr
    if eq_constr is True:
        nineq = 0
    problem.addConstraintCallback(nconstr, nineq, constr_callback_fun,
                                  constr_callback_grad)

    # Use Quasi-Newton Update Correction if specified

    # Compliance only
    if qn_correction_comp and not qn_correction_freq:
        problem.addQnCorrectionCallback(nconstr, obj_callback.qn_correction)

    # Frequency only
    if qn_correction_freq and not qn_correction_comp:
        if has_freq_constr:
            problem.addQnCorrectionCallback(nconstr, freqconstr.qn_correction)
        else:
            print("[Warning]doesn't have frequency constraint, \
                discard frequency correction!")

    # Compliance and frequency
    if qn_correction_freq and qn_correction_comp:
        if not has_freq_constr:
            print("[Warning]doesn't have frequency constraint, \
                discard frequency correction!")
        def qn_callback_fun(x, z, zw, s, y):
            obj_callback.qn_correction(x, z, zw, s, y)
            freqconstr.qn_correction(x, z, zw, s, y)
            return
        problem.addQnCorrectionCallback(nconstr, qn_callback_fun)

    # Set output callback
    cb = OutputCallback(assembler, iter_offset=iter_offset)
    problem.setOutputCallback(cb.write_output)

    return problem, obj_callback

class GEP_solver:
    """
    This utility class solves the generalized eigenvalue problem and
    computes the smallest eigenvalue for optimal designs
    """

    def __init__(self, fltr, mg, max_jd_size, max_gmres_size, num_eigenvalues=2):

        self.fltr = fltr
        self.mg = mg
        self.max_jd_size = max_jd_size
        self.max_gmres_size = max_gmres_size
        self.num_eigenvalues = num_eigenvalues

        self.assembler = fltr.getAssembler()
        self.kmat = self.assembler.createMat()
        self.mmat = self.assembler.createMat()
        self.comm = self.assembler.getMPIComm()
        self.eig = np.zeros(self.num_eigenvalues)
        self.eigv = []
        for i in range(self.num_eigenvalues):
            self.eigv.append(self.assembler.createVec())

        # Create the Jacobi-Davidson operator
        self.oper = TACS.JDFrequencyOperator(self.assembler, self.kmat,
            self.mmat, self.mg.getMat(), self.mg)

        # Create the eigenvalue solver and set the number of recycling eigenvectors
        self.jd = TACS.JacobiDavidson(self.oper, self.num_eigenvalues,
            self.max_jd_size, self.max_gmres_size)
        self.jd.setTolerances(eig_rtol=1e-6, eig_atol=1e-8, rtol=1e-6, atol=1e-12)
        self.jd.setRecycle(self.num_eigenvalues)

        return

    def solve(self, dv):
        """
        Solve
        """

        # Assemble matrices
        self.assembler.getDesignVars(dv)
        self.assembler.assembleMatType(TACS.STIFFNESS_MATRIX, self.kmat)
        self.assembler.assembleMatType(TACS.MASS_MATRIX, self.mmat)

        # Assemble the multigrid preconditioner
        self.mg.assembleMatType(TACS.STIFFNESS_MATRIX)
        self.mg.factor()

        # Solve
        self.jd.solve(print_flag=True, print_level=0)

        # Check success
        if self.jd.getNumConvergedEigenvalues() < self.num_eigenvalues:
            if self.comm.rank == 0:
                print("[Warning] Jacobi-Davidson failed to converge for the first run.")

            # Extract the eigenvalues
            for i in range(self.num_eigenvalues):
                self.eig[i], error = self.jd.extractEigenvalue(i)

            # Update preconditioner
            theta = 0.9*np.min(self.eig)
            self.mg.assembleMatCombo(TACS.STIFFNESS_MATRIX, 1.0, TACS.MASS_MATRIX, -theta)
            self.mg.factor()

            # Rerun the solver
            self.jd.solve(print_flag=True, print_level=1)
            nconvd = self.jd.getNumConvergedEigenvalues()

            # If it still fails, raise error and exit
            if nconvd < self.num_eigenvalues:
                msg = "No enough eigenvalues converged! ({:d}/{:d})".format(
                    nconvd, self.num_eigenvalues)
                raise ValueError(msg)

        # Extract eigenvalues and eigenvectors
        for i in range(self.num_eigenvalues):
            self.eig[i], error = self.jd.extractEigenvector(i, self.eigv[i])

        # Get smallest eigenvalue
        min_eigv = np.min(self.eig)

        return min_eigv

class CompFreqOmAnalysis(OmAnalysis):

    def __init__(self, ncon, comm, problem, obj_callback, sizes, offsets):

        OmAnalysis.__init__(self, comm, problem, obj_callback, sizes, offsets)

        if ncon != 1 and ncon != 2:
            raise ValueError("Only support ncon = 1 or 2")
        else:
            self.ncon = ncon

        if ncon == 2:
            self.A2_PVec = self.problem.createDesignVec()
            self.A2_Vec = TMR.convertPVecToVec(self.A2_PVec)
            self.A2_vals = self.A2_Vec.getArray()

        return

    def setup(self):
        self.add_input('x', shape=(self.global_size,))
        self.add_output('obj', shape=1)
        self.add_output('con', shape=self.ncon)

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

        self.x_vals[:] = x[self.start:self.end]
        fail, fobj, cons = self.problem.evalObjCon(self.ncon, self.x_PVec)

        if fail:
            raise RuntimeError("Failed to evaluate objective and constraints!")
        else:
            outputs['obj'] = fobj
            outputs['con'] = cons

        # Barrier here because we don't do block communication
        self.comm.Barrier()

        return

    def compute_partials(self, inputs, partials):
        if self.comm.rank == 0:
            x = inputs['x']
        else:
            x = None
        x =self.comm.bcast(x, root=0)
        self.x_vals[:] = x[self.start:self.end]

        if self.ncon == 1:
            A_list = [self.A_PVec]
        else:
            A_list = [self.A_PVec, self.A2_PVec]

        fail = self.problem.evalObjConGradient(self.x_PVec, self.g_PVec, A_list)

        if fail:
            raise RuntimeError("Failed to evaluate objective and constraints!")
        else:
            global_g = self.comm.allgather(self.g_vals)
            global_g = np.concatenate(global_g)
            partials['obj', 'x'] = global_g

            global_A = self.comm.allgather(self.A_vals)
            global_A = np.concatenate(global_A)
            if self.ncon == 1:
                partials['con', 'x'] = [global_A]
            else:
                global_A2 = self.comm.allgather(self.A2_vals)
                global_A2 = np.concatenate(global_A2)
                partials['con', 'x'] = [global_A, global_A2]

        return