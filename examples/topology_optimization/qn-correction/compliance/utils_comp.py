from tmr import TMR, TopOptUtils
from tacs import TACS, elements, functions
from egads4py import egads
import numpy as np
import openmdao.api as om
import os
import sys
from mpi4py import MPI

sys.path.append("../eigenvalue")
from utils import OctCreator, CreatorCallback, MFilterCreator, OutputCallback
from utils import MassConstr

# Print colored text in terminal
try:
    from termcolor import colored
except:
    pass


class CompObj:
    """
    A class that evaluates the smallest eigenvalue, the objective is evaluated
    using an objective callback. We also add non-design mass to loaded nodes in
    order to form a well-posed frequency maximization problem
    """

    itr = 0

    def __init__(
        self,
        prefix,
        domain,
        forest,
        len0,
        AR,
        ratio,
        iter_offset,
        comp_scale=1.0,
        con_scale=1.0,
        nodal_force=1.0,
        gmres_subspace=50,
        rtol=1e-10,
        atol=1e-30,
        nrestart=10,
        is_flexible=0,
    ):
        """
        Args:
            eig_scale: scale the eigenvalues internally in order to acquire better
                       KS approximation with smaller skrho
            con_scale: scale the mass constraint
            num_eigenvalues: number of smallest eigenvalues to compute
            ksrho: KS parameter
        """

        # Set objects
        self.forest = forest

        # Set up parameters
        self.prefix = prefix
        self.domain = domain
        self.iter_offset = iter_offset
        self.lx = len0 * AR
        self.ly = len0
        self.lz = len0
        if domain == "lbracket":
            self.ly = len0 * ratio
        self.ratio = ratio
        self.comp_scale = comp_scale
        self.con_scale = con_scale
        self.nodal_force = nodal_force
        self.gmres_subspace = gmres_subspace
        self.rtol = rtol
        self.atol = atol
        self.nrestart = nrestart
        self.is_flexible = is_flexible

        self.fltr = None
        self.mg = None
        self.assembler = None
        self.comm = None
        self.oper = None

        self.kmat = None
        self.force = None
        self.u = None

        self.n_fail_qn_corr = 0
        self.pos_curvs = []
        self.neg_curvs = []
        self.qn_time = []

        # Save snapshots throughout optimization iterations
        self.num_obj_evals = 0
        self.save_snapshot_every = 1
        self.snapshot = {
            "iter": [],
            "obj": [],
            "discreteness": [],
            "discreteness_rho": [],
        }
        self.snapshot_x = None
        self.snapshot_rho = None

        return

    def objective(self, fltr, mg):
        """
        Evaluate the compliance
        """

        if self.fltr is None:
            self.mg = mg
            self.fltr = fltr
            self.assembler = self.fltr.getAssembler()
            self.svec = self.assembler.createDesignVec()
            self.comm = self.assembler.getMPIComm()
            self.rank = self.comm.rank

            self.snapshot_x = self.assembler.createDesignVec()
            self.snapshot_rho = self.assembler.createDesignVec()

            # Allocate vectors for qn correction
            self.rho = self.assembler.createDesignVec()
            self.rho_original = self.assembler.createDesignVec()
            self.update = self.assembler.createDesignVec()

            # Initialize space for matrices and vectors
            self.kmat = self.assembler.createMat()
            self.force = self.assembler.createVec()
            self.u = self.assembler.createVec()

            # Set the solver
            self.gmres = TACS.KSM(
                self.mg.getMat(),
                self.mg,
                self.gmres_subspace,
                self.nrestart,
                self.is_flexible,
            )
            self.gmres.setMonitor(self.comm)
            self.gmres.setTolerances(self.rtol, self.atol)

            """
            Create force vector
            """
            self.force_vals = self.force.getArray()

            # Get nodal locations
            Xpts = self.forest.getPoints()

            # Note: the local nodes are organized as follows:
            # |--- dependent nodes -- | ext_pre | -- owned local -- | - ext_post -|

            # Get number of local nodes in the current processor
            n_local_nodes = Xpts.shape[0]

            # Get number of ext_pre nodes
            n_ext_pre = self.forest.getExtPreOffset()

            # Get numbder of own nodes:
            offset = n_ext_pre

            # Create force vector
            # if self.domain == '3dcantilever':
            # self.force = TopOptUtils.computeVertexLoad('pt1',
            #     self.forest, self.assembler, [0, 0, -self.nodal_force])
            # temp = TopOptUtils.computeVertexLoad('pt2',
            #     self.forest, self.assembler, [0, self.nodal_force, 0])
            # self.force.axpy(1.0, temp)
            # self.n_loaded_nodes = 2

            tol = 1e-6
            if self.domain == "cantilever":
                xmin = self.lx - tol
                xmax = self.lx + tol
                ymin = 0.25 * self.ly - tol
                ymax = 0.75 * self.ly + tol
                zmin = 0.0 * self.lz - tol
                zmax = 0.2 * self.lz + tol

            elif self.domain == "3dcantilever":
                xmin = self.lx - tol
                xmax = self.lx + tol
                ymin = 0.0 * self.ly - tol
                ymax = 0.2 * self.ly + tol
                zmin = 0.0 * self.lz - tol
                zmax = 0.2 * self.lz + tol
                y2min = 0.8 * self.ly - tol
                y2max = 1.0 * self.ly + tol
                z2min = 0.8 * self.lz - tol
                z2max = 1.0 * self.lz + tol

            elif self.domain == "michell":
                xmin = self.lx - tol
                xmax = self.lx + tol
                ymin = 0.25 * self.ly - tol
                ymax = 0.75 * self.ly + tol
                zmin = 0.4 * self.lz - tol
                zmax = 0.6 * self.lz + tol

            elif self.domain == "mbb":
                xmin = 0.0 * self.lx - tol
                xmax = 0.2 * self.lx + tol
                ymin = 0.25 * self.ly - tol
                ymax = 0.75 * self.ly + tol
                zmin = self.lz - tol
                zmax = self.lz + tol

            elif self.domain == "lbracket":
                RATIO = self.ratio
                xmin = self.lx - tol
                xmax = self.lx + tol
                ymin = 0.25 * self.ly - tol
                ymax = 0.75 * self.ly + tol
                zmin = 0.5 * RATIO * self.lz - tol
                zmax = 1.0 * RATIO * self.lz + tol

            else:
                raise ValueError("Unsupported domain type!")

            self.n_loaded_nodes = 0
            for i in range(offset, n_local_nodes):
                x, y, z = Xpts[i]
                if xmin < x < xmax:
                    if ymin < y < ymax:
                        if zmin < z < zmax:
                            self.force_vals[3 * (i - offset) + 2] = -self.nodal_force
                            self.n_loaded_nodes += 1

            if self.domain == "3dcantilever":
                for i in range(offset, n_local_nodes):
                    x, y, z = Xpts[i]
                    if xmin < x < xmax:
                        if y2min < y < y2max:
                            if z2min < z < z2max:
                                self.force_vals[3 * (i - offset) + 1] = self.nodal_force
                                self.n_loaded_nodes += 1

            self.n_loaded_nodes = self.comm.allreduce(self.n_loaded_nodes)

            # Match the reordering of the vector
            self.assembler.reorderVec(self.force)

        """
        Export force vectors to f5 file for verification
        """
        # itr = CompObj.itr + self.iter_offset
        # if itr % 10 == 0:
        #     # Set up flags for data output
        #     flag = (TACS.OUTPUT_CONNECTIVITY |
        #             TACS.OUTPUT_NODES |
        #             TACS.OUTPUT_EXTRAS |
        #             TACS.OUTPUT_DISPLACEMENTS)

        #     # Create f5 file writer
        #     f5 = TACS.ToFH5(self.assembler, TACS.SOLID_ELEMENT, flag)

        #     self.assembler.setVariables(self.force)

        #     f5.writeToFile(os.path.join(self.prefix, 'rhs{:d}.f5'.format(
        #         self.iter_offset+CompObj.itr)))

        # CompObj.itr += 1
        """
        End export non-design mass vectors to f5 file for verification
        """

        # Assemble the multigrid preconditioner
        self.mg.assembleMatType(TACS.STIFFNESS_MATRIX)
        self.mg.factor()

        # Solve
        self.gmres.solve(self.force, self.u)

        # Objective
        obj = self.comp_scale / self.n_loaded_nodes * self.force.dot(self.u)

        # Print values
        if self.comm.rank == 0:
            print("{:30s}{:20.10e}".format("[Obj] scaled compliance:", obj))

        # Save a snapshot of the result
        if self.num_obj_evals % self.save_snapshot_every == 0:
            self.snapshot["iter"].append(self.num_obj_evals)
            self.snapshot["obj"].append(obj)

            # Get x and rho
            fltr.getDesignVars(self.snapshot_x)
            self.assembler.getDesignVars(self.snapshot_rho)

            # Compute discreteness
            snapshot_x_global = self.comm.allgather(self.snapshot_x.getArray())
            snapshot_x_global = np.concatenate(snapshot_x_global)
            discreteness = np.dot(snapshot_x_global, 1.0 - snapshot_x_global) / len(
                snapshot_x_global
            )

            # Compute discreteness for rho
            snapshot_rho_global = self.comm.allgather(self.snapshot_rho.getArray())
            snapshot_rho_global = np.concatenate(snapshot_rho_global)
            discreteness_rho = np.dot(
                snapshot_rho_global, 1.0 - snapshot_rho_global
            ) / len(snapshot_rho_global)

            self.snapshot["discreteness"].append(discreteness)
            self.snapshot["discreteness_rho"].append(discreteness_rho)

        self.num_obj_evals += 1

        return obj

    def objective_gradient(self, fltr, mg, dfdrho):
        """
        Evaluate compliance gradient

        dcdx = -u^T dK/dx u
        """

        # Zero out the gradient vector
        dfdrho.zeroEntries()

        # scale the objective
        scale = -self.comp_scale / self.n_loaded_nodes

        # Compute gradient
        self.assembler.addMatDVSensInnerProduct(
            scale, TACS.STIFFNESS_MATRIX, self.u, self.u, dfdrho
        )

        # Make sure the vector is properly distributed over all processors
        dfdrho.beginSetValues(op=TACS.ADD_VALUES)
        dfdrho.endSetValues(op=TACS.ADD_VALUES)

        # Compute gradient norm
        norm = dfdrho.norm()
        if self.comm.rank == 0:
            print("{:30s}{:20.10e}".format("[Obj] gradient norm:", norm))
        return

    def qn_correction(self, x, z, zw, s, y):
        """
        Update y:
        y <- y + F^T P Fs

        where:
        F: filter matrix
        P: Positive definite part of the objective Hessian

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
            z, zw: dummy variable for qn correction for constraints, not used here.
        """
        # Timer
        t_start = MPI.Wtime()

        # Finite difference step length for computing second order
        # derivative of stiffness matrix
        h = 1e-8

        # Get current nodal density
        self.assembler.getDesignVars(self.rho)
        self.rho_original.copyValues(self.rho)

        # Apply filter to step vector
        self.svec.zeroEntries()
        self.fltr.applyFilter(TMR.convertPVecToVec(s), self.svec)
        self.rho.axpy(h, self.svec)

        """
        Compute the following second derivative using finite difference
        P += phi^T d2Kdx2 phi
        """

        # Compute coefficient
        coeff = self.comp_scale / self.n_loaded_nodes

        # Zero out update vector
        self.update.zeroEntries()

        # Compute g(rho + h*s)
        self.assembler.setDesignVars(self.rho)
        self.assembler.addMatDVSensInnerProduct(
            coeff, TACS.STIFFNESS_MATRIX, self.u, self.u, self.update
        )

        # Compute dg = g(rho + h*s) -g(rho)
        self.assembler.setDesignVars(self.rho_original)
        self.assembler.addMatDVSensInnerProduct(
            -coeff, TACS.STIFFNESS_MATRIX, self.u, self.u, self.update
        )

        # Distribute the vector
        self.update.beginSetValues(op=TACS.ADD_VALUES)
        self.update.endSetValues(op=TACS.ADD_VALUES)

        # Compute dg/h
        self.update.scale(1 / h)

        # Compute curvature and check the norm of the update
        # to see if the magnitude makes sense
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
        if curvature > 0:
            self.pos_curvs.append(curvature)
            y_wrap = TMR.convertPVecToVec(y)
            self.fltr.applyTranspose(self.update, self.update)
            y_wrap.axpy(1.0, self.update)

        else:
            self.n_fail_qn_corr += 1
            self.neg_curvs.append(curvature)

        # Timer
        self.qn_time.append(MPI.Wtime() - t_start)
        return

    def getFailQnCorr(self):
        return self.n_fail_qn_corr, self.neg_curvs, self.pos_curvs

    def getAveragedQnTime(self):
        return np.average(self.qn_time)

    def get_snapshot(self):
        return self.snapshot


def create_problem(
    prefix,
    domain,
    forest,
    bcs,
    props,
    nlevels,
    vol_frac=0.25,
    r0_frac=0.05,
    len0=1.0,
    AR=1.0,
    ratio=0.4,
    density=2600.0,
    iter_offset=0,
    qn_correction=True,
    eq_constr=False,
    comp_scale=1.0,
):
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
    problem = TopOptUtils.createTopoProblem(
        forest, obj.creator_callback, filter_type, use_galerkin=True, nlevels=nlevels
    )

    # Get the assembler object we just created
    assembler = problem.getAssembler()

    # Compute the fixed mass target
    lx = len0 * AR  # mm
    ly = len0  # mm
    lz = len0  # mm
    if domain == "lbracket":
        ly = len0 * ratio
    vol = lx * ly * lz
    if domain == "lbracket":
        S1 = lx * lz
        S2 = lx * lz * (1 - ratio) ** 2
        vol = (S1 - S2) * ly
    m_fixed = vol_frac * (vol * density)

    # Add objective callback
    obj_callback = CompObj(
        prefix, domain, forest, len0, AR, ratio, iter_offset, comp_scale=comp_scale
    )
    problem.addObjectiveCallback(
        obj_callback.objective, obj_callback.objective_gradient
    )

    # Add constraint callback
    constr_callback = MassConstr(m_fixed, assembler.getMPIComm())
    nineq = 1
    if eq_constr is True:
        nineq = 0
    problem.addConstraintCallback(
        1, nineq, constr_callback.constraint, constr_callback.constraint_gradient
    )

    # Use Quasi-Newton Update Correction if specified
    if qn_correction:
        problem.addQnCorrectionCallback(1, obj_callback.qn_correction)

    # Set output callback
    cb = OutputCallback(assembler, iter_offset=iter_offset)
    problem.setOutputCallback(cb.write_output)

    return problem, obj_callback, constr_callback
