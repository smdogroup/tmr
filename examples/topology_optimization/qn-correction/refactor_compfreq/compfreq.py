import argparse
import sys
import os
from mpi4py import MPI
from tmr import TMR, TopOptUtils
from tacs import TACS, constitutive
from mma4py import Problem as MMAProblemBase
from mma4py import Optimizer as MMAOptimizer
import numpy as np
import openmdao.api as om
from paropt import ParOpt
import pickle

sys.path.extend(["../eigenvalue", "../refactor_frequency", "../compliance"])
from refactor_utils_freq import (
    FrequencyConstr,
    ReducedProblem,
    ReduOmAnalysis,
    GeneralEigSolver,
)
from utils_comp import CompObj
from utils import (
    create_forest,
    MFilterCreator,
    CreatorCallback,
    MassConstr,
    OutputCallback,
    getNSkipUpdate,
)


def parse_args():
    # Create the argument parser
    p = argparse.ArgumentParser()

    # os
    p.add_argument("--prefix", type=str, default="./results")

    # Analysis
    p.add_argument(
        "--domain",
        type=str,
        default="cantilever",
        choices=["cantilever", "michell", "mbb", "lbracket"],
    )
    p.add_argument("--AR", type=float, default=1.0)
    p.add_argument("--ratio", type=float, default=0.4)
    p.add_argument("--len0", type=float, default=1.0)
    p.add_argument("--vol-frac", type=float, default=0.4)
    p.add_argument("--r0-frac", type=float, default=0.05)
    p.add_argument("--htarget", type=float, default=1.0)
    p.add_argument("--mg-levels", type=int, default=4)
    p.add_argument("--qval", type=float, default=5.0)
    p.add_argument("--lambda0", type=float, default=0.1)
    p.add_argument("--ksrho", type=float, default=1000)
    p.add_argument("--max-jd-size", type=int, default=200)
    p.add_argument("--max-gmres-size", type=int, default=30)
    p.add_argument("--num-eigenvalues", type=int, default=10)

    # Optimization
    p.add_argument(
        "--optimizer",
        type=str,
        default="paropt",
        choices=["paropt", "snopt", "ipopt", "mma", "mma4py"],
    )
    p.add_argument("--hessian", default="bfgs", choices=["bfgs", "sr1"])
    p.add_argument("--n-mesh-refine", type=int, default=1)
    p.add_argument("--niter-finest", type=int, default=100)
    p.add_argument("--max-iter", type=int, default=100)
    p.add_argument("--qn-correction", action="store_true")
    p.add_argument("--comp-scale", type=float, default=1.0)
    p.add_argument("--eig-scale", type=float, default=1.0)
    p.add_argument("--output-level", type=int, default=0)
    p.add_argument("--simple-filter", action="store_false")
    p.add_argument("--tr-eta", type=float, default=0.25)
    p.add_argument("--tr-min", type=float, default=1e-3)
    p.add_argument("--qn-subspace", type=int, default=10)
    p.add_argument(
        "--paropt-type",
        type=str,
        default="penalty_method",
        choices=["penalty_method", "filter_method"],
    )

    # Tests
    p.add_argument("--test-gradient-check", action="store_true")

    # Parse arguments
    args = p.parse_args()
    return args


def create_paropt_options(args):
    # Set up ParOpt parameters
    if args.paropt_type == "penalty_method":
        adaptive_gamma_update = True
    else:
        adaptive_gamma_update = False
    options = {
        "algorithm": "tr",
        "output_level": args.output_level,
        "norm_type": "l1",
        "tr_init_size": 0.05,
        "tr_min_size": args.tr_min,
        "tr_max_size": 1.0,
        "tr_eta": args.tr_eta,
        "tr_infeas_tol": 1e-6,
        "tr_l1_tol": 0.0,
        "tr_linfty_tol": 0.0,
        "tr_adaptive_gamma_update": adaptive_gamma_update,
        "tr_accept_step_strategy": args.paropt_type,
        "filter_sufficient_reduction": args.simple_filter,
        "filter_has_feas_restore_phase": True,
        "tr_use_soc": False,
        "tr_max_iterations": args.max_iter,
        "penalty_gamma": 50.0,
        "qn_subspace_size": args.qn_subspace,  # try 5 or 10
        "qn_type": args.hessian,
        "qn_diag_type": "yty_over_yts",
        "abs_res_tol": 1e-8,
        "starting_point_strategy": "affine_step",
        "barrier_strategy": "mehrotra_predictor_corrector",
        "tr_steering_barrier_strategy": "mehrotra_predictor_corrector",
        "tr_steering_starting_point_strategy": "affine_step",
        "use_line_search": False,  # subproblem
        "max_major_iters": 200,
    }
    return options


def create_paropt_mma_options(args):
    mma_options = {
        "algorithm": "mma",
        "mma_asymptote_contract": 0.7,
        "mma_asymptote_relax": 1.2,
        "mma_bound_relax": 0,
        "mma_delta_regularization": 1e-05,
        "mma_eps_regularization": 0.001,
        "mma_infeas_tol": 1e-05,
        "mma_init_asymptote_offset": 0.25,
        "mma_l1_tol": 1e-06,
        "mma_linfty_tol": 1e-06,
        "mma_max_asymptote_offset": 10,
        "mma_max_iterations": args.max_iter,
        "mma_min_asymptote_offset": 0.01,
        "mma_use_constraint_linearization": True,
    }
    return mma_options


class Constraints:
    def __init__(self, freq_constr: FrequencyConstr, mass_constr: MassConstr):
        self.freq_constr = freq_constr
        self.mass_constr = mass_constr
        return

    def constraint(self, fltr, mg):
        c_mass = self.mass_constr.constraint(fltr, mg)[0]
        c_freq = self.freq_constr.constraint(fltr, mg)[0]
        return [c_mass, c_freq]

    def constraint_gradient(self, fltr, mg, vecs):
        self.mass_constr.constraint_gradient(fltr, mg, vecs, index=0)
        self.freq_constr.constraint_gradient(fltr, mg, vecs, index=1)
        return


class QnCorr:
    def __init__(self, comp_obj: CompObj, freq_constr: FrequencyConstr):
        self.comp_obj = comp_obj
        self.freq_constr = freq_constr
        return

    def qn_correction_func(self, zero_idx, z, s, y):
        self.comp_obj.qn_correction(None, None, None, s, y, zero_idx)
        self.freq_constr.qn_correction(zero_idx, z, s, y)
        return


def create_problem(
    args, step, bcs, stiffness_props, forest, lx, ly, lz, density, iter_offset
):
    # Create the problem and filter object
    mfilter = MFilterCreator(args.r0_frac, N=20, a=args.len0)
    filter_type = mfilter.filter_callback
    obj = CreatorCallback(bcs, stiffness_props)
    problem = TopOptUtils.createTopoProblem(
        forest,
        obj.creator_callback,
        filter_type,
        use_galerkin=True,
        nlevels=args.mg_levels + step,
    )
    assembler = problem.getAssembler()

    # Compute mass target
    vol = lx * ly * lz
    if args.domain == "lbracket":
        S1 = lx * lz
        S2 = lx * lz * (1.0 - args.ratio) ** 2
        vol = (S1 - S2) * ly
    m_fixed = args.vol_frac * (vol * density)

    # Add objective callback
    comp_obj = CompObj(
        args.prefix,
        args.domain,
        forest,
        args.len0,
        args.AR,
        args.ratio,
        iter_offset,
        comp_scale=args.comp_scale,
    )
    problem.addObjectiveCallback(comp_obj.objective, comp_obj.objective_gradient)

    # Add constraints callback
    freq_constr = FrequencyConstr(
        args.prefix,
        args.domain,
        forest,
        args.len0,
        args.AR,
        args.ratio,
        iter_offset,
        args.lambda0,
        ksrho=args.ksrho,
        eig_scale=args.eig_scale,
        num_eigenvalues=args.num_eigenvalues,
        max_jd_size=args.max_jd_size,
        max_gmres_size=args.max_gmres_size,
        add_non_design_mass=False,
    )
    mass_constr = MassConstr(m_fixed, assembler.getMPIComm())
    constr_callback = Constraints(freq_constr, mass_constr)
    problem.addConstraintCallback(
        2, 2, constr_callback.constraint, constr_callback.constraint_gradient
    )

    # Set output callback
    cb = OutputCallback(assembler, iter_offset=iter_offset)
    problem.setOutputCallback(cb.write_output)

    # Set up the qn correction function

    return problem, comp_obj, freq_constr


class MMAProblem(MMAProblemBase):
    def __init__(self, comm, nvars, nvars_l, prob: ReducedProblem):
        self.ncon = 2
        super().__init__(comm, nvars, nvars_l, self.ncon)
        self.prob = prob

        self.xvec = prob.createDesignVec()
        self.gvec = prob.createDesignVec()
        self.gcvec = prob.createDesignVec()
        return

    def getVarsAndBounds(self, x, lb, ub):
        self.prob.getVarsAndBounds(x, lb, ub)
        return

    def evalObjCon(self, x, cons) -> float:
        fail, obj, con = self.prob.evalObjCon(x)
        cons[:] = -con[:]
        return obj

    def evalObjConGrad(self, x, g, gcon):
        self.prob.evalObjConGradient(x, g, gcon)
        gcon[:] = -gcon[:]
        return


if __name__ == "__main__":
    # Set arguments and create result directory
    args = parse_args()
    comm = MPI.COMM_WORLD

    # Create prefix directory if not exist
    if comm.rank == 0 and not os.path.isdir(args.prefix):
        os.mkdir(args.prefix)

    # Save the command and arguments that executed this script
    cmd = None
    if comm.rank == 0:
        cmd = "python " + " ".join(sys.argv)
        with open(os.path.join(args.prefix, "exe.sh"), "w") as f:
            f.write(cmd + "\n")

    # Allow root processor some time to save execution command
    comm.Barrier()

    # Compute geometry parameters
    lx = args.len0 * args.AR
    ly = args.len0
    lz = args.len0
    if args.domain == "lbracket":
        ly = args.len0 * args.ratio

    # Set up material properties, stiffness properties and initial forest
    material_props = constitutive.MaterialProperties(
        rho=2600.0, E=70e3, nu=0.3, ys=100.0
    )
    stiffness_props = TMR.StiffnessProperties(
        material_props, k0=1e-3, eps=0.2, q=args.qval, qmass=args.qval
    )
    forest = create_forest(
        comm, lx, ly, lz, args.ratio, args.htarget, args.mg_levels - 1, args.domain
    )

    # Set boundary conditions
    bcs = TMR.BoundaryConditions()
    if args.domain == "mbb":
        bcs.addBoundaryCondition("symmetry", [0], [0.0])
        bcs.addBoundaryCondition("support", [1, 2], [0.0, 0.0])
    else:
        bcs.addBoundaryCondition("fixed", [0, 1, 2], [0.0, 0.0, 0.0])

    # Set the original filter to NULL
    old_filter = None

    # This is the mesh refinement loop: each step corresponds to a refined mesh
    # except for the first step where the original (coarse) mesh is used
    for step in range(args.n_mesh_refine):
        iter_offset = step * args.max_iter

        # Create the optimization problem instance
        problem, comp_obj, freq_constr = create_problem(
            args=args,
            step=step,
            bcs=bcs,
            stiffness_props=stiffness_props,
            forest=forest,
            lx=lx,
            ly=ly,
            lz=lz,
            density=2600.0,
            iter_offset=iter_offset,
        )

        # Function handle for qn correction, if specified
        if args.qn_correction:
            qn_corr_func_cls = QnCorr(comp_obj, freq_constr)
            qn_corr_func = qn_corr_func_cls.qn_correction_func
        else:
            qn_corr_func = None

        # Set the prefix
        problem.setPrefix(args.prefix)

        # Initialize the problem, set conuter and the prefix
        problem.initialize()
        problem.setIterationCounter(iter_offset)

        # Extract the filter to interpolate design variables
        new_filter = problem.getFilter()

        # Create a reduced problem by fixing design variables where a constant
        # mass need to be applied to
        fixed_dv_idx = []
        redu_prob = ReducedProblem(problem, [], qn_correction_func=qn_corr_func, ncon=2)

        # Check gradient and exit
        if args.test_gradient_check:
            redu_prob.checkGradients(1e-6)
            exit(0)

        # Create a reduced x0
        redu_x0 = redu_prob.createDesignVec()

        # Helper variable: a full-sized design vector
        _x = problem.createDesignVec()

        if step != 0:  # This is a refined step - interpolation needed
            # Interpolate x from xopt from last refinement step
            TopOptUtils.interpolateDesignVec(old_filter, _xopt, new_filter, _x)

            # copy values from _x to redu_x0
            redu_prob.DVtoreduDV(_x, redu_x0)
            redu_prob.setInitDesignVars(redu_x0)

        else:  # This is the first step, we just set x0 to 0.95
            redu_x0[:] = 0.95

        # Update filter
        old_filter = new_filter

        # Set optimizer options
        paropt_options = create_paropt_options(args)
        mma_options = create_paropt_mma_options(args)
        paropt_options["output_file"] = os.path.join(
            args.prefix, "output_file%d.dat" % (step)
        )
        paropt_options["tr_output_file"] = os.path.join(
            args.prefix, "tr_output_file%d.dat" % (step)
        )
        mma_options["mma_output_file"] = os.path.join(
            args.prefix, "mma_output_file%d.dat" % (step)
        )
        if args.n_mesh_refine > 1:
            if step == args.n_mesh_refine - 1:
                paropt_options["tr_max_iterations"] = args.niter_finest
                mma_options["mma_max_iterations"] = args.niter_finest

        # Allocate space to store reduced/full optimal x/rho
        redu_xopt = redu_prob.createDesignVec()
        redu_rhoopt = redu_prob.createDesignVec()
        xopt = problem.getAssembler().createDesignVec()
        rhoopt = problem.getAssembler().createDesignVec()

        if args.optimizer == "mma4py":
            nvars_l = len(redu_x0)
            nvars = np.zeros(1, dtype=type(nvars_l))
            comm.Allreduce(np.array([nvars_l]), nvars)
            mmaprob = MMAProblem(comm, nvars[0], nvars_l, redu_prob)
            out_file = os.path.join(args.prefix, "mma4py_output_file%d.dat" % (step))
            mmaopt = MMAOptimizer(mmaprob, out_file)
            mmaopt.optimize(args.max_iter)

            # Get optimal objective and constraint
            redu_xopt_vals = mmaopt.getOptimizedDesign()
            redu_xopt = redu_prob.createDesignVec()
            redu_xopt[:] = redu_xopt_vals[:]

            fail, obj, cons = redu_prob.evalObjCon(redu_xopt)
            con = cons[0]

            # Compute discreteness
            redu_xopt_g = comm.allgather(np.array(redu_xopt))
            redu_xopt_g = np.concatenate(redu_xopt_g)
            discreteness = np.dot(redu_xopt_g, 1.0 - redu_xopt_g) / len(redu_xopt_g)

            # Compute discreteness for rho
            redu_prob.reduDVtoDV(redu_xopt, xopt.getArray())
            problem.getTopoFilter().applyFilter(xopt, rhoopt)
            redu_prob.DVtoreduDV(rhoopt.getArray(), redu_rhoopt)
            redu_rhoopt_g = comm.allgather(np.array(redu_rhoopt))
            redu_rhoopt_g = np.concatenate(redu_rhoopt_g)
            discreteness_rho = np.dot(redu_rhoopt_g, 1.0 - redu_rhoopt_g) / len(
                redu_rhoopt_g
            )

        elif args.optimizer == "snopt" or args.optimizer == "ipopt":
            # Create distributed openMDAO component
            omprob = om.Problem()
            analysis = ReduOmAnalysis(comm, redu_prob, redu_x0)
            indeps = omprob.model.add_subsystem("indeps", om.IndepVarComp())

            # Create global design vector
            redu_x0_g = comm.allgather(np.array(redu_x0))
            redu_x0_g = np.concatenate(redu_x0_g)
            indeps.add_output("x", redu_x0_g)
            omprob.model.add_subsystem("topo", analysis)
            omprob.model.connect("indeps.x", "topo.x")
            omprob.model.add_design_var("indeps.x", lower=0.0, upper=1.0)
            omprob.model.add_objective("topo.obj")
            omprob.model.add_constraint("topo.con", lower=0.0)

            # Set up optimizer and options
            if args.optimizer == "snopt":
                omprob.driver = om.pyOptSparseDriver()
                omprob.driver.options["optimizer"] = "SNOPT"
                omprob.driver.opt_settings["Iterations limit"] = 9999999999999
                omprob.driver.opt_settings["Major feasibility tolerance"] = 1e-10
                omprob.driver.opt_settings["Major optimality tolerance"] = 1e-10
                omprob.driver.opt_settings["Summary file"] = os.path.join(
                    args.prefix, "snopt_output_file%d.dat" % (step)
                )
                omprob.driver.opt_settings["Print file"] = os.path.join(
                    args.prefix, "print_output_file%d.dat" % (step)
                )
                omprob.driver.opt_settings["Major print level"] = 1
                omprob.driver.opt_settings["Minor print level"] = 0

                if args.n_mesh_refine > 1 and step == args.n_mesh_refine - 1:
                    omprob.driver.opt_settings["Major iterations limit"] = (
                        args.niter_finest
                    )
                else:
                    omprob.driver.opt_settings["Major iterations limit"] = args.max_iter

            elif args.optimizer == "ipopt":
                omprob.driver = om.pyOptSparseDriver()
                omprob.driver.options["optimizer"] = "IPOPT"
                omprob.driver.opt_settings["tol"] = 1e-10
                omprob.driver.opt_settings["limited_memory_update_type"] = args.hessian
                omprob.driver.opt_settings["constr_viol_tol"] = 1e-10
                omprob.driver.opt_settings["dual_inf_tol"] = 1e-10
                omprob.driver.opt_settings["print_info_string"] = "yes"
                omprob.driver.opt_settings["output_file"] = os.path.join(
                    args.prefix, "ipopt_output_file%d.dat" % (step)
                )

                if args.n_mesh_refine > 1 and step == args.n_mesh_refine - 1:
                    omprob.driver.opt_settings["max_iter"] = args.niter_finest
                else:
                    omprob.driver.opt_settings["max_iter"] = args.max_iter

            omprob.setup()

            # Optimize and write result to f5 file
            omprob.run_model()
            omprob.run_driver()

            # Get optimal result from root processor and broadcast
            if comm.rank == 0:
                redu_xopt_g = omprob.get_val("indeps.x")  # Global vector
            else:
                redu_xopt_g = None
            redu_xopt_g = comm.bcast(redu_xopt_g, root=0)

            # Create a distributed vector and store the optimal solution
            # to hot-start the optimization on finer mesh
            analysis.globalVecToLocalvec(redu_xopt_g, redu_xopt)

            # Compute data of interest
            discreteness = np.dot(redu_xopt_g, 1.0 - redu_xopt_g) / len(redu_xopt_g)
            obj = omprob.get_val("topo.obj")[0]
            con = omprob.get_val("topo.con")[0]

            # Compute discreteness for rho
            redu_prob.reduDVtoDV(redu_xopt, xopt.getArray())
            problem.getTopoFilter().applyFilter(xopt, rhoopt)
            redu_prob.DVtoreduDV(rhoopt.getArray(), redu_rhoopt)
            redu_rhoopt_g = comm.allgather(np.array(redu_rhoopt))
            redu_rhoopt_g = np.concatenate(redu_rhoopt_g)
            discreteness_rho = np.dot(redu_rhoopt_g, 1.0 - redu_rhoopt_g) / len(
                redu_rhoopt_g
            )

        # Otherwise, use ParOpt.Optimizer to optimize
        else:
            if args.optimizer == "mma":
                # opt = ParOpt.Optimizer(problem, mma_options)
                opt = ParOpt.Optimizer(redu_prob, mma_options)
            else:
                # opt = ParOpt.Optimizer(problem, optimization_options)
                opt = ParOpt.Optimizer(redu_prob, paropt_options)
            opt.optimize()
            redu_xopt, z, zw, zl, zu = opt.getOptimizedPoint()

            # Get optimal objective and constraint
            fail, obj, cons = redu_prob.evalObjCon(redu_xopt)
            con = cons[0]

            # Compute discreteness
            redu_xopt_g = comm.allgather(np.array(redu_xopt))
            redu_xopt_g = np.concatenate(redu_xopt_g)
            discreteness = np.dot(redu_xopt_g, 1.0 - redu_xopt_g) / len(redu_xopt_g)

            # Compute discreteness for rho
            redu_prob.reduDVtoDV(redu_xopt, xopt.getArray())
            problem.getTopoFilter().applyFilter(xopt, rhoopt)
            redu_prob.DVtoreduDV(rhoopt.getArray(), redu_rhoopt)
            redu_rhoopt_g = comm.allgather(np.array(redu_rhoopt))
            redu_rhoopt_g = np.concatenate(redu_rhoopt_g)
            discreteness_rho = np.dot(redu_rhoopt_g, 1.0 - redu_rhoopt_g) / len(
                redu_rhoopt_g
            )

        # Manually create the f5 file
        flag = (
            TACS.OUTPUT_CONNECTIVITY
            | TACS.OUTPUT_NODES
            | TACS.OUTPUT_DISPLACEMENTS
            | TACS.OUTPUT_EXTRAS
        )
        f5 = TACS.ToFH5(problem.getAssembler(), TACS.SOLID_ELEMENT, flag)
        f5.writeToFile(os.path.join(args.prefix, "output_refine{:d}.f5".format(step)))

        # Compute infeasibility
        infeas = np.max([-con, 0])

        # Populate _xopt
        _xopt = problem.createDesignVec()
        redu_prob.reduDVtoDV(redu_xopt, _xopt)

        # Solve the generalized eigenvalue problem once to cross-check the feasibility
        ges = GeneralEigSolver(
            problem, max_jd_size=args.max_jd_size, max_gmres_size=args.max_gmres_size
        )
        try:
            evals, evecs, res = ges.compute(_xopt)
        except:
            evals, evecs, res = None, None, None
        del ges

        # print('%15s%15s'%('eigenvalue', 'residual'))
        # for ii, e in enumerate(evals): print('[%2d]%15.5e%15.5e'%(ii, e, res[ii]))

        # Export data to python pickle file
        if comm.rank == 0:
            # Check data
            print("[Optimum] discreteness:{:20.10e}".format(discreteness))
            print("[Optimum] discrete_rho:{:20.10e}".format(discreteness_rho))
            print("[Optimum] obj:         {:20.10e}".format(obj))
            print("[Optimum] con:         {:20.10e}".format(con))
            print("[Optimum] infeas:      {:20.10e}".format(infeas))

            pkl = dict()
            pkl["discreteness"] = discreteness
            pkl["discreteness_rho"] = discreteness_rho
            pkl["obj"] = obj
            pkl["con"] = con
            pkl["infeas"] = infeas
            pkl["domain"] = args.domain
            pkl["AR"] = args.AR
            pkl["ratio"] = args.ratio
            pkl["len0"] = args.len0
            pkl["r0-frac"] = args.r0_frac
            pkl["htarget"] = args.htarget
            pkl["mg-levels"] = args.mg_levels
            pkl["qval"] = args.qval
            pkl["max-jd-size"] = args.max_jd_size
            pkl["max-gmres-size"] = args.max_gmres_size
            pkl["optimizer"] = args.optimizer
            pkl["n-mesh-refine"] = args.n_mesh_refine
            pkl["max-iter"] = args.max_iter
            pkl["qn-correction"] = args.qn_correction
            pkl["eig-scale"] = args.eig_scale
            pkl["qn-subspace"] = args.qn_subspace
            pkl["cmd"] = cmd
            pkl["lambda0"] = args.lambda0
            pkl["problem"] = "frequency"
            pkl["paropt-type"] = args.paropt_type
            pkl["gep-evals"] = evals
            pkl["gep-res"] = res
            pkl["snapshot"] = redu_prob.get_snapshot()
            pkl["qn-time"] = None

            if args.qn_correction:
                pkl["qn-time"] = freq_constr.getAveragedQnTime()

            if args.optimizer == "paropt":
                pkl["curvs"] = freq_constr.getQnUpdateCurvs()
                pkl["n_skipH"] = getNSkipUpdate(
                    os.path.join(args.prefix, "tr_output_file%d.dat" % (step))
                )

            with open(
                os.path.join(args.prefix, "output_refine%d.pkl" % (step)), "wb"
            ) as f:
                pickle.dump(pkl, f)

        # Output for visualization
        assembler = problem.getAssembler()
        forest = forest.duplicate()

        # If not the final step, refine and repartition the mesh
        density_based_refine = False  # Hard-coded - don't use density based refine
        if step != args.n_mesh_refine - 1:
            if density_based_refine:
                # Refine based solely on the value of the density variable
                TopOptUtils.densityBasedRefine(forest, assembler, lower=0.05, upper=0.5)
            else:
                # Perform refinement based on distance
                dist_file = os.path.join(args.prefix, "distance_solution%d.f5" % (step))

                # Compute the characteristic domain length
                vol = lx * ly * lz
                domain_length = vol ** (1.0 / 3.0)
                refine_distance = 0.025 * domain_length
                TopOptUtils.approxDistanceRefine(
                    forest,
                    new_filter,
                    assembler,
                    refine_distance,
                    domain_length=domain_length,
                    filename=dist_file,
                )

            # Repartition the mesh
            forest.balance(1)
            forest.repartition()
