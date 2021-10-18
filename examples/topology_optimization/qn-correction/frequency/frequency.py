"""
This script performs mass minimization with natural frequency constraint
"""

# Import analysis-related libraries
from tmr import TMR, TopOptUtils
from paropt import ParOpt
from tacs import TACS, elements, constitutive, functions
from egads4py import egads

# Import general-purpose libraries
import openmdao.api as om
import numpy as np
from mpi4py import MPI
import argparse
import os
import sys
import pickle

# Import optimization libraries
from paropt.paropt_driver import ParOptDriver

# Import utility classes and functions
from utils_freq import create_problem

sys.path.append('../eigenvalue')
from utils import OctCreator, CreatorCallback, MFilterCreator, OutputCallback
from utils import create_forest, OmAnalysis, getNSkipUpdate

if __name__ == '__main__':

    # Create the argument parser
    p = argparse.ArgumentParser()

    # os
    p.add_argument('--prefix', type=str, default='./results')

    # Analysis
    p.add_argument('--domain', type=str, default='cantilever',
        choices=['cantilever', 'michell', 'mbb', 'lbracket'])
    p.add_argument('--AR', type=float, default=1.0)
    p.add_argument('--ratio', type=float, default=0.4)
    p.add_argument('--len0', type=float, default=1.0)
    p.add_argument('--vol-frac', type=float, default=0.4)
    p.add_argument('--r0-frac', type=float, default=0.05)
    p.add_argument('--htarget', type=float, default=1.0)
    p.add_argument('--mg-levels', type=int, default=4)
    p.add_argument('--qval', type=float, default=5.0)
    p.add_argument('--lambda0', type=float, default=0.1)
    p.add_argument('--ksrho', type=float, default=1000)

    p.add_argument('--eig-method', type=str, default='jd', choices=['jd', 'lanczos'])
    p.add_argument('--max-jd-size', type=int, default=100)
    p.add_argument('--max-gmres-size', type=int, default=30)
    p.add_argument('--num-eigenvalues', type=int, default=10)
    p.add_argument('--max-lanczos', type=int, default=30)
    p.add_argument('--jd-recycle', type=str, default='on', choices=['on', 'off'])
    p.add_argument('--lanczos-shift', type=float, default=-10.0)
    p.add_argument('--jd-use-Amat-shift', type=str, default='on', choices=['on', 'off'])

    # Optimization
    p.add_argument('--optimizer', type=str, default='paropt',
        choices=['paropt', 'paropt-pyoptsparse', 'snopt', 'ipopt'])
    p.add_argument('--n-mesh-refine', type=int, default=3)
    p.add_argument('--max-iter', type=int, default=100)
    p.add_argument('--qn-correction', action='store_true')
    p.add_argument('--non-design-mass', type=float, default=10.0)
    p.add_argument('--eig-scale', type=float, default=1.0)
    p.add_argument('--output-level', type=int, default=0)
    p.add_argument('--simple-filter', action='store_false')
    p.add_argument('--tr-eta', type=float, default=0.25)
    p.add_argument('--tr-min', type=float, default=1e-3)
    p.add_argument('--eq-constr', action='store_true')
    p.add_argument('--qn-subspace', type=int, default=2) # Try 2, 5, 10
    p.add_argument('--qn-type', type=str, choices=['bfgs', 'scaled_bfgs'])
    p.add_argument('--paropt-type', type=str, choices=['penalty_method', 'filter_method'])

    # Test
    p.add_argument('--gradient-check', action='store_true')

    # Parse arguments
    args = p.parse_args()

    mg_levels = args.mg_levels
    prefix = args.prefix

    if args.jd_recycle == 'on':
        jd_use_recycle = True
    else:
        jd_use_recycle = False
    
    if args.jd_use_Amat_shift == 'on':
        jd_use_Amat_shift = True
    else:
        jd_use_Amat_shift = False

    # Set the communicator
    comm = MPI.COMM_WORLD

    # Create prefix directory if not exist
    if comm.rank == 0 and not os.path.isdir(prefix):
        os.mkdir(prefix)

    # Save the command and arguments that executed this script
    if comm.rank == 0:
        cmd = 'python ' + ' '.join(sys.argv)
        with open(os.path.join(prefix,'exe.sh'), 'w') as f:
            f.write(cmd + '\n')

    # Barrier here
    comm.Barrier()

    # Geometry parameters
    lx = args.len0*args.AR
    ly = args.len0
    lz = args.len0
    if args.domain == 'lbracket':
        ly = args.len0*args.ratio

    # Set up material properties
    material_props = constitutive.MaterialProperties(rho=2600.0, E=70e3, nu=0.3, ys=100.0)

    # Create stiffness properties
    stiffness_props = TMR.StiffnessProperties(material_props, k0=1e-3, eps=0.2,
        q=args.qval, qmass=args.qval) # Try larger q val: 8, 10, 20, using same penalty for mass

    # Create initial forest
    forest = create_forest(comm, lx, ly, lz, args.ratio, args.htarget, mg_levels-1, args.domain)

    # Set boundary conditions
    bcs = TMR.BoundaryConditions()
    if args.domain == 'mbb':
        bcs.addBoundaryCondition('symmetry', [0], [0.0])
        bcs.addBoundaryCondition('support', [1,2], [0.0, 0.0])
    else:
        bcs.addBoundaryCondition('fixed', [0,1,2], [0.0, 0.0, 0.0])

    # Set up ParOpt parameters
    if args.paropt_type == 'penalty_method':
        adaptive_gamma_update = True
    else:
        adaptive_gamma_update = False
    optimization_options = {
        'algorithm': 'tr',
        'output_level':args.output_level,
        'norm_type': 'l1',
        'tr_init_size': 0.05,
        'tr_min_size': args.tr_min,
        'tr_max_size': 1.0,
        'tr_eta': args.tr_eta,
        'tr_infeas_tol': 1e-6,
        'tr_l1_tol': 0.0,
        'tr_linfty_tol': 0.0,
        'tr_adaptive_gamma_update': adaptive_gamma_update,
        'tr_accept_step_strategy': args.paropt_type,
        'filter_sufficient_reduction': args.simple_filter,
        'filter_has_feas_restore_phase': True,
        'tr_use_soc': False,
        'tr_max_iterations': args.max_iter,
        'penalty_gamma': 50.0,
        'qn_subspace_size': args.qn_subspace, # try 5 or 10
        'qn_type': 'bfgs',
        'qn_diag_type': 'yty_over_yts',
        'abs_res_tol': 1e-8,
        'starting_point_strategy': 'affine_step',
        'barrier_strategy': 'mehrotra_predictor_corrector',
        'tr_steering_barrier_strategy': 'mehrotra_predictor_corrector',
        'tr_steering_starting_point_strategy': 'affine_step',
        'use_line_search': False,  # subproblem
        'max_major_iters': 200}


    # Set the original filter to NULL
    orig_filter = None
    xopt = None

    # Do not use density-based refinement. Use an approximate distance based refinement.
    density_based_refine = False

    count = 0
    max_iterations = args.n_mesh_refine
    for step in range(max_iterations):
        # Create the problem
        iter_offset = step*optimization_options['tr_max_iterations']

        # Create the optimization problem
        problem, obj_callback, constr_callback = create_problem(prefix=args.prefix, domain=args.domain,
                                                 forest=forest, bcs=bcs,
                                                 props=stiffness_props, nlevels=mg_levels+step,
                                                 lambda0=args.lambda0, ksrho=args.ksrho,
                                                 vol_frac=args.vol_frac,
                                                 r0_frac=args.r0_frac,
                                                 len0=args.len0, AR=args.AR, ratio=args.ratio,
                                                 iter_offset=iter_offset,
                                                 qn_correction=args.qn_correction,
                                                 non_design_mass=args.non_design_mass,
                                                 eig_scale=args.eig_scale,
                                                 eq_constr=args.eq_constr,
                                                 num_eigenvalues=args.num_eigenvalues,
                                                 eig_method=args.eig_method,
                                                 max_jd_size=args.max_jd_size,
                                                 jd_use_recycle=jd_use_recycle,
                                                 max_gmres_size=args.max_gmres_size,
                                                 max_lanczos=args.max_lanczos,
                                                 lanczos_shift=args.lanczos_shift,
                                                 jd_use_Amat_shift=jd_use_Amat_shift)

        # Set the prefix
        problem.setPrefix(prefix)

        # Initialize the problem and set the prefix
        problem.initialize()
        problem.setIterationCounter(count)

        if args.gradient_check:
            for i in range(3):
                # xt = problem.createDesignVec()
                # xt_vals = TMR.convertPVecToVec(xt).getArray()
                # # xt_vals[:] = np.random.rand(len(xt_vals))
                # xt_vals[:] = 0.95
                # problem.setInitDesignVars(xt)
                problem.checkGradients(1e-6)
            exit(0)

        # Extract the filter to interpolate design variables
        filtr = problem.getFilter()

        if args.optimizer == 'paropt':
            if orig_filter is not None:
                # Create one of the new design vectors
                x = problem.createDesignVec()
                TopOptUtils.interpolateDesignVec(orig_filter, xopt, filtr, x)
                problem.setInitDesignVars(x)
        else:
            if orig_filter is not None:
                x = problem.createDesignVec()
                TopOptUtils.interpolateDesignVec(orig_filter, xopt, filtr, x)
                x_init = TMR.convertPVecToVec(x).getArray()
            else:
                x = problem.createDesignVec()
                x_init = TMR.convertPVecToVec(x).getArray()
                x_init[:] = 0.95

        orig_filter = filtr

        if max_iterations > 1:
            if step == max_iterations-1:
                optimization_options['tr_max_iterations'] = 15
        count += optimization_options['tr_max_iterations']

        optimization_options['output_file'] = os.path.join(prefix, 'output_file%d.dat'%(step))
        optimization_options['tr_output_file'] = os.path.join(prefix, 'tr_output_file%d.dat'%(step))

        # Optimize with openmdao/pyoptsparse wrapper if specified
        if args.optimizer != 'paropt':
            # Broadcast local size to all processor
            local_size = len(x_init)
            sizes = [ 0 ]*comm.size
            offsets = [ 0 ]*comm.size
            sizes = comm.allgather(local_size)
            if comm.size > 1:
                for i in range(1,comm.size):
                    offsets[i] = offsets[i-1] + sizes[i-1]
            start = offsets[comm.rank]
            end = start + local_size
            src_indices = np.arange(start, end, dtype=int)

            # Create distributed openMDAO component
            prob = om.Problem()
            analysis = OmAnalysis(comm, problem, obj_callback, sizes, offsets)
            indeps = prob.model.add_subsystem('indeps', om.IndepVarComp())

            # Create global design vector
            x_init_global = comm.allgather(x_init)
            x_init_global = np.concatenate(x_init_global)
            indeps.add_output('x', x_init_global)
            prob.model.add_subsystem('topo', analysis)
            prob.model.connect('indeps.x', 'topo.x')
            prob.model.add_design_var('indeps.x', lower=0.0, upper=1.0)
            prob.model.add_objective('topo.obj')
            if args.eq_constr:
                prob.model.add_constraint('topo.con', lower=0.0, upper=0.0)
            else:
                prob.model.add_constraint('topo.con', lower=0.0)

            # Set up optimizer and options
            if args.optimizer == 'paropt-pyoptsparse':
                prob.driver = ParOptDriver()
                for key in optimization_options:
                    prob.driver.options[key] = optimization_options[key]

                if args.qn_correction:
                    prob.driver.use_qn_correction(analysis.qn_correction)

            elif args.optimizer == 'snopt':
                prob.driver = om.pyOptSparseDriver()
                prob.driver.options['optimizer'] = 'SNOPT'
                prob.driver.opt_settings['Iterations limit'] = 9999999999999
                prob.driver.opt_settings['Major feasibility tolerance'] = 1e-10
                prob.driver.opt_settings['Major optimality tolerance'] = 1e-10
                prob.driver.opt_settings['Summary file'] = os.path.join(prefix, 'snopt_output_file%d.dat'%(step))
                prob.driver.opt_settings['Print file'] = os.path.join(prefix, 'print_output_file%d.dat'%(step))
                prob.driver.opt_settings['Major print level'] = 1
                prob.driver.opt_settings['Minor print level'] = 0

                if max_iterations > 1 and step == max_iterations - 1:
                    prob.driver.opt_settings['Major iterations limit'] = 15
                else:
                    prob.driver.opt_settings['Major iterations limit'] = args.max_iter

            elif args.optimizer == 'ipopt':
                prob.driver = om.pyOptSparseDriver()
                prob.driver.options['optimizer'] = 'IPOPT'
                prob.driver.opt_settings['tol'] = 1e-10
                prob.driver.opt_settings['constr_viol_tol'] = 1e-10
                prob.driver.opt_settings['dual_inf_tol'] = 1e-10
                prob.driver.opt_settings['output_file'] = os.path.join(prefix, 'ipopt_output_file%d.dat'%(step))

                if max_iterations > 1 and step == max_iterations - 1:
                    prob.driver.opt_settings['max_iter'] = 15
                else:
                    prob.driver.opt_settings['max_iter'] = args.max_iter

            prob.setup()

            # Optimize and write result to f5 file
            prob.run_model()
            prob.run_driver()
            analysis.write_output(prefix, step)

            # Get optimal result from root processor and broadcast
            if comm.rank == 0:
                xopt_global = prob.get_val('indeps.x')
            else:
                xopt_global = None
            xopt_global = comm.bcast(xopt_global, root=0)

            # Create a distributed vector and store the optimal solution
            # to hot-start the optimization on finer mesh
            xopt = problem.createDesignVec()
            xopt_vals = TMR.convertPVecToVec(xopt).getArray()
            xopt_vals[:] = xopt_global[start:end]

            # Compute data of interest
            discreteness = np.dot(xopt_global, 1.0-xopt_global) / len(xopt_global)
            obj = prob.get_val('topo.obj')[0]
            con = prob.get_val('topo.con')[0]

        # Otherwise, use ParOpt.Optimizer to optimize
        else:
            opt = ParOpt.Optimizer(problem, optimization_options)
            opt.optimize()
            xopt, z, zw, zl, zu = opt.getOptimizedPoint()

            # Get optimal objective and constraint
            fail, obj, cons = problem.evalObjCon(1, xopt)
            con = cons[0]

            # Compute discreteness
            xopt_vals = TMR.convertPVecToVec(xopt).getArray()
            xopt_global = comm.allgather(xopt_vals)
            xopt_global = np.concatenate(xopt_global)
            discreteness = np.dot(xopt_global, 1.0-xopt_global) / len(xopt_global)

        # Compute infeasibility
        if args.eq_constr:
            infeas = np.abs(con)
        else:
            infeas = np.max([-con, 0])

        # Export data to python pickle file
        if comm.rank == 0:

            # Check data
            print('[Optimum] discreteness:{:20.10e}'.format(discreteness))
            print('[Optimum] obj:         {:20.10e}'.format(obj))
            print('[Optimum] con:         {:20.10e}'.format(con))
            print('[Optimum] infeas:      {:20.10e}'.format(infeas))

            pkl = dict()
            pkl['discreteness'] = discreteness
            pkl['obj'] = obj
            pkl['con'] = con
            pkl['infeas'] = infeas
            pkl['domain'] = args.domain
            pkl['AR'] = args.AR
            pkl['ratio'] = args.ratio
            pkl['len0'] = args.len0
            pkl['r0-frac'] = args.r0_frac
            pkl['htarget'] = args.htarget
            pkl['mg-levels'] = args.mg_levels
            pkl['qval'] = args.qval
            pkl['max-jd-size'] = args.max_jd_size
            pkl['max-gmres-size'] = args.max_gmres_size
            pkl['optimizer'] = args.optimizer
            pkl['n-mesh-refine'] = args.n_mesh_refine
            pkl['max-iter'] = args.max_iter
            pkl['qn-correction'] = args.qn_correction
            pkl['non-design-mass'] = args.non_design_mass
            pkl['eig-scale'] = args.eig_scale
            pkl['eq-constr'] = args.eq_constr
            pkl['qn-subspace'] = args.qn_subspace
            pkl['cmd'] = cmd
            pkl['problem'] = 'frequency'

            if args.optimizer == 'paropt' or args.optimizer == 'paropt-pyoptsparse':
                pkl['curvs'] = constr_callback.getQnUpdateCurvs()
                pkl['n_skipH'] = getNSkipUpdate(os.path.join(prefix, 'tr_output_file%d.dat'%(step)))

            with open(os.path.join(prefix, 'output_refine%d.pkl'%(step)), 'wb') as f:
                pickle.dump(pkl, f)

        # Output for visualization
        assembler = problem.getAssembler()
        forest = forest.duplicate()

        # If not the final step, refine and repartition the mesh
        if step != max_iterations-1:
            if density_based_refine:
                # Refine based solely on the value of the density variable
                TopOptUtils.densityBasedRefine(forest, assembler, lower=0.05, upper=0.5)
            else:
                # Perform refinement based on distance
                dist_file = os.path.join(prefix, 'distance_solution%d.f5'%(step))

                # Compute the characteristic domain length
                vol = lx*ly*lz
                domain_length = vol**(1.0/3.0)
                refine_distance = 0.025*domain_length
                TopOptUtils.approxDistanceRefine(forest, filtr, assembler, refine_distance,
                                                domain_length=domain_length,
                                                filename=dist_file)

            # Repartition the mesh
            forest.balance(1)
            forest.repartition()
