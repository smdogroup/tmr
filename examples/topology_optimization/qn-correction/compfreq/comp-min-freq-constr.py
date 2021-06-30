"""
This script performs compliance minimization with mass and frequency constraint
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
from utils_compfreq import GEP_solver, create_problem
sys.path.append('../eigenvalue')
from utils import create_forest, OmAnalysis, getNSkipUpdate


if __name__ == '__main__':

    # Create the argument parser
    p = argparse.ArgumentParser()

    # os
    p.add_argument('--prefix', type=str, default='./results')

    # Analysis
    p.add_argument('--domain', type=str, default='cantilever',
        choices=['cantilever', '3dcantilever', 'michell', 'mbb', 'lbracket'])
    p.add_argument('--AR', type=float, default=1.0)
    p.add_argument('--ratio', type=float, default=0.4)
    p.add_argument('--len0', type=float, default=1.0)
    p.add_argument('--vol-frac', type=float, default=0.4)
    p.add_argument('--r0-frac', type=float, default=0.05)
    p.add_argument('--htarget', type=float, default=1.0)
    p.add_argument('--mg-levels', type=int, default=4)
    p.add_argument('--qval', type=float, default=5.0)
    p.add_argument('--max-jd-size', type=int, default=100)
    p.add_argument('--max-gmres-size', type=int, default=30)
    p.add_argument('--lambda0', type=float, default=0.1)
    p.add_argument('--ksrho', type=float, default=1000)

    # Optimization
    p.add_argument('--optimizer', type=str, default='paropt',
        choices=['paropt', 'paropt-pyoptsparse', 'snopt', 'ipopt'])
    p.add_argument('--n-mesh-refine', type=int, default=3)
    p.add_argument('--max-iter', type=int, default=100)
    p.add_argument('--constr', type=str, default='mass',
        choices=['mass', 'massfreq'])
    p.add_argument('--qn-correction-comp', action='store_true')
    p.add_argument('--qn-correction-freq', action='store_true')
    p.add_argument('--comp-scale', type=float, default=1.0)
    p.add_argument('--eig-scale', type=float, default=1.0)
    p.add_argument('--output-level', type=int, default=0)
    p.add_argument('--simple-filter', action='store_false')
    p.add_argument('--tr-eta', type=float, default=0.25)
    p.add_argument('--tr-min', type=float, default=1e-3)
    p.add_argument('--eq-constr', action='store_true')
    p.add_argument('--qn-subspace', type=int, default=2)

    # Test
    p.add_argument('--gradient-check', action='store_true')
    p.add_argument('--test-om', action='store_true')

    # Parse arguments
    args = p.parse_args()

    mg_levels = args.mg_levels
    prefix = args.prefix

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
    stiffness_props = TMR.StiffnessProperties(material_props, k0=1e-3, eps=0.2, q=args.qval) # Try larger q val: 8, 10, 20

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
        'tr_adaptive_gamma_update': False,
        'tr_accept_step_strategy': 'filter_method',
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
        if args.constr == 'mass':
            has_freq_constr = False
        elif args.constr == 'massfreq':
            has_freq_constr = True
        problem, obj_callback = create_problem(prefix=args.prefix, domain=args.domain,
                                    forest=forest, bcs=bcs,
                                    props=stiffness_props, nlevels=mg_levels+step,
                                    lambda0=args.lambda0, ksrho=args.ksrho,
                                    has_freq_constr=has_freq_constr,
                                    vol_frac=args.vol_frac, r0_frac=args.r0_frac,
                                    len0=args.len0, AR=args.AR, ratio=args.ratio,
                                    iter_offset=iter_offset,
                                    qn_correction_comp=args.qn_correction_comp,
                                    qn_correction_freq=args.qn_correction_freq,
                                    comp_scale=args.comp_scale,
                                    eig_scale=args.eig_scale,
                                    eq_constr=args.eq_constr,
                                    max_jd_size=args.max_jd_size,
                                    max_gmres_size=args.max_gmres_size)

        # Set the prefix
        problem.setPrefix(prefix)

        # Initialize the problem and set the prefix
        problem.initialize()
        problem.setIterationCounter(count)

        if args.gradient_check:
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

            # Optimize
            prob.setup()
            prob.run_model()
            prob.run_driver()

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

            # Write result to f5 file
            analysis.write_output(prefix, step)

            # Compute data of interest
            discreteness = np.dot(xopt_global, 1.0-xopt_global) / len(xopt_global)
            obj = prob.get_val('topo.obj')[0]
            cons = prob.get_val('topo.con')

        # Otherwise, use ParOpt.Optimizer to optimize
        else:
            opt = ParOpt.Optimizer(problem, optimization_options)
            opt.optimize()
            xopt, z, zw, zl, zu = opt.getOptimizedPoint()

            # Get optimal objective and constraint
            fail, obj, cons = problem.evalObjCon(1, xopt)

            # Compute discreteness
            xopt_vals = TMR.convertPVecToVec(xopt).getArray()
            xopt_global = comm.allgather(xopt_vals)
            xopt_global = np.concatenate(xopt_global)
            discreteness = np.dot(xopt_global, 1.0-xopt_global) / len(xopt_global)

        # Compute infeasibility
        infeas = 0.0
        if args.eq_constr:
            for con in cons:
                infeas += np.abs(con)
        else:
            for con in cons:
                infeas += np.max([-con, 0])

        # Compute the smallest eigenvalue for optimal design
        gep_solver = GEP_solver(problem.getTopoFilter(), problem.getMg(),
            args.max_jd_size, args.max_gmres_size)
        min_eig = gep_solver.solve(TMR.convertPVecToVec(xopt))

        # Export data to python pickle file
        if comm.rank == 0:

            # Check data
            print('[Optimum] discreteness:{:20.10e}'.format(discreteness))
            print('[Optimum] obj:         {:20.10e}'.format(obj))
            print('[Optimum] con[0]:      {:20.10e}'.format(cons[0]))
            try:
                print('[Optimum] con[1]:      {:20.10e}'.format(cons[1]))
            except:
                pass
            print('[Optimum] infeas:      {:20.10e}'.format(infeas))

            pkl = dict()
            pkl['discreteness'] = discreteness
            pkl['obj'] = obj
            pkl['con'] = None
            pkl['cons'] = cons
            pkl['infeas'] = infeas
            pkl['domain'] = args.domain
            pkl['AR'] = args.AR
            pkl['ratio'] = args.ratio
            pkl['len0'] = args.len0
            pkl['vol-frac'] = args.vol_frac
            pkl['r0-frac'] = args.r0_frac
            pkl['htarget'] = args.htarget
            pkl['mg-levels'] = args.mg_levels
            pkl['qval'] = args.qval
            pkl['optimizer'] = args.optimizer
            pkl['n-mesh-refine'] = args.n_mesh_refine
            pkl['max-iter'] = args.max_iter
            pkl['qn-correction-comp'] = args.qn_correction_comp
            pkl['qn-correction-freq'] = args.qn_correction_freq
            pkl['comp-scale'] = args.comp_scale
            pkl['eq-constr'] = args.eq_constr
            pkl['qn-subspace'] = args.qn_subspace
            pkl['cmd'] = cmd
            pkl['problem'] = 'comp-min-freq-constr'
            pkl['min_eig'] = min_eig

            if args.optimizer == 'paropt' or args.optimizer == 'paropt-pyoptsparse':
                pkl['n_fail_qn_corr'], pkl['neg_curvs'], pkl['pos_curvs'] = \
                    obj_callback.getFailQnCorr()
                pkl['n_skipH'] = getNSkipUpdate(os.path.join(prefix, 'tr_output_file%d.dat'%(step)))

            with open(os.path.join(prefix, 'output_refine%d.pkl'%(step)), 'wb') as f:
                pickle.dump(pkl, f)

        # Output for visualization (Are these two lines needed?)
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
