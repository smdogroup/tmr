"""
This script performs eigenvalue minimization with mass constraint
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

# Import optimization libraries
from paropt.paropt_driver import ParOptDriver

# Import utility classes and functions
from utils import OctCreator, CreatorCallback, MFilterCreator, OutputCallback
from utils import FrequencyObj, MassConstr
from utils import create_forest, create_problem, OmAnalysis

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
    p.add_argument('--max-jd-size', type=int, default=100)
    p.add_argument('--max-gmres-size', type=int, default=30)

    # Optimization
    p.add_argument('--optimizer', type=str, default='paropt',
        choices=['paropt', 'paropt-pyoptsparse', 'snopt', 'ipopt'])
    p.add_argument('--n-mesh-refine', type=int, default=3)
    p.add_argument('--tr-max-iter', type=int, default=100)
    p.add_argument('--qn-correction', action='store_true')
    p.add_argument('--non-design-mass', type=float, default=10.0)
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
        'tr_max_iterations': args.tr_max_iter,
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
        problem, obj_callback = create_problem(prefix=args.prefix, domain=args.domain,
                                    forest=forest, bcs=bcs,
                                    props=stiffness_props, nlevels=mg_levels+step,
                                    vol_frac=args.vol_frac, r0_frac=args.r0_frac,
                                    len0=args.len0, AR=args.AR, ratio=args.ratio,
                                    iter_offset=iter_offset,
                                    qn_correction=args.qn_correction,
                                    non_design_mass=args.non_design_mass,
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
        if args.optimizer == 'paropt-pyoptsparse':
            prob = om.Problem()
            analysis = OmAnalysis(problem, obj_callback)
            indeps = prob.model.add_subsystem('indeps', om.IndepVarComp())
            indeps.add_output('x', x_init)
            prob.model.add_subsystem('topo', analysis)
            prob.model.connect('indeps.x', 'topo.x')
            prob.model.add_design_var('indeps.x', lower=0.0, upper=1.0)
            prob.model.add_objective('topo.obj')
            if args.eq_constr:
                prob.model.add_constraint('topo.con', lower=0.0, upper=0.0)
            else:
                prob.model.add_constraint('topo.con', lower=0.0)

            prob.driver = ParOptDriver()
            for key in optimization_options:
                prob.driver.options[key] = optimization_options[key]

            if args.qn_correction:
                prob.driver.use_qn_correction(analysis.qn_correction)

            prob.setup()
            prob.run_driver()

            xopt = problem.createDesignVec()
            xopt_vals = TMR.convertPVecToVec(xopt).getArray()
            xopt_vals[:] = prob.get_val('indeps.x')[:]

            # Write result to f5 file
            analysis.write_output(prefix, step)

        # Otherwise, use ParOpt.Optimizer to optimize
        else:
            opt = ParOpt.Optimizer(problem, optimization_options)
            opt.optimize()
            xopt, z, zw, zl, zu = opt.getOptimizedPoint()

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
