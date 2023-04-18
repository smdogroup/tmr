"""
This script performs eigenvalue minimization with mass constraint
"""
import numpy as np
from mpi4py import MPI
import argparse
import os
import sys
from tmr import TMR, TopOptUtils
from paropt import ParOpt
from tacs import constitutive, TACS
from utils import create_forest, create_problem


def get_mma_options(prefix, step, maxit):
    mma_options = {
        "algorithm": "mma",
        "mma_asymptote_contract": 0.7,
        "mma_asymptote_relax": 1.2,
        # "mma_bound_relax": 0,
        # "mma_max_asymptote_offset": 10,
        "mma_delta_regularization": 1e-05,
        "mma_eps_regularization": 0.001,
        "mma_infeas_tol": 1e-05,
        "mma_init_asymptote_offset": 0.25,
        "mma_l1_tol": 1e-06,
        "mma_linfty_tol": 1e-06,
        "mma_max_iterations": maxit,
        "mma_min_asymptote_offset": 0.01,
        "mma_use_constraint_linearization": True,
        "output_file": os.path.join(prefix, "output_file%d.dat" % (step)),
        "mma_output_file": os.path.join(prefix, "mma_output_file%d.dat" % (step)),
    }
    return mma_options


def get_args():
    # Create the argument parser
    p = argparse.ArgumentParser()

    # os
    p.add_argument("--prefix", type=str, default="./results")

    # Analysis
    p.add_argument(
        "--domain",
        type=str,
        default="cantilever",
        choices=[
            "cantilever",
            "michell",
            "mbb",
            "lbracket",
            "8pts",
            "4edges",
        ],
    )
    p.add_argument("--AR", type=float, default=1.0)
    p.add_argument("--ratio", type=float, default=0.4)
    p.add_argument("--len0", type=float, default=1.0, help="characteristic length")
    p.add_argument("--vol-frac", type=float, default=0.4, help="volume fraction")
    p.add_argument("--r0-frac", type=float, default=0.05)
    p.add_argument("--htarget", type=float, default=0.5)
    p.add_argument("--mg-levels", type=int, default=4)
    p.add_argument("--write-f5-every", type=int, default=None)

    # Optimization
    p.add_argument("--optimizer", type=str, default="mma", choices=["mma"])
    p.add_argument("--max-iters", type=int, nargs="+", default=[100])
    p.add_argument("--non-design-mass", type=float, default=None)
    p.add_argument("--eig-scale", type=float, default=1.0)

    # Test
    p.add_argument("--gradient-check", action="store_true")

    # Parse arguments
    args = p.parse_args()

    return args


if __name__ == "__main__":
    # Parse cmd arguments
    args = get_args()

    # Set the communicator
    comm = MPI.COMM_WORLD

    # Create prefix directory if not exist
    prefix = args.prefix
    if comm.rank == 0 and not os.path.isdir(prefix):
        os.mkdir(prefix)

    # Save the command and arguments that executed this script
    if comm.rank == 0:
        cmd = "python " + " ".join(sys.argv)
        with open(os.path.join(prefix, "exe.sh"), "w") as f:
            f.write(cmd + "\n")

    # Barrier here
    comm.Barrier()

    # Geometry parameters
    lx = args.len0 * args.AR
    ly = args.len0
    lz = args.len0
    if args.domain == "lbracket":
        ly = args.len0 * args.ratio

    # Set up material properties
    material_props = constitutive.MaterialProperties(
        rho=2600.0, E=70e3, nu=0.3, ys=100.0
    )

    # Create stiffness properties
    stiffness_props = TMR.StiffnessProperties(
        material_props, k0=1e-3, eps=0.2, q=5.0, qmass=5.0
    )

    # Create initial forest
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
    orig_filter = None
    xopt = None

    count = 0
    for step, maxit in enumerate(args.max_iters):
        # Create the optimization problem
        problem, obj_callback = create_problem(
            prefix=args.prefix,
            domain=args.domain,
            forest=forest,
            bcs=bcs,
            props=stiffness_props,
            nlevels=args.mg_levels + step,
            vol_frac=args.vol_frac,
            r0_frac=args.r0_frac,
            len0=args.len0,
            AR=args.AR,
            ratio=args.ratio,
            iter_offset=step * maxit,
            qn_correction=False,
            non_design_mass=args.non_design_mass,
            eig_scale=args.eig_scale,
            eq_constr=False,
            max_jd_size=200,
            max_gmres_size=30,
            write_f5_every=args.write_f5_every,
            f5_dir=os.path.join(args.prefix, "f5_refine%d" % step),
        )

        # Set the prefix
        problem.setPrefix(prefix)

        # Initialize the problem and set the prefix
        problem.initialize()
        problem.setIterationCounter(count)

        # Check gradient and exit, if specified
        if args.gradient_check:
            for i in range(3):
                xt = problem.createDesignVec()
                xt_vals = TMR.convertPVecToVec(xt).getArray()
                xt_vals[:] = np.random.rand(len(xt_vals))
                problem.setInitDesignVars(xt)
                problem.checkGradients(1e-6)
            exit(0)

        # Extract the filter to interpolate design variables
        filtr = problem.getFilter()

        if orig_filter is not None:
            # Create one of the new design vectors
            x = problem.createDesignVec()
            TopOptUtils.interpolateDesignVec(orig_filter, xopt, filtr, x)
            problem.setInitDesignVars(x)

        orig_filter = filtr

        # Set up options
        mma_options = get_mma_options(args.prefix, step, maxit)
        count += maxit

        opt = ParOpt.Optimizer(problem, mma_options)
        opt.optimize()
        xopt, z, zw, zl, zu = opt.getOptimizedPoint()

        # Get optimal objective and constraint
        fail, obj, cons = problem.evalObjCon(1, xopt)

        # Manually create the f5 file
        flag = (
            TACS.OUTPUT_CONNECTIVITY
            | TACS.OUTPUT_NODES
            | TACS.OUTPUT_DISPLACEMENTS
            | TACS.OUTPUT_EXTRAS
        )
        f5 = TACS.ToFH5(problem.getAssembler(), TACS.SOLID_ELEMENT, flag)
        f5.writeToFile(os.path.join(args.prefix, "output_refine{:d}.f5".format(step)))

        # Output for visualization
        assembler = problem.getAssembler()
        forest = forest.duplicate()

        # If not the final step, refine and repartition the mesh
        if step != len(args.max_iters) - 1:
            # Perform refinement based on distance
            dist_file = os.path.join(prefix, "distance_solution%d.f5" % (step))

            # Compute the characteristic domain length
            vol = lx * ly * lz
            domain_length = vol ** (1.0 / 3.0)
            refine_distance = 0.025 * domain_length
            TopOptUtils.approxDistanceRefine(
                forest,
                filtr,
                assembler,
                refine_distance,
                domain_length=domain_length,
                filename=dist_file,
            )

            # Repartition the mesh
            forest.balance(1)
            forest.repartition()
