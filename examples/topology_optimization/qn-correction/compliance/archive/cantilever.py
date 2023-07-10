"""
Cantilever example using mass-constrained compliance minimization.

This example demonstrates:

1) Creating meshes using the TMR.Creator classes
2) Lagrange-type filters
3) Design-feature based adaptive refinement

Recommended arguments:

mpirun -np n python cantilever.py

This code performs a minimum compliance optimization with a fixed mass
constraint. The design domain is a prismatic beam with a square cross-section.
The domain has an aspect ratio of 5. After one cycle of optimization, the
domain is refined, refining elements where there is material, and coarsening
where there is void.
"""

from mpi4py import MPI
from tmr import TMR, TopOptUtils
from paropt import ParOpt
from tacs import TACS, elements, constitutive, functions
from egads4py import egads
import numpy as np
import argparse
import os


class OctCreator(TMR.OctConformTopoCreator):
    """
    An instance of an OctCreator class.

    This creates discretization for a Largange type filter, where the density is
    interpolated from the nodes of a coarser finite-element mesh. In this type of
    creator, the filter element mesh and the octree element mesh need be the same.
    (In a conformal filter, they must have the same element mesh but may have
    different degree of approximation.)
    """

    def __init__(self, bcs, filt, props=None):
        TMR.OctConformTopoCreator.__init__(bcs, filt)
        self.props = props

        # Create the constitutive object - one for the entire mesh
        self.con = TMR.OctConstitutive(props=props, forest=filt)

        # Create the model (the type of physics we're using)
        self.model = elements.LinearElasticity3D(self.con)

        # Set the basis functions and create the element
        self.basis = elements.LinearHexaBasis()
        self.element = elements.Element3D(self.model, self.basis)

        return

    def createElement(self, order, octant, index, weights):
        """
        Create the element for the given octant.

        This callback provides the global indices for the filter mesh and the weights
        applied to each nodal density value to obtain the element density. The
        local octant is also provided (but not used here).

        Args:
            order (int): Order of the underlying mesh
            octant (Octant): The TMR.Octant class
            index (list): List of the global node numbers referenced by the element
            weights (list): List of weights to compute the element density

        Returns:
            TACS.Element: Element for the given octant
        """
        return self.element


class CreatorCallback:
    def __init__(self, bcs, props):
        self.bcs = bcs
        self.props = props

    def creator_callback(self, forest):
        """
        Create the creator class and filter for the provided OctForest object.

        This is called for every mesh level when the topology optimization
        problem is created.

        Args:
            forest (OctForest): The OctForest for this mesh level

        Returns:
            OctTopoCreator, OctForest: The creator and filter for this forest
        """
        creator = OctCreator(self.bcs, forest, props=self.props)
        return creator, forest


def create_geo(comm, AR, ly=10.0):
    """
    Create a TMR.Model geometry object given aspect ratio of design domain
    """

    rank = comm.Get_rank()

    ctx = egads.context()

    # Dimensions
    Lx = ly * AR
    Ly = ly
    Lz = ly

    x0 = [0.0, 0.0, 0.0]
    x1 = [Lx, Ly, Lz]
    b1 = ctx.makeSolidBody(egads.BOX, rdata=[x0, x1])
    m1 = ctx.makeTopology(egads.MODEL, children=[b1])
    if rank == 0:
        m1.saveModel("geo.egads", overwrite=True)
    comm.Barrier()

    geo = TMR.LoadModel("geo.egads", print_lev=0)
    verts = []
    edges = []
    faces = []
    vols = []
    verts.extend(geo.getVertices())
    edges.extend(geo.getEdges())
    faces.extend(geo.getFaces())
    vols.extend(geo.getVolumes())

    # Set all of the matching faces
    TMR.setMatchingFaces(geo)

    # Create the geometry
    geo = TMR.Model(verts, edges, faces, vols)

    return geo


def create_forest(comm, depth, geo, htarget=5.0):
    """
    Create an initial forest for analysis and optimization

    This code loads in the model, sets names, meshes the geometry and creates
    a QuadForest from the mesh. The forest is populated with quadtrees with
    the specified depth.

    Args:
        comm (MPI_Comm): MPI communicator
        depth (int): Depth of the initial trees
        htarget (float): Target global element mesh size

    Returns:
        OctForest: Initial forest for topology optimization
    """

    # Mark the boundary condition faces
    verts = geo.getVertices()
    faces = geo.getFaces()
    volumes = geo.getVolumes()

    # Set source and target faces
    faces[0].setName("fixed")
    faces[0].setSource(volumes[0], faces[1])
    verts[4].setName("pt4")
    verts[5].setName("pt5")
    verts[6].setName("pt6")
    verts[7].setName("pt7")

    # Create the mesh
    mesh = TMR.Mesh(comm, geo)

    # Set the meshing options
    opts = TMR.MeshOptions()

    # Create the surface mesh
    mesh.mesh(htarget, opts)

    # Create a model from the mesh
    model = mesh.createModelFromMesh()

    # Create the corresponding mesh topology from the mesh-model
    topo = TMR.Topology(comm, model)

    # Create the quad forest and set the topology of the forest
    forest = TMR.OctForest(comm)
    forest.setTopology(topo)

    # Create the trees, rebalance the elements and repartition
    forest.createTrees(depth)

    return forest


class OutputCallback:
    def __init__(self, assembler, iter_offset=0):
        self.fig = None
        self.assembler = assembler
        self.xt = self.assembler.createDesignVec()

        # Set the output file name
        flag = TACS.OUTPUT_CONNECTIVITY | TACS.OUTPUT_NODES | TACS.OUTPUT_EXTRAS
        self.f5 = TACS.ToFH5(self.assembler, TACS.SOLID_ELEMENT, flag)
        self.iter_offset = iter_offset

        return

    def write_output(self, prefix, itr, oct_forest, quad_forest, x):
        self.f5.writeToFile(
            os.path.join(prefix, "output%d.f5" % (itr + self.iter_offset))
        )

        self.assembler.getDesignVars(self.xt)
        TMR.writeSTLToBin(
            os.path.join(prefix, "level_set_output%d.bstl" % (itr + self.iter_offset)),
            oct_forest,
            self.xt,
        )

        return


class MFilterCreator:
    def __init__(self, r0_frac, N, a=0.1):
        self.a = a
        self.r0_frac = r0_frac
        self.N = N

    def filter_callback(self, assemblers, filters):
        """
        Create and initialize a filter with the specified parameters
        """
        # Find the characteristic length of the domain and set the filter length scale
        r0 = self.r0_frac * self.a
        mfilter = TopOptUtils.Mfilter(self.N, assemblers, filters, dim=3, r=r0)
        mfilter.initialize()
        return mfilter


def create_problem(
    forest, bcs, props, nlevels, vol_frac=0.25, density=2600.0, iter_offset=0, AR=2.0
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

    # Characteristic length of the domain
    len0 = 10.0
    r0_frac = 0.05
    N = 20

    # Create the problem and filter object
    mfilter = MFilterCreator(r0_frac, N, a=len0)
    filter_type = mfilter.filter_callback
    obj = CreatorCallback(bcs, props)
    problem = TopOptUtils.createTopoProblem(
        forest, obj.creator_callback, filter_type, use_galerkin=True, nlevels=nlevels
    )
    # problem = TopOptUtils.createTopoProblem(forest, obj.creator_callback,
    #                                        'helmholtz', nlevels=nlevels,
    #                                        r0=r0_frac)

    # Get the assembler object we just created
    assembler = problem.getAssembler()

    # Set the load
    P = 1.0e3
    force = TopOptUtils.computeVertexLoad("pt5", forest, assembler, [0, 0, -P])
    temp = TopOptUtils.computeVertexLoad("pt6", forest, assembler, [0, P, 0])
    force.axpy(1.0, temp)

    # Set the load cases into the topology optimization problem
    problem.setLoadCases([force])

    # Compute the fixed mass target
    lx = 10.0 * AR  # mm
    ly = 10.0  # mm
    lz = 10.0  # mm
    vol = lx * ly * lz
    m_fixed = vol_frac * (vol * density)

    # Set the mass constraint
    funcs = [functions.StructuralMass(assembler)]
    problem.addConstraints(0, funcs, [-m_fixed], [-1.0 / m_fixed])

    # Set the objective (scale the compliance objective)
    problem.setObjective([1.0e3])

    cb = OutputCallback(assembler, iter_offset=iter_offset)
    problem.setOutputCallback(cb.write_output)

    return problem


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--AR", type=float, default=2.0)
    p.add_argument("--prefix", type=str, default="results")
    p.add_argument("--paropt-use-filter", action="store_true")
    p.add_argument("--qn-correction", action="store_true")
    p.add_argument("--n-mesh-refine", type=int, default=1)
    p.add_argument("--tr-max-iter", type=int, default=50)
    p.add_argument("--gradient-check", action="store_true")
    args = p.parse_args()

    # Set up ParOpt parameters
    strategy = "penalty_method"
    if args.paropt_use_filter:
        strategy = "filter_method"
    optimization_options = {
        "algorithm": "tr",
        "output_level": 0,
        "norm_type": "l1",
        "tr_init_size": 0.05,
        "tr_min_size": 1e-3,
        "tr_max_size": 10.0,
        "tr_eta": 0.25,
        "tr_infeas_tol": 1e-6,
        "tr_l1_tol": 0.0,
        "tr_linfty_tol": 0.0,
        "tr_adaptive_gamma_update": False,
        "tr_accept_step_strategy": strategy,
        "filter_sufficient_reduction": True,
        "filter_has_feas_restore_phase": True,
        "tr_use_soc": False,
        "tr_max_iterations": args.tr_max_iter,
        "penalty_gamma": 50.0,
        "qn_subspace_size": 5,
        "qn_type": "bfgs",
        "qn_diag_type": "yty_over_yts",
        "abs_res_tol": 1e-8,
        "starting_point_strategy": "affine_step",
        "barrier_strategy": "mehrotra_predictor_corrector",
        "tr_steering_barrier_strategy": "mehrotra_predictor_corrector",
        "tr_steering_starting_point_strategy": "affine_step",
        "use_line_search": False,
        "max_major_iters": 200,
    }

    prefix = args.prefix

    # Set the communicator
    comm = MPI.COMM_WORLD

    # Create prefix directory if not exist
    if comm.rank == 0 and not os.path.isdir(prefix):
        os.mkdir(prefix)

    # Barrier here
    comm.Barrier()

    nlevels = 4  # Number of multigrid levels
    geo = create_geo(comm, args.AR)
    forest = create_forest(comm, nlevels - 1, geo)

    # Set the boundary conditions for the problem
    bcs = TMR.BoundaryConditions()
    bcs.addBoundaryCondition("fixed")

    # Create the material properties
    material_properties = constitutive.MaterialProperties(
        rho=2600.0, E=70e9, nu=0.3, ys=350e6
    )
    props = TMR.StiffnessProperties(material_properties, q=8.0)

    # Set the original filter to NULL
    orig_filter = None
    xopt = None

    # Do not use density-based refinement. Use an approximate distance based refinement.
    density_based_refine = False

    count = 0
    max_iterations = args.n_mesh_refine
    for step in range(max_iterations):
        # Create the problem
        iter_offset = step * optimization_options["tr_max_iterations"]
        problem = create_problem(
            forest, bcs, props, nlevels + step, iter_offset=iter_offset
        )

        # Set the prefix
        problem.setPrefix(prefix)

        # Use Quasi-Newton Update Correction if specified
        if args.qn_correction:
            problem.useQnCorrectionComplianceObj()

        # Initialize the problem and set the prefix
        problem.initialize()
        problem.setIterationCounter(count)

        if args.gradient_check:
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

        if max_iterations > 1:
            if step == max_iterations - 1:
                optimization_options["tr_max_iterations"] = 10
        count += optimization_options["tr_max_iterations"]

        optimization_options["output_file"] = os.path.join(
            prefix, "output_file%d.dat" % (step)
        )
        optimization_options["tr_output_file"] = os.path.join(
            prefix, "tr_output_file%d.dat" % (step)
        )

        # Optimize
        opt = ParOpt.Optimizer(problem, optimization_options)
        opt.optimize()
        xopt, z, zw, zl, zu = opt.getOptimizedPoint()

        # Output for visualization
        assembler = problem.getAssembler()
        forest = forest.duplicate()

        if density_based_refine:
            # Refine based solely on the value of the density variable
            TopOptUtils.densityBasedRefine(forest, assembler, lower=0.05, upper=0.5)
        else:
            # Perform refinement based on distance
            dist_file = os.path.join(prefix, "distance_solution%d.f5" % (step))

            # Compute the characteristic domain length
            lx = 10.0 * args.AR  # mm
            ly = 10.0  # mm
            lz = 10.0  # mm
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
