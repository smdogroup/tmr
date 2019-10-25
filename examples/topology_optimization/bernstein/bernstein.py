"""
Bernstein-parametrized density field with thermo-elastic analysis.

This example demonstrates:

1) Creating meshes using the TMR.Creator classes
2) Analysis and optimization using higher-order elements
3) Design-feature based adaptive refinement
4) Use of a conforming filter

Recommended arguments:

mpirun -np n python bernstein.py --max_opt_iters 25 25 25 --mg_levels 2 3 4

This code performs a minimum compliance optimization with a fixed mass
constraint using high-order plane stress elements with thermo-elastic
analysis. The design parametrization using a one-degree lower density
parametrization. This makes use of the "TMRConformFilter" in which the
density parametrization is represented on the same element mesh as the
analysis, but with a different degree polynomial. Note that there are
multiple filters that inherit from the TMRConformFilter.

The mesh is updated based on whether there is material within the element.
If material has formed, then the element is refined, if there is almost
not material, the element is coarsened.
"""

from __future__ import print_function
from mpi4py import MPI
from tmr import TMR, TopOptUtils
from paropt import ParOpt
from tacs import TACS, elements, constitutive, functions
import numpy as np
import argparse
import os

class QuadConformCreator(TMR.QuadConformTopoCreator):
    """
    This class is called to create a TACSAssembler class with an underlying
    filter mesh that conforms with Assembler mesh. The input to the class
    consists of the boundary conditions, the forest of quadtrees for the
    analysis mesh, the order of the conforming filter QuadForest mesh, and
    the type of interpolant to be used.
    """
    def __init__(self, bcs, forest, order=2, interp=TMR.BERNSTEIN_POINTS,
                 design_vars_per_node=1, props=None):
        # Set the interpolation for the new filter
        super(TMR.QuadConformTopoCreator, self).__init__(bcs, forest,
                                                         order=order, interp=interp,
                                                         design_vars_per_node=design_vars_per_node)

        # Store the properties
        self.props = props
        self.element = None

    def createTopoElement(self, order, filtr):
        """
        Create an element for the entire mesh
        """

        # Create the constitutive object - one for the entire mesh
        self.con = TMR.QuadConstitutive(props=self.props, forest=filtr)

        # Create the model (the type of physics we're using)
        self.model = elements.LinearThermoelasticity2D(self.con)

        # Set the basis functions and create the element
        if order == 2:
            self.basis = elements.LinearQuadBasis()
        elif order == 3:
            self.basis = elements.QuadraticQuadBasis()
        elif order == 4:
            self.basis = elements.CubicQuadBasis()
        elif order == 5:
            self.basis = elements.QuarticQuadBasis()

        elements.TestElementBasis(self.basis)

        # Create the elemtn type
        self.element = elements.Element2D(self.model, self.basis)

        return

    def createElement(self, order, quadrant, index, filtr):
        """
        Create the element for the specified quadrant

        Args:
            order (int): The order of the element
            quadrant (TMR.Quadrant): The quadrant to be build for this element
            index (list): The global numbers for the quadrant nodes
            filtr (QuadForest): The QuadForest for the filter

        Returns:
            TACS.Element: Element to place within the Assembler
        """
        if self.element is None:
            self.createTopoElement(order, filtr)
        return self.element

class CreatorCallback:
    def __init__(self, bcs, props):
        self.bcs = bcs
        self.props = props

    def creator_callback(self, forest):
        """
        Given the forest, instantiate a creator class that will populate a
        TACSAssembler object with the elements. This allocates the
        QuadConformCreator class above and also returns the QuadForest object
        associated with the filter. This is needed in the createTopoProblem
        call.
        """
        order = forest.getMeshOrder()-1
        interp = TMR.BERNSTEIN_POINTS
        dvs_per_node = self.props.getDesignVarsPerNode()
        creator = QuadConformCreator(self.bcs, forest, order=order,
                                     design_vars_per_node=dvs_per_node,
                                     interp=interp, props=self.props)
        return creator, creator.getFilter()

def create_forest(comm, depth, htarget):
    """
    Create an initial forest for analysis. and optimization

    This code loads in the model, sets names, meshes the geometry and creates
    a QuadForest from the mesh. The forest is populated with quadtrees with
    the specified depth.

    Args:
        comm (MPI_Comm): MPI communicator
        depth (int): Depth of the initial trees
        htarget (float): Target global element mesh size

    Returns:
        QuadForest: Initial forest for topology optimization
    """
    # Load the geometry model
    geo = TMR.LoadModel('biclamped_traction.stp')

    # Mark the boundary condition faces
    verts = geo.getVertices()
    edges = geo.getEdges()
    faces = geo.getFaces()

    edges[1].setName('fixed')
    edges[9].setName('fixed')
    edges[4].setName('traction')

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
    forest = TMR.QuadForest(comm)
    forest.setTopology(topo)

    # Set parameters for later usage
    forest.createTrees(depth)

    return forest

def create_problem(forest, bcs, props, nlevels):
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

    Returns:
        TopoProblem: Topology optimization problem instance
    """
    # Allocate the creator callback function
    obj = CreatorCallback(bcs, props)

    # Create a conforming filter
    filter_type = 'conform'

    # Create the problem and filter object
    problem = TopOptUtils.createTopoProblem(forest, obj.creator_callback, filter_type,
                                            nlevels=nlevels, lowest_order=3,
                                            design_vars_per_node=design_vars_per_node)

    # Get the assembler object we just created
    assembler = problem.getAssembler()

    # Set the constraint type
    funcs = [functions.StructuralMass(assembler)]

    # Create the traction objects that will be used later..
    T = 2.5e6
    tx = np.zeros(args.order)
    ty = T*np.ones(args.order)

    # thermal_tractions = []
    # for findex in range(4):
    #     thermal_tractions.append(elements.PSThermoQuadTraction(findex, tx, ty))

    # Allocate a thermal traction boundary condition
    # force1 = TopOptUtils.computeTractionLoad('traction', forest, assembler,
    #                                          thermal_tractions)
    force1 = assembler.createVec()
    force1.getArray()[:] = 1.0
    assembler.applyBCs(force1)

    # Set the load cases
    problem.setLoadCases([force1])

    # Set the mass constraint
    # (m_fixed - m(x))/m_fixed >= 0.0
    problem.addConstraints(0, funcs, [-m_fixed], [-1.0/m_fixed])
    problem.setObjective(obj_array) # , [functions.Compliance(assembler)])

    # Initialize the problem and set the prefix
    problem.initialize()

    # Set the output file name
    flag = (TACS.OUTPUT_CONNECTIVITY |
            TACS.OUTPUT_NODES |
            TACS.OUTPUT_DISPLACEMENTS |
            TACS.OUTPUT_STRAINS |
            TACS.OUTPUT_EXTRAS)
    problem.setF5OutputFlags(1, TACS.PLANE_STRESS_ELEMENT, flag)

    return problem

# Set the optimization parameters
optimization_options = {
    'optimizer': 'Interior Point',

    # Parameters for the trust region method
    'tr_init_size': 0.01,
    'tr_max_size': 0.05,
    'tr_min_size': 1e-6,
    'tr_eta': 0.25,
    'tr_penalty_gamma': 20.0,

    # Parameters for the interior point method (used to solve the
    # trust region subproblem)
    'max_qn_subspace': 2,
    'tol': 1e-8,
    'maxiter': 500,
    'norm_type': 'L1',
    'barrier_strategy': 'Complementarity fraction',
    'start_strategy': 'Affine step'}

# Create an argument parser to read in arguments from the command line
p = argparse.ArgumentParser()
p.add_argument('--prefix', type=str, default='./results')
p.add_argument('--vol_frac', type=float, default=0.25)
p.add_argument('--htarget', type=float, default=2.5e-3)
p.add_argument('--max_opt_iters', type=int, nargs='+',
               default=[5])
p.add_argument('--init_depth', type=int, default=1)
p.add_argument('--mg_levels', type=int, nargs='+', default=[2])
p.add_argument('--use_L1', action='store_true', default=True)
p.add_argument('--use_Linf', action='store_true', default=False)
p.add_argument('--order', type=int, default=4)
p.add_argument('--q_penalty', type=float, default=8.0)
p.add_argument('--omega', type=float, default=1.0)
args = p.parse_args()

# Set the communicator
comm = MPI.COMM_WORLD

# Print out all of the arguments to the command line
if comm.rank == 0:
    for arg in vars(args):
        print('%-20s'%(arg), getattr(args, arg))

# Set the max nmber of iterations
mg_levels = args.mg_levels
max_iterations = len(mg_levels)

# Input-dependent parameters
optimization_options['output_file'] = os.path.join(args.prefix, 'output_file.dat')
optimization_options['tr_output_file'] = os.path.join(args.prefix, 'tr_output_file.dat')

# Compute the volume of the bracket
r = 0.06
a = 0.04
vol = r*a
vol_frac = args.vol_frac

# Set the values of the objective array
obj_array = [ 1e0 ]

# The thickness value
t = 0.01

# Create the first material properties object
rho = 2600.0*t
E = 70e9*t
nu = 0.3
alpha = 23.5e-6
kcond = 130.0
ys = 450e6
mat1 = constitutive.MaterialProperties(rho=rho, E=E,
                                       nu=nu, alpha=alpha,
                                       kappa=kcond, ys=ys)

# Create the second material properties object
rho = 1300.0*t
E = 35e9*t
nu = 0.3
alpha = 0.5*23.5e-6
kcond = 65.0
ys = 275e6
mat2 = constitutive.MaterialProperties(rho=rho, E=E,
                                       nu=nu, alpha=alpha,
                                       kappa=kcond, ys=ys)

prop_list = [mat1, mat2]

# Set the fixed mass
max_density = 0.5*t*(2600.0 + 1300.0)
initial_mass = vol*max_density
m_fixed = vol_frac*initial_mass

# Set the number of variables per node
design_vars_per_node = 1
if (len(prop_list) > 1):
    design_vars_per_node = 1+len(prop_list)

# Create the stiffness properties object
props = TMR.StiffnessProperties(prop_list, q=args.q_penalty, qtemp=0.0,
                                qcond=0.0, eps=0.2, k0=1e-3)

# Set the boundary conditions for the problem
bcs = TMR.BoundaryConditions()
bcs.addBoundaryCondition('fixed', [0, 1, 2], [0.0, 0.0, 0.0])

time_array = np.zeros(sum(args.max_opt_iters[:]))
t0 = MPI.Wtime()

# Create the initial forest
forest = create_forest(comm, args.init_depth, args.htarget)
forest.setMeshOrder(args.order, TMR.GAUSS_LOBATTO_POINTS)

# Set the original filter to NULL
orig_filter = None
xopt = None

for step in range(max_iterations):
    # Create the TMRTopoProblem instance
    nlevels = mg_levels[step]
    problem = create_problem(forest, bcs, props, nlevels)
    problem.setPrefix(args.prefix)

    # Extract the filter to interpolate design variables
    filtr = problem.getFilter()

    if orig_filter is not None:
        # Create one of the new design vectors
        x = problem.createDesignVec()
        TopOptUtils.interpolateDesignVec(orig_filter, xopt, filtr, x)
        problem.setInitDesignVars(x)

    # Set the new original filter
    orig_filter = filtr

    # Set the max number of iterations
    optimization_options['maxiter'] = args.max_opt_iters[step]

    # Optimize the problem
    opt = TopOptUtils.TopologyOptimizer(problem, optimization_options)
    opt.opt.checkGradients(1e-6)
    xopt = opt.optimize()

    # Refine based solely on the value of the density variable
    assembler = problem.getAssembler()
    if design_vars_per_node == 1:
        lower = 0.05
        upper = 0.5
        TopOptUtils.densityBasedRefine(forest, assembler, lower=lower, upper=upper)
    else:
        # Refine based on the topology variable = 1 - t. Since the topology
        # variable is stored as 1-t, we reverse the refinement criterion.
        index = 2
        lower = 0.5
        upper = 0.95
        TopOptUtils.densityBasedRefine(forest, assembler, index=index,
            lower=lower, upper=upper, reverse=True)

    # Repartition the mesh
    forest.balance(1)
    forest.repartition()

# Do an averaging of all the values and write to text file
new_array=np.zeros(len(time_array))
comm.Allreduce(time_array, new_array,op=MPI.SUM)
new_array /= comm.size
if comm.rank == 0:
    np.savetxt(os.path.join(args.prefix,'time.out'), new_array, fmt='%1.6e')
