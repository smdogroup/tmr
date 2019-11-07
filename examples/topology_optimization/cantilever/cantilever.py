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

def create_forest(comm, depth, htarget=5.0, filename='cantilever.stp'):
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
    # Load the geometry model
    geo = TMR.LoadModel(filename)

    # Mark the boundary condition faces
    verts = geo.getVertices()
    faces = geo.getFaces()
    volumes = geo.getVolumes()

    # Set source and target faces
    faces[3].setName('fixed')
    faces[4].setSource(volumes[0], faces[5])
    verts[4].setName('pt1')
    verts[3].setName('pt2')

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

        # Set the output file name
        flag = (TACS.OUTPUT_CONNECTIVITY |
                TACS.OUTPUT_NODES |
                TACS.OUTPUT_DISPLACEMENTS |
                TACS.OUTPUT_STRAINS |
                TACS.OUTPUT_EXTRAS)
        self.f5 = TACS.ToFH5(assembler, TACS.SOLID_ELEMENT, flag)
        self.iter_offset = iter_offset

        self.has_plotly = True
        try:
            import plotly.graph_objects as go
            self.go = go
        except:
            self.has_plotly = False
        return

    def write_output(self, prefix, itr, oct_forest, quad_forest, x):
        self.f5.writeToFile('results/output%d.f5'%(itr + self.iter_offset))

        if self.has_plotly and oct_forest is not None:
            points = TMR.getSTLTriangles(oct_forest, x)

            n = points.shape[0]
            if n > 0 and iter % 10 == 0:
                mesh = self.go.Mesh3d(x=points[:,0], y=points[:,1], z=points[:,2],
                                      i=np.arange(0, n, 3, dtype=np.intc),
                                      j=np.arange(1, n, 3, dtype=np.intc),
                                      k=np.arange(2, n, 3, dtype=np.intc))
                self.fig = self.go.Figure(data=[mesh])
                self.fig.update_layout(scene_aspectmode='data')
                self.fig.show()

def create_problem(forest, bcs, props, nlevels, iter_offset=0):
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

    # Create the problem and filter object
    filter_type = 'matrix'
    obj = CreatorCallback(bcs, props)
    problem = TopOptUtils.createTopoProblem(forest,
        obj.creator_callback, filter_type, nlevels=nlevels)

    # Get the assembler object we just created
    assembler = problem.getAssembler()

    # Set the load
    P = 1.0e3
    force = TopOptUtils.computeVertexLoad('pt1', forest, assembler, [0, P, 0])
    temp = TopOptUtils.computeVertexLoad('pt2', forest, assembler, [0, 0, P])
    force.axpy(1.0, temp)

    # Set the load cases into the topology optimization problem
    problem.setLoadCases([force])

    # Compute the fixed mass target
    lx = 50.0 # mm
    ly = 10.0 # mm
    lz = 10.0 # mm
    vol = lx*ly*lz
    vol_frac = 0.25
    density = 2600.0
    m_fixed = vol_frac*(vol*density)

    # Set the mass constraint
    funcs = [functions.StructuralMass(assembler)]
    problem.addConstraints(0, funcs, [-m_fixed], [-1.0/m_fixed])

    # Set the objective (scale the compliance objective)
    problem.setObjective([1.0e3])

    cb = OutputCallback(assembler, iter_offset=iter_offset)
    problem.setOutputCallback(cb.write_output)

    return problem

if __name__ == '__main__':
    # Set the optimization parameters
    optimization_options = {
        # Parameters for the trust region method
        'tr_init_size': 0.01,
        'tr_max_size': 0.1,
        'tr_min_size': 0.01,
        'tr_eta': 0.25,
        'tr_penalty_gamma': 20.0,
        'tr_write_output_freq': 1,

        # Parameters for the interior point method (used to solve the
        # trust region subproblem)
        'max_qn_subspace': 2,
        'bfgs_update_type': 'Damped',
        'tol': 1e-8,
        'maxiter': 25,
        'norm_type': 'L1',
        'barrier_strategy': 'Complementarity fraction',
        'start_strategy': 'Affine step'}

    prefix = 'results'

    # Set the communicator
    comm = MPI.COMM_WORLD

    order = 2 # Order of the mesh
    nlevels = 4 # Number of multigrid levels
    forest = create_forest(comm, nlevels-1)

    # Set the boundary conditions for the problem
    bcs = TMR.BoundaryConditions()
    bcs.addBoundaryCondition('fixed')

    # Create the material properties
    material_properties = constitutive.MaterialProperties(rho=2600.0, E=70e9,
        nu=0.3, ys=350e6)
    props = TMR.StiffnessProperties(material_properties, q=8.0)

    # Set the original filter to NULL
    orig_filter = None
    xopt = None

    count = 0
    max_iterations = 3
    for step in range(max_iterations):
        # Create the problem
        iter_offset = step*optimization_options['maxiter']
        problem = create_problem(forest, bcs, props, nlevels, iter_offset=iter_offset)

        # Initialize the problem and set the prefix
        problem.initialize()
        problem.setPrefix(prefix)
        problem.setIterationCounter(count)

        # Extract the filter to interpolate design variables
        filtr = problem.getFilter()

        if orig_filter is not None:
            # Create one of the new design vectors
            x = problem.createDesignVec()
            TopOptUtils.interpolateDesignVec(orig_filter, xopt, filtr, x)
            problem.setInitDesignVars(x)

        orig_filter = filtr

        if step == max_iterations-1:
            optimization_options['maxiter'] = 10
        count += optimization_options['maxiter']

        optimization_options['output_file'] = os.path.join(prefix, 'output_file%d.dat'%(step))
        optimization_options['tr_output_file'] = os.path.join(prefix, 'tr_output_file%d.dat'%(step))

        # Optimize
        opt = TopOptUtils.TopologyOptimizer(problem, optimization_options)
        xopt = opt.optimize()

        # Output for visualization
        flag = (TACS.OUTPUT_CONNECTIVITY |
                TACS.OUTPUT_NODES |
                TACS.OUTPUT_DISPLACEMENTS |
                TACS.OUTPUT_STRAINS |
                TACS.OUTPUT_EXTRAS)
        assembler = problem.getAssembler()
        f5 = TACS.ToFH5(assembler, TACS.SOLID_ELEMENT, flag)
        f5.writeToFile(os.path.join(prefix, 'cantilever%d.f5'%(step)))

        # Refine based solely on the value of the density variable
        assembler = problem.getAssembler()
        forest = forest.duplicate()
        TopOptUtils.densityBasedRefine(forest, assembler, lower=0.05, upper=0.5)

        # Repartition the mesh
        forest.balance(1)
        forest.repartition()
