Cantilever topology optimization
================================

.. code-block:: python

    from mpi4py import MPI
    from tmr import TMR, TopOptUtils
    from paropt import ParOpt
    from tacs import TACS, elements, constitutive, functions
    import numpy as np
    import argparse
    import os

    class OctCreator(TMR.OctTopoCreator):
        def __init__(self, bcs, filt, props):
            TMR.OctTopoCreator.__init__(bcs, filt)
            self.props = props

        def createElement(self, order, octant, index, weights):
            '''Create the element'''
            stiff = TMR.OctStiffness(self.props, index, weights)
            elem = elements.Solid(2, stiff)
            return elem

    class CreatorCallback:
        def __init__(self, bcs, props):
            self.bcs = bcs
            self.props = props

        def creator_callback(self, forest):
            filter = forest.duplicate()
            filter.coarsen()
            creator = OctCreator(self.bcs, filter, self.props)
            return creator, filter 

    def create_forest(comm, depth):
        # Load the geometry model
        geo = TMR.LoadModel('beam.stp')

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
        htarget = 4.0
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

    # Set the optimization parameters
    optimization_options = {
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

    prefix = 'results'
    optimization_options['output_file'] = os.path.join(prefix, 'output_file.dat')
    optimization_options['tr_output_file'] = os.path.join(prefix, 'tr_output_file.dat')

    # Set the communicator
    comm = MPI.COMM_WORLD

    # Set the type of filter to use
    filter_type = 'lagrange'

    # Compute the volume of the bracket from the known geometry
    r = 7.5
    a = 50.0
    t = 25.0
    vol = (3*a*a*t - 3*np.pi*r*r*t)*2600
    vol_frac = 0.2

    order = 2 # Order of the mesh
    nlevels = 2 # Number of multigrid levels
    forest = create_forest(comm, nlevels-1)

    # Set the boundary conditions for the problem
    bcs = TMR.BoundaryConditions()
    bcs.addBoundaryCondition('fixed')

    # Create the material properties
    rho = [2600.0]
    E = [70e9]
    nu = [0.3]
    props = TMR.StiffnessProperties(rho, E, nu)

    # Create the problem and filter object
    obj = CreatorCallback(bcs, props)
    problem, filter_obj = TopOptUtils.createTopoProblem(forest,
        obj.creator_callback, filter_type, nlevels=nlevels)

    # Get the assembler object we just created
    assembler = filter_obj.getAssembler()

    force1 = TopOptUtils.compute3DTractionLoad('surface', forest, assembler, [1, 1, 1])
    force2 = TopOptUtils.compute3DTractionLoad('surface', forest, assembler, [0, 0, 1])

    # Set the load cases into the topology optimization problem
    forces = [force1, force2]
    problem.setLoadCases(forces)

    # Set the constraints
    funcs = [functions.StructuralMass(assembler)]
    initial_mass = assembler.evalFunctions(funcs)
    m_fixed =  vol_frac*vol

    problem.addConstraints(0, funcs, [-m_fixed], [-1.0/m_fixed])
    problem.addConstraints(1, [], [], [])
    problem.setObjective([1.0, 1.0])

    # Initialize the problem and set the prefix
    problem.initialize()
    problem.setPrefix(prefix)

    opt = TopOptUtils.TopologyOptimizer(problem, optimization_options)
    xopt = opt.optimize()

    print(assembler.evalFunctions(funcs)/initial_mass)

    # Output for visualization
    flag = (TACS.ToFH5.NODES |
            TACS.ToFH5.DISPLACEMENTS |
            TACS.ToFH5.STRAINS |
            TACS.ToFH5.STRESSES |
            TACS.ToFH5.EXTRAS)
    f5 = TACS.ToFH5(assembler, TACS.PY_SOLID, flag)
    f5.writeToFile('bracket.f5')
