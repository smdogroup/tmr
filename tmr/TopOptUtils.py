from mpi4py import MPI
from tacs import TACS, elements
from tmr import TMR
from paropt import ParOpt
import numpy as np
from six import iteritems

# Options for the TopologyOptimizer class
_optimizers = ['Interior Point', 'Trust Region']
_qn_types = ['BFGS', 'SR1', 'No Hessian approx']
_norm_types = ['Infinity', 'L1', 'L2']
_barrier_types = ['Monotone', 'Mehrotra', 'Complementarity fraction']
_start_types = ['None', 'Least squares multipliers', 'Affine step']
_bfgs_updates = ['Skip negative', 'Damped']

def createTopoProblem(forest, callback, filter_type, nlevels=2,
                      repartition=True, design_vars_per_node=1,
                      r0=0.05, N=10, lowest_order=2,
                      ordering=TACS.MULTICOLOR_ORDER,
                      use_galerkin=False,
                      scale_coordinate_factor=1.0):
    """
    Create a topology optimization problem instance and a hierarchy of meshes.
    This code takes in the OctForest or QuadForest on the finest mesh level
    and creates a series of coarser meshes for analysis and optimization.
    The discretization at each level is created via a callback function that
    generates the appropriate TACSCreator object and its associated filter (the
    QuadForest or OctForest on which the design parametrization is defined.)
    The code then creates a TMRTopoFilter class which stores information about
    the design parametrization and hierarchy. It creates a multigrid object and
    finally a TMRTopoProblem instance for optimization.

    The callback function takes in a forest object, corresponding to the finite-
    element discretization and returns a creator object and a filter object in
    the following form:

    creator, filter = callback(forest)

    Args:
        callback: A callback function that takes in the forest and
                  returns the filter and the associated creator class
        filter_type (str): Type of filter to create
        forest (TMROctForest or TMRQuadForest): Forest type
        repartition (bool): Repartition the mesh
        design_vars_per_node (int): number of design variables for each node
        r0 (float): Helmholtz/matrix filter radius
        N (int): Matrix filter approximation parameter
        lowest_order (int): Lowest order mesh to create
        ordering: TACS Assembler ordering type
        use_galerkin: Use Galerkin projection to obtain coarse grid operators
        scale_coordinate_factor (float): Scale all coordinates by this factor

    Returns:
        problem (TopoProblem): The allocated topology optimization problem
    """

    # Store data
    forests = []
    filters = []
    assemblers = []

    # Balance the forest and repartition across processors
    forest.balance(1)
    if repartition:
        forest.repartition()

    # Create the forest object
    creator, filtr = callback(forest)
    forests.append(forest)
    filters.append(filtr)
    assemblers.append(creator.createTACS(forest, ordering))

    for i in range(nlevels-1):
        order = forests[-1].getMeshOrder()
        interp = forests[-1].getInterpType()
        if order > lowest_order:
            forest = forests[-1].duplicate()
            order = order-1
            forest.setMeshOrder(order, interp)
        else:
            forest = forests[-1].coarsen()
            forest.setMeshOrder(order, interp)

            # Balance and repartition if needed
            forest.balance(1)
            if repartition:
                forest.repartition()

        # Create the forest object
        creator, filtr = callback(forest)
        forests.append(forest)
        filters.append(filtr)
        assemblers.append(creator.createTACS(forest, ordering))

    # Scale the coordinates by scale_coordinates factor if it is != 1.0
    if scale_coordinate_factor != 1.0:
        for assembler in assemblers:
            X = assembler.createNodeVec()
            assembler.getNodes(X)
            X.scale(scale_coordinate_factor)
            assembler.setNodes(X)

    # Create the multigrid object
    mg = TMR.createMg(assemblers, forests, use_galerkin=use_galerkin)

    # Create the TMRTopoFilter object
    filter_obj = None
    if callable(filter_type):
        filter_obj = filter_type(assemblers, filters)
    elif isinstance(filter_type, str):
        if filter_type == 'lagrange':
            filter_obj = TMR.LagrangeFilter(assemblers, filters)
        elif filter_type == 'matrix':
            filter_obj = TMR.MatrixFilter(r0, N, assemblers, filters)
        elif filter_type == 'conform':
            filter_obj = TMR.ConformFilter(assemblers, filters)
        elif filter_type == 'helmholtz':
            filter_obj = TMR.HelmholtzFiler(r0, assemblers, filters)

    problem = TMR.TopoProblem(filter_obj, mg)

    return problem

def computeVertexLoad(name, forest, assembler, point_force):
    """
    Add a load at vertices with the given name value. The assembler object must
    be created from the forest. The point_force must be equal to the number of
    variables per node in the assembler object.

    Args:
        name (str): Name of the surface where the traction will be added
        forest (QuadForest or OctForest): Forest for the finite-element mesh
        assembler (Assembler): TACSAssembler object for the finite-element problem
        point_force (list): List of point forces to apply at the vertices

    Returns:
        Vec: A force vector containing the point load
    """

    # Get the number of variable per node from the assembler
    vars_per_node = assembler.getVarsPerNode()
    if vars_per_node != len(point_force):
        raise ValueError('Point force length must be equal to vars_per_node')

    # Create the force vector and extract the array
    force = assembler.createVec()
    force_array = force.getArray()

    # Retrieve the node numbers from the forest
    nodes = forest.getNodesWithName(name)

    comm = assembler.getMPIComm()
    node_range = forest.getNodeRange()

    # Add the point force into the force arrays
    for node in nodes:
        if ((node >= node_range[comm.rank]) and (node < node_range[comm.rank+1])):
            index = node - node_range[comm.rank]
            force_array[vars_per_node*index:vars_per_node*(index+1)] += point_force[:]

    # Match the ordering of the vector
    assembler.reorderVec(force)

    return force

def computeTractionLoad(names, forest, assembler, trac):
    """
    Add a surface traction to all quadrants or octants that touch a face or edge with
    the given name. The assembler must be created from the provided forest. The list
    trac must have a traction for each face (6) for octants or each edge (4) for
    quadrants.

    Note: This code uses the fact that the getOctsWithName or getQuadsWithName returns
    the local face or edge index touching the surface or edge in the info member.

    Args:
        names (str) or list[(str)]: Name or list of names of the surface(s) where the traction will be added
        forest (QuadForest or OctForest): Forest for the finite-element mesh
        assembler (Assembler): TACSAssembler object for the finite-element problem
        trac (list): List of tractions, one for each possible face/edge orientation

    Returns:
        Vec: A force vector containing the traction
    """

    if isinstance(forest, TMR.OctForest):
        octants = forest.getOctants()
        if isinstance(names, str):
            face_octs = forest.getOctsWithName(names)
        else:
            face_octs = []
            for name in names:
                face_octs.extend(forest.getOctsWithName(name))
    elif isinstance(forest, TMR.QuadForest):
        octants = forest.getQuadrants()
        if isinstance(names, str):
            face_octs = forest.getQuadsWithName(names)
        else:
            face_octs = []
            for name in names:
                face_octs.extend(forest.getQuadsWithName(name))

    # Create the force vector and zero the variables in the assembler
    force = assembler.createVec()
    assembler.zeroVariables()

    # Create the auxiliary element class
    aux = TACS.AuxElements()

    for i in range(len(face_octs)):
        index = face_octs[i].tag
        if index is not None:
            aux.addElement(index, trac[face_octs[i].info])

    # Keep auxiliary elements already set in the assembler
    # aux_tmp = assembler.getAuxElements()
    assembler.setAuxElements(aux)

    # Compute the residual where force = -residual
    assembler.assembleRes(force)
    force.scale(-1.0)

    # Reset the auxiliary elements
    assembler.setAuxElements(None) # (aux_tmp)

    return force

def compute3DTractionLoad(name, forest, assembler, tr):
    """
    Add a constant surface traction to all octants that touch a face or edge with
    the given name.

    Args:
        forest (QuadForest or OctForest): Forest for the finite-element mesh
        name (str): Name of the surface where the traction will be added
        assembler (Assembler): TACSAssembler object for the finite-element problem
        tr (list): The 3D components of the traction.

    Returns:
        Vec: A force vector containing the traction
    """

    # Get the basis
    element = assembler.getElements()[0]
    basis = element.getElementBasis()

    # Get the number of variables per node
    vars_per_node = assembler.getVarsPerNode()

    trac = []
    for findex in range(6):
        trac.append(elements.Traction3D(vars_per_node, findex, basis, tr))

    return computeTractionLoad(name, forest, assembler, trac)

def interpolateDesignVec(orig_filter, orig_vec, new_filter, new_vec):
    """
    This function interpolates a design vector from the original design space defined
    on an OctForest or QuadForest and interpolates it to a new OctForest or QuadForest.

    This function is used after a mesh adaptation step to get the new design space.

    Args:
        orig_filter (OctForest or QuadForest): Original filter Oct or QuadForest object
        orig_vec (PVec): Design variables on the original mesh in a ParOpt.PVec
        new_filter (OctForest or QuadForest): New filter Oct or QuadForest object
        new_vec (PVec): Design variables on the new mesh in a ParOpt.PVec (set on ouput)
    """

    # Convert the PVec class to TACSBVec
    orig_x = TMR.convertPVecToVec(orig_vec)
    if orig_x is None:
        raise ValueError('Original vector must be generated by TMR.TopoProblem')
    new_x = TMR.convertPVecToVec(new_vec)
    if new_x is None:
        raise ValueError('New vector must be generated by TMR.TopoProblem')

    if orig_x.getVarsPerNode() != new_x.getVarsPerNode():
        raise ValueError('Number of variables per node must be consistent')

    orig_map = orig_x.getNodeMap()
    new_map = new_x.getNodeMap()
    vars_per_node = orig_x.getVarsPerNode()

    # Create the interpolation class
    interp = TACS.VecInterp(orig_map, new_map, vars_per_node)
    new_filter.createInterpolation(orig_filter, interp)
    interp.initialize()

    # Perform the interpolation
    interp.mult(orig_x, new_x)

    return

def addNaturalFrequencyConstraint(problem, omega_min, **kwargs):
    """
    Add a natural frequency constraint to a TopoProblem optimization problem

    This function automatically sets good default arguments that can be
    overridden with keyword arguments passed in through kwargs.

    Args:
        problem (TopoProblem): TopoProblem optimization problem
        omega_min (float): Minimum natural frequency, Hz
        **kwargs: Frequency constraint parameters; check
                TMR documentation for more detail
    """
    # Convert the provided minimum natural frequency from
    # Hz to rad/s, square it, and make it negative to fit the
    # constraint form: omega^2 - offset >= 0.0
    offset = -(2.0*np.pi*omega_min)**2

    # Define all the possible arguments and set defaults
    opts = {'use_jd':True,
            'num_eigs':10,
            'ks_weight':50.0,
            'offset':offset,
            'sigma':-offset,
            'scale':-0.75/offset,
            'max_lanczos':100,
            'tol':1e-30,
            'eig_tol':5e-7,
            'eig_rtol':1e-6,
            'eig_atol':1e-12,
            'num_recycle':10,
            'fgmres_size':8,
            'max_jd_size':50,
            'recycle_type':'num_recycling',
            'track_eigen_iters':2}

    # Apply the user defined parameters
    for key, value in kwargs.items():
        if key in opts:
            opts[key] = value
        else:
            raise ValueError('%s is not a valid option'%(key))

    if opts['use_jd']:
        # Set the recycling strategy
        if opts['recycle_type'] == 'num_recycling':
            recycle_type = TACS.NUM_RECYCLE
        else:
            recycle_type = TACS.SUM_TWO

        problem.addFrequencyConstraint(opts['sigma'], opts['num_eigs'],
                                       opts['ks_weight'], opts['offset'],
                                       opts['scale'], opts['max_jd_size'],
                                       opts['eig_tol'], opts['use_jd'],
                                       opts['fgmres_size'], opts['eig_rtol'],
                                       opts['eig_atol'], opts['num_recycle'],
                                       opts['recycle_type'],
                                       opts['track_eigen_iters'])
    else: # use the Lanczos method
        problem.addFrequencyConstraint(opts['sigma'], opts['num_eigs'],
                                       opts['ks_weight'], opts['offset'],
                                       opts['scale'],
                                       opts['max_lanczos'], opts['tol'], 0,
                                       0, 0, 0, 0, TACS.SUM_TWO,
                                       opts['track_eigen_iters'])

    return

def densityBasedRefine(forest, assembler, index=0,
                       lower=0.05, upper=0.5, reverse=False,
                       min_lev=0, max_lev=TMR.MAX_LEVEL):
    """
    Apply a density-based refinement criteria.

    This function takes in a Quad or OctForest that has been used for analysis and its
    corresponding Assembler object. It then uses the data set in the constitutive object
    to extract the density within each element. If the density falls below the the bound
    *lower* the element is coarsened, if the density exceeds *upper* the element is
    refined. If *reverse* is set, this scheme is reversed so low design values are
    refined. The refinement is applied directly to the forest.

    Args:
        forest (QuadForest or OctForest): OctForest or QuadForest to refine
        assembler (Assembler): The TACS.Assembler object associated with forest
        index (int): The component index of the design vector used to indicate material
        lower (float): the lower limit used for coarsening
        upper (float): the upper limit used for refinement
        reverse (bool): Reverse the refinement scheme
        min_lev (int): Minimum refinement level
        max_lev (int): Maximum refinement level
    """

    # Create refinement array
    num_elems = assembler.getNumElements()
    refine = np.zeros(num_elems, dtype=np.int32)

    # Get the elements from the Assembler object
    elems = assembler.getElements()

    for i in range(num_elems):
        # Extract the design variables from the element
        dvs_per_node = elems[i].getDesignVarsPerNode()
        dvs = elems[i].getDesignVars(i)

        # Apply the refinement criteria
        if reverse:
            value = np.min(dvs[index::dvs_per_node])
            if value >= upper:
                refine[i] = -1
            elif value <= lower:
                refine[i] = 1
        else:
            value = np.max(dvs[index::dvs_per_node])
            if value >= upper:
                refine[i] = 1
            elif value <= lower:
                refine[i] = -1

    # Refine the forest
    forest.refine(refine, min_lev=min_lev, max_lev=max_lev)

    return

def approxDistanceRefine(forest, fltr, assembler, refine_distance, index=0,
                         domain_length=1.0, tfactor=0.05, cutoff=0.15,
                         filename=None, min_lev=0, max_lev=TMR.MAX_LEVEL):
    """
    Apply a distance-based refinement criteria.

    This function takes in a forest associated with the analysis, a filter associated
    with the design variables and the corresponding assembler object. An approximate
    distance function is computed using TMR which gives an approximation of the distance
    to the closest point on the domain boundary. In this case, the domain boundary is
    approximated as those points that are intermediate in [cutoff, 1-cutoff]. Since these
    are applied to the filtered (not projected) states, there will be intermediate density
    values. Finally, all elements that contain values that are within refine_distance to
    the approximate boundary are refined, while all other elements are coarseend.

    Notes: The index controls which component of the design variable is used to estimate
    the distance (useful for multimaterial cases). The tfactor controls the approximation,
    larger values of tfactor lead to more diffusive approximations, but small values may
    lead to numerical issues. The actual factor value is determined baesd on the domain
    length parameter which gives the characteristic length of the domain.

    Args:
        forest (QuadForest or OctForest): OctForest or QuadForest to refine
        filtr (QuadForest or OctForest): OctForest or QuadForest for the filter object
        assembler (Assembler): The TACS.Assembler object associated with forest
        refine_distance (float): Refine all elements within this distance
        index (int): The design variable component index (!= 0 for multimaterial cases)
        tfactor (float): Factor applied to the domain_length for computing the approx dist.
        cutoff (float): Cutoff to indicate structural interface
        min_lev (int): Minimum refinement level
        max_lev (int): Maximum refinement level
    """

    # Set up and solve for an approximate level set function
    x = assembler.createDesignVec()
    assembler.getDesignVars(x)

    # Approximate the distance to the boundary
    dist = TMR.ApproximateDistance(fltr, x, index=index, cutoff=cutoff,
                                   t=tfactor*domain_length, filename=filename)

    # Create refinement array
    num_elems = assembler.getNumElements()
    refine = np.zeros(num_elems, dtype=np.int32)

    for i in range(num_elems):
        # Apply the refinement criteria
        if dist[i] <= refine_distance:
            refine[i] = 1
        else:
            refine[i] = -1

    # Refine the forest
    forest.refine(refine, min_lev=min_lev, max_lev=max_lev)

    return

def targetRefine(forest, fltr, assembler, refine_distance,
                 interface_lev=2, interior_lev=1,
                 interface_index=-1, interior_index=0, reverse=False,
                 domain_length=1.0, tfactor=0.05, cutoff=0.15,
                 filename=None, min_lev=0, max_lev=TMR.MAX_LEVEL):
    """
    Apply a target-based refinement strategy.

    This refinement strategy employs a targeted refinement strategy. The goal is to
    refine the interface elements, defined from an approximate distance calculation,
    and the interior elements, defined as those elements with a given threshold of
    the density field that are not close to the interface, to a prescribed level at
    the first iteration. All other elements are coarsened aggressively.

    Note: The interface and interior can be computed using different indices in
    multimaterial optimization. When the interface index is negative, all materials are
    considered during the interface distance calculation.

    Args:
        forest (QuadForest or OctForest): OctForest or QuadForest to refine
        filtr (QuadForest or OctForest): OctForest or QuadForest for the filter object
        assembler (Assembler): The TACS.Assembler object associated with forest
        refine_distance (float): Refine all elements within this distance
        interface_lev (int): Target interface refinement level
        interior_lev (int): Target interior refinement level
        interface_index (int): Design variable component index for the interface problem
        interior_index (int): Design variable component index for the interior
        reverse (boolean): Reverse the sense of the interior refinement
        tfactor (float): Factor applied to the domain_length for computing the approx dist.
        cutoff (float): Cutoff to indicate structural interface
        filename (str): File name for the approximate distance calculation
        min_lev (int): Minimum refinement level
        max_lev (int): Maximum refinement level
    """

    # Set up and solve for an approximate level set function
    x = assembler.createDesignVec()
    assembler.getDesignVars(x)

    # Approximate the distance to the boundary
    dist = TMR.ApproximateDistance(fltr, x, index=interface_index, cutoff=cutoff,
                                   t=tfactor*domain_length, filename=filename)

    # Create refinement array
    num_elems = assembler.getNumElements()
    refine = np.zeros(num_elems, dtype=np.int32)

    # Compute the levels
    if isinstance(forest, TMR.OctForest):
        octants = forest.getOctants()
        lev = np.zeros(len(octants))
        for i, oc in enumerate(octants):
            lev[i] = oc.level
    elif isinstance(forest, TMR.QuadForest):
        quads = forest.getQuadrants()
        lev = np.zeros(len(quads))
        for i, quad in enumerate(quads):
            lev[i] = quad.level

    # Get the elements from the Assembler object
    elems = assembler.getElements()

    for i in range(num_elems):
        # Apply the refinement criteria
        if dist[i] <= refine_distance:
            refine[i] = interface_lev - lev[i]
        else:
            # Now check whether this is in the interior or exterior of
            # the domain
            dvs_per_node = elems[i].getDesignVarsPerNode()
            dvs = elems[i].getDesignVars(i)

            # Apply the refinement criteria
            if reverse:
                value = np.min(dvs[interior_index::dvs_per_node])
                if value >= 1.0 - cutoff:
                    refine[i] = -1
                elif value <= cutoff:
                    refine[i] = interior_lev - lev[i]
            else:
                value = np.max(dvs[interior_index::dvs_per_node])
                if value >= 1.0 - cutoff:
                    refine[i] = interior_lev - lev[i]
                elif value <= cutoff:
                    refine[i] = -1

    # Refine the forest
    forest.refine(refine, min_lev=min_lev, max_lev=max_lev)

    return

class OptionData:
    def __init__(self):
        self.options = {}
        self.types = {}
        self.values = {}
        self.desc = {}
        self.lower = {}
        self.upper = {}

        return

    def add_option(self, name, default=None, types=None,
                   values=None, desc=None, lower=None, upper=None):
        '''
        Add an option
        '''
        self.options[name] = default
        self.types[name] = types
        self.values[name] = values
        self.desc[name] = desc
        self.lower[name] = lower
        self.upper[name] = upper

        return

    def __getitem__(self, name):
        if not name in self.options:
            raise KeyError('Key %s not in OptionData'%(name))

        return self.options[name]

    def __setitem__(self, name, value):
        '''Set the item into the options dictionary'''
        if not name in self.options:
            desc = 'Key %s not in OptionData. '%(name)
            desc += 'Set new item through add_option()'
            raise KeyError(desc)
        if (self.types[name] is not None and
            not isinstance(value, self.types[name])):
            raise ValueError('Value type does not match')
        if (self.lower[name] is not None and
            value < self.lower[name]):
            raise ValueError('Value violates lower bound')
        if (self.upper[name] is not None and
            value > self.upper[name]):
            raise ValueError('Value violates upper bound')
        if (self.values[name] is not None and
            not value in self.values[name]):
            raise ValueError('Value not in value set %s'%(str(self.values[name])))

        # Set the value
        self.options[name] = value

        return

class TopologyOptimizer:
    """
    Optimizer wrapper for topology optimization problems
    """
    def __init__(self, problem, options={}):
        # Set the basic problem class
        self.options = OptionData()
        self._init_all_options()
        self.opt = None

        # Set the option names
        for name, data in iteritems(options):
            try:
                self.options[name] = data
            except ValueError as err:
                print(err)

        # Initialize the problem
        self._initialize(problem)

        return

    def _init_all_options(self):
        """
        Declare options before kwargs are processed in the init method.
        """

        self.options.add_option('optimizer', 'Trust Region', values=_optimizers,
                                desc='Type of optimization algorithm')
        self.options.add_option('tol', 1.0e-6, lower=0.0,
                                desc='Tolerance for termination')
        self.options.add_option('maxiter', 200, lower=0, types=int,
                                desc='Maximum number of iterations')

        desc = 'Finite difference step size. If no gradient check will be performed.'
        self.options.add_option('dh', desc=desc)
        self.options.add_option('norm_type', values=_norm_types,
                                desc='Norm type')
        self.options.add_option('barrier_strategy', values=_barrier_types,
                                desc='Barrier strategy')
        self.options.add_option('start_strategy', values=_start_types,
                                desc='Starting point strategy')
        self.options.add_option('penalty_gamma',
                                desc='Value of penalty parameter gamma')
        self.options.add_option('barrier_fraction',
                                desc='Barrier fraction')
        self.options.add_option('barrier_power',
                                desc='Barrier power')
        self.options.add_option('hessian_reset_freq', types=int,
                                desc='Hessian reset frequency')
        self.options.add_option('qn_type', default='BFGS', values=_qn_types,
                                desc='Type of Hessian approximation to use')
        self.options.add_option('max_qn_subspace', default=10, types=int,
                                desc='Size of the QN subspace')
        self.options.add_option('qn_diag_factor',
                                desc='QN diagonal factor')
        self.options.add_option('bfgs_update_type', values=_bfgs_updates,
                                desc='Type of BFGS update to apply')

        desc = 'Boolean to indicate if a sequential linear method should be used'
        self.options.add_option('use_sequential_linear', types=bool,
                                desc=desc)
        self.options.add_option('affine_step_multiplier_min',
                                desc='Minimum multiplier for affine step')
        self.options.add_option('init_barrier_parameter',
                                desc='Initial barrier parameter')
        self.options.add_option('relative_barrier',
                                desc='Relative barrier parameter')
        self.options.add_option('set_qn',
                                desc='Quasi-Newton')
        self.options.add_option('qn_updates', types=bool,
                                desc='Update the Quasi-Newton')

        # Line-search parameters
        self.options.add_option('use_line_search', types=bool,
                                desc='Use line search')
        self.options.add_option('max_ls_iters', types=int,
                                desc='Max number of line search iterations')
        self.options.add_option('backtrack_ls', types=bool,
                                desc='Use backtracking line search')
        self.options.add_option('armijo_param',
                                desc='Armijo parameter for line search')
        self.options.add_option('penalty_descent_frac',
                                desc='Descent fraction penalty')
        self.options.add_option('min_penalty_param',
                                desc='Minimum line search penalty')

        # GMRES parameters
        self.options.add_option('use_hvec_prod', types=bool,
                                desc='Use Hvec product with GMRES')
        self.options.add_option('use_diag_hessian', types=bool,
                                desc='Use a diagonal Hessian')
        self.options.add_option('use_qn_gmres_precon', types=bool,
                                desc='Use QN GMRES preconditioner')
        self.options.add_option('set_nk_switch_tol',
                                desc='NK switch tolerance')
        self.options.add_option('eisenstat_walker_param',
                                desc='Eisenstat Walker parameters: array([gamma, alpha])')
        self.options.add_option('gmres_tol',
                                desc='GMRES tolerances: array([rtol, atol])')
        self.options.add_option('gmres_subspace_size', types=int,
                                desc='GMRES subspace size')

        # Output options
        self.options.add_option('output_freq', types=int,
                                desc='Output frequency')
        self.options.add_option('output_file', desc='Output file name')
        self.options.add_option('major_iter_step_check', types=int,
                                desc='Major iter step check')
        self.options.add_option('output_level', types=int,
                                desc='Output level')
        self.options.add_option('grad_check_freq',
                                desc='Gradient check frequency: array([freq, step_size])')

        # Set options for the trust region method
        self.options.add_option('tr_adaptive_gamma_update', default=True, types=bool,
                                desc='Use the adaptive penalty algorithm')
        self.options.add_option('tr_min_size', default=1e-6, lower=0.0,
                                desc='Minimum trust region radius size')
        self.options.add_option('tr_max_size', default=10.0, lower=0.0,
                                desc='Maximum trust region radius size')
        self.options.add_option('tr_init_size', default=1.0, lower=0.0,
                                desc='Initial trust region radius size')
        self.options.add_option('tr_eta', default=0.25, lower=0.0, upper=1.0,
                                desc='Trust region radius acceptance ratio')
        self.options.add_option('tr_penalty_gamma', default=10.0, lower=0.0,
                                desc='Trust region penalty parameter value')
        self.options.add_option('tr_penalty_gamma_max', default=1e4, lower=0.0,
                                desc='Trust region maximum penalty parameter value')

        # Trust region convergence tolerances
        self.options.add_option('tr_infeas_tol', default=1e-5, lower=0.0,
                                desc='Trust region infeasibility tolerance (l1 norm)')
        self.options.add_option('tr_l1_tol', default=1e-5, lower=0.0,
                                desc='Trust region optimality tolerance (l1 norm)')
        self.options.add_option('tr_linfty_tol', default=1e-5, lower=0.0,
                                desc='Trust region optimality tolerance (l-infinity norm)')

        # Trust region output file name
        self.options.add_option('tr_output_file',
                                desc='Trust region output file name')
        self.options.add_option('tr_write_output_freq', default=10, types=int,
                                desc='Trust region output frequency')

        return

    def _initialize(self, problem):
        """
        Prepare the driver for execution.

        This is the final thing to run during setup.

        Args:
            problem (ParOpt.Problem): ParOpt.Problem optimization problem
        """
        # TODO:
        # - logic for different opt algorithms
        # - treat equality constraints

        opt_type = self.options['optimizer']

        # Set the limited-memory options
        max_qn_subspace = self.options['max_qn_subspace']
        if self.options['qn_type'] == 'BFGS':
            qn_type = ParOpt.BFGS
        elif self.options['qn_type'] == 'SR1':
            qn_type = ParOpt.SR1
        elif self.options['qn_type'] == 'No Hessian approx':
            qn_type = ParOpt.NO_HESSIAN_APPROX
        else:
            qn_type = ParOpt.BFGS

        # Create the problem
        if opt_type == 'Trust Region':
            # For the trust region method, you have to use a Hessian
            # approximation
            if qn_type == ParOpt.NO_HESSIAN_APPROX:
                qn_type = ParOpt.BFGS
            if max_qn_subspace < 1:
                max_qn_subspace = 1

            # Create the quasi-Newton method
            if qn_type == ParOpt.SR1:
                qn = ParOpt.LSR1(problem, subspace=max_qn_subspace)
            else:
                qn = ParOpt.LBFGS(problem, subspace=max_qn_subspace)

            # Retrieve the options for the trust region problem
            tr_min_size = self.options['tr_min_size']
            tr_max_size = self.options['tr_max_size']
            tr_eta = self.options['tr_eta']
            tr_penalty_gamma = self.options['tr_penalty_gamma']
            tr_init_size = self.options['tr_init_size']

            # Create the trust region sub-problem
            tr_init_size = min(tr_max_size, max(tr_init_size, tr_min_size))
            tr = ParOpt.TrustRegion(problem, qn, tr_init_size,
                                    tr_min_size, tr_max_size,
                                    tr_eta, tr_penalty_gamma)

            # Set the penalty parameter
            tr.setAdaptiveGammaUpdate(self.options['tr_adaptive_gamma_update'])
            tr.setPenaltyGammaMax(self.options['tr_penalty_gamma_max'])
            tr.setMaxTrustRegionIterations(self.options['maxiter'])

            # Trust region convergence tolerances
            infeas_tol = self.options['tr_infeas_tol']
            l1_tol = self.options['tr_l1_tol']
            linfty_tol = self.options['tr_linfty_tol']
            tr.setTrustRegionTolerances(infeas_tol, l1_tol, linfty_tol)

            # Trust region output file name
            if self.options['tr_output_file'] is not None:
                tr.setOutputFile(self.options['tr_output_file'])
                tr.setOutputFrequency(self.options['tr_write_output_freq'])

            # Create the interior-point optimizer for the trust region sub-problem
            opt = ParOpt.InteriorPoint(tr, 0, ParOpt.NO_HESSIAN_APPROX)
            self.tr = tr
        else:
            # Create the ParOpt object with the interior point method
            opt = ParOpt.InteriorPoint(problem, max_qn_subspace, qn_type)
            opt.setMaxMajorIterations(self.options['maxiter'])

        # Apply the options to ParOpt
        opt.setAbsOptimalityTol(self.options['tol'])
        if self.options['dh']:
            opt.checkGradients(self.options['dh'])
        if self.options['norm_type']:
            if self.options['norm_type'] == 'Infinity':
                opt.setNormType(ParOpt.INFTY_NORM)
            elif self.options['norm_type'] == 'L1':
                opt.setNormType(ParOpt.L1_NORM)
            elif self.options['norm_type'] == 'L2':
                opt.setNormType(ParOpt.L2_NORM)

        # Set barrier strategy
        if self.options['barrier_strategy']:
            if self.options['barrier_strategy'] == 'Monotone':
                barrier_strategy = ParOpt.MONOTONE
            elif self.options['barrier_strategy'] == 'Mehrotra':
                barrier_strategy = ParOpt.MEHROTRA
            elif self.options['barrier_strategy'] == 'Complementarity fraction':
                barrier_strategy = ParOpt.COMPLEMENTARITY_FRACTION
            opt.setBarrierStrategy(barrier_strategy)

        # Set starting point strategy
        if self.options['start_strategy']:
            if self.options['start_strategy'] == 'None':
                start_strategy = ParOpt.NO_START_STRATEGY
            elif self.options['start_strategy'] == 'Least squares multipliers':
                start_strategy = ParOpt.LEAST_SQUARES_MULTIPLIERS
            elif self.options['start_strategy'] == 'Affine step':
                start_strategy = ParOpt.AFFINE_STEP
            opt.setStartingPointStrategy(start_strategy)

        # Set norm type
        if self.options['norm_type']:
            if self.options['norm_type'] == 'Infinity':
                norm_type = ParOpt.INFTY_NORM
            elif self.options['norm_type'] == 'L1':
                norm_type = ParOpt.L1_NORM
            elif self.options['norm_type'] == 'L2':
                norm_type = ParOpt.L2_NORM
            opt.setBarrierStrategy(norm_type)

        # Set BFGS update strategy
        if self.options['bfgs_update_type']:
            if self.options['bfgs_update_type'] == 'Skip negative':
                bfgs_update_type = ParOpt.SKIP_NEGATIVE_CURVATURE
            elif self.options['bfgs_update_type'] == 'Damped':
                bfgs_update_type = ParOpt.DAMPED_UPDATE
            opt.setBFGSUpdateType(bfgs_update_type)

        if self.options['penalty_gamma']:
            opt.setPenaltyGamma(self.options['penalty_gamma'])

        if self.options['barrier_fraction']:
            opt.setBarrierFraction(self.options['barrier_fraction'])

        if self.options['barrier_power']:
            opt.setBarrierPower(self.options['barrier_power'])

        if self.options['hessian_reset_freq']:
            opt.setHessianResetFrequency(self.options['hessian_reset_freq'])

        if self.options['qn_diag_factor']:
            opt.setQNDiagonalFactor(self.options['qn_diag_factor'])

        if self.options['use_sequential_linear']:
            opt.setSequentialLinearMethod(self.options['use_sequential_linear'])

        if self.options['affine_step_multiplier_min']:
            opt.setStartAffineStepMultiplierMin(self.options['affine_step_multiplier_min'])

        if self.options['init_barrier_parameter']:
            opt.setInitBarrierParameter(self.options['init_barrier_parameter'])

        if self.options['relative_barrier']:
            opt.setRelativeBarrier(self.options['relative_barrier'])

        if self.options['set_qn']:
            opt.setQuasiNewton(self.options['set_qn'])

        if self.options['qn_updates']:
            opt.setUseQuasiNewtonUpdates(self.options['qn_updates'])

        if self.options['use_line_search']:
            opt.setUseLineSearch(self.options['use_line_search'])

        if self.options['max_ls_iters']:
            opt.setMaxLineSearchIters(self.options['max_ls_iters'])

        if self.options['backtrack_ls']:
            opt.setBacktrackingLineSearch(self.options['backtrack_ls'])

        if self.options['armijo_param']:
            opt.setArmijoParam(self.options['armijo_param'])

        if self.options['penalty_descent_frac']:
            opt.setPenaltyDescentFraction(self.options['penalty_descent_frac'])

        if self.options['min_penalty_param']:
            opt.setMinPenaltyParameter(self.options['min_penalty_param'])

        if self.options['use_hvec_prod']:
            opt.setUseHvecProduct(self.options['use_hvec_prod'])

        if self.options['use_diag_hessian']:
            opt.setUseDiagHessian(self.options['use_diag_hessian'])

        if self.options['use_qn_gmres_precon']:
            opt.setUseQNGMRESPreCon(self.options['use_qn_gmres_precon'])

        if self.options['set_nk_switch_tol']:
            opt.setNKSwitchTolerance(self.options['set_nk_switch_tol'])

        if self.options['eisenstat_walker_param']:
            opt.setEisenstatWalkerParameters(self.options['eisenstat_walker_param'][0],
                                             self.options['eisenstat_walker_param'][1])

        if self.options['gmres_tol']:
            opt.setGMRESTolerances(self.options['gmres_tol'][0],
                                   self.options['gmres_tol'][1])

        if self.options['gmres_subspace_size']:
            opt.setGMRESSubspaceSize(self.options['gmres_subspace_size'])

        if self.options['output_freq']:
            opt.setOutputFrequency(self.options['output_freq'])

        if self.options['output_file']:
            opt.setOutputFile(self.options['output_file'])

        if self.options['major_iter_step_check']:
            opt.setMajorIterStepCheck(self.options['major_iter_step_check'])

        if self.options['output_level']:
            opt.setOutputLevel(self.options['output_level'])

        if self.options['grad_check_freq']:
            opt.setGradCheckFrequency(self.options['grad_check_freq'])

        # This opt object will be used again when 'run' is executed
        self.opt = opt

        return

    def optimize(self):
        """
        Run the optimization problem using the settings. Runs either InteriorPoint
        or TrustRegion optimization problems

        Returns:
            PVec: Optimized design point x
        """
        # Run the optimization, everything else has been setup
        if self.options['optimizer'] == 'Trust Region':
            self.tr.optimize(self.opt)
        else:
            self.opt.optimize()

        # Keep the values of the design variables/multipliers
        x, z, zw, zl, zu = self.opt.getOptimizedPoint()

        return x
