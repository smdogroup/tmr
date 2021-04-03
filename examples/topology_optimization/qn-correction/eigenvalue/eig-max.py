"""
Smallest eigenvalue (for the generalized eigenvalue problem) with mass constraint
"""

# Import analysis-related libraries
from tmr import TMR, TopOptUtils
from paropt import ParOpt
from tacs import TACS, elements, constitutive, functions
from egads4py import egads

# Import general-purpose libraries
import numpy as np
from mpi4py import MPI
import argparse
import os
import sys

# Import the cantilever example setup
sys.path.append('../../cantilever/')
from cantilever import OctCreator, CreatorCallback, OutputCallback, MFilterCreator

class FrequencyObj:
    """
    A class that evaluates the smallest eigenvalue, the objective is evaluated
    using an objective callback. We also add non-design mass to loaded nodes in
    order to form a well-posed frequency maximization problem
    """

    def __init__(self, forest, eig_scale=1e2, con_scale=1.0,
                 num_eigenvalues=10, max_jd_size=50, max_gmres_size=15, ksrho=50):
        """
        Args:
            eig_scale: scale the eigenvalues internally in order to acquire better
                       KS approximation with smaller skrho
            con_scale: scale the mass constraint
            num_eigenvalues: number of smallest eigenvalues to compute
            ksrho: KS parameter
        """

        # Set objects
        self.forest = forest

        # Set up parameters
        self.eig_scale = eig_scale
        self.con_scale = con_scale
        self.num_eigenvalues = num_eigenvalues
        self.max_jd_size = max_jd_size
        self.max_gmres_size = max_gmres_size
        self.ksrho = ksrho

        self.fltr = None
        self.mg = None
        self.assembler = None
        self.comm = None
        self.oper = None
        self.jd = None
        self.kmat = None
        self.mmat = None
        self.temp = None
        self.eigs = None
        self.eigvs = None

        return

    def objective(self, fltr, mg):
        """
        Evaluate the KS aggregation of the smallest eigenvalue for the generalized
        eigenvalue problem:

        K \Phi = M \Phi \Lambda

        And the aggregation is:
        ks = \lambda_m - \frac{1}{p} \ln{\sum_{i=1}^N{e^{-p(\lambda_i - \lambda_m)}}}

        where:
          - lambda_m = min(lambda_i)
          - p: KS parameter, in code it's called ksrho
          - N: number of eigenvalues computed
        """

        if self.fltr is None:
            self.mg = mg
            self.fltr = fltr
            self.assembler = self.fltr.getAssembler()
            self.comm = self.assembler.getMPIComm()

            # Initialize space for matrices and vectors
            self.mmat = self.assembler.createMat()
            self.kmat = self.assembler.createMat()
            self.temp = self.assembler.createVec()
            self.eigs = np.zeros(self.num_eigenvalues)
            self.eigvs = []
            for i in range(self.num_eigenvalues):
                self.eigvs.append(self.assembler.createVec())
            self.deigs = []
            for i in range(self.num_eigenvalues):
                self.deigs.append(self.assembler.createDesignVec())

            # Create the Jacobi-Davidson operator
            self.oper = TACS.JDFrequencyOperator(self.assembler, self.kmat, self.mmat,
                                                 self.mg.getMat(), self.mg)

            # Create the eigenvalue solver and set the number of recycling eigenvectors
            self.jd = TACS.JacobiDavidson(self.oper, self.num_eigenvalues,
                                          self.max_jd_size, self.max_gmres_size)
            self.jd.setTolerances(1e-6, 1e-8)
            self.jd.setRecycle(self.num_eigenvalues)

        # Assemble the stiffness matrix for the generalized eigenvalue problem
        self.assembler.assembleMatType(TACS.STIFFNESS_MATRIX, self.kmat)

        # Get current design variable
        dv = self.assembler.createDesignVec()
        self.assembler.getDesignVars(dv)

        # Add non-design mass to the mass matrix by modifying the design variable
        mvec = TopOptUtils.createNonDesignMassVec(self.assembler, ['pt5'],
                                                  self.comm, self.forest, m0=200.0)
        self.fltr.applyFilter(mvec, mvec)

        # Update dv <- dv + mvec and update design variable
        dv.axpy(1.0, mvec)
        self.assembler.setDesignVars(dv)

        # Construct mass matrix
        self.assembler.assembleMatType(TACS.MASS_MATRIX, self.mmat)

        # Reset design variable
        dv.axpy(-1.0, mvec)
        self.assembler.setDesignVars(dv)
        self.mg.assembleMatType(TACS.STIFFNESS_MATRIX)
        self.mg.factor()

        # Solve
        self.jd.solve(print_flag=True, print_level=0)

        # Check if succeeded, otherwise try again
        if self.jd.getNumConvergedEigenvalues() < self.num_eigenvalues:
            self.jd.solve(print_flag=True, print_level=1)

        # Extract eigenvalues and eigenvectors
        for i in range(self.num_eigenvalues):
            self.eigs[i], _ = self.jd.extractEigenvector(i, self.eigvs[i])

        # Scale eigenvalues for a better KS approximation
        self.eigs[:] *= self.eig_scale

        # Compute the minimal eigenvalue
        eig_min = np.min(self.eigs)

        # Compute KS aggregation
        self.eta = np.exp(-self.ksrho*(self.eigs - eig_min))
        self.beta = np.sum(self.eta)
        ks = (eig_min - np.log(self.beta)/self.ksrho)
        self.eta = self.eta/self.beta

        # Scale eigenvalue back
        self.eigs[:] /= self.eig_scale

        # Objective
        obj = -ks

        # Print values
        if self.comm.rank == 0:
            print('{:30s}{:20.10e}'.format('[Obj] KS eigenvalue:', ks))
            print('{:30s}{:20.10e}'.format('[Obj] min eigenvalue:', eig_min))
            print('{:30s}{:20.10e}'.format('[Obj] objective:', obj))

        return obj


    def objective_gradient(self, fltr, mg, dfdrho):
        """
        Compute the gradient of frequency ks aggragation w.r.t. \rho: \frac{\partial ks}{\partial \rho}

        Note that if the filter takes the following matrix form:
        \rho = Fx
        where \rho is the filtered nodal density, x is raw design vector, then
        we have:
        \frac{\partial ks}{\partial x} = F^T \frac{\partial ks}{\partial \rho}
        But the the filter transpose will be applied in the c++ function that calls this callback,
        so here we only compute dksdrho

        The derivative of ks w.r.t. \rho is:
        \frac{\partial ks}{\partial \rho}_i = \sum_j{\eta_j \frac{\partial \lambda_j}{\partial \rho_i}}
        where:
        \eta_j = \frac{e^{-p(\lambda_j - \lambda_m)}}{\sum_k{e^{-p(\lambda_k - \lambda_m)}}}
        and:
        \frac{\partial \lambda_j}{\partial \rho_i} = \Phi_j^T (\frac{\partial K}{\partial \rho_i} -
        \lambda_j \frac{\partial M}{\rho_i}) \Phi_j
        """

        # Zero out the gradient vector
        dfdrho.zeroEntries()

        temp = self.assembler.createDesignVec()

        for i in range(self.num_eigenvalues):

            # This is a maximization problem
            scale = -self.eta[i]*self.eig_scale

            # M-orthogonalize eigenvectors
            self.temp.zeroEntries()
            self.mmat.mult(self.eigvs[i], self.temp)
            ortho = self.temp.dot(self.eigvs[i])
            if self.comm.rank == 0:
                print("eig[{:2d}] M-orthogonality:{:20.10e}".format(i, ortho))
            scale /= ortho  # Always 1.0

            # temp.zeroEntries()
            # self.assembler.addMatDVSensInnerProduct(
            #     -scale*self.eigs[i], TACS.MASS_MATRIX,
            #     self.eigvs[i], self.eigvs[i], temp)
            # self.assembler.addMatDVSensInnerProduct(
            #     scale, TACS.STIFFNESS_MATRIX,
            #     self.eigvs[i], self.eigvs[i], temp)

            # dfdrho.axpy(1.0, temp)

            self.assembler.addMatDVSensInnerProduct(
                -scale*self.eigs[i], TACS.MASS_MATRIX,
                self.eigvs[i], self.eigvs[i], dfdrho)
            self.assembler.addMatDVSensInnerProduct(
                scale, TACS.STIFFNESS_MATRIX,
                self.eigvs[i], self.eigvs[i], dfdrho)

        dfdrho.beginSetValues(op=TACS.ADD_VALUES)
        dfdrho.endSetValues(op=TACS.ADD_VALUES)

        rand = self.assembler.createDesignVec()
        rand.setRand()
        rnorm = rand.norm()
        norm = dfdrho.norm()
        proj = dfdrho.dot(rand)
        if self.comm.rank == 0:
            print("{:30s}{:20.10e}".format('[Obj] gradient norm:', norm))
            print("{:30s}{:20.10e}".format('[Obj] random vec norm:', rnorm))
            print("{:30s}{:20.10e}".format('[Obj] random proj:', proj))

        return

    def qn_correction(self, x, z, zw, s, y):
        """
        Update y:
        y <- y + F^T P Fs

        where:
        F: filter matrix
        P: Positive definite part of the objective Hessian

        Note:
        x is raw design variable (unfiltered) and it's NOT equal to
        the design variable in the assembler:

        if:
        self.assembler.getDesignVars(dv)

        Then:
        dv == Fx
        """

        update = self.assembler.createDesignVec()
        temp = self.assembler.createDesignVec()
        rand = self.assembler.createDesignVec()
        rand.setRand()

        h = 1e-6

        # Create PVec type wrappers
        x_wrap = TMR.convertPVecToVec(x)
        s_wrap = TMR.convertPVecToVec(s)

        # Get current nodal density
        rho = self.assembler.createDesignVec()
        self.assembler.getDesignVars(rho)

        self.fltr.applyFilter(s_wrap, s_wrap)

        if self.comm.rank == 0:
            # Just want to make sure the callback is implemented correctly
            print("calling qn_correction!...")

        for i in range(self.num_eigenvalues):
            """
            Compute the first part using finite difference
            P += phi^T d2Kdx2 phi
            """

            # Compute g(rho + h*s)
            rho.axpy(h, s_wrap)
            self.assembler.setDesignVars(rho)
            temp.zeroEntries()
            self.assembler.addMatDVSensInnerProduct(self.eta[i], TACS.STIFFNESS_MATRIX,
                self.eigvs[i], self.eigvs[i], temp)
            # print("g+dg:", temp.norm())

            # Compute g(rho)
            rho.axpy(-h, s_wrap)
            self.assembler.setDesignVars(rho)
            self.assembler.addMatDVSensInnerProduct(-self.eta[i], TACS.STIFFNESS_MATRIX,
                self.eigvs[i], self.eigvs[i], temp)
            # print("dg:", temp.norm())

            # Compute dg/h
            temp.scale(1/h)

            # Add to the update
            update.axpy(1.0, temp)

            """
            Compute the second part:
            P -= phi^T (deig*dMdx) phi^T
            """
            temp.zeroEntries()
            coeff = self.eta[i]*self.deigs[i].dot(s_wrap)
            self.assembler.addMatDVSensInnerProduct(-coeff,
                TACS.MASS_MATRIX, self.eigvs[i], self.eigvs[i], temp)

            # Get update
            update.axpy(1.0, temp)

        # Compute curvature and check the norm of the update
        # to see if the magnitude makes sense
        x_norm = x.norm()
        s_norm = s.norm()
        y_norm = y.norm()
        dy_norm = update.norm()
        curvature = s_wrap.dot(update)
        if self.comm.rank == 0:
            print("curvature: {:20.10e}".format(curvature))
            print("norm(x):   {:20.10e}".format(x_norm))
            print("norm(s):   {:20.10e}".format(s_norm))
            print("norm(y):   {:20.10e}".format(y_norm))
            print("norm(dy):  {:20.10e}".format(dy_norm))

        # Update y
        if curvature > 0:
            y_wrap = TMR.convertPVecToVec(y)
            y_wrap.axpy(1.0, temp)

        return

class MassConstr:
    """
    Mass constraint takes the form of c = -m/m_fixed + 1 >= 0
    """

    def __init__(self, m_fixed, comm):

        self.m_fixed = m_fixed
        self.comm = comm
        self.rank = self.comm.Get_rank()

        self.assembler = None
        self.fltr = None
        self.mass_func = None

        return

    def constraint(self, fltr, mg):
        if self.fltr is None:
            self.fltr = fltr
            self.assembler = self.fltr.getAssembler()
            self.mass_func = functions.StructuralMass(self.assembler)

        # Eval mass
        mass = self.assembler.evalFunctions([self.mass_func])[0]
        mass_constr = -mass/self.m_fixed + 1
        if self.rank == 0:
            print("{:30s}{:20.10e}".format('[Con] mass constraint:',mass_constr))

        return [mass_constr]

    def constraint_gradient(self, fltr, mg, vecs):
        # We only have one constraint
        dcdrho = vecs[0]
        dcdrho.zeroEntries()

        # Evaluate the mass gradient
        self.assembler.addDVSens([self.mass_func], [dcdrho], alpha=-1/self.m_fixed)

        # Compute norm
        norm = dcdrho.norm()
        if self.rank == 0:
            print("{:30s}{:20.10e}".format('[Con] gradient norm:', norm))
        return

def create_geo(comm, AR, ly=10.0):
    """
    Create a TMR.Model geometry object given aspect ratio of design domain
    """

    rank = comm.Get_rank()

    ctx = egads.context()

    # Dimensions
    Lx = ly*AR
    Ly = ly
    Lz = ly

    x0 = [0.0, 0.0, 0.0]
    x1 = [Lx, Ly, Lz]
    b1 = ctx.makeSolidBody(egads.BOX, rdata=[x0, x1])
    m1 = ctx.makeTopology(egads.MODEL, children=[b1])
    if rank == 0:
        m1.saveModel('geo.egads', overwrite=True)
    comm.Barrier()

    geo = TMR.LoadModel('geo.egads', print_lev=0)
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
    # Get MPI rank
    rank = comm.Get_rank()

    # Mark the boundary condition faces
    verts = geo.getVertices()
    faces = geo.getFaces()
    volumes = geo.getVolumes()

    if rank == 0:
        for index, vert in enumerate(verts):
            x,y,z = vert.evalPoint()
            print("verts[{:d}]: ({:.2f}, {:.2f}, {:.2f})".format(index, x, y, z))

        for index, face in enumerate(faces):
            umin, vmin, umax, vmax = face.getRange()
            umid = 0.5*(umin + umax)
            vmid = 0.5*(vmin + vmax)
            x,y,z = face.evalPoint(umid, vmid)
            print("faces[{:d}]: ({:.2f}, {:.2f}, {:.2f})".format(index, x, y, z))



    # Set source and target faces
    faces[0].setName('fixed')
    faces[0].setSource(volumes[0], faces[1])
    verts[4].setName('pt4')
    verts[5].setName('pt5')
    verts[6].setName('pt6')
    verts[7].setName('pt7')

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

def create_problem(forest, bcs, props, nlevels, vol_frac=0.25,
                   density=2600.0, iter_offset=0, AR=2.0, qn_correction=True):
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
    problem = TopOptUtils.createTopoProblem(forest, obj.creator_callback,
                                            filter_type, use_galerkin=True,
                                            nlevels=nlevels)

    # Get the assembler object we just created
    assembler = problem.getAssembler()

    # Compute the fixed mass target
    lx = 10.0*AR # mm
    ly = 10.0 # mm
    lz = 10.0 # mm
    vol = lx*ly*lz
    m_fixed = vol_frac*(vol*density)

    # Add objective callback
    obj_callback = FrequencyObj(forest)
    problem.addObjectiveCallback(obj_callback.objective, obj_callback.objective_gradient)

    # Add constraint callback
    constr_callback = MassConstr(m_fixed, assembler.getMPIComm())
    problem.addConstraintCallback(1, constr_callback.constraint, constr_callback.constraint_gradient)

    # Use Quasi-Newton Update Correction if specified
    if qn_correction:
        problem.addQnCorrectionCallback(1, obj_callback.qn_correction)

    # Set output callback
    cb = OutputCallback(assembler, iter_offset=iter_offset)
    problem.setOutputCallback(cb.write_output)

    return problem

def bind(instance, func, as_name=None):
    """
    Bind the function *func* to *instance*, with either provided name *as_name*
    or the existing name of *func*. The provided *func* should accept the
    instance as the first argument, i.e. "self".
    """
    if as_name is None:
        as_name = func.__name__
    bound_method = func.__get__(instance, instance.__class__)
    setattr(instance, as_name, bound_method)
    return bound_method

if __name__ == '__main__':

    # Create the argument parser
    p = argparse.ArgumentParser()
    p.add_argument('--AR', type=float, default=2.0)
    p.add_argument('--prefix', type=str, default='./results')
    p.add_argument('--n-mesh-refine', type=int, default=1)
    p.add_argument('--tr-max-iter', type=int, default=50)
    p.add_argument('--paropt-use-filter', action='store_true')
    p.add_argument('--qn-correction', action='store_true')
    p.add_argument('--gradient-check', action='store_true')
    args = p.parse_args()

    htarget = 10.0
    mg_levels = 4
    prefix = args.prefix

    # Set the communicator
    comm = MPI.COMM_WORLD

    # Create prefix directory if not exist
    if comm.rank == 0 and not os.path.isdir(prefix):
        os.mkdir(prefix)

    # Barrier here
    comm.Barrier()

    # Set up material properties
    material_props = constitutive.MaterialProperties(rho=2600.0, E=70e3, nu=0.3, ys=100.0)

    # Create stiffness properties
    stiffness_props = TMR.StiffnessProperties(material_props, k0=1e-3, eps=0.2, q=5.0)

    # Set boundary conditions
    bcs = TMR.BoundaryConditions()
    bcs.addBoundaryCondition('fixed', [0,1,2], [0.0, 0.0, 0.0])

    # Create initial forest
    geo = create_geo(comm, args.AR)
    forest = create_forest(comm, mg_levels-1, geo, htarget)

    # Set up ParOpt parameters
    strategy = 'penalty_method'
    if args.paropt_use_filter:

        strategy = 'filter_method'
    optimization_options = {
        'algorithm': 'tr',
        'output_level':0,
        'norm_type': 'l1',
        'tr_init_size': 0.05,
        'tr_min_size': 1e-3,
        'tr_max_size': 10.0,
        'tr_eta': 0.25,
        'tr_infeas_tol': 1e-6,
        'tr_l1_tol': 0.0,
        'tr_linfty_tol': 0.0,
        'tr_adaptive_gamma_update': False,
        'tr_accept_step_strategy': strategy,
        'filter_sufficient_reduction': True,
        'filter_has_feas_restore_phase': True,
        'tr_use_soc': False,
        'tr_max_iterations': args.tr_max_iter,
        'penalty_gamma': 50.0,
        'qn_subspace_size': 5,
        'qn_type': 'bfgs',
        'qn_diag_type': 'yty_over_yts',
        'abs_res_tol': 1e-8,
        'starting_point_strategy': 'affine_step',
        'barrier_strategy': 'mehrotra_predictor_corrector',
        'tr_steering_barrier_strategy': 'mehrotra_predictor_corrector',
        'tr_steering_starting_point_strategy': 'affine_step',
        'use_line_search': False,
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
        problem = create_problem(forest, bcs, stiffness_props, mg_levels+step,
                                 vol_frac=0.25, iter_offset=iter_offset,
                                 AR=args.AR, qn_correction=args.qn_correction)

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

        if orig_filter is not None:
            # Create one of the new design vectors
            x = problem.createDesignVec()
            TopOptUtils.interpolateDesignVec(orig_filter, xopt, filtr, x)
            problem.setInitDesignVars(x)

        orig_filter = filtr

        if max_iterations > 1:
            if step == max_iterations-1:
                optimization_options['tr_max_iterations'] = 15
        count += optimization_options['tr_max_iterations']

        optimization_options['output_file'] = os.path.join(prefix, 'output_file%d.dat'%(step))
        optimization_options['tr_output_file'] = os.path.join(prefix, 'tr_output_file%d.dat'%(step))

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
            dist_file = os.path.join(prefix, 'distance_solution%d.f5'%(step))

            # Compute the characteristic domain length
            lx = 10.0*args.AR # mm
            ly = 10.0 # mm
            lz = 10.0 # mm
            vol = lx*ly*lz
            domain_length = vol**(1.0/3.0)
            refine_distance = 0.025*domain_length
            TopOptUtils.approxDistanceRefine(forest, filtr, assembler, refine_distance,
                                             domain_length=domain_length,
                                             filename=dist_file)

        # Repartition the mesh
        forest.balance(1)
        forest.repartition()
