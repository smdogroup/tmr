from tmr import TMR, TopOptUtils
from tacs import TACS, elements, functions
from egads4py import egads
import numpy as np
import openmdao.api as om
import os

# Print colored text in terminal
try:
    from termcolor import colored
except:
    print("[Unavailable module] termcolor is not installed!")

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
        r0 = self.r0_frac*self.a
        mfilter = TopOptUtils.Mfilter(self.N, assemblers, filters, dim=3, r=r0)
        mfilter.initialize()
        return mfilter

class OutputCallback:
    def __init__(self, assembler, iter_offset=0):
        self.fig = None
        self.assembler = assembler
        self.xt = self.assembler.createDesignVec()

        # Set the output file name
        flag = (TACS.OUTPUT_CONNECTIVITY |
                TACS.OUTPUT_NODES |
                TACS.OUTPUT_EXTRAS)
        self.f5 = TACS.ToFH5(self.assembler, TACS.SOLID_ELEMENT, flag)
        self.iter_offset = iter_offset

        return

    def write_output(self, prefix, itr, oct_forest, quad_forest, x):

        self.f5.writeToFile(os.path.join(prefix, 'output%d.f5'%(itr + self.iter_offset)))

        self.assembler.getDesignVars(self.xt)
        TMR.writeSTLToBin(os.path.join(prefix, 'level_set_output%d.bstl'%(itr + self.iter_offset)),
                          oct_forest, self.xt)

        return

class FrequencyObj:
    """
    A class that evaluates the smallest eigenvalue, the objective is evaluated
    using an objective callback. We also add non-design mass to loaded nodes in
    order to form a well-posed frequency maximization problem
    """

    itr = 0

    def __init__(self, prefix, domain, forest, len0, AR, ratio, iter_offset,
                 eig_scale=1.0, con_scale=1.0, num_eigenvalues=10,
                 max_jd_size=100, max_gmres_size=30,
                 ksrho=50, non_design_mass=5.0):
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
        self.prefix = prefix
        self.domain = domain
        self.iter_offset = iter_offset
        self.lx = len0*AR
        self.ly = len0
        self.lz = len0
        if domain == 'lbracket':
            self.ly = len0*ratio
        self.ratio = ratio
        self.eig_scale = eig_scale
        self.con_scale = con_scale
        self.num_eigenvalues = num_eigenvalues
        self.max_jd_size = max_jd_size
        self.max_gmres_size = max_gmres_size
        self.ksrho = ksrho
        self.non_design_mass = non_design_mass

        self.fltr = None
        self.mg = None
        self.assembler = None
        self.comm = None
        self.oper = None
        self.jd = None
        self.kmat = None
        self.mmat = None
        self.temp = None
        self.eig = None
        self.eigv = None

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
            self.rho_original = self.assembler.createDesignVec()
            self.svec = self.assembler.createDesignVec()
            self.comm = self.assembler.getMPIComm()
            self.rank = self.comm.rank

            # Initialize space for matrices and vectors
            self.mmat = self.assembler.createMat()
            self.kmat = self.assembler.createMat()
            self.temp = self.assembler.createVec()
            self.eig = np.zeros(self.num_eigenvalues)
            self.eigv = []
            for i in range(self.num_eigenvalues):
                self.eigv.append(self.assembler.createVec())
            self.deig = []
            for i in range(self.num_eigenvalues):
                self.deig.append(self.assembler.createDesignVec())

            # Create the Jacobi-Davidson operator
            self.oper = TACS.JDFrequencyOperator(self.assembler, self.kmat, self.mmat,
                                                 self.mg.getMat(), self.mg)

            # Create the eigenvalue solver and set the number of recycling eigenvectors
            self.jd = TACS.JacobiDavidson(self.oper, self.num_eigenvalues,
                                          self.max_jd_size, self.max_gmres_size)
            self.jd.setTolerances(eig_rtol=1e-6, eig_atol=1e-8, rtol=1e-6, atol=1e-12)
            # self.jd.setTolerances(eig_rtol=1e-6, eig_atol=1e-8, rtol=1e-12, atol=1e-15)
            self.jd.setRecycle(self.num_eigenvalues)

            '''
            Create a non-design mass vector
            '''
            self.mvec = self.assembler.createDesignVec()
            mvals = self.mvec.getArray()

            # Get nodal locations
            Xpts = self.forest.getPoints()

            # Note: the local nodes are organized as follows:
            # |--- dependent nodes -- | ext_pre | -- owned local -- | - ext_post -|

            # Get number of local nodes in the current processor
            n_local_nodes = Xpts.shape[0]

            # Get numbder of dependent nodes
            _ptr, _conn, _weights = self.forest.getDepNodeConn()

            # Get number of ext_pre nodes
            n_ext_pre = self.forest.getExtPreOffset()

            # Get numbder of own nodes:
            offset = n_ext_pre

            # # Loop over all owned nodes and set non-design mass values
            tol = 1e-6
            if self.domain == 'cantilever':
                xmin = self.lx - tol
                xmax = self.lx + tol
                ymin = 0.25*self.ly - tol
                ymax = 0.75*self.ly + tol
                zmin = 0.0*self.lz - tol
                zmax = 0.2*self.lz + tol

            elif self.domain == 'michell':
                xmin = self.lx - tol
                xmax = self.lx + tol
                ymin = 0.25*self.ly - tol
                ymax = 0.75*self.ly + tol
                zmin = 0.4*self.lz - tol
                zmax = 0.6*self.lz + tol

            elif self.domain == 'mbb':
                xmin = 0.0*self.lx - tol
                xmax = 0.2*self.lx + tol
                ymin = 0.25*self.ly - tol
                ymax = 0.75*self.ly + tol
                zmin = self.lz - tol
                zmax = self.lz + tol

            elif self.domain == 'lbracket':
                RATIO = self.ratio
                xmin = self.lx - tol
                xmax = self.lx + tol
                ymin = 0.25*self.ly - tol
                ymax = 0.75*self.ly + tol
                zmin = 0.5*RATIO*self.lz - tol
                zmax = 1.0*RATIO*self.lz + tol

            else:
                raise ValueError("Unsupported domain type!")

            for i in range(offset, n_local_nodes):
                x, y, z = Xpts[i]
                if xmin < x < xmax:
                    if ymin < y < ymax:
                        if zmin < z < zmax:
                            mvals[i-offset] = self.non_design_mass

        # Assemble the stiffness matrix for the generalized eigenvalue problem
        self.assembler.assembleMatType(TACS.STIFFNESS_MATRIX, self.kmat)

        # Get current design variable
        dv = self.assembler.createDesignVec()
        self.assembler.getDesignVars(dv)

        # Update dv <- dv + mvec and update design variable
        dv.axpy(1.0, self.mvec)
        self.assembler.setDesignVars(dv)

        # Construct mass matrix
        self.assembler.assembleMatType(TACS.MASS_MATRIX, self.mmat)

        # Reset design variable
        dv.axpy(-1.0, self.mvec)
        self.assembler.setDesignVars(dv)

        '''
        Export non-design mass vectors to f5 file for verification
        '''
        # itr = FrequencyObj.itr + self.iter_offset
        # if itr % 10 == 0:
        #     # Set up flags for data output
        #     flag = (TACS.OUTPUT_CONNECTIVITY |
        #             TACS.OUTPUT_NODES |
        #             TACS.OUTPUT_EXTRAS |
        #             TACS.OUTPUT_DISPLACEMENTS)

        #     # Create f5 file writer
        #     f5 = TACS.ToFH5(self.assembler, TACS.SOLID_ELEMENT, flag)

        #     # Set non-design mass vector as design variable
        #     self.assembler.setDesignVars(self.mvec)

        #     # Set dM*e as state variable
        #     dMe = self.assembler.createVec()
        #     dMe_vals = dMe.getArray()
        #     dMe_vals[:] = 1.0
        #     self.mmat.mult(dMe, dMe)
        #     mmat_temp = self.assembler.createMat()
        #     self.assembler.assembleMatType(TACS.MASS_MATRIX, mmat_temp)
        #     e = self.assembler.createVec()
        #     e_vals = e.getArray()
        #     e_vals[:] = 1.0
        #     mmat_temp.mult(e, e)
        #     dMe.axpy(-1.0, e)
        #     sv = self.assembler.createVec()
        #     self.assembler.getVariables(sv)
        #     self.assembler.setVariables(dMe)

        #     f5.writeToFile(os.path.join(self.prefix, 'non-design-mass{:d}.f5'.format(
        #         self.iter_offset+FrequencyObj.itr)))

        #     # Set dv and sv back
        #     self.assembler.setDesignVars(dv)
        #     self.assembler.setVariables(sv)
        # FrequencyObj.itr += 1
        '''
        End export non-design mass vectors to f5 file for verification
        '''

        # Assemble the multigrid preconditioner
        self.mg.assembleMatType(TACS.STIFFNESS_MATRIX)
        self.mg.factor()

        # Check the whole mass
        e = self.assembler.createVec()
        t = self.assembler.createVec()
        evals = e.getArray()
        evals[:] = 1.0
        self.mmat.mult(e, t)
        eTMe = e.dot(t)
        if self.comm.rank == 0:
            print('[Mmat] eTMe = {:20.10e}'.format(eTMe))

        # Solve
        self.jd.solve(print_flag=True, print_level=0)

        # Check if succeeded, otherwise try again
        if self.jd.getNumConvergedEigenvalues() < self.num_eigenvalues:
            if self.comm.rank == 0:
                print("[Warning] Jacobi-Davidson failed to converge for the first run.")

            # Extract the eigenvalues
            for i in range(self.num_eigenvalues):
                self.eig[i], error = self.jd.extractEigenvalue(i)

            # Update preconditioner
            theta = 0.9*np.min(self.eig)
            self.mg.assembleMatCombo(TACS.STIFFNESS_MATRIX, 1.0, TACS.MASS_MATRIX, -theta)
            self.mg.factor()

            # Rerun the solver
            self.jd.solve(print_flag=True, print_level=1)
            nconvd = self.jd.getNumConvergedEigenvalues()

            # If it still fails, raise error and exit
            if nconvd < self.num_eigenvalues:
                msg = "No enough eigenvalues converged! ({:d}/{:d})".format(
                    nconvd, self.num_eigenvalues)
                raise ValueError(msg)

        # Extract eigenvalues and eigenvectors
        for i in range(self.num_eigenvalues):
            self.eig[i], error = self.jd.extractEigenvector(i, self.eigv[i])

        # Scale eigenvalues for a better KS approximation
        self.eig[:] *= self.eig_scale

        # Compute the minimal eigenvalue
        eig_min = np.min(self.eig)

        # Compute KS aggregation
        self.eta = np.exp(-self.ksrho*(self.eig - eig_min))
        self.beta = np.sum(self.eta)
        ks = (eig_min - np.log(self.beta)/self.ksrho)
        self.eta = self.eta/self.beta

        # Scale eigenvalue back
        self.eig[:] /= self.eig_scale

        # Objective
        obj = -ks

        # Print values
        if self.comm.rank == 0:
            print('{:30s}{:20.10e}'.format('[Obj] KS eigenvalue:', ks))
            print('{:30s}{:20.10e}'.format('[Obj] min eigenvalue:', eig_min))

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

        for i in range(self.num_eigenvalues):

            # This is a maximization problem
            scale = -self.eta[i]*self.eig_scale

            # Compute gradient of eigenvalue
            self.deig[i].zeroEntries()
            self.assembler.addMatDVSensInnerProduct(
                scale, TACS.STIFFNESS_MATRIX,
                self.eigv[i], self.eigv[i], self.deig[i])

            self.assembler.addMatDVSensInnerProduct(
                -scale*self.eig[i], TACS.MASS_MATRIX,
                self.eigv[i], self.eigv[i], self.deig[i])

            # Make sure the vector is properly distributed over all processors
            self.deig[i].beginSetValues(op=TACS.ADD_VALUES)
            self.deig[i].endSetValues(op=TACS.ADD_VALUES)

            dfdrho.axpy(1.0, self.deig[i])

        # Compute gradient norm
        norm = dfdrho.norm()
        if self.comm.rank == 0:
            print("{:30s}{:20.10e}".format('[Obj] gradient norm:', norm))
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

        Args:
            x (PVec): unfiltered design vector (not used because we get x directly from assembler)
            s (PVec): unfiltered update step
            y (PVec): y = Bs
            z, zw: dummy variable for qn correction for constraints, not used here.
        """

        update = self.assembler.createDesignVec()
        temp = self.assembler.createDesignVec()

        h = 1e-8

        # Get current nodal density
        rho = self.assembler.createDesignVec()
        self.assembler.getDesignVars(rho)

        self.rho_original.copyValues(rho)

        self.svec.zeroEntries()
        self.fltr.applyFilter(TMR.convertPVecToVec(s), self.svec)
        rho.axpy(h, self.svec)

        for i in range(self.num_eigenvalues):
            """
            Compute the first part using finite difference
            P += phi^T d2Kdx2 phi
            """

            # Zero out temp vector
            temp.zeroEntries()

            # Compute g(rho + h*s)
            self.assembler.setDesignVars(rho)
            self.assembler.addMatDVSensInnerProduct(self.eta[i], TACS.STIFFNESS_MATRIX,
                self.eigv[i], self.eigv[i], temp)

            # Compute g(rho)
            self.assembler.setDesignVars(self.rho_original)
            self.assembler.addMatDVSensInnerProduct(-self.eta[i], TACS.STIFFNESS_MATRIX,
                self.eigv[i], self.eigv[i], temp)

            # Distribute the vector
            temp.beginSetValues(op=TACS.ADD_VALUES)
            temp.endSetValues(op=TACS.ADD_VALUES)

            # Compute dg/h
            temp.scale(1/h)

            # Add to the update
            update.axpy(1.0, temp)

            """
            Compute the second part:
            P -= phi^T (deig*dMdx) phi^T
            """
            temp.zeroEntries()
            coeff = self.eta[i]*self.deig[i].dot(self.svec)
            self.assembler.addMatDVSensInnerProduct(-coeff,
                TACS.MASS_MATRIX, self.eigv[i], self.eigv[i], temp)

            temp.beginSetValues(op=TACS.ADD_VALUES)
            temp.endSetValues(op=TACS.ADD_VALUES)

            # Get update
            update.axpy(1.0, temp)

        # Compute curvature and check the norm of the update
        # to see if the magnitude makes sense
        x_norm = x.norm()
        s_norm = s.norm()
        y_norm = y.norm()
        dy_norm = update.norm()
        curvature = self.svec.dot(update)
        if self.comm.rank == 0:
            if curvature < 0:
                try:
                    print(colored("curvature: {:20.10e}".format(curvature), "red"))
                except:
                    print("curvature: {:20.10e}".format(curvature))
            else:
                try:
                    print(colored("curvature: {:20.10e}".format(curvature), "green"))
                except:
                    print("curvature: {:20.10e}".format(curvature))
            print("norm(x):   {:20.10e}".format(x_norm))
            print("norm(s):   {:20.10e}".format(s_norm))
            print("norm(y):   {:20.10e}".format(y_norm))
            print("norm(dy):  {:20.10e}".format(dy_norm))

        # Update y
        if curvature > 0:
            y_wrap = TMR.convertPVecToVec(y)
            # is it ok to have the same input and output? yes
            self.fltr.applyTranspose(update, update)
            y_wrap.axpy(1.0, update)

        return

class MassConstr:
    """
    Mass constraint takes the following form:
        m <= m_fixed

    or equivalently:
        c = -m/m_fixed + 1 >= 0
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
        mass_constr = -mass/self.m_fixed + 1.0
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

def cantilever_egads(comm, lx, ly, lz):
    '''
    Create egads model file
    '''
    prefix = './models'
    name = 'cantilever_{:.1f}_{:.1f}_{:.1f}.egads'.format(lx, ly, lz)

    if comm.rank == 0 and not os.path.isdir(prefix):
        os.mkdir(prefix)

    if os.path.isfile(os.path.join(prefix, name)):
        return

    # Create an EGADS context
    ctx = egads.context()

    # Create the domain geometry
    x0 = [0.0, 0.0, 0.0]
    x1 = [lx, ly, lz]
    b1 = ctx.makeSolidBody(egads.BOX, rdata=[x0, x1])
    m1 = ctx.makeTopology(egads.MODEL, children=[b1])
    if comm.rank == 0:
        m1.saveModel(os.path.join(prefix, name), overwrite=True)
    comm.Barrier()

    return

def cantilever_geo(comm, lx, ly, lz):

    prefix = './models'
    name = 'cantilever_{:.1f}_{:.1f}_{:.1f}.egads'.format(lx, ly, lz)

    try:
        geo = TMR.LoadModel(os.path.join(prefix, name), print_lev=0)
    except:
        cantilever_egads(comm, lx, ly, lz)
        geo = TMR.LoadModel(os.path.join(prefix, name), print_lev=0)

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

def lbracket_egads(comm, lx, ly, lz, ratio):

    prefix = './models'
    base_name = 'lbracket_base_{:.1f}_{:.1f}_{:.1f}_r{:.1f}.egads'.format(lx, ly, lz, ratio)
    arm1_name = 'lbracket_arm1_{:.1f}_{:.1f}_{:.1f}_r{:.1f}.egads'.format(lx, ly, lz, ratio)
    arm2_name = 'lbracket_arm2_{:.1f}_{:.1f}_{:.1f}_r{:.1f}.egads'.format(lx, ly, lz, ratio)

    if comm.rank == 0 and not os.path.isdir(prefix):
        os.mkdir(prefix)

    if os.path.isfile(os.path.join(prefix, base_name)) and \
       os.path.isfile(os.path.join(prefix, base_name)) and \
       os.path.isfile(os.path.join(prefix, base_name)):
       return

    RATIO = ratio

    # Create an EGADS context
    ctx = egads.context()

    # Create base
    x0 = [0.0, 0.0, 0.0]
    x1 = [lx*RATIO, ly, lz*RATIO]
    B0 = ctx.makeSolidBody(egads.BOX, rdata=[x0, x1])
    m1 = ctx.makeTopology(egads.MODEL, children=[B0])
    if comm.rank == 0:
        m1.saveModel(os.path.join(prefix, base_name), overwrite=True)
    comm.Barrier()

    # Create arm 1
    x0 = [lx*RATIO, 0.0, 0.0]
    x1 = [lx*(1-RATIO), ly, lz*RATIO]
    B1 = ctx.makeSolidBody(egads.BOX, rdata=[x0, x1])
    m2 = ctx.makeTopology(egads.MODEL, children=[B1])
    if comm.rank == 0:
        m2.saveModel(os.path.join(prefix, arm1_name), overwrite=True)
    comm.Barrier()

    # Create arm 2
    x0 = [0.0, 0.0, lz*RATIO]
    x1 = [lx*RATIO, ly, lz*(1-RATIO)]
    B2 = ctx.makeSolidBody(egads.BOX, rdata=[x0, x1])
    m3 = ctx.makeTopology(egads.MODEL, children=[B2])
    if comm.rank == 0:
        m3.saveModel(os.path.join(prefix, arm2_name), overwrite=True)
    comm.Barrier()

    return

def lbracket_geo(comm, lx, ly, lz, ratio):

    prefix = './models'
    base_name = 'lbracket_base_{:.1f}_{:.1f}_{:.1f}_r{:.1f}.egads'.format(lx, ly, lz, ratio)
    arm1_name = 'lbracket_arm1_{:.1f}_{:.1f}_{:.1f}_r{:.1f}.egads'.format(lx, ly, lz, ratio)
    arm2_name = 'lbracket_arm2_{:.1f}_{:.1f}_{:.1f}_r{:.1f}.egads'.format(lx, ly, lz, ratio)

    geos = []
    names = [base_name, arm1_name, arm2_name]

    for name in names:
        try:
            geos.append(TMR.LoadModel(os.path.join(prefix, name), print_lev=0))
        except:
            lbracket_egads(comm, lx, ly, lz, ratio)
            geos.append(TMR.LoadModel(os.path.join(prefix, name), print_lev=0))

    # Create the full list of vertices, edges, faces and volumes
    verts = []
    edges = []
    faces = []
    vols = []
    for geo in geos:
        verts.extend(geo.getVertices())
        edges.extend(geo.getEdges())
        faces.extend(geo.getFaces())
        vols.extend(geo.getVolumes())

    # Set all of the matching faces
    TMR.setMatchingFaces(geos)

    # Create the geometry
    geo = TMR.Model(verts, edges, faces, vols)

    return geo

def create_forest(comm, lx, ly, lz, ratio, htarget, depth, domain_type):
    """
    Create an initial forest for analysis and optimization

    This code loads in the model, sets names, meshes the geometry and creates
    a QuadForest from the mesh. The forest is populated with quadtrees with
    the specified depth.

    Args:
        comm (MPI_Comm): MPI communicator
        lx, ly, lz (float): sizes of each edge
        htarget (float): Target global element mesh size
        depth (int): Depth of the initial trees

    Returns:
        OctForest: Initial forest for topology optimization
    """

    if domain_type == 'cantilever' or domain_type == 'michell':
        # Create geo
        geo = cantilever_geo(comm, lx, ly, lz)

        # Mark the boundary condition faces
        verts = geo.getVertices()
        edges = geo.getEdges()
        faces = geo.getFaces()
        volumes = geo.getVolumes()
        faces[0].setName('fixed')

        # Set source and target faces
        faces[0].setSource(volumes[0], faces[1])

    elif domain_type == 'mbb':
        # Create geo
        geo = cantilever_geo(comm, lx, ly, lz)

        # Mark the boundary condition faces
        verts = geo.getVertices()
        edges = geo.getEdges()
        faces = geo.getFaces()
        volumes = geo.getVolumes()
        faces[0].setName('symmetry')
        edges[3].setName('support')

        # Set source and target faces
        faces[0].setSource(volumes[0], faces[1])

    elif domain_type == 'lbracket':
        # Create geo
        geo = lbracket_geo(comm, lx, ly, lz, ratio)

        # Mark the boundary condition faces
        verts = geo.getVertices()
        edges = geo.getEdges()
        faces = geo.getFaces()
        volumes = geo.getVolumes()
        faces[17].setName('fixed')

        # Set source and target faces
        faces[0].setSource(volumes[0], faces[1])
        faces[6].setSource(volumes[1], faces[7])
        faces[12].setSource(volumes[2], faces[13])

    else:
        raise ValueError("Unsupported domain type!")

    # # Print geometry element locations
    # if comm.rank == 0:
    #     for index, vert in enumerate(verts):
    #         x,y,z = vert.evalPoint()
    #         print("verts[{:d}]: ({:.2f}, {:.2f}, {:.2f})".format(index, x, y, z))

    #     for index, edge in enumerate(edges):
    #         x,y,z = edge.evalPoint(0.5)
    #         print("edge[{:d}]: ({:.2f}, {:.2f}, {:.2f})".format(index, x, y, z))

    #     for index, face in enumerate(faces):
    #         umin, vmin, umax, vmax = face.getRange()
    #         umid = 0.5*(umin + umax)
    #         vmid = 0.5*(vmin + vmax)
    #         x,y,z = face.evalPoint(umid, vmid)
    #         print("faces[{:d}]: ({:.2f}, {:.2f}, {:.2f})".format(index, x, y, z))

    #     for index, vol in enumerate(volumes):
    #         print("vol[{:d}]:".format(index), vol)
    # exit(0)

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

def create_problem(prefix, domain, forest, bcs, props, nlevels, vol_frac=0.25, r0_frac=0.05,
                   len0=1.0, AR=1.0, ratio=0.4, density=2600.0, iter_offset=0,
                   qn_correction=True, non_design_mass=5.0, eig_scale=1.0, eq_constr=False,
                   max_jd_size=100, max_gmres_size=30):
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

    # Create the problem and filter object
    N = 20
    mfilter = MFilterCreator(r0_frac, N, a=len0)
    filter_type = mfilter.filter_callback
    obj = CreatorCallback(bcs, props)
    problem = TopOptUtils.createTopoProblem(forest, obj.creator_callback,
                                            filter_type, use_galerkin=True,
                                            nlevels=nlevels)

    # Get the assembler object we just created
    assembler = problem.getAssembler()

    # Compute the fixed mass target
    lx = len0*AR # mm
    ly = len0 # mm
    lz = len0 # mm
    if domain == 'lbracket':
        ly = len0*ratio
    vol = lx*ly*lz
    if domain == 'lbracket':
        S1 = lx*lz
        S2 = lx*lz*(1-ratio)**2
        vol = (S1-S2)*ly
    m_fixed = vol_frac*(vol*density)

    # Add objective callback
    obj_callback = FrequencyObj(prefix, domain, forest, len0, AR, ratio, iter_offset,
                                non_design_mass=non_design_mass,
                                eig_scale=eig_scale,
                                max_jd_size=max_jd_size,
                                max_gmres_size=max_gmres_size)
    problem.addObjectiveCallback(obj_callback.objective,
                                 obj_callback.objective_gradient)

    # Add constraint callback
    constr_callback = MassConstr(m_fixed, assembler.getMPIComm())
    nineq = 1
    if eq_constr is True:
        nineq = 0
    problem.addConstraintCallback(1, nineq, constr_callback.constraint,
                                  constr_callback.constraint_gradient)

    # Use Quasi-Newton Update Correction if specified
    if qn_correction:
        problem.addQnCorrectionCallback(1, obj_callback.qn_correction)

    # Set output callback
    cb = OutputCallback(assembler, iter_offset=iter_offset)
    problem.setOutputCallback(cb.write_output)

    return problem, obj_callback

class OmAnalysis(om.ExplicitComponent):
    '''
    This class wraps the analyses with openmdao interface such that
    the optimization can be run with different optimizers such as
    SNOPT and IPOPT. Note that this does not support MPI yet.
    '''

    def __init__(self, problem, obj_callback):
        '''
        Args:
            problem (TMR.TopoProblem)
        '''

        super().__init__()
        self.problem = problem
        self.obj_callback = obj_callback

        self.problem.initialize()

        # Initialize ParOpt.PVec vectors for later use
        self.x_PVec = self.problem.createDesignVec()
        self.x_Vec = TMR.convertPVecToVec(self.x_PVec)
        self.x_vals = self.x_Vec.getArray()

        self.g_PVec = self.problem.createDesignVec()
        self.g_Vec = TMR.convertPVecToVec(self.g_PVec)
        self.g_vals = self.g_Vec.getArray()

        self.A_PVec = self.problem.createDesignVec()
        self.A_Vec = TMR.convertPVecToVec(self.A_PVec)
        self.A_vals = self.A_Vec.getArray()

        # Get size of design variable
        self.xlen = len(self.x_vals)

        # Initialize f5 converter
        self.assembler = self.problem.getAssembler()
        flag = (TACS.OUTPUT_CONNECTIVITY |
                TACS.OUTPUT_NODES |
                TACS.OUTPUT_EXTRAS)
        self.f5 = TACS.ToFH5(self.assembler, TACS.SOLID_ELEMENT, flag)

        return

    def setup(self):
        self.add_input('x', shape=(self.xlen,))

        self.add_output('obj')
        self.add_output('con')

        self.declare_partials(of='obj', wrt='x')
        self.declare_partials(of='con', wrt='x')

        return

    def compute(self, inputs, outputs):
        x = inputs['x']
        self.x_vals[:] = x[:]
        fail, fobj, cons = self.problem.evalObjCon(1, self.x_PVec)

        if fail:
            raise RuntimeError("Failed to evaluate objective and constraints!")
        else:
            outputs['obj'] = fobj
            outputs['con'] = cons[0]

        return

    def compute_partials(self, inputs, partials):
        x = inputs['x']
        self.x_vals[:] = x[:]
        fail = self.problem.evalObjConGradient(self.x_PVec, self.g_PVec, [self.A_PVec])

        if fail:
            raise RuntimeError("Failed to evaluate objective and constraints!")
        else:
            partials['obj', 'x'] = self.g_vals
            partials['con', 'x'] = self.A_vals

        return

    def qn_correction(self, x, z, zw, s, y):
        x_PVec = self.problem.createDesignVec()
        s_PVec = self.problem.createDesignVec()
        y_PVec = self.problem.createDesignVec()

        x_vals = TMR.convertPVecToVec(x_PVec).getArray()
        s_vals = TMR.convertPVecToVec(s_PVec).getArray()
        y_vals = TMR.convertPVecToVec(y_PVec).getArray()

        x_vals[:] = x[:]
        s_vals[:] = s[:]
        y_vals[:] = y[:]

        self.obj_callback.qn_correction(x_PVec, z, zw, s_PVec, y_PVec)

        x[:] = x_vals[:]
        s[:] = s_vals[:]
        y[:] = y_vals[:]

        return

    def write_output(self, prefix, refine_step):
        """
        Args:
            x (TACS.Vec)
        """
        self.f5.writeToFile(os.path.join(prefix, 'output_refine{:d}.f5'.format(refine_step)))
        return
