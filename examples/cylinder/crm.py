from __future__ import print_function
from mpi4py import MPI
from tmr import TMR
from tacs import TACS, elements, constitutive, functions
from paropt import ParOpt
import numpy as np
import argparse
import os
import ksFSDT

# Import pyoptsparse
from pyoptsparse import Optimization, OPT
from scipy import sparse

# Guess the spar, rib and skin indices
def get_face_indices(faces):
    ucrm_skins = []
    ucrm_spars = []
    ucrm_ribs = []

    # Find the edge loops
    for index, face in enumerate(faces):
        nx = np.zeros(3)

        for k in range(face.getNumEdgeLoops()):
            loop = face.getEdgeLoop(k)

            elist, dirs = loop.getEdgeLoop()
            pts = np.zeros((len(elist)+1, 3))
            tx = np.zeros((len(elist)+1, 3))
            for i, e in enumerate(elist):
                v1, v2 = e.getVertices()
                if dirs[i] > 0:
                    pt1 = v1.evalPoint()
                    pt2 = v2.evalPoint()
                else:
                    pt1 = v2.evalPoint()
                    pt2 = v1.evalPoint()
                if i == 0:
                    pts[0,:] = pt1[:]
                pts[i+1,:] = pt2[:]

            for i in range(len(elist)):
                tx[i,:] = pts[i+1,:] - pts[i,:]
            tx[-1,:] = pts[1,:] - pts[0,:]

            # Find the average normal direction
            for i in range(len(elist)):
                nx += np.cross(tx[i], tx[i+1])
            nx /= np.sqrt(np.dot(nx, nx))

            # Quit after the first face loop
            break

        # Normal directions
        theta = 30*np.pi/180.0
        face_dir = np.array([0.0, 0.0, 1.0])
        rib_dir = np.array([np.sin(theta), np.cos(theta), 0.0])
        spar_dir = np.array([-np.cos(theta), np.sin(theta), 0.0])

        fd = np.fabs(np.dot(face_dir, nx))
        rd = np.fabs(np.dot(rib_dir, nx))
        sd = np.fabs(np.dot(spar_dir, nx))

        if fd > rd and fd > sd:
            ucrm_skins.append(index)
        elif rd > fd and rd > sd:
            ucrm_ribs.append(index)
        else:
            ucrm_spars.append(index)

    return ucrm_skins, ucrm_spars, ucrm_ribs

def uniquify(seq):
    set = {}
    for e in seq:
        set[e] = 1
    return sorted(set.keys())

def get_face_to_face_index(faces, edges):
    # Set the face indices
    face_index = {}
    for i, f in enumerate(faces):
        face_index[f.getName()] = i

    # Set the edge to face index
    edge_to_face = {}
    for e in edges:
        edge_to_face[e.getName()] = {}

    # Create an edge to face dictionary
    for i, face in enumerate(faces):
        # Get the edge loop
        loop = face.getEdgeLoop(0)
        elist, dirs = loop.getEdgeLoop()

        for e in elist:
            name = e.getName()
            if name in edge_to_face:
                edge_to_face[name][face.getName()] = True

    # Now loop over the edges
    face_to_face = []
    for i, face in enumerate(faces):
        loop = face.getEdgeLoop(0)
        elist, dirs = loop.getEdgeLoop()

        flist = []
        for e in elist:
            ename = e.getName()
            if ename in edge_to_face:
                for name in edge_to_face[e.getName()]:
                    if i != face_index[name]:
                        flist.append(face_index[name])

        # Sort the edge list
        face_to_face.append(uniquify(flist))

    return face_to_face

def get_constraint_pairs(face_to_face, face_set):
    pairs = []
    for index in face_set:
        for adj in face_to_face[index]:
            if adj in face_set:
                pairs.append((index, adj))
    return pairs

class uCRM_VonMisesMassMin:
    '''
    Mass minimization with a von Mises stress constraint
    '''

    def __init__(self, comm, num_components):
        # Set the communicator
        self.comm = comm

        # Scale the mass objective so that it is O(10)
        self.mass_scale = 1e-2

        # Scale the thickness variables so that they are measured in
        # 1/10-ths of inches
        self.thickness_scale = 10.0

        # The number of thickness variables in the problem
        self.nvars = num_components

        # The number of constraints (1 global stress constraint that
        # will use the KS function)
        self.ncon = 1

        # Set the iteration counter
        self.iter_count = 0

        return

    def setAssembler(self, assembler, mg, gmres, ksfunc):

        # Create tacs assembler object from mesh loader
        self.assembler = assembler

        # Create the list of functions
        self.funcs = [functions.StructuralMass(self.assembler), ksfunc]

        # Set up the solver
        self.ans = self.assembler.createVec()
        self.res = self.assembler.createVec()
        self.adjoint = self.assembler.createVec()
        self.dfdu = self.assembler.createVec()

        self.use_multigrid = False
        if self.use_multigrid:
            self.mg = mg
            self.mat = self.mg.getMat()
            self.gmres = gmres
        else:
            self.mat = self.assembler.createFEMat()
            self.pc = TACS.Pc(self.mat)
            self.gmres = TACS.KSM(self.mat, self.pc, 10)

        # For visualization
        flag = (TACS.ToFH5.NODES |
                TACS.ToFH5.DISPLACEMENTS |
                TACS.ToFH5.STRAINS |
                TACS.ToFH5.EXTRAS)
        self.f5 = TACS.ToFH5(self.assembler, TACS.PY_SHELL, flag)

        return

    def getVarsAndBounds(self, x, lb, ub):
        '''Set the values of the bounds'''
        xvals = np.zeros(self.nvars, TACS.dtype)
        self.assembler.getDesignVars(xvals)
        xvals[:] = self.thickness_scale*xvals

        xlb = 1e20*np.ones(self.nvars, TACS.dtype)
        xub = -1e20*np.ones(self.nvars, TACS.dtype)
        self.assembler.getDesignVarRange(xlb, xub)
        xlb[:] = self.thickness_scale*xlb
        xub[:] = self.thickness_scale*xub

        self.comm.Allreduce(xvals, x, op=MPI.MAX)
        self.comm.Allreduce(xlb, lb, op=MPI.MIN)
        self.comm.Allreduce(xub, ub, op=MPI.MAX)

        return

    def objcon(self, xdict):
        '''Evaluate the objective and constraint'''

        # Extract the values of x
        x = xdict['x']

        # Evaluate the objective and constraints
        fail = 0
        con = np.zeros(1)

        # Set the new design variable values
        self.assembler.setDesignVars(x[:]/self.thickness_scale)

        # Assemble the Jacobian and factor the matrix
        alpha = 1.0
        beta = 0.0
        gamma = 0.0
        self.assembler.zeroVariables()
        if self.use_multigrid:
            self.mg.assembleJacobian(alpha, beta, gamma, self.res)
            self.mg.factor()
        else:
            self.assembler.assembleJacobian(alpha, beta, gamma, self.res, self.mat)
            self.pc.factor()

        # Solve the linear system and set the varaibles into TACS
        self.gmres.solve(self.res, self.ans)
        self.ans.scale(-1.0)
        self.assembler.setVariables(self.ans)

        # Evaluate the function
        fvals = self.assembler.evalFunctions(self.funcs)

        # Set the mass as the objective
        fobj = self.mass_scale*fvals[0]

        # Set the KS function (the approximate maximum ratio of the
        # von Mises stress to the design stress) so that it is less
        # than or equal to 1.0
        con[0] = 1.0 - fvals[1] # ~= 1.0 - max (sigma/design) >= 0

        # Create the dictionary of functions
        funcs = {'objective': fobj, 'con': con}

        return funcs, fail

    def gobjcon(self, xdict, funcs):
        '''Evaluate the objective and constraint gradient'''
        fail = 0

        # Extract the values of x
        x = xdict['x']

        # Evaluate the derivative of the mass and place it in the
        # objective gradient
        gx = np.zeros(self.nvars, TACS.dtype)
        self.assembler.evalDVSens(self.funcs[0], gx)
        gx[:] *= self.mass_scale/self.thickness_scale

        # Compute the total derivative w.r.t. material design variables
        dfdx = np.zeros(self.nvars, TACS.dtype)
        product = np.zeros(self.nvars, TACS.dtype)

        # Compute the derivative of the function w.r.t. the state
        # variables
        self.assembler.evalDVSens(self.funcs[1], dfdx)
        self.assembler.evalSVSens(self.funcs[1], self.dfdu)
        self.gmres.solve(self.dfdu, self.adjoint)

        # Compute the product of the adjoint with the derivative of the
        # residuals
        self.assembler.evalAdjointResProduct(self.adjoint, product)

        # Set the constraint gradient
        dfdx[:] = -(dfdx - product)/self.thickness_scale

        # Write out the solution file every 10 iterations
        if self.iter_count % 10 == 0:
            self.f5.writeToFile('results/ucrm_iter%d.f5'%(self.iter_count))
        self.iter_count += 1

        # Create the sensitivity dictionary
        sens = {'objective':{'x': gx}, 'con':{'x': dfdx}}

        return sens, fail

class CreateMe(TMR.QuadCreator):
    def __init__(self, bcs, topo, elem_dict):
        TMR.QuadCreator.__init__(bcs)
        self.topo = topo
        self.elem_dict = elem_dict
        return

    def createElement(self, order, quad):
        '''Create the element'''
        # Get the model name and set the face
        attr = topo.getFace(quad.face).getName()
        return self.elem_dict[order-2][attr]

def addFaceTraction(order, assembler):
    # Create the surface traction
    aux = TACS.AuxElements()

    # Get the element node locations
    nelems = assembler.getNumElements()

    # Loop over the nodes and create the traction forces in the x/y/z
    # directions
    nnodes = order*order
    tx = np.zeros(nnodes, dtype=TACS.dtype)
    ty = np.zeros(nnodes, dtype=TACS.dtype)
    tz = np.zeros(nnodes, dtype=TACS.dtype)
    tz[:] = 5.0

    # Create the shell traction
    trac = elements.ShellTraction(order, tx, ty, tz)
    for i in range(nelems):
        aux.addElement(i, trac)

    return aux

def createRefined(topo, elem_dict, forest, bcs, pttype=TMR.UNIFORM_POINTS):
    new_forest = forest.duplicate()
    new_forest.setMeshOrder(forest.getMeshOrder()+1, pttype)
    creator = CreateMe(bcs, topo, elem_dict)
    return new_forest, creator.createTACS(new_forest)

def createProblem(topo, elem_dict, forest, bcs, ordering,
                  order=2, nlevels=2, pttype=TMR.UNIFORM_POINTS):
    # Create the forest
    forests = []
    assemblers = []

    # Create the trees, rebalance the elements and repartition
    forest.balance(1)
    forest.setMeshOrder(order, pttype)
    forest.repartition()
    forests.append(forest)

    # Make the creator class
    creator = CreateMe(bcs, topo, elem_dict)
    assemblers.append(creator.createTACS(forest, ordering))

    while order > 2:
        order = order-1
        forest = forests[-1].duplicate()
        forest.setMeshOrder(order, pttype)
        forest.balance(1)
        forests.append(forest)

        # Make the creator class
        creator = CreateMe(bcs, topo, elem_dict)
        assemblers.append(creator.createTACS(forest, ordering))

    for i in range(nlevels-1):
        forest = forests[-1].coarsen()
        forest.setMeshOrder(2, pttype)
        forest.balance(1)
        forests.append(forest)

        # Make the creator class
        creator = CreateMe(bcs, topo, elem_dict)
        assemblers.append(creator.createTACS(forest, ordering))

    # Create the multigrid object
    mg = TMR.createMg(assemblers, forests, omega=0.5)

    return assemblers[0], mg

# Set the communicator
comm = MPI.COMM_WORLD

# Create an argument parser to read in arguments from the commnad line
p = argparse.ArgumentParser()
p.add_argument('--steps', type=int, default=5)
p.add_argument('--htarget', type=float, default=4.0)
p.add_argument('--order', type=int, default=2)
p.add_argument('--ordering', type=str, default='multicolor')
p.add_argument('--ksweight', type=float, default=10.0)
p.add_argument('--uniform_refinement', action='store_true', default=False)
p.add_argument('--structured', action='store_true', default=False)
p.add_argument('--energy_error', action='store_true', default=False)
p.add_argument('--compute_solution_error', action='store_true', default=False)
p.add_argument('--exact_refined_adjoint', action='store_true', default=False)
p.add_argument('--remesh_domain', action='store_true', default=False)
p.add_argument('--remesh_strategy', type=str, default='fraction')
p.add_argument('--optimizer', type=str, default='snopt')
args = p.parse_args()

# Set the KS parameter
ksweight = args.ksweight

# This flag indicates whether to solve the adjoint exactly on the
# next-finest mesh or not
exact_refined_adjoint = args.exact_refined_adjoint

# Set the number of AMR steps to use
steps = args.steps

# Set the order of the mesh
order = args.order

# Set the type of ordering to use for this problem
ordering = args.ordering
ordering = ordering.lower()

tol = 1e-4
fname = 'results/crm_opt.out'
options = {}
if args.optimizer == 'snopt':
    options['Print file'] = fname
    options['Summary file'] = fname + '_summary'
    options['Major optimality tolerance'] = tol
    options['Major iterations limit'] = 50
elif optimizer == 'ipopt':
    options['print_user_options'] = 'yes'
    options['tol'] = tol
    options['nlp_scaling_method'] = 'none'
    options['limited_memory_max_history'] = 25
    options['bound_relax_factor'] = 0.0
    options['linear_solver'] = 'ma27'
    options['output_file'] = fname
    options['max_iter'] = 10000

# Create the CRM wingbox model and set the names
geo = TMR.LoadModel('ucrm_9_model.step')
verts = geo.getVertices()
edges = geo.getEdges()
faces = geo.getFaces()
geo = TMR.Model(verts, edges, faces)

# Set the edge names
for i, e in enumerate(edges):
    attr = 'E%d'%(i)
    e.setName(attr)

# Set the face names and create the constitutive objects
elem_dict = [{}, {}, {}, {}]
for i, f in enumerate(faces):
    attr = 'F%d'%(i)
    f.setName(attr)

# The root rib for boundary conditions
ucrm_root_edges = [200, 380, 270, 220]

ucrm_skins, ucrm_spars, ucrm_ribs = get_face_indices(faces)

# Get the face to face adjacency data structure
face_to_face = get_face_to_face_index(faces, edges)

# Create the constraint pairs for each group
pairs = get_constraint_pairs(face_to_face, ucrm_skins)
pairs.extend(get_constraint_pairs(face_to_face, ucrm_spars))
pairs.extend(get_constraint_pairs(face_to_face, ucrm_ribs))

# Set the number of design variables
num_design_vars = len(faces)

# Set the material properties
rho = 97.5e-3 # 0.0975 lb/in^3
E = 10000e3 # 10,000 ksi: Young's modulus
nu = 0.3 # Poisson ratio
kcorr = 5.0/6.0 # Shear correction factor
ys = 20e3 # psi
thickness = 1.0
min_thickness = 0.125 # 1/8th of an inch
max_thickness = 100.0

# Set the different directions
skin_spar_dir = np.array([np.sin(30*np.pi/180.0), np.cos(30*np.pi/180), 0.0])
rib_dir = np.array([1.0, 0.0, 0.0])

# Set the face names and create the constitutive objects
elem_dict = [{}, {}, {}, {}]
for i, f in enumerate(faces):
    attr = 'F%d'%(i)
    stiff = ksFSDT.ksFSDT(ksweight, rho, E, nu, kcorr, ys, thickness,
                          i, min_thickness, max_thickness)

    # Set the reference direction
    if i in ucrm_ribs:
        stiff.setRefAxis(rib_dir)
    else:
        stiff.setRefAxis(skin_spar_dir)

    # Set the component number for visualization purposes
    comp = 0
    if i in ucrm_skins:
        comp = 1
    elif i in ucrm_ribs:
        comp = 2

    # Create the elements of different orders
    for j in range(4):
        elem_dict[j][attr] = elements.MITCShell(j+2, stiff,
                                                component_num=comp)

# Initial target mesh spacing
htarget = args.htarget

# Create the new mesh
mesh = TMR.Mesh(comm, geo)

# Set the meshing options
opts = TMR.MeshOptions()

# Set the mesh type
# opts.mesh_type_default = TMR.TRIANGLE
opts.num_smoothing_steps = 10
opts.write_mesh_quality_histogram = 1
opts.triangularize_print_iter = 50000

# Create the surface mesh
mesh.mesh(htarget, opts)
mesh.writeToVTK('results/mesh.vtk')

# The boundary condition object
bcs = TMR.BoundaryConditions()
for edge_num in ucrm_root_edges:
    bcs.addBoundaryCondition('E%d'%(edge_num), [0, 1, 2, 3, 4, 5])

# Set the feature size object
feature_size = None

# Create the corresponding mesh topology from the mesh-model
model = mesh.createModelFromMesh()
topo = TMR.Topology(comm, model)

# Create the optimization problem
opt_problem = uCRM_VonMisesMassMin(comm, num_design_vars)

# Create the quad forest and set the topology of the forest
depth = 0
# if order == 2:
#     depth = 1
forest = TMR.QuadForest(comm)
forest.setTopology(topo)
forest.setMeshOrder(order, TMR.UNIFORM_POINTS)
forest.createTrees(depth)

# Set the ordering to use
if ordering == 'rcm':
    ordering = TACS.PY_RCM_ORDER
elif ordering == 'multicolor':
    ordering = TACS.PY_MULTICOLOR_ORDER
else:
    ordering = TACS.PY_NATURAL_ORDER

# Null pointer to the optimizer
opt = None

for k in range(steps):
    # Create the topology problem
    if args.remesh_domain:
        if order == 2:
            nlevs = 2
        else:
            nlevs = 1
    else:
        nlevs = min(5, depth+k+1)

    # Create the assembler object
    assembler, mg = createProblem(topo, elem_dict, forest, bcs, ordering,
                                  order=order, nlevels=nlevs)
    aux = addFaceTraction(order, assembler)
    assembler.setAuxElements(aux)

    # Create the KS functional
    func = functions.KSFailure(assembler, ksweight)
    func.setKSFailureType('continuous')

    # Create the GMRES object
    gmres = TACS.KSM(mg.getMat(), mg, 100, isFlexible=1)
    gmres.setMonitor(comm, freq=10)
    gmres.setTolerances(1e-10, 1e-30)

    # Set the new assembler object
    opt_problem.setAssembler(assembler, mg, gmres, func)

    # if k == 0:
    #     # Set up the optimization problem
    #     max_lbfgs = 5
    #     opt = ParOpt.pyParOpt(opt_problem, max_lbfgs, ParOpt.BFGS)
    #     opt.setOutputFile('results/crm_opt.out')

    #     # Set the optimality tolerance
    #     opt.setAbsOptimalityTol(1e-4)

    #     # Set optimization parameters
    #     opt.setArmijoParam(1e-5)

    #     # Get the optimized point
    #     x, z, zw, zl, zu = opt.getOptimizedPoint()

    #     # Set the starting point strategy
    #     opt.setStartingPointStrategy(ParOpt.AFFINE_STEP)

    #     # Set the max oiterations
    #     opt.setMaxMajorIterations(5)

    #     # Set the output level to understand what is going on
    #     opt.setOutputLevel(2)
    # else:
    #     beta = 1e-4
    #     opt.resetDesignAndBounds()
    #     opt.setStartAffineStepMultiplierMin(beta)

    # # Optimize the new point
    # opt.optimize()

    # Create the optimization problem
    prob = Optimization('topo', opt_problem.objcon)

    # Add the variable group
    n = opt_problem.nvars
    x0 = np.zeros(n)
    lb = np.zeros(n)
    ub = np.zeros(n)
    opt_problem.getVarsAndBounds(x0, lb, ub)
    prob.addVarGroup('x', n, value=x0, lower=lb, upper=ub)

    # Add the constraints
    prob.addConGroup('con', opt_problem.ncon, lower=0.0, upper=None)

    # Add the objective
    prob.addObj('objective')

    # Set up the sparse constraint Jacobian
    rowp = [0]
    cols = []
    data = []
    for i, pair in enumerate(pairs):
        data.extend([1.0, -1.0])
        cols.extend([pair[0], pair[1]])
        rowp.append(rowp[i] + 2)
    jacobian = {'csr':[rowp, cols, data], 'shape':[len(pairs), n]}

    # Set the linear constraint Jacobian
    tscale = opt_problem.thickness_scale
    prob.addConGroup('lincon', len(pairs),
                     lower=-0.05*tscale, upper=0.05*tscale,
                     wrt='x', jac={'x': jacobian}, linear=True)

    # Set the output file name
    fname = 'results/crm_opt%02d.out'%(k)
    if args.optimizer == 'snopt':
        options['Print file'] = fname
        options['Summary file'] = fname + '_summary'
    elif args.optimizer == 'ipopt':
        options['output_file'] = fname

    # Create the optimizer and optimize it!
    opt = OPT(args.optimizer, options=options)
    sol = opt(prob, sens=opt_problem.gobjcon)

    # Create and compute the function
    fval = assembler.evalFunctions([func])[0]

    # Create the refined mesh
    if exact_refined_adjoint:
        forest_refined = forest.duplicate()
        assembler_refined, mg = createProblem(topo, elem_dict,
                                              forest_refined, bcs, ordering,
                                              order=order+1, nlevels=nlevs+1)
    else:
        forest_refined, assembler_refined = createRefined(topo, elem_dict,
                                                          forest_refined, bcs)
    aux = addFaceTraction(order+1, assembler_refined)
    assembler_refined.setAuxElements(aux)

    # Set the design variables in the refined assembler object
    xvals = np.zeros(opt_problem.nvars, TACS.dtype)
    assembler.getDesignVars(xvals)
    assembler_refined.setDesignVars(xvals)

    # Extract the answer
    ans = opt_problem.ans

    if args.energy_error:
        # Compute the strain energy error estimate
        fval_corr = 0.0
        adjoint_corr = 0.0
        err_est, error = TMR.strainEnergyError(forest, assembler,
            forest_refined, assembler_refined)

        TMR.computeReconSolution(forest, assembler,
            forest_refined, assembler_refined)
    else:
        if exact_refined_adjoint:
            # Compute the reconstructed solution on the refined mesh
            ans_interp = assembler_refined.createVec()
            TMR.computeInterpSolution(forest, assembler,
                forest_refined, assembler_refined, ans, ans_interp)

            # Set the interpolated solution on the fine mesh
            assembler_refined.setVariables(ans_interp)

            # Assemble the Jacobian matrix on the refined mesh
            res_refined = assembler_refined.createVec()
            mg.assembleJacobian(1.0, 0.0, 0.0, res_refined)
            mg.factor()
            pc = mg
            mat = mg.getMat()

            # Compute the functional and the right-hand-side for the
            # adjoint on the refined mesh
            adjoint_rhs = assembler_refined.createVec()
            func_refined = functions.KSFailure(assembler_refined, ksweight)
            func_refined.setKSFailureType('continuous')

            # Evaluate the functional on the refined mesh
            fval_refined = assembler_refined.evalFunctions([func_refined])[0]
            assembler_refined.evalSVSens(func_refined, adjoint_rhs)

            # Create the GMRES object on the fine mesh
            gmres = TACS.KSM(mat, pc, 100, isFlexible=1)
            gmres.setMonitor(comm, freq=10)
            gmres.setTolerances(1e-14, 1e-30)

            # Solve the linear system
            adjoint_refined = assembler_refined.createVec()
            gmres.solve(adjoint_rhs, adjoint_refined)
            adjoint_refined.scale(-1.0)

            # Compute the adjoint correction on the fine mesh
            adjoint_corr = adjoint_refined.dot(res_refined)

            # Compute the reconstructed adjoint solution on the refined mesh
            adjoint = assembler.createVec()
            adjoint_interp = assembler_refined.createVec()
            TMR.computeInterpSolution(forest_refined, assembler_refined,
                forest, assembler, adjoint_refined, adjoint)
            TMR.computeInterpSolution(forest, assembler,
                forest_refined, assembler_refined, adjoint, adjoint_interp)
            adjoint_refined.axpy(-1.0, adjoint_interp)

            err_est, __, error = TMR.adjointError(forest, assembler,
                forest_refined, assembler_refined, ans_interp, adjoint_refined)
        else:
            # Compute the adjoint on the original mesh
            res.zeroEntries()
            assembler.evalSVSens(func, res)
            adjoint = assembler.createVec()
            gmres.solve(res, adjoint)
            adjoint.scale(-1.0)

            # Compute the solution on the refined mesh
            ans_refined = assembler_refined.createVec()
            TMR.computeReconSolution(forest, assembler,
                forest_refined, assembler_refined, ans, ans_refined)

            # Apply Dirichlet boundary conditions
            assembler_refined.setVariables(ans_refined)

            # Compute the functional on the refined mesh
            func_refined = functions.KSFailure(assembler_refined, ksweight)
            func_refined.setKSFailureType('continuous')

            # Evaluate the functional on the refined mesh
            fval_refined = assembler_refined.evalFunctions([func_refined])[0]

            # Approximate the difference between the refined adjoint
            # and the adjoint on the current mesh
            adjoint_refined = assembler_refined.createVec()
            TMR.computeReconSolution(forest, assembler,
                forest_refined, assembler_refined, adjoint,
                adjoint_refined)

            # Compute the adjoint and use adjoint-based refinement
            err_est, adjoint_corr, error = TMR.adjointError(forest, assembler,
                forest_refined, assembler_refined, ans_refined, adjoint_refined)

            # Compute the adjoint and use adjoint-based refinement
            err_est, __, error = TMR.adjointError(forest, assembler,
                forest_refined, assembler_refined, ans_refined, adjoint_refined)

            TMR.computeReconSolution(forest, assembler,
                forest_refined, assembler_refined, adjoint, adjoint_refined)
            assembler_refined.setVariables(adjoint_refined)

        # Compute the refined function value
        fval_corr = fval_refined + adjoint_corr

    flag = (TACS.ToFH5.NODES |
            TACS.ToFH5.DISPLACEMENTS |
            TACS.ToFH5.STRAINS |
            TACS.ToFH5.EXTRAS)
    f5_refine = TACS.ToFH5(assembler_refined, TACS.PY_SHELL, flag)
    f5_refine.writeToFile('results/solution_refined%02d.f5'%(k))

    # Compute the refinement from the error estimate
    low = -16
    high = 4
    bins_per_decade = 10
    nbins = bins_per_decade*(high - low)
    bounds = 10**np.linspace(high, low, nbins+1)
    bins = np.zeros(nbins+2, dtype=np.int)

    # Compute the mean and standard deviations of the log(error)
    ntotal = comm.allreduce(assembler.getNumElements(), op=MPI.SUM)
    mean = comm.allreduce(np.sum(np.log(error)), op=MPI.SUM)
    mean /= ntotal

    # Compute the standard deviation
    stddev = comm.allreduce(np.sum((np.log(error) - mean)**2), op=MPI.SUM)
    stddev = np.sqrt(stddev/(ntotal-1))

    # Get the total number of nodes
    nnodes = comm.allreduce(assembler.getNumOwnedNodes(), op=MPI.SUM)

    # Compute the bins
    for i in range(len(error)):
        if error[i] > bounds[0]:
            bins[0] += 1
        elif error[i] < bounds[-1]:
            bins[-1] += 1
        else:
            for j in range(len(bounds)-1):
                if (error[i] <= bounds[j] and
                    error[i] > bounds[j+1]):
                    bins[j+1] += 1

    # Compute the number of bins
    bins = comm.allreduce(bins, MPI.SUM)

    # Compute the sum of the bins
    total = np.sum(bins)

    # Print out the result
    if comm.rank == 0:
        print('fval      = ', fval)
        print('fval corr = ', fval_corr)
        print('estimate  = ', err_est)
        print('mean      = ', mean)
        print('stddev    = ', stddev)

        # Set the data
        data = np.zeros((nbins, 4))
        for i in range(nbins-1, -1, -1):
            data[i,0] = bounds[i]
            data[i,1] = bounds[i+1]
            data[i,2] = bins[i+1]
            data[i,3] = 100.0*bins[i+1]/total
        np.savetxt('results/crm_data%d.txt'%(k), data)

    # Perform the refinement
    if args.uniform_refinement:
        forest.refine()
    elif k < steps-1:
        if args.remesh_domain:
            # Ensure that we're using an unstructured mesh
            opts.mesh_type_default = TMR.UNSTRUCTURED

            # Find the positions of the center points of each node
            nelems = assembler.getNumElements()

            # Allocate the positions
            Xp = np.zeros((nelems, 3))
            for i in range(nelems):
                # Get the information about the given element
                elem, Xpt, vrs, dvars, ddvars = assembler.getElementData(i)

                # Get the approximate element centroid
                Xp[i,:] = np.average(Xpt.reshape((-1, 3)), axis=0)

            # Prepare to collect things to the root processor (only
            # one where it is required)
            root = 0

            # Get the element counts
            if comm.rank == root:
                size = error.shape[0]
                count = comm.gather(size, root=root)
                count = np.array(count, dtype=np.int)
                ntotal = np.sum(count)

                errors = np.zeros(np.sum(count))
                Xpt = np.zeros(3*np.sum(count))
                comm.Gatherv(error, [errors, count])
                comm.Gatherv(Xp.flatten(), [Xpt, 3*count])

                # Reshape the point array
                Xpt = Xpt.reshape(-1,3)

                # Asymptotic order of accuracy on per-element basis
                s = 1.0
                if args.order >= 3:
                    s = args.order - 1.0

                # Dimension of the problem
                d = 2.0

                # Set the exponent
                exponent = d/(d + s)

                # Compute the target error as a fixed fraction of the
                # error estimate. This will result in the size of the
                # mesh increasing at each iteration.
                if args.remesh_strategy == 'fraction':
                    err_target = 0.1*err_est
                else:
                    # Set a fixed target error
                    err_target = 1e-4   # Target error estimate

                # Set the error estimate
                cval = 1.0
                if args.remesh_strategy == 'fixed_mesh':
                    # Set the target error
                    count_target = 20e3 # Target element count

                    # Compute the constant for element count
                    cval = (count_target)**(-1.0/d)
                    cval *= (np.sum(errors**(exponent)))**(1.0/d)
                else:
                    # Compute the constant for target error
                    cval = (err_target/np.sum(errors**(exponent)))**(1.0/s)

                # Compute the element-wise target error
                hvals = cval*errors**(-(1.0/(d + s)))

                # Set the values of
                if feature_size is not None:
                    for i, hp in enumerate(hvals):
                        hlocal = feature_size.getFeatureSize(Xpt[i,:])
                        hvals[i] = np.min(
                            (np.max((hp*hlocal, 0.25*hlocal)), 2*hlocal))
                else:
                    for i, hp in enumerate(hvals):
                        hvals[i] = np.min(
                            (np.max((hp*htarget, 0.25*htarget)), 2*htarget))

                # Allocate the feature size object
                hmax = args.htarget
                hmin = 0.05*args.htarget
                feature_size = TMR.PointFeatureSize(Xpt, hvals, hmin, hmax)
            else:
                size = error.shape[0]
                comm.gather(size, root=root)
                comm.Gatherv(error, None)
                comm.Gatherv(Xp.flatten(), None)

                # Create a dummy feature size object...
                feature_size = TMR.ConstElementSize(0.5*htarget)

            # Create the surface mesh
            mesh.mesh(fs=feature_size, opts=opts)

            # Create the corresponding mesh topology from the mesh-model
            model = mesh.createModelFromMesh()
            topo = TMR.Topology(comm, model)

            # Create the quad forest and set the topology of the forest
            depth = 0
            # if order == 2:
            #     depth = 1
            forest = TMR.QuadForest(comm)
            forest.setTopology(topo)
            forest.setMeshOrder(order, TMR.UNIFORM_POINTS)
            forest.createTrees(depth)
        else:
            # The refinement array
            refine = np.zeros(len(error), dtype=np.intc)

            # Determine the cutoff values
            cutoff = bins[-1]
            bin_sum = 0
            for i in range(len(bins)+1):
                bin_sum += bins[i]
                if bin_sum > 0.3*ntotal:
                    cutoff = bounds[i]
                    break

            log_cutoff = np.log(cutoff)

            # Element target error is still too high. Adapt based solely
            # on decreasing the overall error
            nrefine = 0
            for i, err in enumerate(error):
                # Compute the log of the error
                logerr = np.log(err)

                if logerr > log_cutoff:
                    refine[i] = 1
                    nrefine += 1

            # Refine the forest
            forest.refine(refine)
