from __future__ import print_function
from mpi4py import MPI
from tmr import TMR
from tacs import TACS, elements, constitutive, functions
from paropt import ParOpt
import numpy as np
import argparse
import os
import topo as const_topo

# Import pyoptsparse
from pyoptsparse import Optimization, OPT

# Import the smoother
import sys
sys.path.append('../smooth')
import locate

# Optimize the mesh spacing to meet criteria
def optimize_mesh(hvals, errors, mass_errors, lower, upper, Ntarget, p, d,
                  beta=0.25):
    def objfunc(xdict):
        '''Evaluate the objective/constraint'''
        x = xdict['x']

        # Minimize predicted error
        fobj = 0.0
        ratio = 0.0
        for i, h in enumerate(hvals):
            fobj += ((1.0 - beta)*errors[i]*(x[i]/h)**s +
                     beta*mass_errors[i]*(x[i]/h)**2)
            ratio += (x[i]/h)**(-d)

        # Set the objective and constraint
        funcs = {}
        funcs['fobj'] = fobj
        funcs['ratio'] = ratio/Ntarget

        fail = 0
        return funcs, fail

    def gobjfunc(xdict, funcs):
        '''Evaluate the objective/constraint derivatives'''
        x = xdict['x']

        # Minimize predicted error
        gobj = np.zeros(x.shape)
        gratio = np.zeros(x.shape)
        for i, h in enumerate(hvals):
            gobj[i] = (s*(errors[i]/h)*(x[i]/h)**(s-1) +
                       2*beta*(mass_errors[i]/h)*(x[i]/h))
            gratio[i] = -(d/h)*(x[i]/h)**(-d-1.0)
        gratio[:] /= Ntarget

        # Set the objective and constraint
        sens = {'fobj': {'x': gobj}, 'ratio': {'x': gratio}}

        fail = 0.0
        return sens, fail

    print('hvals = ', hvals)
    print('errors = ', errors)

    # Create the optimization problem
    mesh_prob = Optimization('mesh', objfunc, comm=MPI.COMM_SELF)

    # Add the variable group
    x0 = 0.5*(lower + upper)
    mesh_prob.addVarGroup('x', len(hvals), value=x0, lower=lower, upper=upper)

    # Add the constraints
    mesh_prob.addConGroup('ratio', 1, lower=1.0, upper=1.0)

    # Add the objective
    mesh_prob.addObj('fobj')

    # The Optimizer is IPOPT
    options = {}
    options['print_user_options'] = 'yes'
    options['tol'] = 1e-6
    options['bound_relax_factor'] = 0.0
    options['linear_solver'] = 'ma27'
    options['output_file'] = 'results/mesh_opt.out'
    options['max_iter'] = 500

    # Create the optimizer and optimize it!
    opt = OPT('ipopt', options=options)
    sol = opt(mesh_prob, sens=gobjfunc)

    return sol.xStar['x']

class MassMin:
    '''
    Mass minimization with a von Mises stress constraint
    '''
    def __init__(self, comm, nvars, locator):
        # Set the communicator
        self.comm = comm

        # Set the locator
        self.locator = locator

        # The number of topo variables
        self.nvars = nvars
        self.x = np.zeros(self.nvars)
        self.x[:] = 0.95
        self.x_scale = 100.0

        # Set the scaling on the mass objective
        xlen = 0.1
        self.mass_scale = 100.0/(0.64*xlen**2)
        self.con_scale = 1.0

        # The number of constraints (1 global stress constraint that
        # will use the KS function)
        self.ncon = 1

        self.iter_count = 0
        return

    def setAssembler(self, assembler, ksfunc):

        # Create tacs assembler object from mesh loader
        self.assembler = assembler

        # Create the list of functions
        self.funcs = [functions.StructuralMass(self.assembler), ksfunc]

        # Set up the solver
        self.ans = self.assembler.createVec()
        self.res = self.assembler.createVec()
        self.adjoint = self.assembler.createVec()
        self.dfdu = self.assembler.createVec()

        self.mat = self.assembler.createFEMat()
        self.pc = TACS.Pc(self.mat)
        self.gmres = TACS.KSM(self.mat, self.pc, 10)

        # For visualization
        flag = (TACS.ToFH5.NODES |
                TACS.ToFH5.DISPLACEMENTS |
                TACS.ToFH5.STRAINS |
                TACS.ToFH5.EXTRAS)
        self.f5 = TACS.ToFH5(self.assembler, TACS.PY_PLANE_STRESS, flag)

        return

    def getPoint(self, pt, xpts):
        '''Get the design point within the element'''
        order = int(np.sqrt(xpts.shape[0]/3))

        if order == 2:
            na = [0.5*(1.0 - pt[0]), 0.5*(1.0 + pt[0])]
            nb = [0.5*(1.0 - pt[1]), 0.5*(1.0 + pt[1])]
            naa = [-0.5, 0.5]
            nbb = [-0.5, 0.5]
        elif order == 3:
            na = [-0.5*(1.0 - pt[0])*pt[0],
                  1.0 - pt[0]*pt[0],
                  0.5*(1.0 + pt[0])*pt[0]]
            nb = [-0.5*(1.0 - pt[1])*pt[1],
                  1.0 - pt[1]*pt[1],
                  0.5*(1.0 + pt[1])*pt[1]]
            naa = [-0.5 + pt[0], -2*pt[0], 0.5 + pt[0]]
            nbb = [-0.5 + pt[1], -2*pt[1], 0.5 + pt[1]]

        # Compute the shape functions
        N = np.outer(nb, na).flatten()
        Na = np.outer(nb, naa).flatten()
        Nb = np.outer(nbb, na).flatten()

        # Compute the derivatives
        xpt = np.dot(N, xpts.reshape(-1, 3))
        xa = np.dot(Na, xpts.reshape(-1, 3))
        xb = np.dot(Nb, xpts.reshape(-1, 3))
        normal = np.cross(xa, xb)

        return xpt, np.sqrt(np.dot(normal, normal))

    def getDesignValue(self, xpt):
        '''Get the exact value of the design variable'''

        # Get the points that are the closest
        m = 100
        index = np.zeros(m, dtype=np.intc)
        dist = np.zeros(m)
        self.locator.locateKClosest(xpt, index, dist)
        dist = np.sqrt(dist)

        # Set the length as a fraction of the overall distance
        xlen = 0.1
        r0 = 0.05*xlen

        nodes = []
        weights = []
        for i in range(m):
            if dist[i] < r0:
                nodes.append(index[i])
                weights.append(((r0 - dist[i])/r0))
            else:
                break

        # Set the nodes/weights
        nodes = np.array(nodes, dtype=np.intc)
        weights = np.array(weights)
        weights[:] = weights/np.sum(weights)
        
        return np.dot(weights, self.x[nodes])

    def estimateMassError(self):
        '''Estimate the error in the mass estimate'''

        nelems = self.assembler.getNumElements()
        errors = np.zeros(nelems)

        # Set the length as a fraction of the overall distance
        xlen = 0.1
        r0 = 0.05*xlen

        # Set the quadrature estimates
        qwts2 = [1.0, 1.0]
        qpts2 = [-0.577350269189626, 0.577350269189626]

        qwts3 = [5.0/9.0, 8.0/9.0, 5.0/9.0]
        qpts3 = [-0.774596669241483, 0.0, 0.774596669241483]

        for i in range(nelems):
            # Get the node locations
            elem, xpts, v, dv, ddv = self.assembler.getElementData(i)

            # # Get the center point location
            # xpt, detJ = self.getPoint([0.0, 0.0], xpts)
            # xconst = self.getDesignValue(xpt)

            # Compute the difference in the area interval
            const = 0.0
            exact = 0.0

            # Compute the mass based on a second-order element
            for m in range(2):
                for n in range(2):
                    pt = [qpts2[n], qpts2[m]]
                    wval = qwts2[n]*qwts2[m]

                    # Get the point
                    xpt, detJ = self.getPoint(pt, xpts)

                    # Get the exact value of the density
                    xexact = self.getDesignValue(xpt)

                    # Compute the constant value/exact
                    const += xexact*wval*detJ

            for m in range(3):
                for n in range(3):
                    pt = [qpts3[n], qpts3[m]]
                    wval = qwts3[n]*qwts3[m]

                    # Get the point
                    xpt, detJ = self.getPoint(pt, xpts)

                    # Get the exact value of the density
                    xexact = self.getDesignValue(xpt)

                    # Compute the constant value/exact
                    exact += xexact*wval*detJ

            # Compute the absolute error in the mass
            errors[i] = np.fabs(const - exact)

        return errors

    def getVarsAndBounds(self, x, lb, ub):
        '''Set the values of the bounds'''
        x[:] = self.x_scale*self.x[:]
        xlb = 1e-3
        xub = 1.0
        lb[:] = self.x_scale*xlb
        ub[:] = self.x_scale*xub
        return

    def objcon(self, xdict):
        '''Evaluate the objective and constraint'''

        # Extract the values of x
        x = xdict['x']

        # Evaluate the objective
        fail, fobj, con = self.evalObjCon(x)

        # Create the dictionary of functions
        funcs = {'objective': fobj, 'con': con}

        return funcs, fail

    def gobjcon(self, xdict, funcs):
        '''Evaluate the objective and constraint gradient'''
        fail = 0

        # Extract the values of x
        x = xdict['x']

        gx = np.zeros(self.nvars)
        dfdx = np.zeros((1, self.nvars))
        fail = self.evalObjConGradient(x, gx, dfdx)

        # Create the sensitivity dictionary
        sens = {'objective':{'x': gx}, 'con':{'x': dfdx}}

        return sens, fail

    def evalObjCon(self, x):
        self.x[:] = x[:]/self.x_scale

        # Evaluate the objective and constraints
        fail = 0
        con = np.zeros(1)

        # Set the new design variable values
        self.assembler.setDesignVars(self.x)

        # Assemble the Jacobian and factor the matrix
        alpha = 1.0
        beta = 0.0
        gamma = 0.0
        self.assembler.zeroVariables()
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
        con[0] = self.con_scale*(1.0 - fvals[1])

        if self.comm.rank == 0:
            print('fobj = ', fobj)
            print('con = ', con)

        return fail, fobj, con

    def evalObjConGradient(self, x, g, A):
        '''Evaluate the objective and constraint gradient'''
        fail = 0

        # Evaluate the derivative of the mass and place it in the
        # objective gradient
        gx = np.zeros(self.nvars, TACS.dtype)
        self.assembler.evalDVSens(self.funcs[0], gx)
        g[:] = self.mass_scale*gx/self.x_scale

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
        A[0][:] = -self.con_scale*(dfdx - product)/self.x_scale

        # Write out the solution file every 10 iterations
        if self.iter_count % 1 == 0:
            if self.comm.rank == 0:
                print('Write out file %d'%(self.iter_count))
            self.f5.writeToFile('results/topo%04d.f5'%(self.iter_count))
        self.iter_count += 1

        return fail

class CreateMe(TMR.QuadCreator):
    def __init__(self, bcs, topo, locator):
        TMR.QuadCreator.__init__(bcs)
        self.topo = topo
        self.locator = locator
        return

    def createElement(self, order, quad):
        '''Create the element'''
        # Get the model name and set the face
        face = self.topo.getFace(quad.face)
        h = (1 << (TMR.MAX_LEVEL - quad.level))

        # Set the properties
        rho = 1.0
        E = 70e9
        nu = 0.3
        ys = 275e6
        p = 3.0
        eps = 0.2

        if False:
            x = 1.0*(quad.x + 0.5*h)/(1 << TMR.MAX_LEVEL)
            y = 1.0*(quad.y + 0.5*h)/(1 << TMR.MAX_LEVEL)

            # Evaluate the node locations
            pt = face.evalPoint(x, y)

            nodes = []
            weights = []

            m = 100
            index = np.zeros(m, dtype=np.intc)
            dist = np.zeros(m)
            self.locator.locateKClosest(pt, index, dist)
            dist = np.sqrt(dist)

            # Set the length as a fraction of the overall distance
            xlen = 0.1
            r0 = 0.05*xlen
            
            for i in range(m):
                if dist[i] < r0:
                    nodes.append(index[i])
                    weights.append(((r0 - dist[i])/r0))
                else:
                    break

            # Compute the nodes
            nodes = np.array(nodes, dtype=np.intc)
            weights = np.array(weights)
            weights[:] = weights/np.sum(weights)
            
            ps = const_topo.pstopo(rho, E, nu, ys, p, eps, nodes, weights)
        else:
            all_nodes = []
            all_weights = []

            for j in range(2):
                for i in range(2):          
                    x = 1.0*(quad.x + h*i)/(1 << TMR.MAX_LEVEL)
                    y = 1.0*(quad.y + h*j)/(1 << TMR.MAX_LEVEL)

                    # Evaluate the node locations
                    pt = face.evalPoint(x, y)

                    nodes = []
                    weights = []

                    m = 100
                    index = np.zeros(m, dtype=np.intc)
                    dist = np.zeros(m)
                    self.locator.locateKClosest(pt, index, dist)
                    dist = np.sqrt(dist)

                    # Set the length as a fraction of the overall distance
                    xlen = 0.1
                    r0 = 0.05*xlen

                    for i in range(m):
                        if dist[i] < r0:
                            nodes.append(index[i])
                            weights.append(((r0 - dist[i])/r0))
                        else:
                            break

                    # Compute the nodes
                    nodes = np.array(nodes, dtype=np.intc)
                    weights = np.array(weights)
                    weights[:] = weights/np.sum(weights)

                    all_nodes.append(nodes)
                    all_weights.append(weights)
            
            # Set the properties
            ps = const_topo.pstopo4(rho, E, nu, ys, p, eps,
                                    all_nodes[0], all_weights[0],
                                    all_nodes[1], all_weights[1],
                                    all_nodes[2], all_weights[2],
                                    all_nodes[3], all_weights[3])

        # Create the plane stree quadrilateral
        elem = elements.PlaneQuad(order, ps)

        return elem

def addFaceTraction(order, forest, assembler):
    # Get the quadrilaterals of interest
    quads = forest.getQuadsWithName('traction')

    # Create the surface traction
    aux = TACS.AuxElements()

    # Loop over the nodes and create the traction forces in the x/y/z
    # directions
    nnodes = order
    tx = np.zeros(nnodes, dtype=TACS.dtype)
    ty = np.zeros(nnodes, dtype=TACS.dtype)
    thickness = 1e-3
    ty[:] = 10e4/thickness

    # Create the shell traction
    surf = 1
    trac = elements.PSQuadTraction(surf, tx, ty)
    for q in quads:
        aux.addElement(q.tag, trac)

    return aux

def createProblem(topo, forest, bcs, locator, order=2, nlevels=2,
                  ordering=TACS.PY_MULTICOLOR_ORDER,
                  pttype=TMR.UNIFORM_POINTS):
    # Create the forest
    forests = []
    assemblers = []

    # Create the trees, rebalance the elements and repartition
    forest.balance(1)
    forest.setMeshOrder(order, pttype)
    forest.repartition()
    forests.append(forest)

    # Make the creator class
    creator = CreateMe(bcs, topo, locator)
    assemblers.append(creator.createTACS(forest, ordering))

    while order > 2:
        order = order-1
        forest = forests[-1].duplicate()
        forest.setMeshOrder(order, pttype)
        forest.balance(1)
        forests.append(forest)

        # Make the creator class
        creator = CreateMe(bcs, topo, locator)
        assemblers.append(creator.createTACS(forest, ordering))

    for i in range(nlevels-1):
        forest = forests[-1].coarsen()
        forest.setMeshOrder(2, pttype)
        forest.balance(1)
        forests.append(forest)

        # Make the creator class
        creator = CreateMe(bcs, topo, locator)
        assemblers.append(creator.createTACS(forest, ordering))

    # Create the multigrid object
    mg = TMR.createMg(assemblers, forests, omega=0.5)

    return assemblers, mg

# Set the communicator
comm = MPI.COMM_WORLD

# Create an argument parser to read in arguments from the commnad line
p = argparse.ArgumentParser()
p.add_argument('--steps', type=int, default=5)
p.add_argument('--htarget', type=float, default=0.001)
p.add_argument('--order', type=int, default=2)
p.add_argument('--ksweight', type=float, default=30.0)
p.add_argument('--optimizer', type=str, default='snopt')
p.add_argument('--remesh_strategy', type=str, default='fixed_mesh')
p.add_argument('--element_count_target', type=float, default=20e3)
args = p.parse_args()

# Set the KS parameter
ksweight = args.ksweight

# Set the number of AMR steps to use
steps = args.steps

# Set the order of the mesh
order = args.order

# Set the count target
element_count_target = args.element_count_target

fname = 'results/opt.out'
options = {}
if args.optimizer == 'snopt':
    options['Print file'] = fname
    options['Summary file'] = fname + '_summary'
    options['Penalty parameter'] = 10.0
    options['Nonderivative linesearch'] = None
    options['Minor print level'] = 0

    # Set the tolerances
    options['Major optimality tolerance'] = 1e-5
    options['Major feasibility tolerance'] = 1e-6
    options['Minor feasibility tolerance'] = 1e-8

    # Set a large number of iterations
    options['Minor iterations limit'] = 10000
elif args.optimizer == 'ipopt':
    options['print_user_options'] = 'yes'
    options['tol'] = tol
    options['nlp_scaling_method'] = 'none'
    options['limited_memory_max_history'] = 25
    options['bound_relax_factor'] = 0.0
    options['linear_solver'] = 'ma27'
    options['output_file'] = fname
    options['max_iter'] = 10000
elif args.optimizer == 'paropt':
    # options['this is an option'] = value
    options['algorithm'] = 'tr'
    options['abs_optimality_tol'] = 1e-8
    options['max_iterations'] = 200
    options['qn_subspace_size'] = 10

    # Set the trust region size
    options['tr_init_size'] = 0.5
    options['tr_max_size'] = 2.0
    options['tr_min_size'] = 1e-5
    options['tr_penalty_gamma'] = 20.0

# Set up the design variables
xlen = 0.1

# Compute the (approximate) area
Area = 0.64*xlen**2

# Set the bounds on the variables
a = 0.4*xlen + 1e-6

# Set the x/y locations of the density "nodes"
n = 100
x = np.linspace(0, xlen, n)

pts = []
for j in range(n):
    for i in range(n):
        if x[i] <= a or x[j] <= a:
            pts.append([x[i], x[j], 0.0])
locator = locate.locate(np.array(pts))

# Set the number of design variables
num_design_vars = len(pts)

# Load in the L-bracket model
geo = TMR.LoadModel('2d-bracket-fillet.stp')
verts = geo.getVertices()
edges = geo.getEdges()
faces = geo.getFaces()
geo = TMR.Model(verts, edges, faces)

# Set the edges
edges[5].setName('clamped')
edges[1].setName('traction')

# Initial target mesh spacing
htarget = args.htarget

# Create the new mesh
mesh = TMR.Mesh(comm, geo)

# Set the meshing options
opts = TMR.MeshOptions()

# Set the mesh type
opts.num_smoothing_steps = 10
opts.write_mesh_quality_histogram = 1
opts.triangularize_print_iter = 50000

# Create the surface mesh
mesh.mesh(htarget, opts)
mesh.writeToVTK('results/mesh.vtk')

# The boundary condition object
bcs = TMR.BoundaryConditions()
bcs.addBoundaryCondition('clamped')

# Set the feature size object
feature_size = None

# Create the corresponding mesh topology from the mesh-model
model = mesh.createModelFromMesh()
topo = TMR.Topology(comm, model)

# Create the optimization problem
opt_problem = MassMin(comm, num_design_vars, locator)

# Create the quad forest and set the topology of the forest
depth = 0
forest = TMR.QuadForest(comm)
forest.setTopology(topo)
forest.setMeshOrder(order, TMR.UNIFORM_POINTS)
forest.createTrees(depth)

# Null pointer to the optimizer
opt = None

for step in range(2*steps):
    # Balance the forest, create the assembler object
    forest.balance(1)
    forest.setMeshOrder(order, TMR.UNIFORM_POINTS)
    forest.repartition()
    creator = CreateMe(bcs, topo, locator)
    assembler = creator.createTACS(forest, order)

    # Add the face tractions to TACS
    aux = addFaceTraction(order, forest, assembler)
    assembler.setAuxElements(aux)

    # Create the KS functional
    func = functions.KSFailure(assembler, ksweight)
    func.setKSFailureType('continuous')

    # Set the new assembler object
    opt_problem.setAssembler(assembler, func)

    if step % 2 == 1:
        # Solve the analysis problem at the first step
        opt_problem.evalObjCon(opt_problem.x*opt_problem.x_scale)
    else:
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

        fname = 'results/opt%02d.out'%(step)
        if args.optimizer == 'snopt':
            options['Print file'] = fname
            options['Summary file'] = fname + '_summary'
        elif args.optimizer == 'ipopt':
            options['output_file'] = fname
        elif args.optimizer == 'paropt':
            options['tr_penalty_gamma'] = 20.0
            options['filename'] = fname
            if step >= 2:
                options['max_iterations'] = 100

        # Create the optimizer and optimize it!
        opt = OPT(args.optimizer, options=options)
        sol = opt(prob, sens=opt_problem.gobjcon)

    # Compute the aggregation function
    fval = assembler.evalFunctions([func])[0]

    # Write out the current solution
    filename = 'results/solution%02d.f5'%(step)
    opt_problem.f5.writeToFile(filename)

    # Create the refined mesh
    forest_refined = forest.duplicate()
    if order == 2:
        nlevs = 2
    else:
        nlevs = 1
    assembler_list, mg = createProblem(topo,
                                       forest_refined, bcs, locator,
                                       order=order+1, nlevels=nlevs+1)
    assembler_refined = assembler_list[0]
    aux = addFaceTraction(order+1, forest_refined, assembler_refined)
    assembler_refined.setAuxElements(aux)

    # Set the design variables for the refined problem
    for assemb in assembler_list:
        assemb.setDesignVars(opt_problem.x)

    # Compute the mass error estimate
    mass_error = opt_problem.estimateMassError()/Area

    # Extract the answer
    ans = opt_problem.ans

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

    # Compute the functional and the right-hand-side for the adjoint
    # on the refined mesh
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
                              forest_refined, assembler_refined,
                              adjoint, adjoint_interp)
    adjoint_refined.axpy(-1.0, adjoint_interp)

    err_est, __, error = TMR.adjointError(forest, assembler,
                                          forest_refined, assembler_refined,
                                          ans_interp, adjoint_refined)

    # Compute the refined function value
    fval_corr = fval_refined + adjoint_corr

    flag = (TACS.ToFH5.NODES |
            TACS.ToFH5.DISPLACEMENTS |
            TACS.ToFH5.STRAINS |
            TACS.ToFH5.EXTRAS)
    f5_refine = TACS.ToFH5(assembler_refined, TACS.PY_PLANE_STRESS, flag)
    f5_refine.writeToFile('results/solution_refined%02d.f5'%(step))

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
        np.savetxt('results/topo_data%d.txt'%(step), data)

    # Perform the refinement
    # Ensure that we're using an unstructured mesh
    opts.mesh_type_default = TMR.UNSTRUCTURED

    # Find the positions of the center points of each node
    nelems = assembler.getNumElements()

    # Allocate the positions
    Xp = np.zeros((nelems, 3))
    xrho = np.zeros(nelems)
    param = np.zeros(3)
    for i in range(nelems):
        # Get the information about the given element
        elem, Xpt, vrs, dvars, ddvars = assembler.getElementData(i)
        c = elem.getConstitutive()
        xrho[i] = c.getDVOutputValue(0, param)

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
        mass_errors = np.zeros(np.sum(count))
        Xpt = np.zeros(3*np.sum(count))
        comm.Gatherv(error, [errors, count])
        comm.Gatherv(mass_error, [mass_errors, count])
        comm.Gatherv(Xp.flatten(), [Xpt, 3*count])

        xr = np.zeros(np.sum(count))
        comm.Gatherv(xrho, [xr, count])

        # Reshape the point array
        Xpt = Xpt.reshape(-1,3)

        # Asymptotic order of accuracy on per-element basis
        s = args.order-1

        # Dimension of the problem
        d = 2.0

        # Set upper/lower limits on the mesh spacing
        hmax = 2*args.htarget
        hmin = 0.01*args.htarget

        if True:
            # Evaluate the exiting element sizes
            h = np.zeros(Xpt.shape[0])
            if feature_size is not None:
                for i in range(len(h)):
                    h[i] = feature_size.getFeatureSize(Xpt[i,:])
            else:
                h[:] = htarget

            lower = np.zeros(h.shape)
            upper = np.zeros(h.shape)
            for i in range(len(h)):
                lower[i] = max(hmin, 0.25*h[i])
                upper[i] = min(hmax, 2.0*h[i])

            # Optimize the mesh
            hvals = optimize_mesh(h, errors, mass_errors, lower, upper,
                                  element_count_target, s, d)
        else:
            # Set the exponent
            exponent = d/(d + s)

            # Compute the target error as a fixed fraction of the error
            # estimate. This will result in the size of the mesh
            # increasing at each iteration.
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

            # Decide on whether to refine further based on the value of
            # the design variables -- target transition areas
            for i, hp in enumerate(hvals):
                hvals[i] = min(hp, 0.5 + 2*(xr[i] - 0.5)**2)

            # Set the new h value based on the previous feature
            # size value
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
        feature_size = TMR.PointFeatureSize(Xpt, hvals, hmin, hmax)
    else:
        size = error.shape[0]
        comm.gather(size, root=root)
        comm.Gatherv(error, None)
        comm.Gatherv(mass_error, None)
        comm.Gatherv(Xp.flatten(), None)
        comm.Gatherv(xrho, None)

        # Create a dummy feature size object...
        feature_size = TMR.ConstElementSize(0.5*htarget)

    # Create the surface mesh
    mesh.mesh(fs=feature_size, opts=opts)

    # Create the corresponding mesh topology from the mesh-model
    model = mesh.createModelFromMesh()
    topo = TMR.Topology(comm, model)

    # Create the quad forest and set the topology of the forest
    depth = 0
    forest = TMR.QuadForest(comm)
    forest.setTopology(topo)
    forest.setMeshOrder(order, TMR.UNIFORM_POINTS)
    forest.createTrees(depth)
