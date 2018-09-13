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

class MassMin:
    '''
    Mass minimization with a von Mises stress constraint
    '''
    def __init__(self, comm, nvars):
        # Set the communicator
        self.comm = comm

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
        self.f5 = TACS.ToFH5(self.assembler, TACS.PY_PLANE_STRESS, flag)

        return

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
    def __init__(self, bcs, topo):
        TMR.QuadCreator.__init__(bcs)
        self.topo = topo
        return

    def createElement(self, order, quad):
        '''Create the element'''
        # Get the model name and set the face
        face = self.topo.getFace(quad.face)
        h = (1 << (TMR.MAX_LEVEL - quad.level))
        x = 1.0*(quad.x + 0.5*h)/(1 << TMR.MAX_LEVEL)
        y = 1.0*(quad.y + 0.5*h)/(1 << TMR.MAX_LEVEL)

        # Evaluate the node locations
        pt = face.evalPoint(x, y)

        nodes = []
        weights = []

        m = 100
        index = np.zeros(m, dtype=np.intc)
        dist = np.zeros(m)
        locator.locateKClosest(pt, index, dist)
        dist = np.sqrt(dist)

        # Set the length as a fraction of the overall distance
        xlen = 0.1
        r0 = 0.05*xlen

        for i in range(m):
            if dist[i] < r0:
                nodes.append(index[i])
                weights.append(((r0 - dist[i])/r0)**1.5)
            else:
                break
        nodes = np.array(nodes, dtype=np.intc)
        weights = np.array(weights)
        weights[:] = weights/np.sum(weights)

        # Set the properties
        rho = 1.0
        E = 70e9
        nu = 0.3
        ys = 275e6
        q = 5.0
        eps = 0.15
        ps = const_topo.pstopo(rho, E, nu, ys, q, eps, nodes, weights)

        # Create the plane stree quadrilateral
        elem = elements.PlaneQuad(order, ps)

        return elem

def addFaceTraction(order, forest, assembler):
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

def createRefined(topo, forest, bcs, pttype=TMR.UNIFORM_POINTS):
    new_forest = forest.duplicate()
    new_forest.setMeshOrder(forest.getMeshOrder()+1, pttype)
    creator = CreateMe(bcs, topo)
    return new_forest, creator.createTACS(new_forest)

def createProblem(topo, forest, bcs, ordering,
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
    creator = CreateMe(bcs, topo)
    assemblers.append(creator.createTACS(forest, ordering))

    while order > 2:
        order = order-1
        forest = forests[-1].duplicate()
        forest.setMeshOrder(order, pttype)
        forest.balance(1)
        forests.append(forest)

        # Make the creator class
        creator = CreateMe(bcs, topo)
        assemblers.append(creator.createTACS(forest, ordering))

    for i in range(nlevels-1):
        forest = forests[-1].coarsen()
        forest.setMeshOrder(2, pttype)
        forest.balance(1)
        forests.append(forest)

        # Make the creator class
        creator = CreateMe(bcs, topo)
        assemblers.append(creator.createTACS(forest, ordering))

    # Create the multigrid object
    mg = TMR.createMg(assemblers, forests, omega=0.5)

    return assemblers[0], mg

# Set the communicator
comm = MPI.COMM_WORLD

# Create an argument parser to read in arguments from the commnad line
p = argparse.ArgumentParser()
p.add_argument('--steps', type=int, default=5)
p.add_argument('--htarget', type=float, default=0.001)
p.add_argument('--order', type=int, default=2)
p.add_argument('--ordering', type=str, default='multicolor')
p.add_argument('--ksweight', type=float, default=30.0)
p.add_argument('--uniform_refinement', action='store_true', default=False)
p.add_argument('--structured', action='store_true', default=False)
p.add_argument('--energy_error', action='store_true', default=False)
p.add_argument('--compute_solution_error', action='store_true', default=False)
p.add_argument('--exact_refined_adjoint', action='store_true', default=False)
p.add_argument('--remesh_domain', action='store_true', default=False)
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
fname = 'results/opt.out'
options = {}
if args.optimizer == 'snopt':
    options['Print file'] = fname
    options['Summary file'] = fname + '_summary'
    options['Major optimality tolerance'] = tol
    options['Penalty parameter'] = 10.0
    options['Nonderivative linesearch'] = None
    options['Minor print level'] = 0
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

# Set up the design variables
xlen = 0.1

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

# Create the CRM wingbox model and set the names
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
# opts.mesh_type_default = TMR.TRIANGLE
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
opt_problem = MassMin(comm, num_design_vars)

# Create the quad forest and set the topology of the forest
depth = 0
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
    assembler, mg = createProblem(topo, forest, bcs, ordering,
                                  order=order, nlevels=nlevs)
    aux = addFaceTraction(order, forest, assembler)
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

    fname = 'results/opt%02d.out'%(k)
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
        assembler_refined, mg = createProblem(topo,
                                              forest_refined, bcs, ordering,
                                              order=order+1, nlevels=nlevs+1)
    else:
        forest_refined, assembler_refined = createRefined(topo, elem_dict,
                                                          forest_refined, bcs)
    aux = addFaceTraction(order+1, forest_refined, assembler_refined)
    assembler_refined.setAuxElements(aux)

    # Set the design variables for the refined problem
    assembler_refined.setDesignVars(opt_problem.x)

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
    f5_refine = TACS.ToFH5(assembler_refined, TACS.PY_PLANE_STRESS, flag)
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
                Xpt = np.zeros(3*np.sum(count))
                comm.Gatherv(error, [errors, count])
                comm.Gatherv(Xp.flatten(), [Xpt, 3*count])

                xr = np.zeros(np.sum(count))
                comm.Gatherv(xrho, [xr, count])

                # Reshape the point array
                Xpt = Xpt.reshape(-1,3)

                # Compute the target relative error in each element.
                # This is set as a fraction of the error in the
                # current mesh level.
                err_target = 0.1*(err_est/ntotal)

                # Compute the element size factor in each element
                # using the order of accuracy
                hvals = (err_target/errors)**(0.5)

                # Decide on whether to refine based on the value of
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
                hmax = args.htarget
                hmin = 0.05*args.htarget
                feature_size = TMR.PointFeatureSize(Xpt, hvals, hmin, hmax)
            else:
                size = error.shape[0]
                comm.gather(size, root=root)
                comm.Gatherv(error, None)
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
