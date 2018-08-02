from __future__ import print_function
from mpi4py import MPI
from tmr import TMR
from tacs import TACS, elements, constitutive, functions
from paropt import ParOpt
import numpy as np
import argparse
import os
import ksFSDT


class uCRM_VonMisesMassMin(ParOpt.pyParOptProblem):
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

        # Initialize the base class - this will run the same problem
        # on all processors
        super(uCRM_VonMisesMassMin, self).__init__(MPI.COMM_SELF,
                                                   self.nvars, self.ncon)

        # Set the inequality options for this problem in ParOpt:
        # The dense constraints are inequalities c(x) >= 0 and
        # use both the upper/lower bounds
        self.setInequalityOptions(dense_ineq=True,
                                  use_lower=True, use_upper=True)
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
        self.mg = mg
        self.mat = self.mg.getMat()
        self.gmres = gmres

        # For visualization
        flag = (TACS.ToFH5.NODES |
                TACS.ToFH5.DISPLACEMENTS |
                TACS.ToFH5.STRAINS |
                TACS.ToFH5.EXTRAS)
        self.f5 = TACS.ToFH5(self.assembler, TACS.PY_SHELL, flag)
        self.iter_count = 0

        return

    def getVarsAndBounds(self, x, lb, ub):
        '''Set the values of the bounds'''
        xvals = np.zeros(self.nvars, TACS.dtype)
        self.assembler.getDesignVars(xvals)
        x[:] = self.thickness_scale*xvals

        xlb = np.zeros(self.nvars, TACS.dtype)
        xub = np.zeros(self.nvars, TACS.dtype)
        self.assembler.getDesignVarRange(xlb, xub)
        lb[:] = self.thickness_scale*xlb
        ub[:] = self.thickness_scale*xub

        return

    def evalObjCon(self, x):
        '''Evaluate the objective and constraint'''
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
        self.mg.assembleJacobian(alpha, beta, gamma, self.res)
        self.mg.factor()

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

        return fail, fobj, con

    def evalObjConGradient(self, x, g, A):
        '''Evaluate the objective and constraint gradient'''
        fail = 0

        # Evaluate the derivative of the mass and place it in the
        # objective gradient
        gx = np.zeros(self.nvars, TACS.dtype)
        self.assembler.evalDVSens(self.funcs[0], gx)
        g[:] = self.mass_scale*gx/self.thickness_scale

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
        A[0][:] = -(dfdx - product)/self.thickness_scale

        # Write out the solution file every 10 iterations
        if self.iter_count % 10 == 0:
            self.f5.writeToFile('results/ucrm_iter%d.f5'%(self.iter_count))
        self.iter_count += 1

        return fail

class CreateMe(TMR.QuadCreator):
    def __init__(self, bcs, topo, elem_dict):
        TMR.QuadCreator.__init__(bcs)
        self.topo = topo
        self.elem_dict = elem_dict
        return

    def createElement(self, order, quad):
        '''Create the element'''
        # Get the model attribute and set the face
        attr = topo.getFace(quad.face).getAttribute()
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
    tz[:] = 10.0

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

# The root rib for boundary conditions
ucrm_root_rib = 149

# The 47 top and bottom skin segments
ucrm_top_skins = [127, 141,  30, 133, 154,  10,  29,  33, 116, 100,
                  42,  104, 164,  94, 159, 160, 124, 142, 150,  78,
                  87,  122, 130, 153, 123,  58,  52, 157,   0, 152,
                  54,  172,  61, 108, 148,  62, 173,  51, 107, 65,
                  115,  37,  72, 112,  73,  15, 106]
ucrm_bottom_skins = [117,  91,  97,   5, 183,  45, 180,  31,  28, 129,
                     155,  96, 162, 139,  63, 147,  77,  26,  59,  20,
                     80,   53,  88,  23, 151,  84,  89, 126, 114, 121,
                     179, 103, 184, 120, 166, 167,  98,  79,  27, 128,
                     41,  146,  39, 101,  83,  57,  90]

# Create the CRM wingbox model and set the attributes
geo = TMR.LoadModel('ucrm_9_full_model.step')
verts = geo.getVertices()
edges = geo.getEdges()
faces = geo.getFaces()
geo = TMR.Model(verts, edges, faces)

# Set the number of design variables
num_design_vars = len(faces)

# Set the material properties
rho = 97.5e-3 # 0.0975 lb/in^3
E = 10000e3 # 10,000 ksi: Young's modulus
nu = 0.3 # Poisson ratio
kcorr = 5.0/6.0 # Shear correction factor
ys = 40e3 # psi
thickness = 1.0
min_thickness = 0.1
max_thickness = 10.0

# Set the face attributes and create the constitutive objects
elem_dict = [{}, {}, {}, {}]
for i, f in enumerate(faces):
    attr = 'F%d'%(i)
    f.setAttribute(attr)
    stiff = constitutive.isoFSDT(rho, E, nu, kcorr, ys, thickness,
                                 i, min_thickness, max_thickness)

    # Set the component number for visualization purposes
    comp = 0
    if i in ucrm_top_skins:
        comp = 1
    elif i in ucrm_bottom_skins:
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
opts.triangularize_print_level = 2
opts.frontal_quality_factor = 1.25
opts.num_smoothing_steps = 50
opts.write_mesh_quality_histogram = 1
opts.triangularize_print_iter = 50000

opts.write_pre_smooth_triangle = 1
opts.write_post_smooth_triangle = 1
opts.write_dual_recombine = 1
opts.write_pre_smooth_quad = 1
opts.write_post_smooth_quad = 1

# Create the surface mesh
mesh.mesh(htarget, opts)
mesh.writeToVTK('mesh.vtk')

exit(0)

# The boundary condition object
bcs = TMR.BoundaryConditions()
bcs.addBoundaryCondition('F%d'%(ucrm_root_rib), [0, 1, 2, 3, 4, 5])

# Set the feature size object
feature_size = None

# Create the corresponding mesh topology from the mesh-model
model = mesh.createModelFromMesh()
topo = TMR.Topology(comm, model)

# Create the optimization problem
opt_problem = uCRM_VonMisesMassMin(comm, num_design_vars)

# Create the quad forest and set the topology of the forest
depth = 0
if order == 2:
    depth = 1
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

    if k == 0:
        # Set up the optimization problem
        max_lbfgs = 5
        opt = ParOpt.pyParOpt(opt_problem, max_lbfgs, ParOpt.BFGS)
        opt.setOutputFile('results/crm_opt.out')

        # Set the optimality tolerance
        opt.setAbsOptimalityTol(1e-4)

        # Set optimization parameters
        opt.setArmijoParam(1e-5)

        # Get the optimized point
        x, z, zw, zl, zu = opt.getOptimizedPoint()

        # Set the starting point strategy
        opt.setStartingPointStrategy(ParOpt.AFFINE_STEP)

        # Set the max oiterations
        opt.setMaxMajorIterations(100)

        # Set the output level to understand what is going on
        opt.setOutputLevel(2)
    else:
        beta = 1e-4
        opt.resetDesignAndBounds()
        opt.setStartAffineStepMultiplierMin(beta)

    # Optimize the new point
    opt.optimize()

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
    aux = addFaceTraction(order+1, assembler)
    assembler_refined.setAuxElements(aux)

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

                # Compute the target relative error in each element.
                # This is set as a fraction of the error in the
                # current mesh level.
                err_target = 0.1*(err_est/ntotal)

                # Compute the element size factor in each element
                # using the order or accuracy
                hvals = (err_target/errors)**(0.5)

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
                hmax = 10.0
                hmin = 0.25
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
            if order == 2:
                depth = 1
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
