from __future__ import print_function
from mpi4py import MPI
from tmr import TMR
from tacs import TACS, elements, functions
import numpy as np
import argparse
import os

def integrate(integrand):
    '''Integrate over equally spaced data'''
    sigma = [17.0/48.0, 59.0/48.0, 43.0/48, 49/48.9]
    r = len(sigma)

    integral = 0.0
    for i, s in enumerate(sigma):
        integral += s*integrand[i]
        integral += s*integrand[-1-i]

    for i in range(r, len(integrand)-r):
        integral += integrand[i]

    return integral

def poisson_evalf(x):
    R = 100.0
    return (3920.0/363)*(1.0 - (x[0]*x[0] + x[1]*x[1])/R**2)**6

def get_disk_aggregate(functional, rho, R, n=1000):
    '''Evaluate the KS functional on a disk'''
    scale = 1.0/R**2
    r = np.linspace(0.0, R, n)

    # This is the solution for constant order
    # phi = (0.1875*R**2 - 0.25*r**2 + (0.0625/R**2)*r**4)

    # Compute the solution scaled to r/R
    x = np.linspace(0.0, 1.0, n)
    phi = R**2*(-(20.0/363)*x**14 + (490.0/1089)*x**12 - (196.0/121)*x**10 +
                (1225.0/363)*x**8 - (4900.0/1089)*x**6 + (490.0/121)*x**4 -
                (980.0/363)*x**2)
    phi -= phi[-1]

    if functional == 'ks':
        ksmax = np.max(phi)
        integrand = r*np.exp(rho*scale*(phi - ksmax))
        kssum = 2*np.pi*(R/(n-1))*integrate(integrand)
        return scale*ksmax, scale*ksmax + np.log(kssum)/rho
    else:
        maxphi = scale*np.max(phi)
        integrand = r*np.power((scale*phi/maxphi), rho)
        psum = 2*np.pi*(R/(n-1))*integrate(integrand)
        return maxphi, maxphi*np.power(psum, 1.0/rho)

class CreateMe(TMR.QuadCreator):
    def __init__(self, bcs, topo, case='disk'):
        self.case = case
        TMR.QuadCreator.__init__(bcs)
        self.topo = topo
        return

    def createElement(self, order, quad):
        '''Create the element'''
        elem = elements.PoissonQuad(order, evalf=poisson_evalf)
        return elem

def createRefined(case, forest, bcs, pttype=TMR.GAUSS_LOBATTO_POINTS):
    new_forest = forest.duplicate()
    new_forest.setMeshOrder(forest.getMeshOrder()+1, pttype)
    creator = CreateMe(bcs, forest.getTopology(), case=case)
    return new_forest, creator.createTACS(new_forest)

def createProblem(case, forest, bcs, ordering, order=2, nlevels=2,
                  pttype=TMR.GAUSS_LOBATTO_POINTS):
    # Create the forest
    forests = []
    assemblers = []

    # Create the trees, rebalance the elements and repartition
    forest.setMeshOrder(order, pttype)
    forests.append(forest)

    # Make the creator class
    creator = CreateMe(bcs, forest.getTopology(), case=case)
    assemblers.append(creator.createTACS(forest, ordering))

    while order > 2:
        order = order-1
        forest = forests[-1].duplicate()
        forest.setMeshOrder(order, pttype)
        forest.balance(1)
        forests.append(forest)

        # Make the creator class
        creator = CreateMe(bcs, forest.getTopology(), case=case)
        assemblers.append(creator.createTACS(forest, ordering))

    for i in range(nlevels-1):
        forest = forests[-1].coarsen()
        forest.setMeshOrder(2, pttype)
        forest.balance(1)
        forests.append(forest)

        # Make the creator class
        creator = CreateMe(bcs, forest.getTopology(), case=case)
        assemblers.append(creator.createTACS(forest, ordering))

    # Create the multigrid object
    mg = TMR.createMg(assemblers, forests, omega=0.5)

    return assemblers[0], mg

# Set the communicator
comm = MPI.COMM_WORLD

# Create an argument parser to read in arguments from the commnad line
p = argparse.ArgumentParser()
p.add_argument('--case', type=str, default='disk')
p.add_argument('--functional', type=str, default='ks')
p.add_argument('--steps', type=int, default=5)
p.add_argument('--htarget', type=float, default=10.0)
p.add_argument('--ordering', type=str, default='multicolor')
p.add_argument('--order', type=int, default=2)
p.add_argument('--ksweight', type=float, default=100.0)
p.add_argument('--uniform_refinement', action='store_true', default=False)
p.add_argument('--structured', action='store_true', default=False)
p.add_argument('--energy_error', action='store_true', default=False)
p.add_argument('--exact_refined_adjoint', action='store_true', default=False)
p.add_argument('--remesh_domain', action='store_true', default=False)
p.add_argument('--remesh_strategy', type=str, default='fraction')
p.add_argument('--element_count_target', type=float, default=20e3)
args = p.parse_args()

# Set the case type
case = args.case

# Set the type of ordering to use for this problem
ordering = args.ordering
ordering = ordering.lower()

# Set the value of the target length scale in the mesh
htarget = args.htarget

# Set the order
order = args.order

# Set the KS parameter
ksweight = args.ksweight

# Set the number of AMR steps to use
steps = args.steps

# Store whether to use a structured/unstructed mesh
structured = args.structured

# Set what functional to use
functional = args.functional

# This flag indicates whether to solve the adjoint exactly on the
# next-finest mesh or not
exact_refined_adjoint = args.exact_refined_adjoint

# Set the count target
element_count_target = args.element_count_target

# The boundary condition object
bcs = TMR.BoundaryConditions()

# The radius of the disk
R = 100.0

# The side-length of the rectangular plate
L = 200.0

if case == 'disk':
    for n in [10000, 100000, 1000000, 10000000]:
        ks_max, approx = get_disk_aggregate(functional, ksweight, R, n=n)
        if comm.rank == 0:
            print('%10d %25.16e %25.16e'%(n, ks_max, approx))
        exact_functional = approx

    geo = TMR.LoadModel('2d-disk.stp')
    verts = geo.getVertices()
    edges = geo.getEdges()
    faces = geo.getFaces()

    # Set the names
    verts[1].setName('midpoint')
    for index in [0, 2, 3, 4]:
        verts[index].setName('clamped')
    for index in [2, 4, 6, 7]:
        edges[index].setName('clamped')

    # Set the boundary conditions
    bcs.addBoundaryCondition('clamped', [0])

if comm.rank == 0:
    print('Exact functional = %25.15e'%(exact_functional))

# Create the new mesh
mesh = TMR.Mesh(comm, geo)

# Set the meshing options
opts = TMR.MeshOptions()

# Set the mesh type
opts.mesh_type_default = TMR.UNSTRUCTURED
if structured:
    opts.mesh_type_default = TMR.STRUCTURED

opts.frontal_quality_factor = 1.25
opts.num_smoothing_steps = 50
opts.write_mesh_quality_histogram = 1
opts.triangularize_print_iter = 50000

# Create the surface mesh
mesh.mesh(htarget, opts)

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

# Set the ordering to use
if ordering == 'rcm':
    ordering = TACS.PY_RCM_ORDER
elif ordering == 'multicolor':
    ordering = TACS.PY_MULTICOLOR_ORDER
else:
    ordering = TACS.PY_NATURAL_ORDER

# The target relative error
target_rel_err = 1e-4

# Open a log file to write
descript = '%s_order%d_poisson'%(case, order)
if args.uniform_refinement:
    descript += '_uniform'
descript += '_' + functional

# Add a description about the meshing strategy
if args.remesh_domain:
    descript += '_' + args.remesh_strategy
    if args.remesh_strategy == 'fixed_mesh':
        descript += '_%g'%(args.element_count_target)

# Create the log file and write out the header
log_fp = open('results/%s.dat'%(descript), 'w')
s = 'Variables = iter, nelems, nnodes, fval, fcorr, abs_err, adjoint_corr, '
s += 'exact, fval_error, fval_corr_error, '
s += 'fval_effectivity, indicator_effectivity\n'
log_fp.write(s)

# Set the feature size object
feature_size = None

for k in range(steps):
    # Create the topology problem
    nlevs = min(4, depth+k+1)
    assembler, mg = createProblem(case, forest, bcs, ordering,
                                  order=order, nlevels=nlevs)

    # Create the assembler object
    res = assembler.createVec()
    ans = assembler.createVec()

    use_direct_solve = False
    if use_direct_solve:
        mat = assembler.createFEMat()
        assembler.assembleJacobian(1.0, 0.0, 0.0, res, mat)

        # Create the direct solver
        pc = TACS.Pc(mat)
        pc.factor()
    else:
        # Solve the linear system
        mg.assembleJacobian(1.0, 0.0, 0.0, res)
        mg.factor()
        pc = mg
        mat = mg.getMat()

    # Create the GMRES object
    gmres = TACS.KSM(mat, pc, 100, isFlexible=1)
    gmres.setMonitor(comm, freq=10)
    gmres.setTolerances(1e-14, 1e-30)
    gmres.solve(res, ans)
    ans.scale(-1.0)

    # Set the variables
    assembler.setVariables(ans)

    flag = (TACS.ToFH5.NODES |
            TACS.ToFH5.DISPLACEMENTS |
            TACS.ToFH5.STRAINS)
    f5 = TACS.ToFH5(assembler, TACS.PY_POISSON_2D_ELEMENT, flag)
    f5.writeToFile('results/solution%02d.f5'%(k))

    # Create and compute the function
    fval = 0.0
    if functional == 'ks':
        direction = [1.0/R**2, 0.0, 0.0]
        func = functions.KSDisplacement(assembler, ksweight, direction)
        func.setKSDispType('continuous')
        fval = assembler.evalFunctions([func])[0]
    else:
        direction = [1.0/R**2, 0.0, 0.0]
        func = functions.KSDisplacement(assembler, ksweight, direction)
        func.setKSDispType('pnorm-continuous')
        fval = assembler.evalFunctions([func])[0]

    # Create the refined mesh
    if exact_refined_adjoint:
        forest_refined = forest.duplicate()
        assembler_refined, mg = createProblem(case, forest_refined,
                                              bcs, ordering,
                                              order=order+1, nlevels=nlevs)
    else:
        forest_refined, assembler_refined = createRefined(case, forest, bcs)

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
                                      forest_refined, assembler_refined,
                                      ans, ans_interp)

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
            if functional == 'ks':
                direction = [1.0/R**2, 0.0, 0.0]
                func_refined = functions.KSDisplacement(assembler_refined,
                                                        ksweight, direction)
                func_refined.setKSDispType('continuous')
            else:
                direction = [1.0/R**2, 0.0, 0.0]
                func_refined = functions.KSDisplacement(assembler_refined,
                                                        ksweight, direction)
                func_refined.setKSDispType('pnorm-continuous')

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
                                      forest, assembler,
                                      adjoint_refined, adjoint)
            TMR.computeInterpSolution(forest, assembler,
                                      forest_refined, assembler_refined,
                                      adjoint, adjoint_interp)
            adjoint_refined.axpy(-1.0, adjoint_interp)

            # Estimate the element-wise errors
            err_est, __, error = TMR.adjointError(forest, assembler,
                                                  forest_refined, assembler_refined,
                                                  ans_interp, adjoint_refined)
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
            if functional == 'ks':
                direction = [1.0/R**2, 0.0, 0.0]
                func_refined = functions.KSDisplacement(assembler_refined,
                                                        ksweight, direction)
                func_refined.setKSDispType('continuous')
            else:
                direction = [1.0/R**2, 0.0, 0.0]
                func_refined = functions.KSDisplacement(assembler_refined,
                                                        ksweight, direction)
                func_refined.setKSDispType('pnorm-continuous')

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

    f5_refine = TACS.ToFH5(assembler_refined, TACS.PY_POISSON_2D_ELEMENT, flag)
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

    # Compute the error from the exact solution
    fval_error = np.fabs(exact_functional - fval)
    fval_corr_error = np.fabs(exact_functional - fval_corr)

    fval_effectivity = (fval_corr - fval)/(exact_functional - fval)
    indicator_effectivity = err_est/np.fabs(exact_functional - fval)

    # Write the log data to a file
    s = '%6d %6d %6d %20.15e %20.15e %20.15e %20.15e '
    s += '%20.15e %20.15e %20.15e %20.15e %20.15e\n'
    log_fp.write(
        s%(k, ntotal, nnodes, fval, fval_corr, err_est, adjoint_corr,
        exact_functional, fval_error, fval_corr_error,
        fval_effectivity, indicator_effectivity))
    log_fp.flush()

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
        np.savetxt('results/%s_data%d.txt'%(descript, k), data)

    # Print out the error estimate
    assembler.setDesignVars(error)
    f5.writeToFile('results/error%02d.f5'%(k))

    # Perform the refinement
    if args.uniform_refinement:
        forest.refine()
        forest.balance(1)
        forest.repartition()
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
                    # Compute the constant for element count
                    cval = (element_count_target)**(-1.0/d)
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
                hmax = 10.0
                hmin = 0.1
                feature_size = TMR.PointFeatureSize(Xpt, hvals, hmin, hmax,
                                                    num_sample_pts=16)
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
            # Allocate the refinement array
            refine = np.zeros(len(error), dtype=np.intc)

            # Determine the cutoff values
            cutoff = bins[-1]
            bin_sum = 0
            for i in range(len(bins)+1):
                bin_sum += bins[i]
                if bin_sum > 0.15*ntotal:
                    cutoff = bounds[i+1]
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
            forest.balance(1)
            forest.repartition()
