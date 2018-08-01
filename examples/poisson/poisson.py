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

def get_disk_aggregate(rho, R, n=1000):
    '''Evaluate the KS functional on a disk'''
    scale = 1.0/R**2    
    r = np.linspace(0.0, R, n)
    phi = (0.1875*R**2 - 0.25*r**2 + (0.0625/R**2)*r**4)
    ksmax = np.max(phi)
    integrand = r*np.exp(rho*scale*(phi - ksmax))
    kssum = 2*np.pi*(R/(n-1))*integrate(integrand)
    return scale*ksmax + np.log(kssum)/rho

def get_midpoint_vector(comm, forest, assembler, attr):
    vec = assembler.createVec()
    v = vec.getArray()
    nodes = forest.getNodesWithAttribute(attr)
    node_range = forest.getNodeRange()
    for n in nodes:
        if n >= node_range[comm.rank] and n < node_range[comm.rank+1]:
            index = n - node_range[comm.rank]
            v[index] = 1.0
    assembler.reorderVec(vec)
    return vec

class CreateMe(TMR.QuadCreator):
    def __init__(self, bcs, topo, case='disk'):
        self.case = case
        TMR.QuadCreator.__init__(bcs)
        self.topo = topo
        return

    def createElement(self, order, quad):
        '''Create the element'''
        f = np.ones(order*order)
        face = topo.getFace(quad.face)
        h = (1 << (TMR.MAX_LEVEL - quad.level))
        x = 1.0*quad.x/(1 << TMR.MAX_LEVEL)
        y = 1.0*quad.y/(1 << TMR.MAX_LEVEL)
        d = 1.0*h/(1 << TMR.MAX_LEVEL)
        u = 0.5*d*(1.0 - np.cos(np.linspace(0, np.pi, order)))

        if case == 'square':
            for j in range(order):
                for i in range(order):
                    # Evaluate the node locations
                    pt = face.evalPoint(x + u[i], y + u[j])
                    sx = np.sin(np.pi*(pt[0] + 0.5*L)/L)
                    sy = np.sin(np.pi*(pt[1] + 0.5*L)/L)
                    f[i + j*order] = sx*sy
        elif case == 'disk':
            for j in range(order):
                for i in range(order):
                   # Evaluate the node locations
                    pt = face.evalPoint(x + u[i], y + u[j])
                    r = np.sqrt(pt[0]*pt[0] + pt[1]*pt[1])
                    f[i + j*order] = 1.0 - (r/R)**2

        elem = elements.PoissonQuad(order, f)
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
p.add_argument('--functional', type=str, default='integral')
p.add_argument('--steps', type=int, default=5)
p.add_argument('--htarget', type=float, default=10.0)
p.add_argument('--ordering', type=str, default='multicolor')
p.add_argument('--order', type=int, default=2)
p.add_argument('--ksweight', type=float, default=50.0)
p.add_argument('--uniform_refinement', action='store_true', default=False)
p.add_argument('--structured', action='store_true', default=False)
p.add_argument('--energy_error', action='store_true', default=False)
p.add_argument('--exact_refined_adjoint', action='store_true', default=False)
    
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

# The boundary condition object
bcs = TMR.BoundaryConditions()

# The radius of the disk
R = 100.0

# The side-length of the rectangular plate
L = 200.0

if case == 'disk':
    if functional == 'integral':
        exact_functional = (1.0/12)*np.pi*R**4 # for f = 1 - (r/R)**2
    elif functional == 'displacement':
        exact_functional = 0.25*R**2
    elif functional == 'aggregate':
        for n in [10000, 100000, 1000000, 10000000]:
            approx = get_disk_aggregate(ksweight, R, n=n)
            if comm.rank == 0:
                print('%10d %25.16e'%(n, approx))
            exact_functional = approx

    geo = TMR.LoadModel('2d-disk.stp')
    verts = geo.getVertices()
    edges = geo.getEdges()
    faces = geo.getFaces()

    # Set the attributes
    verts[1].setAttribute('midpoint')
    for index in [0, 2, 3, 4]:
        verts[index].setAttribute('clamped')
    for index in [2, 4, 6, 7]:
        edges[index].setAttribute('clamped')

    # Set the boundary conditions
    bcs.addBoundaryCondition('clamped', [0])
elif case == 'square':
    if functional == 'integral':
        exact_functional = -2*(L/np.pi)**4
    elif functional == 'displacement':
        exact_functional = 2.946854131254945e+03
    elif functional == 'aggregate':
        raise ValueError('No exact solution for square aggregate')

    # Load the geometry model
    geo = TMR.LoadModel('2d-square.stp')
    verts = geo.getVertices()
    edges = geo.getEdges()
    faces = geo.getFaces()

    # Set the attributes
    verts[1].setAttribute('midpoint')
    for index in [0, 3, 2, 5, 4, 6, 7, 8]:
        verts[index].setAttribute('clamped')
    for index in [2, 3, 5, 6, 7, 8, 10, 11]:
        edges[index].setAttribute('clamped')

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
mesh.writeToVTK('mesh.vtk')

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
descript = 'poisson_unstructured'
if structured:
    descript = 'poisson_structured'
if args.uniform_refinement:
    descript += '_uniform'
descript += '_' + functional

# Create the log file and write out the header
log_fp = open('results/%s_order%d_%s.dat'%(case, order, descript), 'w')
s = 'Variables = iter, nelems, nnodes, fval, fcorr, abs_err, adjoint_corr, '
s += 'exact, fval_error, fval_corr_error, '
s += 'fval_effectivity, indicator_effectivity\n'
log_fp.write(s)

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
    if functional == 'displacement':
        res = get_midpoint_vector(comm, forest,
                                  assembler, 'midpoint')
        fval = ans.dot(res)
    elif functional == 'integral':
        direction = [1.0, 0.0, 0.0]
        func = functions.DisplacementIntegral(assembler, direction)
        fval = assembler.evalFunctions([func])[0]
    elif functional == 'aggregate':
        direction = [1.0/R**2, 0.0, 0.0]
        func = functions.KSDisplacement(assembler,
                                        ksweight, direction)
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
            if functional == 'displacement':
                adjoint_rhs = get_midpoint_vector(comm, forest_refined,
                                                  assembler_refined, 'midpoint')
                fval_refined = ans_interp.dot(adjoint_rhs)
            else:
                if functional == 'integral':
                    direction = [1.0, 0.0, 0.0]
                    func_refined = functions.DisplacementIntegral(assembler_refined,
                                                                  direction)

                elif functional == 'aggregate':
                    direction = [1.0/R**2, 0.0, 0.0]
                    func_refined = functions.KSDisplacement(assembler_refined,
                                                            ksweight, direction)

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
            if functional != 'displacement':
                assembler.evalSVSens(func, res)
            else:
                res = get_midpoint_vector(comm, forest,
                                          assembler, 'midpoint')
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
            if functional == 'displacement':
                adjoint_rhs = get_midpoint_vector(comm, forest_refined,
                                                  assembler_refined, 'midpoint')
                fval_refined = ans_refined.dot(adjoint_rhs)
            else:
                if functional == 'integral':
                    direction = [1.0, 0.0, 0.0]
                    func_refined = functions.DisplacementIntegral(assembler_refined,
                                                                  direction)

                elif functional == 'aggregate':
                    direction = [1.0/R**2, 0.0, 0.0]
                    func_refined = functions.KSDisplacement(assembler_refined,
                                                            ksweight, direction)

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
        np.savetxt('results/%s_data%d.txt'%(case, k), data)

    # Print out the error estimate
    assembler.setDesignVars(error)
    f5.writeToFile('results/error%02d.f5'%(k))

    # Perform the refinement
    if args.uniform_refinement:
        forest.refine()
        forest.balance(1)
        forest.repartition()
    elif k < steps-1:
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
        forest.balance(1)
        forest.repartition()
