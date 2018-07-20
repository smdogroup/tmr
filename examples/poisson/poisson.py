from __future__ import print_function
from mpi4py import MPI
from tmr import TMR
from tacs import TACS, elements, functions
import numpy as np
import argparse
import os

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
    def __init__(self, bcs, case='disk'):
        self.case = case
        TMR.QuadCreator.__init__(bcs)
        return

    def createElement(self, order, quad):
        '''Create the element'''
        f = np.ones(order*order)
        elem = elements.PoissonQuad(order, f)
        return elem

def createRefined(case, forest, bcs, pttype=TMR.GAUSS_LOBATTO_POINTS):
    new_forest = forest.duplicate()
    new_forest.setMeshOrder(forest.getMeshOrder()+1, pttype)
    creator = CreateMe(bcs, case=case)
    return new_forest, creator.createTACS(new_forest)

def createProblem(case, forest, bcs, ordering, order=2, nlevels=2,
                  pttype=TMR.GAUSS_LOBATTO_POINTS):
    # Create the forest
    forests = []
    assemblers = []

    # Create the trees, rebalance the elements and repartition
    forest.balance(1)
    forest.setMeshOrder(order, pttype)
    forest.repartition()
    forests.append(forest)

    # Make the creator class
    creator = CreateMe(bcs, case=case)
    assemblers.append(creator.createTACS(forest, ordering))

    while order > 2:
        order = order-1
        forest = forests[-1].duplicate()
        forest.setMeshOrder(order, pttype)
        forest.balance(1)
        forests.append(forest)

        # Make the creator class
        creator = CreateMe(bcs, case=case)
        assemblers.append(creator.createTACS(forest, ordering))

    for i in range(nlevels-1):
        forest = forests[-1].coarsen()
        forest.setMeshOrder(2, pttype)
        forest.balance(1)
        forests.append(forest)

        # Make the creator class
        creator = CreateMe(bcs, case=case)
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
p.add_argument('--ksweight', type=float, default=1.0)
p.add_argument('--uniform_refinement', action='store_true', default=False)
p.add_argument('--structured', action='store_true', default=False)
p.add_argument('--energy_error', action='store_true', default=False)
args = p.parse_args()

# Set the case type
case = args.case

# Parameters that define the geometry of the cylinder
t = 1.0 # 1 mm thickness
L = 100.0 # 100 mm length
R = 100.0/np.pi # Radius of the cylinder (in mm)

# Set the load to apply to the cylinder
load = 3.0 # 3 MPa
alpha = 4.0/R # Set the value of alpha
beta = 3.0*np.pi/L # Set the value of beta

# Set the material properties
density = 2700.0e-9 # kg/mm^3
E = 70e3 # 70 GPa
nu = 0.3 # Poisson's ratio
ys = 350.0 # 350 MPa
kcorr = 5.0/6.0 # The shear correction factor

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

# The boundary condition object
bcs = TMR.BoundaryConditions()

# The radius of the disk
R = 100.0

# The side-length of the rectangular plate
L = 200.0

if case == 'disk':
    if functional == 'integral':
        exact_functional = (1.0/8)*np.pi*R**4
    elif functional == 'displacement':
        exact_functional = 0.25*R**2
    elif functional == 'aggregate':
        exact_functional = 2.502531024247042e+03
    
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
    find_exact = False

    if find_exact:
        exact_functional = 0.0
        for n in range(1, 10000, 2):
            for m in range(1, 10000, 2):
                func = 1
                if n % 4 == 3:
                    func *= -1
                if m % 4 == 3:
                    func *= -1
                exact_functional += \
                    func*16.0/((m**2 + n**2)*m*n*np.pi**2)
        exact_functional *= (L/np.pi)**2
    
    if functional == 'integral':
        exact_functional = (1.0/8)*np.pi*R**4
    elif functional == 'displacement':
        exact_functional = 0.25*R**2
    elif functional == 'aggregate':
        exact_functional = 2.502531024247042e+03

    exact_functional = 2.946854131254945e+03
    if disp_aggregation:
        exact_functional = 2.502531024247042e+03

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
log_fp = open('%s_order%d_%s.dat'%(case, order, descript), 'w')
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
        direction = [1.0, 0.0, 0.0]
        func = functions.KSDisplacement(assembler,
                                        ksweight, direction)
        fval = assembler.evalFunctions([func])[0]

    # This flag indicates whether to solve the adjoint exactly on the
    # next-finest mesh or not
    exact_refined_adjoint = False

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

        if exact_refined_adjoint:
            # Compute the reconstructed solution on the refined mesh
            ans_interp = assembler_refined.createVec()
            TMR.computeReconSolution(forest, assembler,
                forest_refined, assembler_refined, ans, ans_interp)

            # Set the interpolated solution on the fine mesh
            assembler_refined.setVariables(ans_interp)        

            # Compute the reconstructed adjoint solution on the refined mesh
            adjoint_interp = assembler_refined.createVec()
            TMR.computeReconSolution(forest, assembler,
                forest_refined, assembler_refined, adjoint, adjoint_interp)

            # Solve the linear system
            res_refined = assembler_refined.createVec()
            adjoint_refined = assembler_refined.createVec()

            # Compute the functional on the refined mesh
            if functional == 'displacement':
                res_refined = get_midpoint_vector(comm, forest_refined,
                                                  assembler_refined, 'midpoint')
                fval_refined = ans_interp.dot(res_refined)
            else:
                if functional == 'integral':
                    direction = [1.0, 0.0, 0.0]
                    func_refined = functions.DisplacementIntegral(assembler_refined,
                                                                  direction)

                elif functional == 'aggregate':
                    direction = [1.0, 0.0, 0.0]
                    func_refined = functions.KSDisplacement(assembler_refined,
                                                            ksweight, direction)

                # Evaluate the functional on the refined mesh
                fval_refined = assembler_refined.evalFunctions([func_refined])[0]
                assembler_refined.evalSVSens(func_refined, res_refined)

            # Assemble the adjoint equations
            mg.assembleJacobian(1.0, 0.0, 0.0, adjoint_refined)
            mg.factor()
            pc = mg
            mat = mg.getMat()

            # Create the GMRES object on the fine mesh
            gmres = TACS.KSM(mat, pc, 100, isFlexible=1)
            gmres.setMonitor(comm, freq=10)
            gmres.setTolerances(1e-14, 1e-30)
            gmres.solve(res_refined, adjoint_refined)
            adjoint_refined.scale(-1.0)

            # Compute the adjoint and use adjoint-based refinement
            err_est, adjoint_corr, error = TMR.adjointError(forest, assembler,
                forest_refined, assembler_refined, ans_interp, adjoint_refined)
        else:
            # Compute the solution on the refined mesh
            ans_refined = assembler_refined.createVec()
            TMR.computeReconSolution(forest, assembler,
                forest_refined, assembler_refined, ans, ans_refined)

            # Apply Dirichlet boundary conditions
            assembler_refined.setVariables(ans_refined)

            # Test orthogonality of the residual of the interpolated
            # solution on the fine mesh
            test_res_orthogonality = False
            if test_res_orthogonality:
                # Assemble the residual on the fine mesh from the
                # interpolated solution
                res_interp = assembler_refined.createVec()
                assembler_refined.assembleRes(res_interp)

                # Create a random vector on the coarse mesh
                vec = assembler.createVec()
                vec_interp = assembler_refined.createVec()
                for i in range(5):
                    # Compute the random displacement vector in the
                    # original space and apply Dirichlet bcs
                    vec.setRand(-1.0, 1.0)
                    assembler.applyBCs(vec)            
                    TMR.computeInterpSolution(forest, assembler,
                                              forest_refined, assembler_refined,
                                              vec, vec_interp)

                    # Compute the relative error
                    dprod = vec_interp.dot(res_interp)
                    norm = res_interp.norm()*vec_interp.norm()
                    if comm.rank == 0:
                        print('Orthogonality check: ', dprod/norm)

            # Compute the functional on the refined mesh
            if functional == 'displacement':
                res_refined = get_midpoint_vector(comm, forest_refined,
                                                  assembler_refined, 'midpoint')
                fval_refined = ans_refined.dot(res_refined)
            else:
                if functional == 'integral':
                    direction = [1.0, 0.0, 0.0]
                    func_refined = functions.DisplacementIntegral(assembler_refined,
                                                                  direction)

                elif functional == 'aggregate':
                    direction = [1.0, 0.0, 0.0]
                    func_refined = functions.KSDisplacement(assembler_refined,
                                                            ksweight, direction)

                # Evaluate the functional on the refined mesh
                fval_refined = assembler_refined.evalFunctions([func_refined])[0]

            # Approximate the difference between the refined adjoint and the
            # adjoint on the current mesh
            adjoint_refined = assembler_refined.createVec()
            TMR.computeReconSolution(forest, assembler,
                forest_refined, assembler_refined, adjoint, adjoint_refined)

            # Compute the adjoint and use adjoint-based refinement
            err_est, adjoint_corr, error = TMR.adjointError(forest, assembler,
                forest_refined, assembler_refined, ans_refined, adjoint_refined)

        # Compute the refined function value
        fval_corr = fval_refined + adjoint_corr

    # f5_refine = TACS.ToFH5(assembler_refined, TACS.PY_SHELL, flag)
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
    elif k < steps-1:
        # The refinement array
        refine = np.zeros(len(error), dtype=np.intc)

        # Compute the target relative error
        element_target_error = target_rel_err*fval/ntotal
        log_elem_target_error = np.log(element_target_error)

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
