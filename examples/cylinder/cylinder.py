from __future__ import print_function
from mpi4py import MPI
from tmr import TMR
from tacs import TACS, elements, constitutive, functions
import numpy as np
import argparse
import os

# Parameters that define the geometry of the cylinder
t = 1.0 # 1 mm thickness
L = 100.0 # 100 mm length
R = 100.0/np.pi # Radius of the cylinder (in mm)

# Set the load to apply to the cylinder
load = 3.0 # 3 MPa
alpha = 4.0/R # Set the value of alpha
beta = 3.0*np.pi/L # Set the value of beta

# Set the material properties
rho = 2700.0e-9 # kg/mm^3
E = 70e3 # 70 GPa
nu = 0.3 # Poisson's ratio
ys = 350.0 # 350 MPa
kcorr = 5.0/6.0 # The shear correction factor

class CreateMe(TMR.QuadCreator):
    def __init__(self, bcs):
        TMR.QuadCreator.__init__(bcs)

    def createElement(self, order, quad):
        '''Create the element'''
        tmin = 0.0
        tmax = 1e6
        tnum = quad.tag
        stiff = constitutive.isoFSDT(rho, E, nu, kcorr, ys, t,
                                     tnum, tmin, tmax)
        stiff.setRefAxis(np.array([0.0, 0.0, 1.0]))
        elem = elements.MITCShell(order, stiff)
        return elem

def addFaceTraction(order, assembler, load):
    # Create the surface traction
    aux = TACS.AuxElements()

    # Get the element node locations
    nelems = assembler.getNumElements()
    for i in range(nelems):
        # Get the information about the given element
        elem, Xpt, vars0, dvars, ddvars = assembler.getElementData(i)

        # Loop over the nodes and create the traction forces in the
        # x/y/z directions
        nnodes = order*order
        tx = np.zeros(nnodes, dtype=TACS.dtype)
        ty = np.zeros(nnodes, dtype=TACS.dtype)
        tz = np.zeros(nnodes, dtype=TACS.dtype)

        # Set the components of the surface traction
        for j in range(nnodes):
            x = Xpt[3*j]
            y = Xpt[3*j+1]
            z = Xpt[3*j+2]
            theta = -R*np.arctan2(y, x)
            p = -load*np.sin(beta*z)*np.sin(alpha*theta)

            tx[j] = p*Xpt[3*j]/R
            ty[j] = p*Xpt[3*j+1]/R
            tz[j] = 0.0

        trac = elements.ShellTraction(order, tx, ty, tz)
        aux.addElement(i, trac)
        
    return aux

def createRefined(forest, bcs, pttype=TMR.UNIFORM_POINTS):
    new_forest = forest.duplicate()
    new_forest.setMeshOrder(forest.getMeshOrder()+1, pttype)
    creator = CreateMe(bcs)
    return new_forest, creator.createTACS(new_forest)

def createProblem(forest, bcs, ordering, order=2, nlevels=2,
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
    creator = CreateMe(bcs)
    assemblers.append(creator.createTACS(forest, ordering))

    while order > 2:
        order = order-1
        forest = forests[-1].duplicate()
        forest.setMeshOrder(order, pttype)
        forest.balance(1)
        forests.append(forest)

        # Make the creator class
        creator = CreateMe(bcs)
        assemblers.append(creator.createTACS(forest, ordering))

    for i in range(nlevels-1):
        forest = forests[-1].coarsen()
        forest.setMeshOrder(2, pttype)
        forest.balance(1)
        forests.append(forest)

        # Make the creator class
        creator = CreateMe(bcs)
        assemblers.append(creator.createTACS(forest, ordering))

    # Create the multigrid object
    mg = TMR.createMg(assemblers, forests, omega=0.5)

    return assemblers[0], mg

# Set the communicator
comm = MPI.COMM_WORLD

# Create an argument parser to read in arguments from the commnad line
p = argparse.ArgumentParser()
p.add_argument('--steps', type=int, default=5)
p.add_argument('--htarget', type=float, default=10.0)
p.add_argument('--ordering', type=str, default='multicolor')
p.add_argument('--order', type=int, default=2)
p.add_argument('--ksweight', type=float, default=100.0)
p.add_argument('--structured', type=int, default=0)
p.add_argument('--adjoint_error_est', action='store_true', default=False)
args = p.parse_args()

# Set the type of ordering to use for this problem
ordering = args.ordering
ordering = ordering.lower()

# Set whether to use the adjoint error estimate or not
adjoint_error_est = args.adjoint_error_est

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

# Load the geometry model
geo = TMR.LoadModel('cylinder.stp')
verts = geo.getVertices()
edges = geo.getEdges()
faces = geo.getFaces()

# Create the simplified geometry with only faces
geo = TMR.Model(verts, edges, [faces[0]])

# Set the boundary conditions
verts[0].setAttribute('Clamped')
verts[1].setAttribute('Restrained')
edges[0].setAttribute('Restrained')
edges[2].setAttribute('Restrained')

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

# Set the boundary conditions for the problem
bcs = TMR.BoundaryConditions()
bcs.addBoundaryCondition('Clamped', [0, 1, 2, 5])
bcs.addBoundaryCondition('Restrained', [0, 1, 5])

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
if order == 2:
    log_fp = open('cylinder_refine2nd.dat', 'w')
elif order == 3:
    log_fp = open('cylinder_refine3rd.dat', 'w')
elif order == 4:
    log_fp = open('cylinder_refine4th.dat', 'w')

# Write the first line to the file
log_fp.write('Variables = iter, nelems, nnodes, fval, fcorr, abs_err\n')

for k in range(steps):
    # Create the topology problem
    nlevs = min(3, depth+k+1)
    assembler, mg = createProblem(forest, bcs, ordering, 
                                  order=order, nlevels=nlevs)

    # Add the surface traction
    aux = addFaceTraction(order, assembler, load)
    assembler.setAuxElements(aux)

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
    
    gmres = TACS.KSM(mat, pc, 100, isFlexible=0)
    gmres.setMonitor(comm, freq=1)
    gmres.solve(res, ans)
    ans.scale(-1.0)

    # Set the variables
    assembler.setVariables(ans)

    # Output for visualization
    flag = (TACS.ToFH5.NODES |
            TACS.ToFH5.DISPLACEMENTS |
            TACS.ToFH5.STRAINS |
            TACS.ToFH5.EXTRAS)
    f5 = TACS.ToFH5(assembler, TACS.PY_SHELL, flag)
    f5.writeToFile('results/solution%02d.f5'%(k))

    # Create the function
    func = functions.KSFailure(assembler, ksweight)
    func.setKSFailureType('continuous')
    fval = assembler.evalFunctions([func])[0]

    # Create the refined mesh
    forest_refined, assembler_refined = createRefined(forest, bcs)

    # Find the refined solution
    TMR.computeReconSolution(forest, assembler,
        forest_refined, assembler_refined)

    f5_refine = TACS.ToFH5(assembler_refined, TACS.PY_SHELL, flag)
    f5_refine.writeToFile('results/solution_refined%02d.f5'%(k))

    if adjoint_error_est:
        # Compute the adjoint
        res.zeroEntries()
        assembler.evalSVSens(func, res)

        # Compute the adjoint solution
        adjoint = assembler.createVec()
        gmres.solve(res, adjoint)
        adjoint.scale(-1.0)

        # Compute the adjoint and use adjoint-based refinement
        err_est, func_corr, error = \
            TMR.adjointError(forest, assembler,
                             forest_refined, assembler_refined, adjoint)
    else:
        # Compute the strain energy error estimate
        func_corr = 0.0
        err_est, error = TMR.strainEnergyError(forest, assembler,
            forest_refined, assembler_refined)

    # Compute the refinement from the error estimate
    nbins = 60
    low = -10
    high = 2
    bounds = 10**np.linspace(low, high, nbins+1)
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

    # Write the log data to a file
    log_fp.write('%d %d %d %e %e %e\n'%(
        k, ntotal, nnodes, fval, fval + func_corr, err_est))
    log_fp.flush()               
    
    # Compute the bins
    for i in range(len(error)):
        if error[i] < bounds[0]:
            bins[0] += 1
        elif error[i] > bounds[-1]:
            bins[-1] += 1
        else:
            for j in range(len(bounds)-1):
                if (error[i] >= bounds[j] and
                    error[i] < bounds[j+1]):
                    bins[j+1] += 1

    # Compute the number of bins
    bins = comm.allreduce(bins, MPI.SUM)

    # Compute the sum of the bins
    total = np.sum(bins)
    
    # Print out the result
    if comm.rank == 0:
        print('fval      = ', fval)
        print('corr func = ', func_corr)
        print('estimate =  ', err_est)
        print('mean =      ', mean)
        print('stddev =    ', stddev)

        print('%10s   %10s   %12s   %12s'%(
            'low', 'high', 'bins', 'percentage'))
        print('%10.2e   %10s   %12d   %12.2f'%(
            bounds[-1], ' ', bins[-1], 100.0*bins[-1]/total))
        for i in range(nbins-1, -1, -1):
            print('%10.2e   %10.2e   %12d   %12.2f'%(
                bounds[i], bounds[i+1], bins[i+1], 100.0*bins[i+1]/total))
        print('%10s   %10.2e   %12d   %12.2f'%(
            ' ', bounds[0], bins[0], 100.0*bins[0]/total))

        # Set the data 
        data = np.zeros((60, 4))
        for i in range(nbins-1, -1, -1):
            data[i,0] = bounds[i]
            data[i,1] = bounds[i+1]
            data[i,2] = bins[i+1]
            data[i,3] = 100.0*bins[i+1]/total
        np.savetxt('results/cylinder_data%d.txt'%(k), data)

    # Print out the error estimate
    assembler.setDesignVars(error)
    f5.writeToFile('results/error%02d.f5'%(k))

    # Perform the refinement
    if k < steps-1:
        # The refinement array
        refine = np.zeros(len(error), dtype=np.intc)

        # Compute the target relative error
        element_target_error = target_rel_err*fval/ntotal
        log_elem_target_error = np.log(element_target_error)

        # Determine the cutoff values
        cutoff = 0.0
        bin_sum = 0
        for i in range(len(bins)-2, -1, -1):
            bin_sum += bins[i]
            if bin_sum > 0.15*ntotal:
                cutoff = bounds[i]
                break

        log_cutoff = np.log(cutoff)

        if comm.rank == 0:
            print('elem_target_error =        %15.3e'%(
                element_target_error))
            print('log10(elem_target_error) = %15.3e'%(
                np.log10(element_target_error)))
            print('log10(mean) =              %15.3e'%(
                np.log10(np.exp(mean))))
            print('log10(stddev) =            %15.3e'%(
                np.log10(np.exp(stddev))))
            print('log_elem_target_error =    %15.3e'%(
                log_elem_target_error))
            print('mean =                     %15.3e'%(mean))
            print('stddev =                   %15.3e'%(stddev))
            print('cutoff =                   %15.3e'%(cutoff))

        if log_elem_target_error < mean - stddev:
            if comm.rank == 0:
                print('-----------------------------------------------')
                print('First refinement phase')
                print('-----------------------------------------------')

            # Element target error is still too high. Adapt based
            # solely on decreasing the overall error
            for i, err in enumerate(error):
                # Compute the log of the error
                logerr = np.log(err)

                if logerr > log_cutoff:
                    refine[i] = 1
                # elif logerr < mean - 2*stddev:
                #     refine[i] = -1
        else:
            if comm.rank == 0:
                print('-----------------------------------------------')
                print('Entering final phase refinement')
                print('-----------------------------------------------')

            # Try to target the relative error
            for i, err in enumerate(error):
                # Compute the log of the error
                logerr = np.log(err)

                if logerr > log_elem_target_error:
                    refine[i] = 1
                elif logerr < log_elem_target_error - 2*stddev:
                    refine[i] = -1
            
        # Refine the forest
        forest.refine(refine)
