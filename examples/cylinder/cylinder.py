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
        stiff = constitutive.isoFSDT(rho, E, nu, kcorr, ys, t, tnum,
                                     tmin, tmax)
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

def createRefined(forest, bcs, order=2):
    new_forest = forest.duplicate()
    new_forest.refine()
    new_forest.balance(1)
    creator = CreateMe(bcs)
    return creator.createTACS(order, new_forest)

def createProblem(forest, bcs, ordering, order=2, nlevels=2):
    # Create the forest
    forests = []
    assemblers = []

    # Create the trees, rebalance the elements and repartition
    forest.balance(1)
    forest.repartition()
    forests.append(forest)

    # Make the creator class
    creator = CreateMe(bcs)
    assemblers.append(creator.createTACS(order, forest, ordering))

    if order == 3:
        forest = forests[-1].duplicate()
        forest.balance(1)
        forests.append(forest)

        # Make the creator class
        creator = CreateMe(bcs)
        assemblers.append(creator.createTACS(2, forest, ordering))

    for i in xrange(nlevels-1):
        forest = forests[-1].coarsen()
        forest.balance(1)
        forest.repartition()
        forests.append(forest)

        # Make the creator class
        creator = CreateMe(bcs)
        assemblers.append(creator.createTACS(2, forest, ordering))

    # Create the multigrid object
    mg = TMR.createMg(assemblers, forests,
                      use_coarse_direct_solve=True, 
                      use_chebyshev_smoother=False)

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
args = p.parse_args()

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

for k in range(steps):
    # Create the topology problem
    assembler, mg = createProblem(forest, bcs, ordering, 
                                  order=order, nlevels=depth+k+1)

    # Add the surface traction
    aux = addFaceTraction(order, assembler, load)
    assembler.setAuxElements(aux)

    # Create the assembler object
    res = assembler.createVec()
    ans = assembler.createVec()
    mg.assembleJacobian(1.0, 0.0, 0.0, res)

    # Factor the matrix
    mg.factor()
    gmres = TACS.KSM(mg.getMat(), mg, 50, isFlexible=1)
    gmres.setMonitor(comm, freq=10)
    gmres.solve(res, ans)
    ans.scale(-1.0)

    # Set the variables
    assembler.setVariables(ans)

    # Output for visualization
    flag = (TACS.ToFH5.NODES |
            TACS.ToFH5.DISPLACEMENTS |
            TACS.ToFH5.EXTRAS)
    f5 = TACS.ToFH5(assembler, TACS.PY_SHELL, flag)
    f5.writeToFile('solution%02d.f5'%(k))

    # Create the function
    func = functions.KSFailure(assembler, ksweight)
    func.setKSFailureType('continuous')
    fval = assembler.evalFunctions([func])

    if False:
        # Compute the strain energy error estimate
        err_est, error = TMR.strainEnergyError(assembler, forest)
    else:
        # Compute the adjoint
        assembler.evalSVSens(func, res)
        
        # Compute the adjoint solution
        adjoint = assembler.createVec()
        gmres.solve(res, adjoint)
        adjoint.scale(-1.0)
        
        # Create the refined mesh
        refined = createRefined(forest, bcs, order=order)
        
        # Compute the adjoint and use adjoint-based refinement
        err_est, func_corr, error = \
            TMR.adjointError(assembler, refined, adjoint, forest)

    # Print the error estimate
    if comm.rank == 0:
        print 'estimate = ', err_est

    # Compute the refinement from the error estimate
    nbins = 30
    low = -15
    high = 0
    bounds = 10**np.linspace(low, high, nbins+1)
    bins = np.zeros(nbins+2, dtype=np.int)

    # Compute the mean and standard deviations of the log(error)
    ntotal = comm.allreduce(assembler.getNumElements(), op=MPI.SUM)
    mean = comm.allreduce(np.sum(np.log(error)), op=MPI.SUM)
    mean /= ntotal

    # Compute the standard deviation
    stddev = comm.allreduce(np.sum((np.log(error) - mean)**2), op=MPI.SUM)
    stddev = np.sqrt(stddev/(ntotal-1))
    
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
        print '%10s  %10s  %12s  %12s'%(
            'low', 'high', 'bins', 'percentage')
        print '%10.2e  %10s  %12d  %12.2f'%(
            bounds[-1], ' ', bins[-1], 100.0*bins[-1]/total)
        for i in range(nbins-1, -1, -1):
            print '%10.2e  %10.2e  %12d  %12.2f'%(
                bounds[i], bounds[i+1], bins[i+1], 100.0*bins[i+1]/total)
        print '%10s  %10.2e  %12d  %12.2f'%(
            ' ', bounds[0], bins[0], 100.0*bins[0]/total)

    # Print out the error estimate
    assembler.setDesignVars(error)
    f5.writeToFile('error%02d.f5'%(k))

    # Perform the refinement
    if k < steps-1:
        # Compute the refinement
        refine = np.zeros(len(error), dtype=np.intc)
        for i in range(len(error)):
            if np.log(error[i]) > mean:
                refine[i] = 1

        # Refine the forest
        forest.refine(refine)
