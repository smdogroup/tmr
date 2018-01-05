from mpi4py import MPI
from tmr import TMR
from tacs import TACS, elements, constitutive, functions
import numpy as np
import argparse
import os

class CreateMe(TMR.OctCreator):
    def __init__(self, bcs):
        TMR.OctTopoCreator.__init__(bcs)

    def createElement(self, order, octant):
        '''Create the element'''
        rho = 0.1015 # lbs/in
        E = 100e6 # 100,000 ksi
        nu = 0.3
        stiff = constitutive.isoSolidStiff(rho, E, nu)
        elem = elements.Solid(order, stiff)
        return elem

def addFaceTraction(order, forest, attr, assembler, tr):
    trac = []
    for findex in range(6):
        trac.append(elements.Traction3D(order, findex, tr[0], tr[1], tr[2]))

    # Retrieve octants from the forest
    octants = forest.getOctants()
    face_octs = forest.getOctsWithAttribute(attr)
    aux = TACS.AuxElements()

    for i in range(len(face_octs)):
        index = octants.findIndex(face_octs[i])
        if index is not None:
            aux.addElement(index, trac[face_octs[i].tag])

    return aux

def createRefined(forest, bcs, order=2):
    refined = forest.refine()
    creator = CreateMe(bcs)
    return creator.createTACS(order, refined)

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
    mg = TMR.createMg(assemblers, forests, use_coarse_direct_solve=True, 
        use_chebyshev_smoother=False)

    return assemblers[0], mg

# Set the communicator
comm = MPI.COMM_WORLD

# Create an argument parser to read in arguments from the commnad line
p = argparse.ArgumentParser()
p.add_argument('--htarget', type=float, default=0.15)
p.add_argument('--ordering', type=str, default='natural')
args = p.parse_args()

# Set the type of ordering to use for this problem
ordering = args.ordering
ordering = ordering.lower()

# Set the filename
filename = 'crank.stp'

# Set the value of the target length scale in the mesh
htarget = args.htarget

# Load the geometry model
geo = TMR.LoadModel(filename)

# Set the source/target meshes
faces = geo.getFaces()
vols = geo.getVolumes()
faces[7].setSource(vols[0], faces[6])

# Set the attributes
v0 = None
faces[4].setAttribute('fixed')
for i in range(faces[4].getNumEdgeLoops()):
    eloop = faces[4].getEdgeLoop(i)
    edges, d = eloop.getEdgeLoop()
    for e in edges:
        e.setAttribute('fixed')
        v1, v2 = e.getVertices()
        if v0 is None:
            v1.setAttribute('fully fixed')
            v0 = v1
        else:
            v1.setAttribute('fixed')
            v2.setAttribute('fixed')

# Set the loads
faces[5].setAttribute('load')

# Create the new mesh
mesh = TMR.Mesh(comm, geo)

# Set the meshing options
opts = TMR.MeshOptions()
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
forest = TMR.OctForest(comm)
forest.setTopology(topo)
forest.createTrees(depth)

# Set the boundary conditions for the problem
bcs = TMR.BoundaryConditions()
bcs.addBoundaryCondition('fixed', [0, 1])
bcs.addBoundaryCondition('fully fixed')

# Set the ordering to use
if ordering == 'rcm':
    ordering = TACS.PY_RCM_ORDER
elif ordering == 'multicolor':
    ordering = TACS.PY_MULTICOLOR_ORDER
else:
    ordering = TACS.PY_NATURAL_ORDER

order = 3
niters = 3
for k in range(niters):
    # Create the topology problem
    assembler, mg = createProblem(forest, bcs, ordering, 
                                  order=order, nlevels=depth+1+k)

    # Computet the surface traction magnitude
    diameter = 1.0
    circ = np.pi*diameter
    Area = 0.25*circ
    ty = 5000.0/Area
    tr = [0.0, ty, 0.0]

    # Add the surface traction
    aux = addFaceTraction(order, forest, 'load', assembler, tr)
    assembler.setAuxElements(aux)

    # Create the assembler object
    res = assembler.createVec()
    ans = assembler.createVec()
    mg.assembleJacobian(1.0, 0.0, 0.0, res)

    # Factor the matrix
    mg.factor()
    gmres = TACS.KSM(mg.getMat(), mg, 50, isFlexible=1)
    gmres.setMonitor(comm, freq=1)
    gmres.solve(res, ans)
    ans.scale(-1.0)

    # Set the variables
    assembler.setVariables(ans)

    # Output for visualization
    flag = (TACS.ToFH5.NODES |
            TACS.ToFH5.DISPLACEMENTS |
            TACS.ToFH5.STRESSES)
    f5 = TACS.ToFH5(assembler, TACS.PY_SOLID, flag)
    f5.writeToFile('beam%d.f5'%(k))

    ksweight = 10
    func = functions.KSFailure(assembler, ksweight)
    func.setKSFailureType('continuous')

    func = functions.Compliance(assembler)
    fval = assembler.evalFunctions([func])

    if k < niters:
        if True:
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

        # Compute the cutoff
        bsum = bins[-1]
        cutoff = bounds[-1]
        for i in range(len(bounds), -1, -1):
            if bsum > 0.25*total:
                cutoff = bounds[i]
                break
            bsum += bins[i]

        # Print out the result
        if comm.rank == 0:
            print '%10s  %10s  %12s  %12s'%(
                'low', 'high', 'bins', 'percentage')
            print '%10.2e  %10s  %12d  %12.2f'%(
                bounds[-1], ' ', bins[-1], 1.0*bins[-1]/total)
            for k in range(nbins-1, -1, -1):
                print '%10.2e  %10.2e  %12d  %12.2f'%(
                    bounds[k], bounds[k+1], bins[k+1], 1.0*bins[k]/total)
            print '%10s  %10.2e  %12d  %12.2f'%(
                ' ', bounds[0], bins[0], 1.0*bins[0]/total)
            print 'cutoff:  %15.2e'%(cutoff)

        # Compute the refinement
        refine = np.zeros(len(error), dtype=np.intc)
        for i in range(len(error)):
            if error[i] > cutoff:
                refine[i] = 1

        # Refine the forest
        forest.refine(refine)
