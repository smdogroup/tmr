from __future__ import print_function
from mpi4py import MPI
from tmr import TMR
from tacs import TACS, elements, constitutive, functions
import numpy as np
import argparse
import os

class CreateMe(TMR.OctCreator):
    def __init__(self, bcs):
        TMR.OctTopoCreator.__init__(bcs)

        # Create the stiffness object
        self.props = constitutive.MaterialProperties(rho=2570.0, E=70e9, nu=0.3, ys=350e6)
        self.stiff = constitutive.SolidConstitutive(self.props)

        # Set up the basis function
        self.model = elements.LinearElasticity3D(self.stiff)

    def createElement(self, order, octant):
        """Create the element"""
        if order == 2:
            basis = elements.LinearHexaBasis()
        elif order == 3:
            basis = elements.QuadraticHexaBasis()
        elif order == 4:
            basis = elements.CubicHexaBasis()

        return elements.Element3D(self.model, basis)

def addFaceTraction(order, forest, attr, assembler, tr):
    vpn = assembler.getVarsPerNode()
    if order == 2:
        basis = elements.LinearHexaBasis()
    elif order == 3:
        basis = elements.QuadraticHexaBasis()
    elif order == 4:
        basis = elements.CubicHexaBasis()

    trac = []
    for findex in range(6):
        trac.append(elements.Traction3D(vpn, findex, basis, tr))

    # Retrieve octants from the forest
    octants = forest.getOctants()
    face_octs = forest.getOctsWithName(attr)
    aux = TACS.AuxElements()

    for i in range(len(face_octs)):
        aux.addElement(face_octs[i].tag, trac[face_octs[i].info])

    return aux

def createRefined(forest, bcs, pttype=TMR.GAUSS_LOBATTO_POINTS):
    refined = forest.duplicate()
    refined.setMeshOrder(forest.getMeshOrder()+1, pttype)
    refined.balance(1)
    creator = CreateMe(bcs)
    return refined, creator.createTACS(refined)

def createProblem(forest, bcs, ordering, mesh_order=2, nlevels=2,
                  pttype=TMR.GAUSS_LOBATTO_POINTS):
    # Create the forest
    forests = []
    assemblers = []

    # Create the trees, rebalance the elements and repartition
    forest.balance(1)
    forest.setMeshOrder(mesh_order, pttype)
    forest.repartition()
    forests.append(forest)

    # Make the creator class
    creator = CreateMe(bcs)
    assemblers.append(creator.createTACS(forest, ordering))

    while mesh_order > 2:
        mesh_order = mesh_order-1
        forest = forests[-1].duplicate()
        forest.balance(1)
        forest.setMeshOrder(mesh_order, pttype)
        forests.append(forest)

        # Make the creator class
        creator = CreateMe(bcs)
        assemblers.append(creator.createTACS(forest, ordering))

    for i in range(nlevels-1):
        forest = forests[-1].coarsen()
        forest.setMeshOrder(2, pttype)
        forest.balance(1)
        forest.repartition()
        forests.append(forest)

        # Make the creator class
        creator = CreateMe(bcs)
        assemblers.append(creator.createTACS(forest, ordering))

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
p.add_argument('--order', type=int, default=3)
args = p.parse_args()

# Set the element order
order = args.order
if order < 2:
    order = 2
elif order > 5:
    order = 5

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

# Set the names
v0 = None
faces[4].setName('fixed')
for i in range(faces[4].getNumEdgeLoops()):
    eloop = faces[4].getEdgeLoop(i)
    edges, d = eloop.getEdgeLoop()
    for e in edges:
        e.setName('fixed')
        v1, v2 = e.getVertices()
        if v0 is None:
            v1.setName('fully fixed')
            v2.setName('fixed')
            v0 = v1
        else:
            v1.setName('fixed')
            v2.setName('fixed')

# Set the loads
faces[5].setName('load')

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
forest = TMR.OctForest(comm)
forest.setTopology(topo)
forest.createTrees(2)

# Set the boundary conditions for the problem
bcs = TMR.BoundaryConditions()
bcs.addBoundaryCondition('fixed', [0, 1, 2])
bcs.addBoundaryCondition('fully fixed')

# Set the ordering to use
if ordering == 'rcm':
    ordering = TACS.RCM_ORDER
elif ordering == 'multicolor':
    ordering = TACS.MULTICOLOR_ORDER
else:
    ordering = TACS.NATURAL_ORDER

niters = 1
for k in range(niters):
    t = MPI.Wtime()
    # Create the topology problem
    assembler, mg = createProblem(forest, bcs, ordering,
                                  mesh_order=order, nlevels=3+k)
    if comm.rank == 0:
        print('Creating TACS assembler objects', MPI.Wtime() - t)

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
    t = MPI.Wtime()
    mg.assembleJacobian(1.0, 0.0, 0.0, res)
    if comm.rank == 0:
        print('Assembling the Jacobians', MPI.Wtime() - t)

    # Factor the matrix
    t = MPI.Wtime()
    mg.factor()
    if comm.rank == 0:
        print('Factoring the Jacobians', MPI.Wtime() - t)
    gmres = TACS.KSM(mg.getMat(), mg, 100, isFlexible=1)
    gmres.setMonitor(comm, freq=1)
    gmres.solve(res, ans)
    ans.scale(-1.0)

    # Set the variables
    assembler.setVariables(ans)

    # Output for visualization
    flag = (TACS.OUTPUT_CONNECTIVITY |
            TACS.OUTPUT_NODES |
            TACS.OUTPUT_DISPLACEMENTS |
            TACS.OUTPUT_STRESSES)
    f5 = TACS.ToFH5(assembler, TACS.SOLID_ELEMENT, flag)
    f5.writeToFile('crank%d.f5'%(k))

    if order >= 4:
        if comm.rank == 0:
            print('Cannot perform adaptive refinement with order >= 4')
        break

    ksweight = 10
    func = functions.KSFailure(assembler, ksweight)
    func.setKSFailureType('continuous')

    func = functions.Compliance(assembler)
    fval = assembler.evalFunctions([func])

    if k < niters:
        # Create the refined mesh
        forest_refined, assembler_refined = createRefined(forest, bcs)

        if True:
            # Compute the strain energy error estimate
            err_est, error = TMR.strainEnergyError(forest, assembler,
                                                   forest_refined,
                                                   assembler_refined)
        else:
            # Compute the adjoint
            assembler.evalSVSens(func, res)

            # Compute the adjoint solution
            adjoint = assembler.createVec()
            gmres.solve(res, adjoint)
            adjoint.scale(-1.0)

            # Compute the adjoint and use adjoint-based refinement
            err_est, func_corr, error = \
                     TMR.adjointError(forest, assembler,
                                      forest_refined, assembler_refined, adjoint)

        # Print the error estimate
        if comm.rank == 0:
            print('estimate = ', err_est)

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
            print('%10s  %10s  %12s  %12s'%(
                'low', 'high', 'bins', 'percentage'))
            print('%10.2e  %10s  %12d  %12.2f'%(
                bounds[-1], ' ', bins[-1], 1.0*bins[-1]/total))
            for k in range(nbins-1, -1, -1):
                print('%10.2e  %10.2e  %12d  %12.2f'%(
                    bounds[k], bounds[k+1], bins[k+1],
                    100.0*bins[k]/total))
            print('%10s  %10.2e  %12d  %12.2f'%(
                ' ', bounds[0], bins[0], 100.0*bins[0]/total))
            print('cutoff:  %15.2e'%(cutoff))

        # Compute the refinement
        refine = np.zeros(len(error), dtype=np.intc)
        for i in range(len(error)):
            if error[i] > cutoff:
                refine[i] = 1

        # Refine the forest
        forest.refine(refine)
