from mpi4py import MPI
from tmr import TMR
from tacs import TACS, elements, constitutive, functions
import numpy as np
import argparse
import os

class CreateMe(TMR.OctTopoCreator):
    def __init__(self, bcs, filt):
        TMR.OctTopoCreator.__init__(bcs, filt)

    def createElement(self, order, octant, index, weights):
        '''Create the element'''
        rho = 0.1015 # lbs/in
        E = 100e6 # 100,000 ksi
        nu = 0.3
        stiff = constitutive.isoSolidStiff(rho, E, nu)
        elem = elements.Solid(2, stiff)
        return elem

def addVertexLoad(comm, order, forest, attr, assembler, F):
    # Retrieve octants from the forest
    octants = forest.getOctants()
    node_octs = forest.getNodesWithAttribute(attr)
    force = assembler.createVec()
    f_array = force.getArray()
    node_range = forest.getNodeRange()
    mpi_rank = comm.Get_rank()
    for i in range(len(node_octs)):
        if (node_octs[i].tag >= node_range[mpi_rank]) and \
               (node_octs[i].tag < node_range[mpi_rank+1]): 
            index = node_octs[i].tag-node_range[mpi_rank]
            
            f_array[3*index] -= F[0]
            f_array[3*index+1] -= F[1]
            f_array[3*index+2] -= F[2]
    return force

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

def createTopoProblem(forest, ordering, order=2, nlevels=2):
    # Create the forest
    forests = []
    filters = []
    assemblers = []
    varmaps = []
    vecindices = []

    # Create the trees, rebalance the elements and repartition
    forest.balance(1)
    forest.repartition()
    forests.append(forest)

    # Create the filter
    filtr = forest.coarsen()
    filtr.balance(1)
    filters.append(filtr)

    # Make the creator class
    creator = CreateMe(bcs, filters[-1])
    assemblers.append(creator.createTACS(order, forest, ordering))
    varmaps.append(creator.getMap())
    vecindices.append(creator.getIndices())

    for i in xrange(nlevels-1):
        forest = forests[-1].coarsen()
        forest.balance(1)
        forest.repartition()
        forests.append(forest)

        # Create the filter
        filtr = forest.coarsen()
        filtr.balance(1)
        filters.append(filtr)

        # Make the creator class
        creator = CreateMe(bcs, filters[-1])
        assemblers.append(creator.createTACS(order, forest, ordering))
        varmaps.append(creator.getMap())
        vecindices.append(creator.getIndices())

    # Create the multigrid object
    mg = TMR.createMg(assemblers, forests)

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
faces[4].setAttribute('fixed')
for i in range(faces[4].getNumEdgeLoops()):
    eloop = faces[4].getEdgeLoop(i)
    edges, d = eloop.getEdgeLoop()
    for e in edges:
        e.setAttribute('fixed')
        v1, v2 = e.getVertices()
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
nlevs = 2
forest = TMR.OctForest(comm)
forest.setTopology(topo)
forest.createTrees(nlevs)

# Set the boundary conditions for the problem
bcs = TMR.BoundaryConditions()
bcs.addBoundaryCondition('fixed')

target_rel_err = 1e-4

# Set the ordering to use
if ordering == 'rcm':
    ordering = TACS.PY_RCM_ORDER
elif ordering == 'multicolor':
    ordering = TACS.PY_MULTICOLOR_ORDER
else:
    ordering = TACS.PY_NATURAL_ORDER

niters = 3
for k in range(niters):
    # Create the topology problem
    assembler, mg = createTopoProblem(forest, ordering, nlevels=nlevs+1)

    # Computet the surface traction magnitude
    diameter = 1.0
    circ = np.pi*diameter
    Area = 0.25*circ
    ty = 5000.0/Area
    tr = [0.0, ty, 0.0]

    # Add the surface traction
    order = 2
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

    if k < niters-1:
        # Get the number of elements
        nelems = assembler.getNumElements()
        nelems = comm.allreduce(nelems, op=MPI.SUM)

        # Compute the strain energy
        obj = -0.5*res.dot(ans)

        # Set the target error
        target_err = target_rel_err*obj/nelems

        # Create the refined mesh
        refined = forest.duplicate()
        refined.refine()
        refined.balance(1)
        
        # Create the filter
        filtr = refined.coarsen()
        filtr.balance(1)

        # Make the creator class
        creator = CreateMe(bcs, filtr)
        assembler_refined = creator.createTACS(order, refined)

        # Do an adjoint-based refinement of the mesh
        TMR.adjointRefine(assembler, assembler_refined, ans, forest, target_err)

        # Code for strain energy based refinement
        # TMR.strainEnergyRefine(assembler, forest, target_err)

