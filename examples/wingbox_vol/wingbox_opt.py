from mpi4py import MPI
from tmr import TMR
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

class CreateMe(TMR.OctTopoCreator):
    def __init__(self, bcs, filt):
        TMR.OctTopoCreator.__init__(bcs, filt)

    def createElement(self, order, octant, index, weights):
        '''Create the element'''
        rho = 2600.0
        E = 70e9
        nu = 0.3
        stiff = TMR.OctStiffness(rho, E, nu, index, weights, q=5.0)
        elem = elements.Solid(2, stiff)
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

def createTopoProblem(forest, order=2, nlevels=2):
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
    assemblers.append(creator.createTACS(order, forest))
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
        assemblers.append(creator.createTACS(order, forest))
        varmaps.append(creator.getMap())
        vecindices.append(creator.getIndices())

    # Create the multigrid object
    mg = TMR.createMg(assemblers, forests)

    # Create the topology optimization problem
    problem = TMR.TopoProblem(assemblers, filters, varmaps, vecindices, mg)

    return assemblers[0], problem

def initGeo():
    # Load the geometry model
    filename = 'wingbox_solid1.stp'
    geo = TMR.LoadModel(filename)

    # Create a model by discarding the volumes
    verts = geo.getVertices()
    edges = geo.getEdges()
    faces = geo.getFaces()
    vols = geo.getVolumes()

    # set top edges to bottom edges
    top_edges = [66, 38, 39, 40, 6, 7, 8, 9, 10, 11, 12,
                 65, 36, 42, 41, 4, 0, 18, 17, 16, 15, 14,
                 67, 37, 5, 13]
    bottom_edges = [73, 61, 60, 59, 33, 32, 31, 30, 29, 28, 27,
                    71, 57, 58, 56, 20, 2, 21, 22, 23, 24, 25,
                    72, 62, 34, 26]
    index = 0
    for i in range(len(top_edges)):
        edges[top_edges[i]].setSource(edges[bottom_edges[i]])
        index += 1
        
    # set left edges to right edges: LE
    vert_edges_LE = [75, 69, 68, 70, 44, 43, 45, 47, 49, 51, 53, 54]
    index = 0
    for i in range(len(vert_edges_LE)-1):
        edges[vert_edges_LE[i]].setSource(edges[vert_edges_LE[i+1]])
        index += 1

    # Set LE inboard edge to TE inboard edge
    vert_edges_TE = [74, 63, 64, 55, 19, 1, 3, 35, 46, 48, 50, 52]
    edges[vert_edges_LE[0]].setSource(edges[vert_edges_TE[0]])

    # Set left edges to right edges: TE
    index = 0
    for i in range(len(vert_edges_TE)-1):
        edges[vert_edges_TE[i]].setSource(edges[vert_edges_TE[i+1]])
        index += 1

    # Create TFI Faces
    edges1 = [edges[62], edges[63], edges[37], edges[69]]
    dirs1 = [1, -1, 1, 1]
    verts1 = [verts[43], verts[39], verts[33], verts[34]]
    faces.append(TMR.TFIFace(edges1, dirs1, verts1))
    edges2 = [edges[34], edges[19], edges[5], edges[44]]
    dirs2 = [1, -1, 1, 1]
    verts2 = [verts[31], verts[18], verts[4], verts[5]]
    faces.append(TMR.TFIFace(edges2, dirs2, verts2))

    # Define new volumes from faces
    faces1 = [faces[22], faces[26], faces[28], faces[30], faces[29], faces[27]]
    dir1 = [1, 1, -1, -1, -1, -1]
    faces2 = [faces[5], faces[19], faces[30], faces[31], faces[23], faces[24],
              faces[25], faces[20], faces[21], faces[18]]
    dir2 = [-1, -1, 1, -1, 1, 1, 1, 1, 1, 1]
    faces3 = [faces[1], faces[3], faces[31], faces[16], faces[6], faces[7], faces[9],
              faces[11], faces[13], faces[15], faces[17], faces[2], faces[0], faces[4],
              faces[8], faces[10], faces[12], faces[14]]
    dir3 = [-1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    
    vols = [TMR.Volume(faces1, dir1),
            TMR.Volume(faces2, dir2),
            TMR.Volume(faces3, dir3)]

    # set top faces to bottom faces
    faces[26].setSource(vols[0], faces[22])
    faces[19].setSource(vols[1], faces[5])
    faces[3].setSource(vols[2], faces[1])

    geo = TMR.Model(verts, edges, faces, vols)

    return geo

# Set the communicator
comm = MPI.COMM_WORLD

# Create the geometry
geo = initGeo()

# Mark the boundary condition faces
faces = geo.getFaces()
volumes = geo.getVolumes()
faces[28].setAttribute('fixed') # fix root
faces[22].setAttribute('surface') # top inboard face
faces[5].setAttribute('surface') # top middle face
faces[1].setAttribute('surface') # top outboard face

# Set the boundary conditions for the problem
bcs = TMR.BoundaryConditions()
bcs.addBoundaryCondition('fixed')

# Create the geometry    
mesh = TMR.Mesh(comm, geo)

# # Mesh the wing box
opts = TMR.MeshOptions()
# opts.num_smoothing_steps = 10

# Taper elements in the spanwise direction
ymax = 914.0
hmin = 1.0
hmax = 10.0
c = hmax
ay = -(hmax - hmin)/ymax
fs = TMR.LinearElementSize(hmin, hmax,
                           c=c, ay=ay)

mesh.mesh(opts=opts, fs=fs)
mesh.writeToVTK('wingbox_vol.vtk', 'hex')
# Create a model from the mesh
model = mesh.createModelFromMesh()

# Create the corresponding mesh topology from the mesh-model 
topo = TMR.Topology(comm, model)

#sys.exit(0)

# Create the quad forest and set the topology of the forest
forest = TMR.OctForest(comm)
forest.setTopology(topo)

# Create the trees, rebalance the elements and repartition
nlevels = 2
order = 2
#forest.createTrees(nlevels-1)

forest.createTrees(1) #nlevels-1)
octs = forest.getOctants()
q, he = mesh.getMeshConnectivity()
Xp = mesh.getMeshPoints()
refine_array = np.zeros(np.shape(octs))

# linear refinement in the y-direction
max_level = 3
min_level = 1
span = 100

index = 0
#for i in range(len(octs)):
    # oc = octs[i]
    # y = ((min_level-max_level)/span)*y + mex_level
    # i += 1
    
sys.exit(0)

assembler, problem = createTopoProblem(forest,
                                       nlevels=nlevels)
aux = addFaceTraction(order, forest, 'surface', assembler,
                      [0.0, 0.0, 1.0])


assembler.zeroVariables()
force = assembler.createVec()
assembler.setAuxElements(aux)
assembler.assembleRes(force)
force.scale(-1.0)
forces = [force]

problem.setLoadCases(forces)
funcs = [functions.StructuralMass(assemblers[0])]
m_fixed = 10.0
problem.addConstraints(0, funcs, [-m_fixed], [1.0])
problem.setObjective([1.0])
problem.initialize()
problem.setPrefix('./')

max_bfgs = 20
opt = ParOpt.pyParOpt(problem, max_bfgs, ParOpt.BFGS)
opt.optimize()

# Output for visualization
flag = (TACS.ToFH5.NODES |
        TACS.ToFH5.DISPLACEMENTS |
        TACS.ToFH5.STRAINS |
        TACS.ToFH5.STRESSES |
        TACS.ToFH5.EXTRAS)
f5 = TACS.ToFH5(assembler, TACS.PY_SOLID, flag)
f5.writeToFile('bracket.f5')
