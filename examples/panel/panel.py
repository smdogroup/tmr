from mpi4py import MPI
from tmr import TMR
from tacs import TACS, elements, constitutive
import numpy as np

class CreateMe(TMR.QuadCreator):
    def __init__(self, bcs):
        TMR.QuadCreator.__init__(bcs)

    def createElement(self, order, quad):
        # Set constitutive properties
        rho = 2500.0 # density, kg/m^3
        E = 70e3 # elastic modulus, Pa
        nu = 0.3 # poisson's ratio
        kcorr = 5.0 / 6.0 # shear correction factor
        ys = 350e3 # yield stress, Pa
        min_thickness = 0.2
        max_thickness = 1.0
        thickness = 2.5
        
        stiff = constitutive.isoFSDT(rho, E, nu, kcorr, ys, 
                                     thickness, quad.face,
                                     min_thickness, max_thickness)
        stiff.setRefAxis(np.array([1.0, 0.0, 0.0]))

        element = None
        if order == 2:
            element = elements.MITCShell(2, stiff)
        elif order == 3:
            element = elements.MITCShell(3, stiff)
        elif order == 4:
            element = elements.MITCShell(4, stiff)
        return element

    def createAuxElements(self, assembler, order=3):
        aux = TACS.AuxElements()
        tx = np.zeros(order*order)
        ty = np.zeros(order*order)
        tz = -np.ones(order*order)
        trac = elements.ShellTraction(order, tx, ty, tz)
        for i in range(assembler.getNumElements()):
            aux.addElement(i, trac)
        return aux

    def createMg(self, forest, nlevels=2, order=3):
        # Create the forest
        forest.balance(1)
        forest.repartition()
        forest.setMeshOrder(order)

        # Create the forests
        forests = [ forest ]
        assemblers = [ self.createTACS(forest) ]

        if order == 3:
            forests.append(forests[-1].duplicate())
            forests[-1].balance(1)
            forests[-1].repartition()
            forests[-1].setMeshOrder(2)
            assemblers.append(self.createTACS(forests[-1]))

        for i in range(nlevels-1):
            forests.append(forests[-1].coarsen())
            forests[-1].balance(1)
            forests[-1].repartition()
            forests[-1].setMeshOrder(2)
            assemblers.append(self.createTACS(forests[-1]))

        # Create the multigrid object
        mg = TMR.createMg(assemblers, forests)

        return assemblers[0], mg

def create_panel(Lx, Ly, use_hole=True):
    '''
    Create a panel with a whole in it
    '''
    # Set the number of load cases
    nu = 2
    nv = 2
    x = np.linspace(-0.5*Lx, 0.5*Lx, nu)
    y = np.linspace(-0.5*Ly, 0.5*Ly, nv)
    pts = np.zeros((nu, nv, 3))
    for j in range(nv):
        for i in range(nu):
            pts[i,j,0] = x[i]
            pts[i,j,1] = y[j]

    # Create the b-spline surface
    surf = TMR.BsplineSurface(pts)
    face = TMR.FaceFromSurface(surf)

    r = 0.2
    c = 0.5
    v1 = TMR.VertexFromFace(face, 0.0, 0.0)
    v2 = TMR.VertexFromFace(face, 1.0, 0.0)
    v3 = TMR.VertexFromFace(face, 1.0, 1.0)
    v4 = TMR.VertexFromFace(face, 0.0, 1.0)
    v5 = TMR.VertexFromFace(face, c-r, c)

    # Set up the first edge
    pcurve1 = TMR.BsplinePcurve(np.array([[0.0, 0.0], [1.0, 0.0]]))
    edge1 = TMR.EdgeFromFace(face, pcurve1)
    edge1.setVertices(v1, v2)
    edge1.setName('y-')

    # Set up the first edge
    pcurve2 = TMR.BsplinePcurve(np.array([[1.0, 0.0], [1.0, 1.0]]))
    edge2 = TMR.EdgeFromFace(face, pcurve2)
    edge2.setVertices(v2, v3)
    edge2.setName('x+')

    # Set up the first edge
    pcurve3 = TMR.BsplinePcurve(np.array([[1.0, 1.0], [0.0, 1.0]]))
    edge3 = TMR.EdgeFromFace(face, pcurve3)
    edge3.setVertices(v3, v4)
    edge3.setName('y+')

    # Set up the first edge
    pcurve4 = TMR.BsplinePcurve(np.array([[0.0, 1.0], [0.0, 0.0]]))
    edge4 = TMR.EdgeFromFace(face, pcurve4)
    edge4.setVertices(v4, v1)
    edge4.setName('x-')

    # Create the inner edge loop
    # (c-r, c+r) -- (c, c+r) -- (c+r, c+r)
    #    |                          |
    #    |                          |
    # (c-r, c)                  (c+r, c)
    #    |                          |
    #    |                          |
    # (c-r, c-r) -- (c, c-r) -- (c+r, c-r)

    pts = [[c-r, c], [c-r, c+r], [c, c+r], [c+r, c+r],
        [c+r, c], [c+r, c-r], [c, c-r], [c-r, c-r], [c-r, c]]
    wts = [1.0, 1.0/np.sqrt(2), 1.0, 1.0/np.sqrt(2), 
           1.0, 1.0/np.sqrt(2), 1.0, 1.0/np.sqrt(2), 1.0]
    Tu = [0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0]
    pcurve5 = TMR.BsplinePcurve(np.array(pts), 
                                tu=np.array(Tu), wts=np.array(wts), k=3)
    edge5 = TMR.EdgeFromFace(face, pcurve5)
    edge5.setVertices(v5, v5)

    # Create the loop
    dirs = [1, 1, 1, 1]
    loop = TMR.EdgeLoop([edge1, edge2, edge3, edge4], dirs)
    face.addEdgeLoop(1, loop)

    if use_hole:
        # Create the second edge loop
        loop = TMR.EdgeLoop([edge5], [1])
        face.addEdgeLoop(-1, loop)

        verts = [v1, v2, v3, v4, v5]
        edges = [edge1, edge2, edge3, edge4, edge5]
    else:
        verts = [v1, v2, v3, v4]
        edges = [edge1, edge2, edge3, edge4]

    faces = [face]

    # Create the TMRModel
    geo = TMR.Model(verts, edges, faces)
    return geo

geo = create_panel(100.0, 100.0, use_hole=True)

# Create the mesh
comm = MPI.COMM_WORLD
mesh = TMR.Mesh(comm, geo)

# Mesh the part
opts = TMR.MeshOptions()
opts.num_smoothing_steps = 20
opts.write_mesh_quality_histogram = 1

# Mesh the geometry with the given target size
htarget = 10.0
mesh.mesh(htarget, opts=opts)
mesh.writeToVTK('surface-mesh.vtk')

# Create a model from the mesh
model = mesh.createModelFromMesh()

# Create the corresponding mesh topology from the mesh-model 
topo = TMR.Topology(comm, model)

# Create the quad forest and set the topology of the forest
forest = TMR.QuadForest(comm)
forest.setTopology(topo)

# Make the creator class
bcs = TMR.BoundaryConditions()
bcs.addBoundaryCondition('x+')
bcs.addBoundaryCondition('x-')
bcs.addBoundaryCondition('y+')
bcs.addBoundaryCondition('y-')

# Allocate the creator class
creator = CreateMe(bcs)

# Create the initial forest
nlevels = 2
forest.createTrees(nlevels-1)

# Target relative error
target_rel_err = 1e-5

order = 3

# Create the problem
assembler, mg = creator.createMg(forest, nlevels=nlevels, order=order)
aux = creator.createAuxElements(assembler, order=order)
assembler.setAuxElements(aux)

# Create a solution vector
ans = assembler.createVec()
res = assembler.createVec()
alpha = 1.0
beta = 0.0
gamma = 0.0
mg.assembleJacobian(alpha, beta, gamma, res)

# Factor the preconditioner
mg.factor()

subspace = 100
gmres = TACS.KSM(mg.getMat(), mg, subspace, isFlexible=1)
gmres.setMonitor(comm, 'GMRES', 1)
gmres.solve(res, ans)
ans.scale(-1.0)
assembler.setVariables(ans)

# Output for visualization 
flag = (TACS.ToFH5.NODES |
        TACS.ToFH5.DISPLACEMENTS |
        TACS.ToFH5.STRAINS)
f5 = TACS.ToFH5(assembler, TACS.PY_SHELL, flag)
f5.writeToFile('visualization.f5')

# Duplicate the forest and refine it uniformly
forest_refined = forest.duplicate()
forest_refined.setMeshOrder(4)
forest_refined.balance(1)
assembler_refined = creator.createTACS(forest_refined)
ans_refined = assembler_refined.createVec()
TMR.computeReconSolution(forest, assembler,
                         forest_refined, assembler_refined, 
                        ans, ans_refined)
assembler_refined.setVariables(ans_refined)

# Output for visualization 
flag = (TACS.ToFH5.NODES |
        TACS.ToFH5.DISPLACEMENTS |
        TACS.ToFH5.STRAINS)
f5 = TACS.ToFH5(assembler_refined, TACS.PY_SHELL, flag)
f5.writeToFile('visualization_refined.f5')

