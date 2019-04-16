from mpi4py import MPI
from tmr import TMR
from egads4py import egads
import numpy as np

ctx = egads.context()

# Create the bodies that are 
x = [0, 0, 0]
d = [0.5, 0.5, 0.5]
b1 = ctx.makeSolidBody(egads.BOX, rdata=[x, d])

x = [0, 0, 0.5]
d = [0.5, 0.5, 0.5]
b2 = ctx.makeSolidBody(egads.BOX, rdata=[x, d])

# Write out to a STEP file
m1 = ctx.makeTopology(egads.MODEL, children=[b1])
m2 = ctx.makeTopology(egads.MODEL, children=[b2])
m1.saveModel('box1.egads', overwrite=True)
m2.saveModel('box2.egads', overwrite=True)

comm = MPI.COMM_WORLD
htarget = 0.05

geo1 = TMR.LoadModel('box1.egads', print_lev=1)
geo2 = TMR.LoadModel('box2.egads', print_lev=1)

TMR.setMatchingFaces([geo1, geo2])

verts = []
edges = []
faces = []
vols = []
sfi = [0, 4]
for i, geo in enumerate([geo1, geo2]):
    vols = geo.getVolumes()
    fail = vols[0].setExtrudeFaces(source_face_index = sfi[i])
    print(fail)

    verts.extend(geo.getVertices())
    edges.extend(geo.getEdges())
    faces.extend(geo.getFaces())
    # for i, f in enumerate(geo.getFaces()):
    #     umin, vmin, umax, vmax = f.getRange()
    #     pt1 = f.evalPoint(umin, vmin)
    #     pt2 = f.evalPoint(umin, vmax)
    #     pt3 = f.evalPoint(umax, vmin)
    #     pt4 = f.evalPoint(umax, vmax)
    #     xmin = np.amin(np.array([pt1[0], pt2[0], pt3[0], pt4[0]]))
    #     ymin = np.amin(np.array([pt1[1], pt2[1], pt3[1], pt4[1]]))
    #     zmin = np.amin(np.array([pt1[2], pt2[2], pt3[2], pt4[2]]))
    #     xmax = np.amax(np.array([pt1[0], pt2[0], pt3[0], pt4[0]]))
    #     ymax = np.amax(np.array([pt1[1], pt2[1], pt3[1], pt4[1]]))
    #     zmax = np.amax(np.array([pt1[2], pt2[2], pt3[2], pt4[2]]))
    #     print("Face {}: x in [{}, {}], y in [{}, {}], z in [{}, {}]".format(i,
    #                                                                         xmin,
    #                                                                         xmax,
    #                                                                         ymin,
    #                                                                         ymax,
    #                                                                         zmin,
    #                                                                         zmax))
    #     print("")
    vols.extend(geo.getVolumes())

geo = TMR.Model(verts, edges, faces)#, vols)

# Create the new mesh
mesh = TMR.Mesh(comm, geo)

# Set the meshing options
opts = TMR.MeshOptions()
opts.write_mesh_quality_histogram = 1
opts.triangularize_print_iter = 50000

# Create the surface mesh
mesh.mesh(htarget, opts)

# Write the surface mesh to a file
mesh.writeToVTK('output.vtk', 'hex')

