import os
from egads4py import egads
from mpi4py import MPI
from tmr import TMR

comm = MPI.COMM_WORLD

# Create the egads context
ctx = egads.context()

# Set the dimensions
h1 = 10.0
h2 = 15.0
Lx = 100.0
Ly = 100.0
r1 = 15.0
r2 = 25.0

cx1 = 35.0
cy1 = 35.0

cx2 = 45.0
cy2 = 45.0

parts = []

# Create the lower box
x0 = [0, 0, 0]
x1 = [Lx, Ly, h1]
B1 = ctx.makeSolidBody(egads.BOX, rdata=[x0, x1])

# Create the cylinder cutout for the bottom box
x0 = [cx1, cy1, 0]
x1 = [cx1, cy1, h1]
C0 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x0, x1, r2])

x0 = [cx2, cy2, 0]
x1 = [cx2, cy2, h1]
C1 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x0, x1, r2])
C12 = C0.solidBoolean(C1, egads.FUSION)
C12b = C12.getChildren()[0]

# x0 = [cx1, cy1, 0]
# x1 = [cx1, cy1, h1]
# C11 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x0, x1, r1])
# parts.append(C12b.solidBoolean(C11, egads.SUBTRACTION))
parts.append(B1.solidBoolean(C12b, egads.SUBTRACTION))

# Create the upper box
x0 = [0, 0, h1]
x1 = [Lx, Ly, h2]
B2 = ctx.makeSolidBody(egads.BOX, rdata=[x0, x1])

# Create the cylinder cutout for the upper box
x0 = [cx1, cy1, h1]
x1 = [cx1, cy1, h1+h2]
C0 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x0, x1, r2])

x0 = [cx2, cy2, h1]
x1 = [cx2, cy2, h1+h2]
C1 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x0, x1, r2])
C21 = C1.solidBoolean(C0, egads.FUSION)
C21b = C21.getChildren()[0]

parts.append(B2.solidBoolean(C21b, egads.SUBTRACTION))

# Create all of the models
geos = []
for p in parts:
    geos.append(TMR.ConvertEGADSModel(p))

# Create the full list of vertices, edges, faces and volumes
verts = []
edges = []
faces = []
vols = []
for geo in geos:
    verts.extend(geo.getVertices())
    edges.extend(geo.getEdges())
    faces.extend(geo.getFaces())
    vols.extend(geo.getVolumes())

# Set all of the matching faces
TMR.setMatchingFaces(geos)

# Set the swept directions
for vol in vols:
    fail = vol.setExtrudeFaces()

# Create the geometry
geo = TMR.Model(verts, edges, faces) #, vols)

# Create the new mesh
mesh = TMR.Mesh(comm, geo)

# Set the meshing options
opts = TMR.MeshOptions()
opts.write_mesh_quality_histogram = 1
opts.triangularize_print_iter = 50000

# Create the surface mesh
htarget = 2.0
mesh.mesh(htarget, opts)

# Write the surface mesh to a file
mesh.writeToVTK('block.vtk', 'quads')
