import os
from egads4py import egads
from mpi4py import MPI
from tmr import TMR

def setMatchingFaces(geo_list, atol=1e-6):
    """
    Take in a list of geometries, find the matching faces,
    and set them as copies
    """

    swept_pairs = []
    for geo in geo_list:
        vols = geo.getVolumes()
        for v in vols:
            pair = v.getSweptFacePairs()
            swept_pairs.append((pair[0], pair[1], v))

    copy_pairs = []
    for i in range(len(geo_list)):
        ifaces = geo_list[i].getFaces()
        for j in range(i+1, len(geo_list)):
            jfaces = geo_list[j].getFaces()
            for fi in ifaces:
                for fj in jfaces:
                    if fi != fj and fi.checkMatching(fj, atol=atol):
                        copy_pairs.append((fi, fj))

    while len(swept_pairs) > 0 or len(copy_pairs) > 0:
        # Find a pair of faces and determine
        edges = []
        if len(swept_pairs) > 0:
            pair = swept_pairs.pop(0)
            edges = [(pair[0], pair[1], pair[2])]
        else:
            pair = copy_pairs.pop(0)
            edges = [(pair[0], pair[1])]

        while True:
            len_edges = len(edges)
            for index, pair in enumerate(swept_pairs):
                vol = pair[2]
                length = len(edges)
                if edges[-1][1].isSameObject(pair[0]):
                    edges.append((pair[0], pair[1], vol))
                elif edges[-1][1].isSameObject(pair[1]):
                    edges.append((pair[1], pair[0], vol))
                elif edges[0][0].isSameObject(pair[0]):
                    edges.insert(0, (pair[1], pair[0], vol))
                elif edges[0][0].isSameObject(pair[1]):
                    edges.insert(0, (pair[0], pair[1], vol))
                if len(edges) > length:
                    swept_pairs.pop(index)
                    break

            for index, pair in enumerate(copy_pairs):
                length = len(edges)
                if edges[-1][1].isSameObject(pair[0]):
                    edges.append((pair[0], pair[1]))
                elif edges[-1][1].isSameObject(pair[1]):
                    edges.append((pair[1], pair[0]))
                elif edges[0][0].isSameObject(pair[0]):
                    edges.insert(0, (pair[1], pair[0]))
                elif edges[0][0].isSameObject(pair[1]):
                    edges.insert(0, (pair[0], pair[1]))
                if len(edges) > length:
                    copy_pairs.pop(index)
                    break

            if len_edges == len(edges):
                break

        # Now set the copy/target edges
        for e in edges:
            if len(e) == 3:
                vol = e[2]
                e[0].setSource(vol, e[1])
            else:
                e[0].setCopyFaces(e[1])

    return


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
C12 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x0, x1, r2])

x0 = [cx1, cy1, 0]
x1 = [cx1, cy1, h1]
C11 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x0, x1, r1])
parts.append(C12.solidBoolean(C11, egads.SUBTRACTION))
parts.append(B1.solidBoolean(C12, egads.SUBTRACTION))

# Create the upper box
x0 = [0, 0, h1]
x1 = [Lx, Ly, h2]
B2 = ctx.makeSolidBody(egads.BOX, rdata=[x0, x1])

# Create the cylinder cutout for the upper box
x0 = [cx1, cy1, h1]
x1 = [cx1, cy1, h1+h2]
C21 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x0, x1, r2])

parts.append(B2.solidBoolean(C21, egads.SUBTRACTION))

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
setMatchingFaces(geos)

# Create the geometry
geo = TMR.Model(verts, edges, faces, vols)

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
