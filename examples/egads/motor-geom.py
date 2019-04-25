import os
from egads4py import egads
from mpi4py import MPI
from tmr import TMR
import numpy as np
import argparse

def jacobian3d(sh, x):
    detJ = 0.0
    invsqrt = 1.0/np.sqrt(2)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                eta = [(i-1)*invsqrt,
                       (j-1)*invsqrt,
                       (k-1)*invsqrt]

                detJ += detJacobian3d(h, x, eta)
    return detJ

def detJacobian3d(h, x, eta):
    N1 = 0.125*np.array([-(1.0 - eta[1])*(1.0 - eta[2]),
                         (1.0 - eta[1])*(1.0 - eta[2]),
                         (1.0 + eta[1])*(1.0 - eta[2]),
                         -(1.0 + eta[1])*(1.0 - eta[2]),
                         -(1.0 - eta[1])*(1.0 + eta[2]),
                         (1.0 - eta[1])*(1.0 + eta[2]),
                         (1.0 + eta[1])*(1.0 + eta[2]),
                         -(1.0 + eta[1])*(1.0 + eta[2])])
    N2 = 0.125*np.array([-(1.0 - eta[0])*(1.0 - eta[2]),
                         -(1.0 + eta[0])*(1.0 - eta[2]),
                         (1.0 + eta[0])*(1.0 - eta[2]),
                         (1.0 - eta[0])*(1.0 - eta[2]),
                         -(1.0 - eta[0])*(1.0 + eta[2]),
                         -(1.0 + eta[0])*(1.0 + eta[2]),
                         (1.0 + eta[0])*(1.0 + eta[2]),
                        (1.0 - eta[0])*(1.0 + eta[2])])
    N3 = 0.125*np.array([-(1.0 - eta[0])*(1.0 - eta[1]),
                         -(1.0 + eta[0])*(1.0 - eta[1]),
                         -(1.0 + eta[0])*(1.0 + eta[1]),
                         -(1.0 - eta[0])*(1.0 + eta[1]),
                         (1.0 - eta[0])*(1.0 - eta[1]),
                         (1.0 + eta[0])*(1.0 - eta[1]),
                         (1.0 + eta[0])*(1.0 + eta[1]),
                         (1.0 - eta[0])*(1.0 + eta[1])])

    Xd = np.zeros((3, 3))
    for i, N in enumerate([N1, N2, N3]):
        Xd[:,i] = np.dot(N, x[h, :])

    return np.linalg.det(Xd)

# Load in the model to TMR
comm = MPI.COMM_WORLD

# Create an argument parser to read in arguments from the commnad line
p = argparse.ArgumentParser()
p.add_argument('--htarget', type=float, default=2.0)
p.add_argument('--extension', type=str, default='egads', help='egads or step')
p.add_argument('--model_type', type=str, default='full',
               help='full or anything else')
args = p.parse_args()

htarget = args.htarget
extension = args.extension
model_type = args.model_type

# Create the egads context
ctx = egads.context()

def getBodyFacesAndDirs(body):
    # takes in an egads body object, and returns it's faces

    children = body.getChildren()

    # Check that there is only one shell
    if len(children) > 1:
        print("More than one shell associated with this body")
        return None
    # Get the shell associated with the body
    else:
        shell = children[0]

    # geo, oclass, mtype, lim, children, sens = shell.getTopology()
    children = shell.getChildren()

    return children

# Create the top shell
x1 = [0, 0, 6.35]
x2 = [0, 0, 67.06]
radius = 80.78
c1 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x1, x2, radius])

x1 = [0, 0, 6.35]
x2 = [0, 0, 67.06]
radius = 70.78
c2 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x1, x2, radius])

shell_model = c1.solidBoolean(c2, egads.SUBTRACTION)
top_shell = shell_model.getChildren()[0]
top_shell.attributeAdd('name', egads.ATTRSTRING, 'top_shell')

# Get the faces associated with the top shell, then compute the range
# over each face to determine which face it is
faces = getBodyFacesAndDirs(top_shell)

# Set the face attributes for the top shell volume
faces[3].attributeAdd('name', egads.ATTRSTRING, 'shell-top-face')
faces[1].attributeAdd('name', egads.ATTRSTRING, 'shell-bottom-face')

# Save the shell model
if comm.rank == 0:
    shell_model.saveModel('shell.%s'%(extension), overwrite=True)

# Create the bottom outer ring
x1 = [0, 0, 6.35]
x2 = [0, 0, 0]
radius = 80.78
c3 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x1, x2, radius])

x1 = [0, 0, 6.35]
x2 = [0, 0, 0]
radius = 70.78
c4 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x1, x2, radius])

ring_model = c3.solidBoolean(c4, egads.SUBTRACTION)
ring = ring_model.getChildren()[0]
ring.attributeAdd('name', egads.ATTRSTRING, 'ring')

faces = getBodyFacesAndDirs(ring)

# Set the face attributes for the ring volume
faces[3].attributeAdd('name', egads.ATTRSTRING, 'ring-top-face')
faces[1].attributeAdd('name', egads.ATTRSTRING, 'ring-bottom-face')
faces[4].attributeAdd('name', egads.ATTRSTRING, 'ring-inner-face1')
faces[5].attributeAdd('name', egads.ATTRSTRING, 'ring-inner-face2')

# Save the ring model
if comm.rank == 0:
    ring_model.saveModel('ring.%s'%(extension), overwrite=True)

# Create the bottom plate
x1 = [0, 0, 0]
x2 = [0, 0, 6.35]
radius = 70.78
plate = ctx.makeSolidBody(egads.CYLINDER, rdata=[x1, x2, radius])
plate.attributeAdd('name', egads.ATTRSTRING, 'plate')

# Subtract the holes from the plate
x1 = [0, 0, 0]
x2 = [0, 0, 6.35]
radius = 9.91
mid_hole = ctx.makeSolidBody(egads.CYLINDER, rdata=[x1, x2, radius])

x1 = [0, 45.98, 0]
x2 = [0, 45.98, 6.35]
radius = 4.08
hole1 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x1, x2, radius])

x1 = [39.82, -22.99, 0]
x2 = [39.82, -22.99, 6.35]
radius = 4.08
hole2 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x1, x2, radius])

x1 = [-39.82, -22.99, 0]
x2 = [-39.82, -22.99, 6.35]
radius = 4.08
hole3 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x1, x2, radius])

for hole in [mid_hole, hole1, hole2, hole3]:
    plate_model = plate.solidBoolean(hole, egads.SUBTRACTION)
    plate = plate_model.getChildren()[0]

faces = getBodyFacesAndDirs(plate)

# Set the face attributes for the plate volume
faces[2].attributeAdd('name', egads.ATTRSTRING, 'plate-top-face')
faces[3].attributeAdd('name', egads.ATTRSTRING, 'plate-bottom-face')
faces[1].attributeAdd('name', egads.ATTRSTRING, 'plate-outer-face1')
faces[0].attributeAdd('name', egads.ATTRSTRING, 'plate-outer-face2')

# Save the plate model
if comm.rank == 0:
    plate_model.saveModel('plate.%s'%(extension), overwrite=True)

# Set the meshing options
opts = TMR.MeshOptions()
opts.write_mesh_quality_histogram = 1
opts.triangularize_print_iter = 50000

# Load the separate geometries and mesh each
if extension == 'egads':
    shell_geo = TMR.ConvertEGADSModel(shell_model)
    ring_geo = TMR.ConvertEGADSModel(ring_model)
    plate_geo = TMR.ConvertEGADSModel(plate_model)
else:
    # Load in the files
    comm.Barrier()
    for rank in range(comm.size):
        if comm.rank == rank:
            shell_geo = TMR.LoadModel('shell.%s'%(extension))
            ring_geo = TMR.LoadModel('ring.%s'%(extension))
            plate_geo = TMR.LoadModel('plate.%s'%(extension))
        comm.Barrier()

# All the model objects
if model_type == 'full':
    all_geos = [ring_geo, plate_geo, shell_geo]
else:
    all_geos = [plate_geo, ring_geo]

# Create the full list of vertices, edges, faces and volumes
verts = []
edges = []
faces = []
vols = []
for geo in all_geos:
    verts.extend(geo.getVertices())
    edges.extend(geo.getEdges())
    faces.extend(geo.getFaces())
    vols.extend(geo.getVolumes())

for vol in vols:
    fail = vol.setExtrudeFaces(reverse_extrude=True)
    if fail:
        print('Setting the swept directions failed')

# Combine the geometries and mesh the assembly
if model_type == 'full':
    TMR.setMatchingFaces([plate_geo, ring_geo])
    TMR.setMatchingFaces([shell_geo, ring_geo])
else:
    TMR.setMatchingFaces([plate_geo, ring_geo])

# Create the geometry
geo = TMR.Model(verts, edges, faces, vols)

# Create the new mesh
mesh = TMR.Mesh(comm, geo)

# Set the meshing options
opts = TMR.MeshOptions()
opts.write_mesh_quality_histogram = 1
opts.triangularize_print_iter = 50000

# Create the surface mesh
mesh.mesh(htarget, opts)

# Write the surface mesh to a file
mesh.writeToVTK('motor.vtk', 'hex')

for index, face in enumerate(faces):
    orient, src = face.getCopySource()
    if src is not None:
        src_mesh = src.getMesh()
        src_mesh.writeToVTK('source_face_mesh%d.vtk'%(index))
        face_mesh = face.getMesh()
        face_mesh.writeToVTK('copied_face_mesh%d.vtk'%(index))

# Get the mesh
x = mesh.getMeshPoints()
x = x.reshape((-1, 3))
quads = mesh.getQuadConnectivity()
hexa = mesh.getHexConnectivity()

if comm.rank == 0:
    count = 0
    vol = 0.0
    for h in hexa:
        j = jacobian3d(h, x)
        if j < 0.0:
            count += 1
        vol += j

    if count > 0:
        print('Warning: %d negative volume Jacobians'%(count))
    print('Volume = %e'%(vol))

# Create the model from the unstructured volume mesh
model = mesh.createModelFromMesh()

# Create the corresponding mesh topology from the mesh-model
topo = TMR.Topology(comm, model)

# Create the quad forest and set the topology of the forest
# forest = TMR.QuadForest(comm)
forest = TMR.OctForest(comm)
forest.setTopology(topo)

# Create random trees and balance the mesh. Print the output file
forest.createRandomTrees(nrand=1, max_lev=6)
forest.balance(1)
filename = 'motor_forest%d.vtk'%(comm.rank)
forest.writeForestToVTK(filename)
