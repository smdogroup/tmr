from mpi4py import MPI
from tmr import TMR
import sys
import numpy as np

# Set the communicator
comm = MPI.COMM_WORLD

# Load the geometry model
filename = 'wingbox_solid.stp'
geo = TMR.LoadModel(filename)

# Create a model by discarding the volumes
verts = geo.getVertices()
edges = geo.getEdges()
faces = geo.getFaces()
vols = geo.getVolumes()

# for i, f in enumerate(faces):
#     f.writeToVTK('face{0}.vtk'.format(i))

# for i, e in enumerate(edges):
#     e.writeToVTK('edge{0}.vtk'.format(i))

# set top edges to bottom edges
edges[77].setSource(edges[69])
edges[60].setSource(edges[44])
edges[61].setSource(edges[43])
edges[62].setSource(edges[42])
edges[58].setSource(edges[41])
edges[20].setSource(edges[18])
edges[2].setSource(edges[0])
edges[22].setSource(edges[17])
edges[23].setSource(edges[16])
edges[24].setSource(edges[15])
edges[25].setSource(edges[14])
edges[26].setSource(edges[13])
edges[27].setSource(edges[12])
edges[28].setSource(edges[11])
edges[29].setSource(edges[10])
edges[30].setSource(edges[9])
edges[31].setSource(edges[8])
edges[32].setSource(edges[7])
edges[33].setSource(edges[6])
edges[34].setSource(edges[5])
edges[63].setSource(edges[40])
edges[64].setSource(edges[39])
edges[65].setSource(edges[38])
edges[66].setSource(edges[37])
edges[79].setSource(edges[67])
edges[78].setSource(edges[68])
edges[59].setSource(edges[36])
edges[21].setSource(edges[4])
# set left edges to right edges: front faces
edges[80].setSource(edges[70])
edges[70].setSource(edges[71])
edges[71].setSource(edges[74])
edges[74].setSource(edges[57])
edges[57].setSource(edges[19])
edges[19].setSource(edges[1])
edges[1].setSource(edges[3])
edges[3].setSource(edges[35])
edges[35].setSource(edges[49])
edges[49].setSource(edges[51])
edges[51].setSource(edges[53])
edges[53].setSource(edges[55])
# set forward inboard edge to aft inboard edge
edges[80].setSource(edges[81])
# set left edges to right edges: aft faces
edges[81].setSource(edges[73])
edges[73].setSource(edges[72])
edges[72].setSource(edges[75])
edges[75].setSource(edges[76])
edges[76].setSource(edges[46])
edges[46].setSource(edges[45])
edges[45].setSource(edges[47])
edges[47].setSource(edges[48])
edges[48].setSource(edges[50])
edges[50].setSource(edges[52])
edges[52].setSource(edges[54])
edges[54].setSource(edges[56])

new_edges_1 = [edges[36], edges[73], edges[59], edges[70]]
dir1 = [1, 1, 1, -1]
new_verts_1 = [verts[32], verts[33], verts[41], verts[42]]
faces.append(TMR.TFIFace(new_edges_1, dir1, new_verts_1))
new_edges_2 = [edges[4], edges[46], edges[21], edges[19]]
dir2 = [1, 1, 1, -1]
new_verts_2 = [verts[4], verts[5], verts[19], verts[18]]
faces.append(TMR.TFIFace(new_edges_2, dir2, new_verts_2))

# # Write out the edge loops to indicate the direction of the faces
# index = 0
# for f in faces:
#     fp = open('face_orient%d.vtk'%(index), 'w')
#     index += 1

#     fp.write('# vtk DataFile Version 3.0\n')
#     fp.write('vtk output\n')
#     fp.write('ASCII\n')
#     fp.write('DATASET UNSTRUCTURED_GRID\n')

#     npts = 0
#     nedges = 0
#     for k in range(f.getNumEdgeLoops()):
#         loop = f.getEdgeLoop(k)
#         e, dirs = loop.getEdgeLoop()
#         nedges += len(e)
#         npts += len(e)+1
#     fp.write('POINTS %d float\n'%(npts))

#     for k in range(f.getNumEdgeLoops()):
#         loop = f.getEdgeLoop(k)
#         e, dirs = loop.getEdgeLoop()
#         for i in range(len(e)):
#             v1, v2 = e[i].getVertices()
#             if dirs[i] > 0:
#                 pt1 = v1.evalPoint()
#                 pt2 = v2.evalPoint()
#             else:
#                 pt1 = v2.evalPoint()
#                 pt2 = v1.evalPoint()
#             if i == 0:
#                 fp.write('%e %e %e\n'%(pt1[0], pt1[1], pt1[2]))
#             fp.write('%e %e %e\n'%(pt2[0], pt2[1], pt2[2]))

#     fp.write('CELLS %d %d\n'%(nedges, 3*nedges))
#     npts = 0
#     for k in range(f.getNumEdgeLoops()):
#         loop = f.getEdgeLoop(k)
#         e, dirs = loop.getEdgeLoop()
#         for i in range(len(e)):
#             fp.write('2 %d %d\n'%(npts+i, npts+i+1))
#         npts += len(e)+1

#     fp.write('CELL_TYPES %d\n'%(nedges))
#     for i in range(nedges):
#         fp.write('3\n')

#     # write vector data for each edge
#     fp.write('CELL_DATA %d\n'%(nedges))
#     fp.write('VECTORS scalars float\n')
#     for k in range(f.getNumEdgeLoops()):
#         loop = f.getEdgeLoop(k)
#         e, dirs = loop.getEdgeLoop()
#         for i in range(len(e)):
#             v1, v2 = e[i].getVertices()
#             if dirs[i] > 0:
#                 pt1 = v1.evalPoint()
#                 pt2 = v2.evalPoint()
#             else:
#                 pt1 = v2.evalPoint()
#                 pt2 = v1.evalPoint()
#             vec = pt2 - pt1
#             fp.write('%e %e %e\n'%(vec[0], vec[1], vec[2]))
            
#     # Close the file
#     fp.close()

# Write out the edge loops to indicate the direction of the faces
fp = open('face_orient.dat', 'w')
fp.write('Variables = x, y, z, tx, ty, tz\n')
index = 0
for f in faces:
    for k in range(f.getNumEdgeLoops()):
        fp.write('Zone T = \"Face %d Loop %d\"\n'%(index, k))

        loop = f.getEdgeLoop(k)
        e, dirs = loop.getEdgeLoop()
        pts = np.zeros((len(e)+1, 3))
        tx = np.zeros((len(e)+1, 3))
        for i in range(len(e)):
            v1, v2 = e[i].getVertices()
            if dirs[i] > 0:
                pt1 = v1.evalPoint()
                pt2 = v2.evalPoint()
            else:
                pt1 = v2.evalPoint()
                pt2 = v1.evalPoint()
            if i == 0:
                pts[0,:] = pt1[:]
            pts[i+1,:] = pt2[:]

        for i in xrange(len(e)):
            tx[i,:] = pts[i+1,:] - pts[i,:]

        for i in xrange(len(e)+1):
            fp.write('%e %e %e %e %e %e\n'%(
                pts[i,0], pts[i,1], pts[i,2], tx[i,0], tx[i,1], tx[i,2]))

    index += 1
    
# Close the file
fp.close()

vol_faces_1 = [faces[20], faces[28], faces[31], faces[32], faces[29], faces[30]]
dir1 = [1, 1, -1, -1, -1, -1]
vol_faces_2 = [faces[5], faces[19], faces[32], faces[33], faces[21], faces[23],
               faces[25], faces[18], faces[22], faces[24], faces[26], faces[27]]
dir2 = [-1, -1, 1, -1, 1, 1, 1, 1, 1, 1, 1, 1]
vol_faces_3 = [faces[1], faces[3], faces[33], faces[17], faces[2], faces[0], faces[4],
               faces[9], faces[11], faces[13], faces[15], faces[6], faces[7], faces[8],
               faces[10], faces[12], faces[14], faces[16]]
dir3 = [-1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

vols = [TMR.Volume(vol_faces_1, dir1),
        TMR.Volume(vol_faces_2, dir2),
        TMR.Volume(vol_faces_3, dir3)]

# set top faces to bottom faces
faces[3].setSource(vols[2], faces[1])
faces[19].setSource(vols[1], faces[5])
faces[28].setSource(vols[0], faces[20])

geo = TMR.Model(verts, edges, faces, vols)

# Create the geometry    
mesh = TMR.Mesh(comm, geo)
 
# Mesh the wing box
opts = TMR.MeshOptions()
opts.num_smoothing_steps = 10
#opts.write_post_smooth_quad = 1

# Mesh the geometry with the given target size
htarget = 5.0

ymax = 914.0
hmin = 2.0
hmax = 20.0
c = hmax
ay = -(hmax - hmin)/ymax
fs = TMR.LinearElementSize(hmin, hmax,
                           c=c, ay=ay)

mesh.mesh(opts=opts, fs=fs)

# Write the mesh to file
mesh.writeToBDF('wingbox_vol.bdf', 'hex')
mesh.writeToVTK('wingbox_vol.vtk', 'hex')

# Create a model from the mesh
model = mesh.createModelFromMesh()

# Create the corresponding mesh topology from the mesh-model 
topo = TMR.Topology(comm, model)

# Create the quad forest and set the topology of the forest
forest = TMR.OctForest(comm)
forest.setTopology(topo)

# Create random trees and balance the mesh. Print the output file
forest.createRandomTrees(nrand=3, min_lev=0, max_lev=3)
forest.balance(1)
filename = 'forest-mesh%d.vtk'%(comm.rank)
forest.writeForestToVTK(filename)
