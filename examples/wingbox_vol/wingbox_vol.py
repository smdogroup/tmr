from mpi4py import MPI
from tmr import TMR
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

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

# Create new surfaces with TFI
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

# Define new volumes from faces
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

# Taper elements in the spanwise direction
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

#topo_mesh = topo.getMesh()

# Create the quad forest and set the topology of the forest
forest = TMR.OctForest(comm)
forest.setTopology(topo)

# Create random trees and balance the mesh. Print the output file
forest.createRandomTrees(nrand=3, min_lev=0, max_lev=3)
forest.createTrees(1)
forest.balance(1)
filename = 'forest-mesh%d.vtk'%(comm.rank)
forest.writeForestToVTK(filename)
'''
Create histogram of hex element mesh quality 
'''
X = mesh.getMeshPoints()
quads, hexes = mesh.getMeshConnectivity()

forest.createNodes(2)
conn = forest.createMeshConn()

def dist(x1, x2):
    d = ((x2[0]-x1[0])**2 + (x2[1]-x1[1])**2 + (x2[2]-x1[2])**2)**0.5
    return d

# compute aspect ratio of each hex
AR = np.zeros(len(hexes))
index = 0
for h in hexes:
    # get side lengths
    sides = np.zeros(12)
    for i in range(3):
        sides[i] = dist(X[h[i]], X[h[i+1]])
    sides[3] = dist(X[h[3]], X[h[0]])
    for i in range(4, 7):
        sides[i] = dist(X[h[i]], X[h[i+1]])
    sides[7] = dist(X[h[7]], X[h[4]])
    for i in range(4):
        sides[i+8] = dist(X[h[i]], X[h[i+4]])
    AR[index] = np.amax(sides)/np.amin(sides)
    index += 1


# Compute min angle of each element
def computeAngle(v1, v2, d1, d2):
    angle = np.arccos(np.dot(v1, v2)/(d1*d2))
    angle *= 180.0/np.pi

    return angle

edge_angles = [[0, 3], [0, 8 ], [3, 8 ],
               [0, 1], [0, 9 ], [1, 9 ],
               [1, 2], [1, 10], [2, 10],
               [2, 3], [2, 11], [3, 11],
               [4, 7], [4, 8 ], [7, 8 ],
               [4, 5], [4, 9 ], [5, 9 ],
               [5, 6], [5, 10], [6, 10],
               [6, 7], [6, 11], [7, 11]]
min_angles = np.zeros(len(hexes))
index = 0
for h in hexes:
    # get side lengths
    sides = np.zeros(12)
    vectors = np.zeros((12, 3))
    angles = np.zeros(24)
    for i in range(3):
        sides[i] = dist(X[h[i]], X[h[i+1]])
        vectors[i, :] = X[h[i+1]] - X[h[i]]
    sides[3] = dist(X[h[3]], X[h[0]])
    vectors[3, :] = X[h[0]] - X[h[3]]
    for i in range(4, 7):
        sides[i] = dist(X[h[i]], X[h[i+1]])
        vectors[i, :] = X[h[i+1]] - X[h[i]]
    sides[7] = dist(X[h[7]], X[h[4]])
    vectors[7, :] = X[h[4]] - X[h[7]]
    for i in range(4):
        sides[i+8] = dist(X[h[i]], X[h[i+4]])
        vectors[i+8, :] = X[h[i+4]] - X[h[i]]
    for i, ea in enumerate(edge_angles):
        v1 = vectors[ea[0], :]
        v2 = vectors[ea[1], :]
        d1 = sides[ea[0]]
        d2 = sides[ea[1]]
        angles[i] = computeAngle(v1, v2, d1, d2)

    min_angles[index] = np.amin(angles)
    index += 1

''' Plot the histogram of aspect ratio '''    
    
hist_min = 0.0
hist_max = np.ceil(max(AR))
nbins = 2*hist_max
cutoffs = np.linspace(hist_min, hist_max, nbins+1)
AR_hist, bin_edges = np.histogram(AR, bins=cutoffs)  

AR_hist_percent = 100.0*AR_hist/(len(AR)+1)
# plot the histogram

fig, ax = plt.subplots()
n, bins, pathches = plt.hist(AR, bins=nbins, range=(hist_min, hist_max), align='mid', color=(0.1176, 0.5647, 1.0, 1.0), normed=False)
# width = 0.2 
# xlocs = np.arrange(range(len(AR_hist)))#+width#/2.0
# rects = ax.bar(xlocs, AR_hist_percent)#, width=width)
#plt.xticks(width*cutoffs+width, cutoffs)
plt.title('Aspect Ratio of Hexahedral Elements')
plt.xlabel('Aspect Ratio')
plt.ylabel('Number of Elements')

majorLocator = MultipleLocator(2)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(1)

ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_major_formatter(majorFormatter)
ax.xaxis.set_minor_locator(minorLocator)
#ax.yaxis.set_major_formatter(plt.FuncFormatter('{:.0f}%'.format))
ax.tick_params(which='both', direction='out', top='off', right='off')
plt.tight_layout()
plt.show()

''' Plot the histogram of minimum angle '''    
hist_min = 0.0
hist_max = 90.0
nbins = 1*hist_max

fig, ax = plt.subplots()
n, bins, pathches = plt.hist(min_angles, bins=nbins, range=(hist_min, hist_max), align='mid', color=(0.1176, 0.5647, 1.0, 1.0), normed=False)
#plt.xticks(width*cutoffs+width, cutoffs)
plt.title('Minimum Angle of Hexahedral Elements')
plt.xlabel('Minimum Angle [degrees]')
plt.ylabel('Number of Elements')

majorLocator = MultipleLocator(5)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(1)

ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_major_formatter(majorFormatter)
ax.xaxis.set_minor_locator(minorLocator)
#ax.yaxis.set_major_formatter(plt.FuncFormatter('{:.0f}%'.format))
ax.tick_params(which='both', direction='out', top='off', right='off')
plt.tight_layout()
plt.show()
