"""
Create the TMR logo
"""

from PIL import Image, ImageDraw, ImageFont
import numpy as np
from mpi4py import MPI
from tmr import TMR

width = 72
height = 32

img = Image.new("RGB", (width, height), color=(255, 255, 255))
fnt = ImageFont.truetype("/Library/Fonts/Arial Bold.ttf", 32)
d = ImageDraw.Draw(img)
d.text((2, -2), "TMR", font=fnt, fill=(0, 0, 0))
img.save("tmr.png")
pixels = list(img.getdata())
pixels = [pixels[i * width : (i + 1) * width] for i in range(height)]

# Set the bounding box for now
xlow = [0.0, 0.0, -5.0]
xhigh = [width, height, 5.0]

# Set hmin/hmax
hmin = 0.15
hmax = 0.5
h = 0.25

# Create a box
box = TMR.BoxFeatureSize(xlow, xhigh, hmin, hmax)

cutoff = 128
for j in range(height):
    for i in range(width):
        (R, G, B) = pixels[j][i]
        gscale = (0.3 * R) + (0.59 * G) + (0.11 * B)
        if gscale < cutoff:
            xlow = [i, height - j, -1]
            xhigh = [i + 1, height - (j + 1), 1]
            box.addBox(xlow, xhigh, h)

# Set the number of load cases
nu = 2
nv = 2
x = [0.0, width]
y = [0.0, height]
pts = np.zeros((nu, nv, 3))
for j in range(nv):
    for i in range(nu):
        pts[i, j, 0] = x[i]
        pts[i, j, 1] = y[j]

# Create the b-spline surface
surf = TMR.BsplineSurface(pts)
face = TMR.FaceFromSurface(surf)

v1 = TMR.VertexFromFace(face, 0.0, 0.0)
v2 = TMR.VertexFromFace(face, 1.0, 0.0)
v3 = TMR.VertexFromFace(face, 1.0, 1.0)
v4 = TMR.VertexFromFace(face, 0.0, 1.0)
verts = [v1, v2, v3, v4]

# Set up the first edge
pcurve1 = TMR.BsplinePcurve(np.array([[0.0, 0.0], [1.0, 0.0]]))
edge1 = TMR.EdgeFromFace(face, pcurve1)
edge1.setVertices(v1, v2)

# Set up the second edge
pcurve2 = TMR.BsplinePcurve(np.array([[1.0, 0.0], [1.0, 1.0]]))
edge2 = TMR.EdgeFromFace(face, pcurve2)
edge2.setVertices(v2, v3)

# Set up the third edge
pcurve3 = TMR.BsplinePcurve(np.array([[1.0, 1.0], [0.0, 1.0]]))
edge3 = TMR.EdgeFromFace(face, pcurve3)
edge3.setVertices(v3, v4)

# Set up the fourth edge
pcurve4 = TMR.BsplinePcurve(np.array([[0.0, 1.0], [0.0, 0.0]]))
edge4 = TMR.EdgeFromFace(face, pcurve4)
edge4.setVertices(v4, v1)
edges = [edge1, edge2, edge3, edge4]

# Create the loop
dirs = [1, 1, 1, 1]
loop = TMR.EdgeLoop([edge1, edge2, edge3, edge4], dirs)
face.addEdgeLoop(1, loop)

# Create the TMRModel
faces = [face]
geo = TMR.Model(verts, edges, faces)

# Create the mesh
comm = MPI.COMM_WORLD
mesh = TMR.Mesh(comm, geo)

# Mesh the part
opts = TMR.MeshOptions()
opts.num_smoothing_steps = 20
opts.write_mesh_quality_histogram = 1
opts.mesh_type_default = TMR.UNSTRUCTURED

# Mesh the geometry with the given target size
htarget = 10.0
mesh.mesh(fs=box, opts=opts)
mesh.writeToVTK("surface-mesh.vtk")
