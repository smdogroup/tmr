from egads4py import egads
from mpi4py import MPI
from tmr import TMR

ctx = egads.context()

def getBodyFacesAndDirs(body):
    # takes in an egads body object, and returns it's faces

    children = body.getChildren()
    print(children)

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

model = c1.solidBoolean(c2, egads.SUBTRACTION)
model.saveModel('shell.step', overwrite=True)
top_shell = model.getChildren()[0]
top_shell.attributeAdd('name', egads.ATTRSTRING, 'top_shell')

# # Get the faces associated with the top shell, then
# # compute the range over each face to determine which face it is
# faces = getBodyFacesAndDirs(top_shell)
# print(faces)
# for i, f in enumerate(faces):
#     r1, r2 = f.getBoundingBox()
#     print("Top shell: face {0}".format(i))
#     print("xyz min = {0}, xyz max = {1}".format(r1, r2))
#     print("")

# Set the face attributes for the top shell volume
faces[3].attributeAdd('name', egads.ATTRSTRING, 'shell-top-face')
faces[1].attributeAdd('name', egads.ATTRSTRING, 'shell-bottom-face')

# Create the bottom outer ring
x1 = [0, 0, 6.35]
x2 = [0, 0, 0]
radius = 80.78
c3 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x1, x2, radius])

x1 = [0, 0, 6.35]
x2 = [0, 0, 0]
radius = 70.78
c4 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x1, x2, radius])

model = c3.solidBoolean(c4, egads.SUBTRACTION)
model.saveModel('ring.step', overwrite=True)
ring = model.getChildren()[0]
ring.attributeAdd('name', egads.ATTRSTRING, 'ring')

# faces = getBodyFacesAndDirs(ring)
# for i, f in enumerate(faces):
#     r1, r2 = f.getBoundingBox()
#     print("ring: face {0}".format(i))
#     print("xyz min = {0}, xyz max = {1}".format(r1, r2))
#     print("")
    
# Set the face attributes for the ring volume
faces[3].attributeAdd('name', egads.ATTRSTRING, 'ring-top-face')
faces[1].attributeAdd('name', egads.ATTRSTRING, 'ring-bottom-face')
faces[4].attributeAdd('name', egads.ATTRSTRING, 'ring-inner-face1')
faces[5].attributeAdd('name', egads.ATTRSTRING, 'ring-inner-face2')

# Create the bottom plate
x1 = [0, 0, 0]
x2 = [0, 0, 6.35]
radius = 70.78
plate = ctx.makeSolidBody(egads.CYLINDER, rdata=[x1, x2, radius])
plate.attributeAdd('name', egads.ATTRSTRING, 'plate')

# faces = getBodyFacesAndDirs(plate)
# for i, f in enumerate(faces):
#     r1, r2 = f.getBoundingBox()
#     print("plate: face {0}".format(i))
#     print("xyz min = {0}, xyz max = {1}".format(r1, r2))
#     print("")

# Set the face attributes for the plate volume
faces[2].attributeAdd('name', egads.ATTRSTRING, 'plate-top-face')
faces[3].attributeAdd('name', egads.ATTRSTRING, 'plate-bottom-face')
faces[1].attributeAdd('name', egads.ATTRSTRING, 'plate-outer-face1')
faces[0].attributeAdd('name', egads.ATTRSTRING, 'plate-outer-face2')

# assemble the components
model = top_shell.solidBoolean(ring, egads.FUSION)
body = model.getChildren()[0]
model = body.solidBoolean(plate, egads.FUSION)

# Get the fully fused motor
motor = model.getChildren()[0]
motor.attributeAdd('name', egads.ATTRSTRING, 'motor')

model.saveModel('model.step', overwrite=True)
model.saveModel('model.egads', overwrite=True)

# Load in the model to TMR
comm = MPI.COMM_WORLD
htarget = 5.0

geo = TMR.LoadModel('model.step', print_lev=1)
verts = geo.getVertices()
edges = geo.getEdges()
faces = geo.getFaces()
vols = geo.getVolumes()

print(len(verts))
print(len(edges))
print(len(faces))

for vol in vols:
    print(vol.getName())

geo = TMR.LoadModel('model.egads', print_lev=1)
verts = geo.getVertices()
edges = geo.getEdges()
faces = geo.getFaces()
vols = geo.getVolumes()

print(len(verts))
print(len(edges))
print(len(faces))

for vol in vols:
    print(vol.getName())

geo_new = TMR.Model(verts, edges, faces)

# Create the new mesh
mesh = TMR.Mesh(comm, geo_new)

# Set the meshing options
opts = TMR.MeshOptions()
opts.write_mesh_quality_histogram = 1
opts.triangularize_print_iter = 50000

# Create the surface mesh
mesh.mesh(htarget, opts)

# Write the surface mesh to a file
mesh.writeToVTK('output.vtk')
