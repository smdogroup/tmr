from mpi4py import MPI
from tmr import TMR
import argparse
import os

# Set the communicator
comm = MPI.COMM_WORLD

# Create an argument parser to read in arguments from the commnad line
p = argparse.ArgumentParser()
p.add_argument('--htarget', type=float, default=4.0)
p.add_argument('--filename', type=str, default=None, help='STEP file name')
args = p.parse_args()

# Get the value of the filename
filename = args.filename
if not os.path.isfile(filename):
    raise ValueError('File %s does not exist'%(filename))

# Set the value of the target length scale in the mesh
htarget = args.htarget

# Load the geometry model
geo = TMR.LoadModel(filename)

# Create a model by discarding the volumes
verts = geo.getVertices()
edges = geo.getEdges()
faces = geo.getFaces()
geo_new = TMR.Model(verts, edges, faces)

# Create the new mesh
mesh = TMR.Mesh(comm, geo_new)

# Set the meshing options
opts = TMR.MeshOptions()
# opts.frontal_quality_factor = 1.25
opts.num_smoothing_steps = 2
opts.write_mesh_quality_histogram = 1
opts.write_post_smooth_triangle = 1
opts.write_pre_smooth_quad = 1
opts.write_post_smooth_quad = 1

# Create the surface mesh
mesh.mesh(htarget, opts)

# Write the surface mesh to a file
mesh.writeToVTK('surface-mesh.vtk')