from mpi4py import MPI
from tmr import TMR
import argparse

# Create the communicator
comm = MPI.COMM_WORLD

# Create an argument parser to read in command-line arguments
p = argparse.ArgumentParser()
p.add_argument('--reverse', default=False, action='store_true')
args = p.parse_args()

# Load the model from the STEP file
geo = TMR.LoadModel('first-section.stp')

# Get the volumes
vols = geo.getVolumes()

# Get the edges/faces from the geometry
faces = geo.getFaces()
edges = geo.getEdges()

# Set the source/target relationships
if args.reverse:
    faces[4].setSource(vols[0], faces[5])
else:
    faces[5].setSource(vols[0], faces[4])
edges[8].setSource(edges[5])

# Create the geometry
mesh = TMR.Mesh(comm, geo)

# Mesh the part
opts = TMR.MeshOptions()
opts.num_smoothing_steps = 10

# Mesh the geometry with the given target size
htarget = 4.0
mesh.mesh(htarget, opts=opts)

# Write the mesh to a bdf file
mesh.writeToBDF('volume-mesh.bdf', 'hex')
