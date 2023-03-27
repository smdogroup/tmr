import os
from egads4py import egads
from mpi4py import MPI
from tmr import TMR
import numpy as np
import argparse

# Load in the model to TMR
comm = MPI.COMM_WORLD

# Create an argument parser to read in arguments from the commnad line
p = argparse.ArgumentParser()
p.add_argument("--htarget", type=float, default=8.0)
p.add_argument("--extension", type=str, default="egads", help="egads or step")
args = p.parse_args()

htarget = args.htarget
extension = args.extension

# Create the egads context
ctx = egads.context()

# Create the top shell
x1 = [0, 0, 0]
x2 = [0, 0, 100.0]
radius = 100.0
c1 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x1, x2, radius])

x1 = [50, 0, 0]
x2 = [50, 0, 100.0]
radius = 75.0
c2 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x1, x2, radius])

model = c1.solidBoolean(c2, egads.SUBTRACTION)
model.saveModel("cylinder_sweep.%s" % (extension), overwrite=True)

# Set the meshing options
opts = TMR.MeshOptions()
opts.write_mesh_quality_histogram = 1
opts.triangularize_print_iter = 50000

# Load the separate geometries and mesh each
geo = TMR.LoadModel("cylinder_sweep.%s" % (extension))

faces = geo.getFaces()
vols = geo.getVolumes()
faces[3].setSource(vols[0], faces[1])

# Create the new mesh
mesh = TMR.Mesh(comm, geo)

# Set the meshing options
opts = TMR.MeshOptions()
opts.write_mesh_quality_histogram = 1
opts.triangularize_print_iter = 50000

# Create the surface mesh
mesh.mesh(htarget, opts)

# Write the surface mesh to a file
mesh.writeToVTK("cylinder_sweep.vtk", "hex")

# Create the model from the unstructured volume mesh
model = mesh.createModelFromMesh()

# Create the corresponding mesh topology from the mesh-model
topo = TMR.Topology(comm, model)

# Create the quad forest and set the topology of the forest
forest = TMR.OctForest(comm)
forest.setTopology(topo)

# Create random trees and balance the mesh. Print the output file
forest.createRandomTrees(nrand=3, max_lev=2)
forest.balance(1)
filename = "cylinder_forest%d.vtk" % (comm.rank)
forest.writeForestToVTK(filename)
