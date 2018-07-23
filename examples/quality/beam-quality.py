from __future__ import print_function
from mpi4py import MPI
from tmr import TMR
from paropt import ParOpt
from tacs import TACS, elements, constitutive, functions
import numpy as np
import argparse
import os
from OctMeshQuality import *

def get_edge_dirs_verts(elist):
    edge_list = elist[:]

    edges = [edge_list.pop()]
    dirs = [1]
    v1, vnext = edges[-1].getVertices()
    verts = [v1, vnext]

    nedges = len(edge_list)
    for k in range(nedges):
        for i, edge in enumerate(edge_list):
            v1, v2 = edge.getVertices()        
            if v1.getEntityId() == vnext.getEntityId():
                dirs.append(1)
                edges.append(edge_list.pop(i))
                vnext = v2
                break
            elif v2.getEntityId() == vnext.getEntityId():
                dirs.append(-1)
                edges.append(edge_list.pop(i))
                vnext = v1
                break
        verts.append(vnext)

    return edges, dirs, verts[:-1]

def load_model():
    geo = TMR.LoadModel('beam.step')

    # Get the faces/volume from the model
    faces = geo.getFaces()
    edges = geo.getEdges()
    verts = geo.getVertices()
    vols = geo.getVolumes()

    # Set the attributes
    faces[6].setAttribute('hole')
    faces[7].setAttribute('hole')
    faces[0].setAttribute('fixed')

    # Set the source/destination faces
    faces[1].setSource(vols[0], faces[3])

    return geo

# The communicator
comm = MPI.COMM_WORLD

# Load the geometry model
geo = load_model()

# Set the boundary conditions for the problem
bcs = TMR.BoundaryConditions()
bcs.addBoundaryCondition('fixed')

# Create the mesh
mesh = TMR.Mesh(comm, geo)

# Set the meshing options
opts = TMR.MeshOptions()
opts.frontal_quality_factor = 1.25
opts.num_smoothing_steps = 50
opts.triangularize_print_iter = 50000
opts.write_mesh_quality_histogram = 1

# Create the surface mesh
mesh.mesh(0.02, opts)

# Create a model from the mesh
model = mesh.createModelFromMesh()

# Create the corresponding mesh topology from the mesh-model 
topo = TMR.Topology(comm, model)

# Create the quad forest and set the topology of the forest
forest = TMR.OctForest(comm)
forest.setTopology(topo)
forest.createTrees(1)
forest.createNodes()

ar = computeAR(forest)
min_ang = computeMinAngle(forest)
fshape = computeShape(forest)

# Wrtie the mesh quality to vtk
writeQualityToVtk(forest, ar, min_ang, fshape)
plotShapeHist(fshape, xmin=np.amin(fshape), fname='beam_shape_hist.pdf')
plotARHist(ar, fname='beam_ar_hist.pdf', xmax=np.ceil(4.0*np.amax(ar))/4.0)
plotMinAngHist(min_ang, xmin=np.amin(min_ang)-1.0, fname='beam_ang_hist.pdf')
