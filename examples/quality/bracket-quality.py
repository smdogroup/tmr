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
    geo = TMR.LoadModel('model.step')

    # Get the faces/volume from the model
    faces = geo.getFaces()
    edges = geo.getEdges()
    verts = geo.getVertices()

    # Create the edge loops
    elist = [edges[3], edges[5], edges[45], edges[32]]
    ex, dx, vx = get_edge_dirs_verts(elist)

    elist = [edges[45], edges[6], edges[7], edges[51]]
    ey, dy, vy = get_edge_dirs_verts(elist)

    elist = [edges[2], edges[32], edges[51], edges[8]]
    ez, dz, vz = get_edge_dirs_verts(elist)

    # Create the faces
    fx = TMR.TFIFace(ex, dx, vx)
    fy = TMR.TFIFace(ey, dy, vy)
    fz = TMR.TFIFace(ez, dz, vz)
    faces.extend([fx, fy, fz])

    # Make the volumes
    s1 = [fx, faces[3], faces[19], faces[6], faces[10],
          faces[9], faces[11], faces[12]]
    d1 = [1, -1, -1, -1, 1, 1, -1, -1]

    s2 = [fy, faces[16], faces[8], faces[5], faces[23],
          faces[15], faces[17], faces[18]]
    d2 = [-1, 1, -1, -1, 1, 1, -1, -1]

    s3 = [fz, faces[4], faces[20], faces[13], faces[7],
          faces[14], faces[21], faces[22]]
    d3 = [1, -1, 1, 1, 1, 1, -1, -1]

    s4 = [fx, fy, fz, faces[0], faces[1], faces[2]]
    d4 = [-1, 1, -1, -1, -1, -1]

    # Set the attributes
    faces[11].setAttribute('fx')
    faces[12].setAttribute('fx')
    faces[17].setAttribute('fy')
    faces[18].setAttribute('fy')
    faces[21].setAttribute('fz')
    faces[22].setAttribute('fz')

    # Form the 4 independent bodies that are connected through
    v1 = TMR.Volume(s1, d1)
    v2 = TMR.Volume(s2, d2)
    v3 = TMR.Volume(s3, d3)
    v4 = TMR.Volume(s4, d4)
    vols = [v1, v2, v3, v4]

    # Set the source/destination faces
    faces[19].setSource(v1, faces[3])
    faces[5].setSource(v2, faces[23])
    faces[7].setSource(v3, faces[13])
    faces[1].setSource(v4, fz)

    # Create a new model
    geo = TMR.Model(verts, edges, faces, vols)

    return geo

# The communicator
comm = MPI.COMM_WORLD

# Load the geometry model
geo = load_model()

# Set the boundary conditions for the problem
bcs = TMR.BoundaryConditions()
bcs.addBoundaryCondition('fz')

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

# Wrtie the mesh quality to vtk
writeQualityToVtk(forest, ar, min_ang)
plotARHist(ar, fname='bracket_ar_hist.pdf', xmax=np.ceil(4.0*np.amax(ar))/4.0)
plotMinAngHist(min_ang, fname='bracket_ang_hist.pdf', xmin=np.amin(min_ang)-1.0)
