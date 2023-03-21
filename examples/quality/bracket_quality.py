from __future__ import print_function
from mpi4py import MPI
from tmr import TMR
from egads4py import egads
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
    # Create the egads context
    ctx = egads.context()

    parts = []

    r0 = 0.05

    # Create the boxes
    x0 = [0, 0, 0]
    x1 = [0.25, 0.25, 0.25]
    B0 = ctx.makeSolidBody(egads.BOX, rdata=[x0, x1])
    parts.append(ctx.makeTopology(egads.MODEL, children=[B0]))

    # Create the x-arm
    x0 = [0.25, 0, 0]
    x1 = [0.75, 0.25, 0.25]
    B1 = ctx.makeSolidBody(egads.BOX, rdata=[x0, x1])

    x0 = [0.85, 0.125, 0]
    x1 = [0.85, 0.125, 0.25]
    C1 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x0, x1, r0])
    parts.append(B1.solidBoolean(C1, egads.SUBTRACTION))

    # Create the y-arm
    x0 = [0, 0.25, 0]
    x1 = [0.25, 0.75, 0.25]
    B2 = ctx.makeSolidBody(egads.BOX, rdata=[x0, x1])

    x0 = [0, 0.85, 0.125]
    x1 = [0.25, 0.85, 0.125]
    C2 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x0, x1, r0])
    parts.append(B2.solidBoolean(C2, egads.SUBTRACTION))

    # Create the z-arm
    x0 = [0, 0, 0.25]
    x1 = [0.25, 0.25, 0.75]
    B3 = ctx.makeSolidBody(egads.BOX, rdata=[x0, x1])

    x0 = [0.125, 0, 0.85]
    x1 = [0.125, 0.25, 0.85]
    C3 = ctx.makeSolidBody(egads.CYLINDER, rdata=[x0, x1, r0])
    parts.append(B3.solidBoolean(C3, egads.SUBTRACTION))

    # Create all of the models
    geos = []
    for p in parts:
        geos.append(TMR.ConvertEGADSModel(p))

    # Create the full list of vertices, edges, faces and volumes
    verts = []
    edges = []
    faces = []
    vols = []
    for geo in geos:
        verts.extend(geo.getVertices())
        edges.extend(geo.getEdges())
        faces.extend(geo.getFaces())
        vols.extend(geo.getVolumes())

    # Set all of the matching faces
    TMR.setMatchingFaces(geos)

    # Create the geometry
    geo = TMR.Model(verts, edges, faces, vols)

    return geo


# The communicator
comm = MPI.COMM_WORLD

# Load the geometry model
geo = load_model()

# Create the mesh
mesh = TMR.Mesh(comm, geo)

# Set the meshing options
opts = TMR.MeshOptions()
opts.frontal_quality_factor = 1.25
opts.num_smoothing_steps = 50
opts.triangularize_print_iter = 5
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
writeQualityToVtk(forest, ar, min_ang, fshape, fname="quality-bracket.vtk")
plotShapeHist(fshape, xmin=np.amin(fshape), fname="bracket_shape_hist.pdf")
plotARHist(ar, fname="bracket_ar_hist.pdf", xmax=np.ceil(4.0 * np.amax(ar)) / 4.0)
plotMinAngHist(min_ang, fname="bracket_ang_hist.pdf", xmin=np.amin(min_ang) - 1.0)
