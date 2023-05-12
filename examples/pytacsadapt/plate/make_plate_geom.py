"""
A function to generate the plate geometry used for testing
"""

# imports
import numpy as np
from egads4py import egads


def makePlateGeom(width=1.0, height=1.0, npanels=1, makeIGES=False):
    """
    Writes a plate.step/.iges geometry file for the planar plate model of the
    following form:
     __________ __________
    |          |          |
    |          |          |    ^ y
    |          |          |    |
    |__________|__________|    ---> x

    Parameters
    ----------
    width : float
        the length of the plate geometry along the x-axis

    height : float
        the length of the plate geometry along the y-axis

    npanels : int
        the number of panels used define the plate geometry along the x-dimension

    makeIGES : bool
        boolean flag on whether to return an IGES/STEP geometry
    """
    # get the model size info
    nverts = 2 * (npanels + 1)
    nedges = (3 * npanels) + 1
    nfaces = npanels

    # get the vertex coordinates
    nx = npanels + 1
    ny = 2
    x = np.linspace(0, width, nx)
    y = np.linspace(0, height, ny)
    X, Y = np.meshgrid(x, y)
    coords = np.zeros([nverts, 3])
    coords[:, 0] = X.ravel()
    coords[:, 1] = Y.ravel()
    # print(coords)

    # make the quad node connectivity - used to make the edge and face connectivity
    quad_conn = []
    for iface in range(nfaces):
        i = iface % 4
        j = iface // 4
        quad_conn.append(
            [i + nx * j, (i + 1) + nx * j, (i + 1) + nx * (j + 1), i + nx * (j + 1)]
        )
    quad_conn = np.array(quad_conn)
    # print(quad_conn)

    # make the edge connectivity
    edge_conn = []
    for iedge in range(nfaces * 4):
        ipanel = iedge // 4
        local_edge = iedge % 4
        if ipanel > 0 and local_edge == 3:
            continue
        v1 = quad_conn[ipanel, local_edge]
        if local_edge == 3:
            v2 = quad_conn[ipanel, 0]
        else:
            v2 = quad_conn[ipanel, local_edge + 1]
        edge_conn.append([v1, v2])
    edge_conn = np.array(edge_conn)
    # print(edge_conn)

    # make the face connectivity
    face_conn = []
    for iface in range(nfaces):
        conn = []
        for iedge in range(4):
            v1 = quad_conn[iface, iedge]
            if iedge == 3:
                v2 = quad_conn[iface, 0]
            else:
                v2 = quad_conn[iface, iedge + 1]
            edge = np.array([v1, v2])
            edge_ind = np.where((edge_conn == edge).all(axis=1))[0]
            if edge_ind.size == 0:
                edge_ind = np.where((edge_conn == np.flip(edge)).all(axis=1))[0]
            conn.append(edge_ind[0])
        face_conn.append(conn)
    face_conn = np.array(face_conn)
    # print(face_conn)

    # create egads
    ctx = egads.context()
    ctx.setOutLevel(0)

    # create the node topology
    nodes = []
    for i in range(nverts):
        nodes.append(ctx.makeTopology(egads.NODE, rdata=coords[i]))

    # create the line geometry
    lines = []
    for i in range(nedges):
        n1_ind = edge_conn[i, 0]
        n2_ind = edge_conn[i, 1]
        delta = coords[n2_ind] - coords[n1_ind]
        lines.append(
            ctx.makeGeometry(
                egads.CURVE, mtype=egads.LINE, rdata=[coords[n1_ind], delta]
            )
        )

    # create the edge topology
    edges = []
    for i in range(nedges):
        n1_ind = edge_conn[i, 0]
        n2_ind = edge_conn[i, 1]
        delta = coords[n2_ind] - coords[n1_ind]
        dist = np.linalg.norm(delta, 2)
        edges.append(
            ctx.makeTopology(
                egads.EDGE,
                mtype=egads.TWONODE,
                geom=lines[i],
                children=[nodes[n1_ind], nodes[n2_ind]],
                rdata=[0, dist],
            )
        )

    # create the edge loops
    edge_loops = []
    edge_loop_nums = []
    for i in range(nfaces):
        e1_ind = face_conn[i, 0]
        e2_ind = face_conn[i, 1]
        e3_ind = face_conn[i, 2]
        e4_ind = face_conn[i, 3]
        eloop, nloop_edges = ctx.makeLoop(
            [edges[e1_ind], edges[e2_ind], edges[e3_ind], edges[e4_ind]]
        )
        edge_loops.append(eloop)
        edge_loop_nums.append(nloop_edges)

    # create the faces
    faces = []
    for i in range(nfaces):
        faces.append(ctx.makeFace(edge_loops[i], egads.SFORWARD))

    # piece it all together and make the model
    shell = ctx.makeTopology(egads.SHELL, egads.OPEN, children=faces)
    body = ctx.makeTopology(egads.BODY, egads.SHEETBODY, children=[shell])
    model = ctx.makeTopology(egads.MODEL, children=[body])
    fname = "plate"
    if makeIGES:
        fname += ".iges"
    else:
        fname += ".step"
    model.saveModel(fname, overwrite=True)
    return

makePlateGeom()
