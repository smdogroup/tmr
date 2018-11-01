'''
Imports a STL file and applies a smoothing to it
'''
from __future__ import print_function

# Import the locate point code
from tmr import TMR
import numpy as np
from itertools import islice
import time
import argparse

def readSTLFile(fname):
    '''
    Reads in the STL file

    Input:
    fname:    STL filename

    Output:
    norm:     Norm of triangle
    P1:       Node 1 of triangle
    P2:       Node 2 of triangle
    P3:       Node 3 of triangle
    '''

    # Normal of triangle
    norm = []

    # Points of triangle
    P1 = []
    P2 = []
    P3 = []

    with open(fname, 'r') as fp:
        # Initialize variables for reading in STL file
        ind = 0
        k = -1

        while True:
            lines = list(islice(fp, 5000))
            if not lines:
                break

            for line in lines:
                # Looping through the STL file
                if ind == 0:
                    Header = line
                else:
                    # Reading in the normal
                    if k % 7 == 0:
                        norm_str = line[13:-1].split()
                        try:
                            norm.extend([float(x) for x in norm_str])
                        except:
                            pass
                    # Reading in the 1st coordinate
                    elif k % 7 == 2:
                        p1_str = line[7:-1].split()
                        P1.extend([float(x) for x in p1_str])
                    # Reading in the 2nd coordinate
                    elif k % 7 == 3:
                        p2_str = line[7:-1].split()
                        P2.extend([float(x) for x in p2_str])
                    # Reading in the 3rd coordinate
                    elif k % 7 == 4:
                        p3_str = line[7:-1].split()
                        P3.extend([float(x) for x in p3_str])

                k = k + 1
                ind = 1

    # Reshape coordinate array to correspond to number of elements
    # row-wise
    norm = np.array(norm).reshape(len(norm)/3,3)
    P1 = np.array(P1).reshape(len(P1)/3,3)
    P2 = np.array(P2).reshape(len(P2)/3,3)
    P3 = np.array(P3).reshape(len(P3)/3,3)

    return norm, P1, P2, P3

def createUniqueList(P1, P2, P3, tol=1e-5):
    '''
    Create unique list of nodes

    Input:
    P1:  Node 1 of triangle
    P2:  Node 2 of triangle
    P3:  Node 3 of triangle

    Output:
    unique_nodes: Unique list of nodes in structure
    conn:         Elemental connectivity
    node_conn:    Adjacency matrix
    '''

    # Tolerance for uniqueness
    Xpts = np.vstack((P1, P2, P3))
    loc = TMR.PointLocator(Xpts)
    node_nums = -np.ones(Xpts.shape[0], dtype='intc')

    # Locate the closest K points
    K = 20
    index = np.zeros(K, dtype='intc')
    dist = np.zeros(K)

    unique_node = 0
    for row in range(Xpts.shape[0]):
        if node_nums[row] < 0:
            # Locate the closest points and label them with
            # the same index
            loc.locateClosest(Xpts[row,:], index, dist)
            for k in range(K):
                if np.sqrt(dist[k]) < tol:
                    node_nums[index[k]] = unique_node
                else:
                    break

            # If we ordered one node, increment the counter
            if dist[0] < tol:
                unique_node += 1

    # Create the unique list of nodes
    unique_nodes = np.zeros((unique_node, 3))
    for row in range(Xpts.shape[0]):
        unique_nodes[node_nums[row], :] = Xpts[row, :]

    # Create the connectivity
    conn = np.zeros((P1.shape[0], 3), dtype='intc')
    for row in range(P1.shape[0]):
        conn[row, 0] = node_nums[row]
        conn[row, 1] = node_nums[row + P1.shape[0]]
        conn[row, 2] = node_nums[row + 2*P1.shape[0]]

    # Return node connectivity (adjacency matrix)
    node_conn = [[] for x in range(unique_node)]

    # Loop over each triangle and add the connectivity
    for k in range(conn.shape[0]):
        u = conn[k,0]
        v = conn[k,1]
        w = conn[k,2]

        if u < v:
            node_conn[u].append(v)
        if v < w:
            node_conn[v].append(w)
        if u < w:
            node_conn[u].append(w)

    return unique_nodes, conn, node_conn

def smoothMesh(unique_nodes, conn, node_conn, w=0.5):
    '''
    Smooth out the mesh using Laplacian and returns the new nodal position
    Input:

    unique_nodes:     Unique list of nodes
    conn:            Elemental connectivity
    node_conn:       Nodal adjacency matrix
    w:               Weighting ratio
    '''
    unique_nodes_new = np.zeros([unique_nodes.shape[0],3])

    # Go through the adjacency matrix
    for k in range(unique_nodes.shape[0]):
        # Node k coordinates of interest
        node_k = unique_nodes[k,:]

        # Adjacent nodes coordinates to node k
        nodes_adj = unique_nodes[node_conn[k],:]

        # Number of adjacent nodes for node k
        N = len(node_conn[k])
        if N > 0:
            # Perform Laplacian smoothing
            x_bar = np.sum(nodes_adj, axis=0)/N
            unique_nodes_new[k,:] = node_k + w*(x_bar - node_k)
        else:
            unique_nodes_new[k,:] = node_k

    return unique_nodes_new

def outputSTLFile(unique_nodes, conn, fname):
    '''
    Output the new STL filename
    '''
    s = 'solid topology\n'
    # Loop over all elements
    for k in range(conn.shape[0]):
        # Calculate the normal
        pt1 = unique_nodes[conn[k,0],:]
        pt2 = unique_nodes[conn[k,1],:]
        pt3 = unique_nodes[conn[k,2],:]

        # Find two edges of the triangle
        u = pt2-pt1
        v = pt3-pt1
        Norm = np.array([u[1]*v[2]-u[2]*v[1],
                         u[2]*v[0]-u[0]*v[2],
                         u[0]*v[1]-u[1]*v[0]])

        # Write the normal and three vertices
        s += 'facet normal %e %e %e\nouter loop\n'%(
            Norm[0], Norm[1], Norm[2])
        s += 'vertex %e %e %e\n'%(pt1[0], pt1[1], pt1[2])
        s += 'vertex %e %e %e\n'%(pt2[0], pt2[1], pt2[2])
        s += 'vertex %e %e %e\n'%(pt3[0], pt3[1], pt3[2])
        s += 'endloop\nendfacet\n'

    s += 'endsolid topology\n'

    fp = open(fname, 'w')
    fp.write(s)
    fp.close()
    return

def smoothSTLFile(infile, outfile):
    # Create the smoothing object
    t1 = time.time()
    norm, P1, P2, P3 = readSTLFile(infile)
    t2 = time.time()

    # Create a unique list of nodes
    unique_nodes, conn, node_conn = createUniqueList(P1, P2, P3)
    t3 = time.time()

    # Perform the smoothing
    x_new = smoothMesh(unique_nodes, conn, node_conn)
    t4 = time.time()

    # Write out the new mesh
    outputSTLFile(x_new, conn, outfile)
    t5 = time.time()

    print('Time to extract nodes:        %10.5e'%(t2-t1))
    print('Time to extract connectivity: %10.5e'%(t3-t2))
    print('Time to smooth mesh:          %10.5e'%(t4-t3))
    print('Time to output STL:           %10.5e'%(t5-t4))

    return

# Define the performance profile objective function
parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, default='input.stl',
                    help='Input stl file name')
parser.add_argument('--output', type=str, default='output.stl',
                    help='Output stl file name')
args = parser.parse_args()

# Assign the input and output names
smoothSTLFile(args.input, args.output)
