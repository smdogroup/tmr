'''
Imports a STL file and applies a smoothing to it 

'''
import numpy
class smoothSTL:
    def __init__(self,fname):
        '''
        Pass in a STL file and returns a smoothed version of the structure
        '''
        self.fname = fname
        
    def readSTL(self, fname):
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
        fp = open(fname,  'rb')

        # Initialize variables for reading in STL file
        ind = 0
        k = -1
        # Normal of triangle
        norm = []
        # Points of triangle
        P1 = []
        P2 = []
        P3 = []

        # Looping through the STL file
        for line in fp:
            if ind == 0:
                Header = line        
            else:
                # Reading in the normal
                if k % 7 == 0:
                    norm_str = line[13:-1].split()
                    norm = numpy.append(norm,[float(x) for x in norm_str])
                # Reading in the 1st coordinate
                elif k % 7 == 2:
                    p1_str = line[7:-1].split()
                    P1 = numpy.append(P1,[float(x) for x in p1_str])
                # Reading in the 2nd coordinate
                elif k % 7 == 3:
                    p2_str = line[7:-1].split()
                    P2 = numpy.append(P2,[float(x) for x in p2_str])
                # Reading in the 3rd coordinate
                elif k % 7 == 4:
                    p3_str = line[7:-1].split()
                    P3 = numpy.append(P3,[float(x) for x in p3_str])

            k = k + 1
            ind = 1

        # Reshape coordinate array to correspond to number of elements row-wise
        norm = norm.reshape(len(norm)/3,3)
        P1 = P1.reshape(len(P1)/3,3)
        P2 = P2.reshape(len(P2)/3,3)
        P3 = P3.reshape(len(P3)/3,3)

        return norm, P1, P2, P3

    def createUniqueList(self, P1, P2, P3):
        '''
        Create unique list of nodes

        Input:
        P1:  Node 1 of triangle
        P2:  Node 2 of triangle
        P3:  Node 3 of triangle

        Output:
        unique_list: Unique list of nodes in structure
        conn:        Elemental connectivity
        node_conn:   Adjacency matrix
        '''
        # Tolerance for uniqueness
        tol = 1e-7
        conn = numpy.zeros([P1.shape[0],3])
        unique_list = []
        # Initialize the connectivity matrix
        conn[0,:] = [0, 1, 2]
        unique_list = P1[0,:]       
        unique_list = numpy.vstack((unique_list,P2[0,:]))
        unique_list = numpy.vstack((unique_list,P3[0,:]))
        # Loop over all coordinates for each element to search for unique nodes
        for row in xrange(1,P1.shape[0]):
            pt1 = P1[row,:]
            pt2 = P2[row,:]
            pt3 = P3[row,:]
            
            # Loop over all the present unique nodes
            for k in xrange(unique_list.shape[0]):
                # Not a unique node
                if numpy.linalg.norm(pt1-unique_list[k,:]) <= tol:
                    conn[row,0] = k
                    break
                # Loop through all nodes and it is unique
                elif k == unique_list.shape[0]-1:
                    conn[row,0] = unique_list.shape[0]
                    unique_list = numpy.vstack((unique_list, pt1))
            for k in xrange(unique_list.shape[0]):
                # Not a unique node
                if numpy.linalg.norm(pt2-unique_list[k,:]) <= tol:
                    conn[row,1] = k
                    break
                # Loop through all nodes and it is unique
                elif k == unique_list.shape[0]-1:
                    conn[row,1] = unique_list.shape[0]
                    unique_list = numpy.vstack((unique_list, pt2))
            for k in xrange(unique_list.shape[0]):
                # Not a unique node
                if numpy.linalg.norm(pt3-unique_list[k,:]) <= tol:
                    conn[row,2] = k
                    break
                # Loop through all nodes and it is unique
                elif k == unique_list.shape[0]-1:
                    conn[row,2] = unique_list.shape[0]
                    unique_list = numpy.vstack((unique_list, pt3))

        # Return node connectivity (adjacency matrix)
        node_conn = numpy.zeros([unique_list.shape[0],unique_list.shape[0]],dtype="intc")
        for k in xrange(conn.shape[0]):
            # Three edges in each triangle
            u = conn[k,0]
            v = conn[k,1]
            w = conn[k,2]
            node_conn[u][v] = 1
            node_conn[v][w] = 1
            node_conn[u][w] = 1
        # Ensure that the matrix is symmetric    
        node_conn = numpy.maximum(node_conn, node_conn.T)
        return unique_list, conn, node_conn

    def smoothMesh(self, unique_list, conn, node_conn, w = 0.5):
        '''
        Smooth out the mesh using Laplacian and returns the new nodal position
        Input:

        unique_list:     Unique list of nodes
        conn:            Elemental connectivity
        node_conn:       Nodal adjacency matrix
        w:               Weighting ratio
        '''
        unique_list_new = numpy.zeros([unique_list.shape[0],3])
        # Go through the adjacency matrix
        for k in xrange(unique_list.shape[0]):
            # Node k coordinates of interest
            node_k = unique_list[k,:]
            # Adjacent nodes coordinates to node k
            nodes_adj = unique_list[numpy.nonzero(node_conn[k,:]),:]
            # Number of adjacent nodes for node k
            N = nodes_adj.shape[1]
            # Perform Laplacian smoothing
            x_bar = numpy.sum(nodes_adj,axis=1)/N
            unique_list_new[k,:] = node_k + w*(x_bar-node_k)
        return unique_list_new

    def outputSTL(self):

        return 
    
fname = "trial.stl"
STL_new = smoothSTL(fname)
norm, P1, P2, P3 = STL_new.readSTL(fname)
unique_list, conn, node_conn = STL_new.createUniqueList(P1, P2, P3)
x_new = STL_new.smoothMesh(unique_list, conn, node_conn)
