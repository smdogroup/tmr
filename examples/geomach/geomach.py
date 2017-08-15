'''
This example demonstrates the conversion of a GeoMACH model to TMR.

GeoMACH is an open source geometry engine developed by John Hwang. 
The model is taken from the examples and was developed by John Hwang
and Davide Ivaldi.
'''

from __future__ import division

from mpi4py import MPI
from tmr import TMR
import numpy as np

from GeoMACH.PGM.core import PGMconfiguration, PGMparameter, PGMdv
from GeoMACH.PGM.components import PGMwing, PGMbody, PGMshell
from GeoMACH.PGM.components import PGMjunction, PGMtip, PGMcone

class Trussbraced(PGMconfiguration):

    def _define_comps(self):
        self.comps['fuse'] = PGMbody(num_x=17, num_y=6, num_z=4)
        self.comps['lwing'] = PGMwing(num_x=7, num_z=7, left_closed=True)
        self.comps['lstrut'] = PGMwing(num_x=4, num_z=4, left_closed=True, right_closed=True)
        self.comps['lv'] = PGMwing(num_z=4)
        self.comps['ltail'] = PGMwing(left_closed=True)
        self.comps['vtail'] = PGMwing(num_x=5, num_z=4, left_closed=True)

        self.comps['fuse_f'] = PGMcone(self, 'fuse', 'front', 18)
        self.comps['fuse_r'] = PGMcone(self, 'fuse', 'rear', 2)
        self.comps['lwing_t'] = PGMtip(self, 'lwing', 'left', 0.1)
        self.comps['ltail_t'] = PGMtip(self, 'ltail', 'left', 0.1)
        self.comps['vtail_t'] = PGMtip(self, 'vtail', 'left', 0.1)
        self.comps['lwing_fuse'] = PGMjunction(self, 'fuse', 'lft', 'E', [0,1], 'lwing', 'right', fweight=4, mweight=2)
        self.comps['lstrut_fuse'] = PGMjunction(self, 'fuse', 'lft', 'E', [4,2], 'lstrut', 'right')
        self.comps['lstrut_lwing'] = PGMjunction(self, 'lwing', 'low', 'S', [4,1], 'lstrut', 'left', fweight=3, mweight=3)
        self.comps['lv_lwing'] = PGMjunction(self, 'lwing', 'low', 'S', [1,3], 'lv', 'left')
        self.comps['lv_lstrut'] = PGMjunction(self, 'lstrut', 'upp', 'S', [1,0], 'lv', 'right')
        self.comps['vtail_fuse'] = PGMjunction(self, 'fuse', 'top', 'E', [1,10], 'vtail', 'right')
        self.comps['ltail_vtail'] = PGMjunction(self, 'vtail', 'low', 'N', [0,1], 'ltail', 'right')

    def _define_params(self):
        fuse = self.comps['fuse'].props
        fuse['nor'].params[''] = PGMparameter(1, 3)
        fuse['pos'].params[''] = PGMparameter(2, 3)
        fuse['pos'].params['nose'] = PGMparameter(2, 3, pos_u=[0,0.12])
        fuse['pos'].params['tail'] = PGMparameter(2, 3, pos_u=[0.76,1.0])
        fuse['scl'].params['rad1'] = PGMparameter(4, 1, pos_u=[0,0.01,0.05,0.12], order_u=4)
        fuse['scl'].params['rad2'] = PGMparameter(2, 1, pos_u=[0.12,0.76])
        fuse['scl'].params['rad3'] = PGMparameter(4, 1, pos_u=[0.76,0.83,0.99,1], order_u=4)
        fuse['scl'].params['tail'] = PGMparameter(2, 3, pos_u=[0.76,1.0])
        fuse['flt'].params['flt1a'] = PGMparameter(4, 2, pos_u=[0.24,0.27,0.33,0.36], pos_v=[0.5,1], order_u=4)
        fuse['flt'].params['flt1b'] = PGMparameter(2, 2, pos_u=[0.36,0.41], pos_v=[0.5,1])
        fuse['flt'].params['flt1c'] = PGMparameter(4, 2, pos_u=[0.41,0.44,0.49,0.52], pos_v=[0.5,1], order_u=4)
        fuse['flt'].params['flt2a'] = PGMparameter(4, 2, pos_u=[0.24,0.27,0.33,0.36], pos_v=[0,0.5], order_u=4)
        fuse['flt'].params['flt2b'] = PGMparameter(2, 2, pos_u=[0.36,0.41], pos_v=[0,0.5])
        fuse['flt'].params['flt2c'] = PGMparameter(4, 2, pos_u=[0.41,0.44,0.49,0.52], pos_v=[0,0.5], order_u=4)

        lwing = self.comps['lwing'].props
        lwing['pos'].params[''] = PGMparameter(1, 3)
        lwing['scl'].params[''] = PGMparameter(2, 1)
        lwing['pos'].params['lin'] = PGMparameter(3, 3, order_u=3)

        lstrut = self.comps['lstrut'].props
        lstrut['pos'].params[''] = PGMparameter(1, 3)
        lstrut['pos'].params['lin'] = PGMparameter(2, 3)
        lstrut['scl'].params[''] = PGMparameter(2, 1)
        lstrut['nor'].params[''] = PGMparameter(1, 1)

        lv = self.comps['lv'].props
        lv['pos'].params[''] = PGMparameter(1, 3)
        lv['pos'].params['lin'] = PGMparameter(2, 3)
        lv['scl'].params[''] = PGMparameter(2, 1)
        lv['nor'].params[''] = PGMparameter(1, 1)
        lv['rot'].params[''] = PGMparameter(2, 3, pos_u=[0,1])

        ltail = self.comps['ltail'].props
        ltail['pos'].params[''] = PGMparameter(1, 3)
        ltail['pos'].params['lin'] = PGMparameter(2, 3)
        ltail['scl'].params[''] = PGMparameter(2, 1)

        vtail = self.comps['vtail'].props
        vtail['pos'].params[''] = PGMparameter(1, 3)
        vtail['pos'].params['lin'] = PGMparameter(2, 3)
        vtail['scl'].params[''] = PGMparameter(2, 1)
        vtail['nor'].params[''] = PGMparameter(1, 3)

    def _compute_params(self):
        fuse = self.comps['fuse'].props
        fuse['nor'].params[''].val([1.0,0.0,1.0])
        fuse['pos'].params[''].val([[0,0,0],[36,0,0]])
        fuse['pos'].params['nose'].val([[0,-0.4,0],[0,0,0]])
        fuse['pos'].params['tail'].val([[0,0,0],[0,1.6,0]])
        fuse['scl'].params['rad1'].val([1,1.2,1.9,2])
        fuse['scl'].params['rad2'].val([2,2])
        fuse['scl'].params['rad3'].val([2,1.7,0.6,0.4])
        fuse['scl'].params['tail'].val([[0,0,0],[-0.3,0,0]])
        fuse['flt'].params['flt1a'].val([[0,0],[0,0],[0.6,0.6],[0.6,0.6]])
        fuse['flt'].params['flt1b'].val([[0.6,0.6],[0.6,0.6]])
        fuse['flt'].params['flt1c'].val([[0.6,0.6],[0.6,0.6],[0,0],[0,0]])
        fuse['flt'].params['flt2a'].val([[0,0],[0,0],[1,1],[1,1]])
        fuse['flt'].params['flt2b'].val([[1,1],[1,1]])
        fuse['flt'].params['flt2c'].val([[1,1],[1,1],[0,0],[0,0]])

        lwing = self.comps['lwing'].props
        lwing['pos'].params[''].val([12,1.7,2.8])
        lwing['scl'].params[''].val([3.6,0.8])
        lwing['pos'].params['lin'].val([[0,0,0],[2.5,-0.1,11],[5,-0.8,22]])

        lstrut = self.comps['lstrut'].props
        lstrut['pos'].params[''].val([13.4,-1.6,2.6])
        lstrut['pos'].params['lin'].val([[0,0,0],[1.6,2.6,11.8]])
        lstrut['scl'].params[''].val([1.8,1.6])
        lstrut['nor'].params[''].val([1.0])

        lv = self.comps['lv'].props
        lv['pos'].params[''].val([14.3,-0.12,8.8])
        lv['pos'].params['lin'].val([[0,0,0],[0,1.58,0]])
        lv['scl'].params[''].val([1.1,1.1])
        lv['nor'].params[''].val([1.0])
        lv['rot'].params[''].val([[0,2,0],[0,-2,0]])

        ltail = self.comps['ltail'].props
        ltail['pos'].params[''].val([35.3,6.6,0.25])
        ltail['pos'].params['lin'].val([[0,0,0],[2.6,0,5]])
        ltail['scl'].params[''].val([3.3,1])

        vtail = self.comps['vtail'].props
        vtail['pos'].params[''].val([30.7,2.1,0])
        vtail['pos'].params['lin'].val([[0,0,0],[4.6,5,0]])
        vtail['scl'].params[''].val([5,4.5])
        vtail['nor'].params[''].val([1.0,0.0,0.0])

        return [], [], []

    def _set_bspline_options(self):
        comps = self.comps

        comps['fuse'].faces['lft'].set_option('num_cp', 'u', [4,4,14,14,4,4])
        comps['fuse'].faces['rgt'].set_option('num_cp', 'v', [85,4,4,4,4,4,4,4,4,4,102,4,4,16,8,4,6])
        comps['vtail'].faces['low'].set_option('num_cp', 'u', [6,4,30,4,4])
        comps['vtail'].faces['low'].set_option('num_cp', 'v', [10,10,10,4])
        comps['lwing'].faces['upp'].set_option('num_cp', 'v', [20,4,4,20,5,4,31])
        comps['lwing'].faces['low'].set_option('num_cp', 'u', [12,12,20,4,4,4,4])
        comps['lstrut'].faces['upp'].set_option('num_cp', 'u', [4,8,12,4])
        comps['lstrut'].faces['upp'].set_option('num_cp', 'v', [4,5,4,4])
        
def geomach_to_tmr(bse):
    '''
    Convert a BSE-GeoMach model to a TMR model
    '''

    # Compute the number of vertices, edges and faces
    num = bse._num
    nverts = num['vert']
    nedges = num['edge']
    nfaces = num['surf']

    # Create the list of edges, faces and 
    verts = nverts*[ None ]
    edges = nedges*[ None ]
    faces = nfaces*[ None ]
    surfs = nfaces*[ None ]

    # Set the topology object
    topo = bse._topo
    size = bse._size
    str_indices = bse._str_indices
    bspline = bse._bspline

    # Extract the control point array
    cp_str = bse.vec['cp_str'].array

    # Point from the edge index on each face to the i/j locations
    # in the surf_ptrs array
    face_to_edge = [[1, 0], [2, 1], [1, 2], [0, 1]]

    # Point from the edge index on each face to the corresponding
    # vertices on each face
    face_to_edge_verts = [
        [[0, 0], [2, 0]],
        [[2, 0], [2, 2]],
        [[0, 2], [2, 2]],
        [[0, 0], [0, 2]]]

    # Create the vertices
    surf_ptrs = topo['surf_ptrs']
    edge_ptrs = topo['edge_ptrs']
    for i in range(nfaces):
        cp_offset = str_indices['cp'][i,0]
        ku = bspline['order'][topo['surf_group'][i,0]-1]
        kv = bspline['order'][topo['surf_group'][i,1]-1]
        nu = bspline['num_cp'][topo['surf_group'][i,0]-1]
        nv = bspline['num_cp'][topo['surf_group'][i,1]-1]

        for jj in range(2):
            for ii in range(2):
                # Get the vertex index
                index = surf_ptrs[i, 2*ii, 2*jj]-1

                # If the vertex has not be allocated, create it now
                if verts[index] is None:
                    # Find the corner corresponding to the vertex location
                    i1 = ii*(nu-1)
                    j1 = jj*(nv-1)
                    cp_index = cp_offset + i1 + j1*nu
                    pt = cp_str[cp_index]
                    verts[index] = TMR.VertexFromPoint(pt)

        # Extract and create the b-spline surfaces
        cp = np.zeros((nu, nv, 3), dtype=np.double)
        for jj in xrange(nv):
            for ii in xrange(nu):
                cp_index = cp_offset + ii + jj*nu
                cp[ii, jj, :] = cp_str[cp_index]

        surfs[i] = TMR.BsplineSurface(cp, ku=ku, kv=kv)
        faces[i] = TMR.FaceFromSurface(surfs[i])

    # Create the edges
    for i in range(nfaces):
        for ii in range(4):
            i1 = face_to_edge[ii][0]
            j1 = face_to_edge[ii][1]

            # Keep track of the direction
            edge_num = surf_ptrs[i, i1, j1]
            index = abs(edge_num)-1

            if edges[index] is None:
                # Find the vertex numbers
                v1 = edge_ptrs[index, 0]-1
                v2 = edge_ptrs[index, 1]-1

                # The start/end vertex location
                vert1 = None
                vert2 = None

                # Get the indices of the vertices within the surf_ptrs array
                i1 = face_to_edge_verts[ii][0][0]
                j1 = face_to_edge_verts[ii][0][1]
                i2 = face_to_edge_verts[ii][1][0]
                j2 = face_to_edge_verts[ii][1][1]

                if (v1 == surf_ptrs[i, i1, j1]-1 and
                    v2 == surf_ptrs[i, i2, j2]-1):
                    vert1 = verts[v1]
                    vert2 = verts[v2]
                    pts = np.array([face_to_edge_verts[ii][0],
                                    face_to_edge_verts[ii][1]], dtype=np.double)
                elif (v2 == surf_ptrs[i, i1, j1]-1 and
                      v1 == surf_ptrs[i, i2, j2]-1):
                    vert1 = verts[v2]
                    vert2 = verts[v1]
                    pts = np.array([face_to_edge_verts[ii][1],
                                    face_to_edge_verts[ii][0]], dtype=np.double)
                pts = pts/2.0

                # Create the parametric curve
                pcurve = TMR.BsplinePcurve(pts)
                edges[index] = TMR.EdgeFromFace(faces[i], pcurve)
                edges[index].setVertices(vert1, vert2)

    for i in range(nedges):
        if edges[i] is None:
            raise ValueError('Edge %d was not initialized in conversion to TMR\n'%(i))

    # After all of the edges are created, create the edge loops
    # for each face. Account for the CCW orientation.
    ccw_order = [1, 1, -1, -1]
    for i in range(nfaces):
        # Find the edges and directions for this edge loop
        e = []
        dirs = []
        for ii in range(4):
            i1 = face_to_edge[ii][0]
            j1 = face_to_edge[ii][1]
            edge_num = surf_ptrs[i, i1, j1]
            e.append(edges[abs(edge_num)-1])
            dirs.append(np.sign(edge_num)*ccw_order[ii])

        # Set the loop
        loop = TMR.EdgeLoop(e, dirs)
        faces[i].addEdgeLoop(loop)

    # Create the model
    geo = TMR.Model(verts, edges, faces)
    return geo

# Create the truss braced wing model
pgm = Trussbraced()
bse = pgm.initialize()

# Convert from GeoMACH to TMR
geo = geomach_to_tmr(bse)

# Create the mesh
comm = MPI.COMM_WORLD
mesh = TMR.Mesh(comm, geo)

# Mesh the part
opts = TMR.MeshOptions()
opts.num_smoothing_steps = 10
opts.write_mesh_quality_histogram = 1

# Mesh the geometry with the given target size
htarget = 0.25
mesh.mesh(htarget, opts=opts)

# # Create a model from the mesh
# model = mesh.createModelFromMesh()

# # Create the corresponding mesh topology from the mesh-model 
# topo = TMR.Topology(comm, model)

# # Create the quad forest and set the topology of the forest
# forest = TMR.QuadForest(comm)
# forest.setTopology(topo)

# # Create random trees and balance the mesh. Print the output file
# forest.createRandomTrees(nrand=3, max_lev=3)
# forest.balance(1)
# forest.writeForestToVTK('surface-mesh%d.vtk'%(comm.rank))
