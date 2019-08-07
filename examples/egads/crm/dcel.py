'''
Implementation of a DECL data structure
'''
import numpy as np

def orient2d(self, ax, ay, bx, by, cx, cy):
    '''Check the relative orientation of the points a, b, and pt'''
    A = np.array([
        [ax - cx, ay - cy],
        [bx - cx, by - cy]])
    # This is NOT robust to precision errors
    return  A[0,0]*A[1,1] - A[0,1]*A[1,0] >= 0.0

class face:
    def __init__(self):
        self.edge = None

    def contains(self, X):
        '''
        Check if this face contains the edge
        '''
        e = self.edge

        while True:
            u = next.u
            v = next.v
            check = oriend2d(X[u][0], X[u][1],
                             X[v][0], X[v][1],
                             px, py)

            # Check if this is false, if so, the
            # point is not in the face
            if not check:
                return False

            # Go to the next edge
            e = e.next

            # If the next edge returns itself, quit
            if e == self.edge:
                break

        return True

class edge:
    def __init__(self, u, v, name):
        self.u = u
        self.v = v
        self.face = None

        # Pointer to the twin (opposite) edge, next
        # and previous edge
        self.twin = None
        self.next = None
        self.prev = None

        # Name the edge
        self.name = name

class dcel:
    def __init__(self, X, conn, names=None):
        '''
        Initialize the DCEL from points and a list of polygons
        '''

        # Create a list of the vertices
        self.verts = []
        for x in X:
            self.verts.append([x[0], x[1]])

        # Create an empty list of the faces
        self.faces = []

        # Dictionary hashed on (u, v)
        self.edges = {}

        if names is None:
            names = [None]*len(conn)

        # Loop over all of the input polygons
        for p, name in zip(conn, names):
            # Look for the edges
            p.append(p[0])

            # Create the face
            f = face()

            # Set the pairs for the first/previous edges
            first = None
            prev = None

            for i in range(len(p)-1):
                u = p[i]
                v = p[i+1]

                # Set the key for the owned edge
                if (u, v) not in self.edges:
                    self.edges[(u, v)] = edge(u, v, name)
                    self.edges[(u, v)].face = f

                # Set the twin if it exists
                if (v, u) in self.edges:
                    self.edges[(v, u)].twin = self.edges[(u, v)]
                    self.edges[(u, v)].twin = self.edges[(v, u)]

                # Set the previous/next pointer
                if prev is not None:
                    self.edges[(u, v)].prev = self.edges[prev]
                    self.edges[prev].next = self.edges[(u, v)]

                # Set the first pair and the previous pairs
                if first is None:
                    first = (u, v)
                prev = (u, v)

            if first is not None:
                self.edges[first].prev = self.edges[prev]
                self.edges[prev].next = self.edges[first]

            # Set the one edge in the face
            f.edge = self.edges[first]
            self.faces.append(f)

        return

    def get_connectivity(self):
        '''
        Retrieve a list of the edges
        '''
        edge_dict = {}

        # Order the edges
        for e in self.edges:
            op = (e[1], e[0])
            if e not in edge_dict and op not in edge_dict:
                index = len(edge_dict)
                edge_dict[e] = index

        # Set the edges
        edges = [None]*len(edge_dict)
        for key in edge_dict:
            edges[edge_dict[key]] = key

        # Create the faces
        faces = []
        sense = []
        for f in self.faces:
            fl = []
            sl = []

            e = f.edge
            ef = (e.u, e.v)
            er = (e.v, e.u)
            if ef in edge_dict:
                fl.append(edge_dict[ef])
                sl.append(1)
            else:
                fl.append(edge_dict[er])
                sl.append(-1)
            e = e.next
            while e != f.edge:
                ef = (e.u, e.v)
                er = (e.v, e.u)
                if ef in edge_dict:
                    fl.append(edge_dict[ef])
                    sl.append(1)
                else:
                    fl.append(edge_dict[er])
                    sl.append(-1)
                e = e.next

            faces.append(fl)
            sense.append(sl)

        return self.verts, edges, faces, sense

    def find_enclosing(self, px, py):
        '''
        Find the face that encloses the point (if any)
        '''
        for f in self.faces:
            if f.contains(self.X, px, py):
                return f

        return None

    def find_closest_edge(self, px, py):
        '''
        Find the closest edge
        '''
        closest = None
        dist = 0.0

        for key in self.edges:
            u = key[0]
            v = key[1]

            # Compute the component normal to the edge
            dx = self.verts[v][0] - self.verts[u][0]
            dy = self.verts[v][1] - self.verts[u][1]
            nx = dy
            ny = -dx
            inv = 1.0/np.sqrt(nx**2 + ny**2)

            ax = px - self.verts[u][0]
            ay = py - self.verts[u][1]

            proj = (ax*dx + ay*dy)*inv*inv
            if proj > 1.0:
                dist = np.sqrt((px - self.verts[v][0])**2 +
                               (py - self.verts[v][1])**2)
            elif proj < 0.0:
                dist = np.sqrt((px - self.verts[u][0])**2 +
                               (py - self.verts[u][1])**2)
            else:
                dist = np.fabs((ax*nx + ay*ny)*inv)

            if closest is None or dist < min_dist:
                closest = (u, v)
                min_dist = 1.0*dist

        return self.edges[closest], min_dist

    def add_vertex(self, e, px, py):
        '''
        Add a vertex along the given edge at the given point
        '''

        # Set the u, v vertex values
        u = e[0]
        v = e[1]

        # Create new half edge objects for
        w = len(self.verts)
        self.verts.append([px, py])

        a1 = None
        a2 = None
        if (u, v) in self.edges:
            # Add the edges (u, w) and (w, v)
            orig = self.edges[(u, v)]
            a1 = edge(u, w, orig.name)
            a2 = edge(w, v, orig.name)
            a1.face = orig.face
            a1.prev = orig.prev
            a1.next = a2
            orig.prev.next = a1

            a2.face = orig.face
            a2.prev = a1
            a2.next = orig.next
            orig.next.prev = a2
            a1.face.edge = a1

            # Delete the entry from the hash table
            del self.edges[(u, v)]

            # Add the edges to the dictionary
            self.edges[(u, w)] = a1
            self.edges[(w, v)] = a2

        b1 = None
        b2 = None
        if (v, u) in self.edges:
            # Add the edges (u, w) and (w, v)
            orig = self.edges[(v, u)]
            b1 = edge(v, w, orig.name)
            b2 = edge(w, u, orig.name)
            b1.face = orig.face
            b1.prev = orig.prev
            b1.next = b2
            orig.prev.next = b1

            b2.face = orig.face
            b2.prev = b1
            b2.next = orig.next
            orig.next.prev = b2
            b2.face.edge = b2

            # Delete the entry from the hash table
            del self.edges[(v, u)]

            # Add the edges to the dictionary
            self.edges[(v, w)] = b1
            self.edges[(w, u)] = b2

        # Set the twins if both edges were created
        if a1 and b1:
            a1.twin = b2
            b2.twin = a1

            a2.twin = b1
            b1.twin = a2

        return w

    def add_edge_from_face(self, f, u, v, name=None):
        '''
        Split a face shared by the
        '''
        e1 = None
        e2 = None

        # Loop around the existing face
        e = f.edge
        reverse = False
        if e.v == u or e.v == v:
            e1 = e

        e = e.next
        while e != f.edge:
            if e.v == u or e.v == v:
                if e1 is None:
                    e1 = e
                elif e2 is None:
                    e2 = e
                    break
            e = e.next

        if e1 is None or e2 is None:
            # errmsg = 'Edge (%d, %d) does not split the given face'%(u, v)
            # raise RuntimeError(errmsg)
            return

        f2 = face()
        a = edge(e1.v, e2.v, name)
        b = edge(e2.v, e1.v, name)

        # Set the edge
        f.edge = a

        # Set the pointers in the first edge
        a.face = f
        a.prev = e1
        a.next = e2.next
        a.twin = b

        # Set the pointers in the twi
        b.face = f2
        b.prev = e2
        b.next = e1.next
        b.twin = a

        # Set the previous and next pointers in the edges
        e2.next.prev = a
        e2.next = b
        e1.next.prev = b
        e1.next = a

        # Relabel all the edges in the second face
        f2.edge = b
        e = f2.edge
        e.face = f2
        e = e.next
        while e != f2.edge:
            e.face = f2
            e = e.next

        # Add the new face to the list
        self.faces.append(f2)

        # Add the edges to the edge list
        self.edges[(e1.v, e2.v)] = a
        self.edges[(e2.v, e1.v)] = b

        return

    def get_intersection(self, u, v, xr, yr, dxr, dyr):
        '''
        Find the intersection between the edge and the line

        Solve the equation
        (x + u*dx, y + u*dy) = (xrb + v*dxrb, yrb + v*dyrb)

        For u and v and return the v component by solving the equations:

        u*dx - v*dxrb = xrb - x
        u*dy - v*dyrb = yrb - y
        '''

        x = self.verts[u][0]
        y = self.verts[u][1]
        dx = self.verts[v][0] - self.verts[u][0]
        dy = self.verts[v][1] - self.verts[u][1]

        # Check if the two lines are parallel
        tol = 1e-9
        if np.fabs(dx*dyr - dy*dxr) < tol*np.fabs(dx*dxr + dy*dyr):
            return None

        # Solve for uv = [u, -v]
        uv = np.linalg.solve(np.array([[dx, dxr],
                                       [dy, dyr]], dtype=np.float),
                             np.array([xr - x, yr - y], dtype=np.float))
        return uv[0]

    def split_face(self, f, x, y, dx, dy):
        '''Split the face'''

        # Set the edges and the edge parameters
        e1 = None
        e2 = None
        p1 = 0.0
        p2 = 0.0

        e = f.edge
        p = self.get_intersection(e.u, e.v, x, y, dx, dy)
        if p is not None and p >= 0.0 and p < 1.0:
            p1 = p
            e1 = e
        e = e.next
        while e != f.edge:
            p = self.get_intersection(e.u, e.v, x, y, dx, dy)
            if p is not None and p >= 0.0 and p < 1.0:
                if e1 is None:
                    p1 = p
                    e1 = e
                elif e2 is None:
                    p2 = p
                    e2 = e
            if e1 is not None and e2 is not None:
                break
            e = e.next

        if e1 is None or e2 is None:
            return

        # Add the first vertex
        px = (1.0 - p1)*self.verts[e1.u][0] + p1*self.verts[e1.v][0]
        py = (1.0 - p1)*self.verts[e1.u][1] + p1*self.verts[e1.v][1]
        w = self.add_vertex((e1.u, e1.v), px, py)

        # Add the first vertex
        px = (1.0 - p2)*self.verts[e2.u][0] + p2*self.verts[e2.v][0]
        py = (1.0 - p2)*self.verts[e2.u][1] + p2*self.verts[e2.v][1]
        x = self.add_vertex((e2.u, e2.v), px, py)

        # Add the vertex
        self.add_edge_from_face(f, w, x)

        return

    def plot(self):
        '''Plot the DCEL'''
        import matplotlib.pylab as plt

        plt.figure()
        for f in self.faces:
            e = f.edge
            x = [self.verts[e.u][0], self.verts[e.v][0]]
            y = [self.verts[e.u][1], self.verts[e.v][1]]
            e = e.next
            while e != f.edge:
                x.append(self.verts[e.v][0])
                y.append(self.verts[e.v][1])
                e = e.next

            plt.plot(x, y, '-o')

        plt.axis('equal')
        return

if __name__ == '__main__':
    import matplotlib.pylab as plt

    X = [[0, 0], [1, 0], [2, 0],
         [0.5, 1], [1.5, 1], [2.5, 1],
         [0.5, 2], [2, 2]]

    conn = [
        [0, 1, 2, 4, 3],
        [4, 2, 5],
        [3, 4, 5, 7, 6]]

    d = dcel(X, conn)
    d.add_vertex((0, 1), 0.5, 0)
    d.add_vertex((3, 4), 1, 1)

    d.add_edge_from_face(d.faces[0], 1, 4)
    d.add_edge_from_face(d.faces[0], 1, 3)

    d.split_face(d.faces[2], 0.25, 0.0, 0.9, 1.0)

    d.plot()
    plt.show()
