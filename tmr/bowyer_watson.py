import numpy as np
import matplotlib.pylab as plt
import sys
import cProfile


class quadnode:
    def __init__(self, low, high, u=0, v=0, level=0):
        # Set the coordinate direction
        self.low = low
        self.high = high

        # Record the lower right-hand parametric point
        self.u = u
        self.v = v

        # Set the level
        self.level = level

        # Set the maximum quad length and the edge length
        hmax = 1 << 30
        self.h = 1 << 30 - (self.level + 1)

        # Set the value of alpha
        ax = (1.0 * (self.u + self.h)) / hmax
        ay = (1.0 * (self.v + self.h)) / hmax

        # Set the new point location
        self.x = (1.0 - ax) * low[0] + ax * high[0]
        self.y = (1.0 - ay) * low[1] + ay * high[1]

        # The left and right nodes
        self.low_left = None
        self.low_right = None
        self.up_left = None
        self.up_right = None

        # Keep track of the points
        self.pts = {}

        return

    def add_point(self, num, pt):
        """Add the point to the left or right node"""

        if self.low_left is not None:
            if pt[0] < self.x and pt[1] < self.y:
                self.low_left.add_point(num, pt)
            elif pt[0] < self.x:
                self.up_left.add_point(num, pt)
            elif pt[1] < self.y:
                self.low_right.add_point(num, pt)
            else:
                self.up_right.add_point(num, pt)
            return

        # Else add the point
        self.pts[num] = np.array(pt)

        # Check if we need to split
        if len(self.pts) > 10:
            # Current (u, v) location for the lower part of the mesh
            u = self.u
            v = self.v

            # Side length for the mesh
            h = self.h
            self.low_left = quadnode(
                self.low, self.high, u=u, v=v, level=self.level + 1
            )
            self.low_right = quadnode(
                self.low, self.high, u=u + h, v=v, level=self.level + 1
            )
            self.up_left = quadnode(
                self.low, self.high, u=u, v=v + h, level=self.level + 1
            )
            self.up_right = quadnode(
                self.low, self.high, u=u + h, v=v + h, level=self.level + 1
            )

            # Add the points to the left or right nodes, respectively
            for key in self.pts.keys():
                item = self.pts.pop(key)
                self.add_point(key, item)

        return

    def delete_point(self, num, pt):
        """Free the point - the number is matched, not the point value"""

        if self.low_left is not None:
            if pt[0] < self.x and pt[1] < self.y:
                return self.low_left.delete_point(num, pt)
            elif pt[0] < self.x:
                return self.up_left.delete_point(num, pt)
            elif pt[1] < self.y:
                return self.low_right.delete_point(num, pt)
            else:
                return self.up_right.delete_point(num, pt)

        # Else add the point
        if num in self.pts:
            item = self.pts.pop(num)
            return True

        return False

    def find_closest(self, pt):
        """Find the closest point"""

        if self.low_left is not None:
            if pt[0] < self.x and pt[1] < self.y:
                num, d = self.low_left.find_closest(pt)
                if self.x - pt[0] <= d:
                    num2, d2 = self.low_right.find_closest(pt)
                    if d2 < d:
                        num, d = num2, d2
                if self.y - pt[1] <= d:
                    num2, d2 = self.up_left.find_closest(pt)
                    if d2 < d:
                        num, d = num2, d2
                if self.x - pt[0] <= d or self.y - pt[1] <= d:
                    num2, d2 = self.up_right.find_closest(pt)
                    if d2 < d:
                        num, d = num2, d2
            elif pt[0] < self.x:
                num, d = self.up_left.find_closest(pt)
                if self.x - pt[0] <= d:
                    num2, d2 = self.up_right.find_closest(pt)
                    if d2 < d:
                        num, d = num2, d2
                if pt[1] - self.y <= d:
                    num2, d2 = self.low_left.find_closest(pt)
                    if d2 < d:
                        num, d = num2, d2
                if self.x - pt[0] <= d or pt[1] - self.y <= d:
                    num2, d2 = self.low_right.find_closest(pt)
                    if d2 < d:
                        num, d = num2, d2
            elif pt[1] < self.y:
                num, d = self.low_right.find_closest(pt)
                if pt[0] - self.x <= d:
                    num2, d2 = self.low_left.find_closest(pt)
                    if d2 < d:
                        num, d = num2, d2
                if self.y - pt[1] <= d:
                    num2, d2 = self.up_right.find_closest(pt)
                    if d2 < d:
                        num, d = num2, d2
                if pt[0] - self.x <= d or self.y - pt[1] <= d:
                    num2, d2 = self.up_left.find_closest(pt)
                    if d2 < d:
                        num, d = num2, d2
            else:
                num, d = self.up_right.find_closest(pt)
                if pt[0] - self.x <= d:
                    num2, d2 = self.up_left.find_closest(pt)
                    if d2 < d:
                        num, d = num2, d2
                if pt[1] - self.y <= d:
                    num2, d2 = self.low_right.find_closest(pt)
                    if d2 < d:
                        num, d = num2, d2
                if pt[0] - self.x <= d or pt[1] - self.y <= d:
                    num2, d2 = self.low_left.find_closest(pt)
                    if d2 < d:
                        num, d = num2, d2
            return num, d

        # This is a leaf, search the points
        dmin = 1e40
        num = -1
        for key in self.pts:
            d = np.sqrt(np.dot(self.pts[key] - pt, self.pts[key] - pt))
            if d < dmin:
                dmin = d
                num = key

        return num, dmin


class triangluate:
    def __init__(self, pts, segs, hole=None, maxd=2):
        self.tri_key = 0
        self.tris = {}
        self.edge_to_tris = {}

        # What triangles have we deleted lately?
        self.deleted_tris = []

        # Keep track of one adjacent triangle to each node
        self.adjacent_tris = {}

        # Set the initial points
        self.pts = [
            np.array((-maxd, -maxd)),
            np.array((maxd, -maxd)),
            np.array((maxd, maxd)),
            np.array((-maxd, maxd)),
        ]

        # Add the initial triangles
        self.add_triangle(0, 1, 3)
        self.add_triangle(3, 1, 2)

        # Create the quadtree root node
        self.root = quadnode([-maxd, -maxd], [maxd, maxd])
        for i in xrange(len(self.pts)):
            self.root.add_point(i, self.pts[i])

        # If there is a hole, add it
        offset = 4
        self.hole_number = -1
        if hole is not None:
            self.hole_number = len(self.pts)
            offset = 5

        # Set the edges that shall not be crossed. These must be initialized
        # before we can call add_vertex()
        self.pslg_edges = []
        for s in segs:
            self.pslg_edges.append((s[0] + offset, s[1] + offset))
            self.pslg_edges.append((s[1] + offset, s[0] + offset))

        if hole is not None:
            self.add_vertex(np.array(hole))

        # Add the vertices in the
        for pt in pts:
            self.add_vertex(np.array(pt))

        # Remove points from the list
        for i in xrange(offset):
            self.root.delete_point(i, self.pts[i])

        # Clean things up, removing the holes/triangles
        for t in self.tris.values():
            if 0 in t or 1 in t or 2 in t or 3 in t:
                self.delete_triangle(t[0], t[1], t[2])

        # Remove all of the holes
        for t in self.tris.values():
            if self.hole_number in t:
                self.remove_hole(t[0], t[1], t[2])

        return

    def remove_hole(self, u, v, w):
        """Remove the specified hole"""

        self.delete_triangle(u, v, w)

        if (u, v) not in self.pslg_edges:
            x = self.adjacent(v, u)
            if x is not None:
                self.remove_hole(v, u, x)
        if (v, w) not in self.pslg_edges:
            x = self.adjacent(w, v)
            if x is not None:
                self.remove_hole(w, v, x)
        if (w, u) not in self.pslg_edges:
            x = self.adjacent(u, w)
            if x is not None:
                self.remove_hole(u, w, x)

        return

    def add_triangle(self, u, v, w):
        """Add the triangle uvw to the triangle list"""
        key = 1 * self.tri_key
        if (u, v) in self.edge_to_tris:
            raise ValueError(
                "Edge 1 (%d, %d) already exists for triangle (%d, %d, %d)"
                % (u, v, t[0], t[1], t[2])
            )
        else:
            self.edge_to_tris[(u, v)] = key

        if (v, w) in self.edge_to_tris:
            raise ValueError(
                "Edge 2 (%d, %d) already exists for triangle (%d, %d, %d)"
                % (v, w, t[0], t[1], t[2])
            )
        else:
            self.edge_to_tris[(v, w)] = key

        if (w, u) in self.edge_to_tris:
            raise ValueError(
                "Edge 3 (%d, %d) already exists for triangle (%d, %d, %d)"
                % (w, u, t[0], t[1], t[2])
            )
        else:
            self.edge_to_tris[(w, u)] = key

        # Set the adjacent triangle to this node
        self.adjacent_tris[u] = key
        self.adjacent_tris[v] = key
        self.adjacent_tris[w] = key

        # Add the triangle itself
        self.tris[key] = (u, v, w)
        self.tri_key += 1
        return

    def delete_triangle(self, u, v, w):
        """Delete the enclosing triangle"""
        if (u, v) in self.edge_to_tris:
            self.edge_to_tris.pop((u, v))
            self.edge_to_tris.pop((v, w))
            key = self.edge_to_tris.pop((w, u))
            self.tris.pop(key)
            self.deleted_tris.append(key)
        return

    def adjacent(self, u, v):
        """Get the triangle that completes the given edge"""
        if (u, v) in self.edge_to_tris:
            t = self.tris[self.edge_to_tris[(u, v)]]
            if t[0] == u and t[1] == v:
                return t[2]
            elif t[1] == u and t[2] == v:
                return t[0]
            elif t[2] == u and t[0] == v:
                return t[1]
        return None

    def get_triangle(self, u, v):
        """Get the unique triangle index associated with the edge"""
        if (u, v) in self.edge_to_tris:
            return self.edge_to_tris[(u, v)]
        return None

    def incircle(self, pt, t):
        """Check if the point is within the circle or not"""
        a = [(self.pts[v][0] - pt[0]) for v in t]
        b = [(self.pts[v][1] - pt[1]) for v in t]
        c = [a[i] ** 2 + b[i] ** 2 for i in range(3)]
        A = np.vstack((a, b, c)).T  # The 3x3 matrix to check
        return np.linalg.det(A) >= 0.0

    def orient2d(self, a, b, pt):
        """Check the relative orientation of the points a, b, and pt"""
        A = np.array(
            [
                [self.pts[a][0] - pt[0], self.pts[a][1] - pt[1]],
                [self.pts[b][0] - pt[0], self.pts[b][1] - pt[1]],
            ]
        )
        # This is NOT robust to precision errors
        return A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0] >= 0.0

    def enclosed(self, u, v, w, pt):
        """Does this triangle enclose this point?"""
        if (
            self.orient2d(u, v, pt)
            and self.orient2d(v, w, pt)
            and self.orient2d(w, u, pt)
        ):
            return True
        return False

    def find_encircled(self, pt):
        """Find the first triangle that encloses the point"""
        for t in self.tris.values():
            if self.incircle(pt, t):
                return t
        return None

    def find_enclosing(self, pt):
        """Find the first triangle that encloses the point"""

        # The closest has node to the given point
        num, dist = self.root.find_closest(pt)

        # Search for nodes that are adjacent to this guy
        tri_num = self.adjacent_tris[num]

        # Access the nodes in the triangle
        t = self.tris[tri_num]

        # Find the orientation of the triangle relative to the
        # node number 'num'
        if t[0] == num:
            u, v, w = t
        elif t[1] == num:
            w, u, v = t
        else:
            v, w, u = t

        winit = w

        # Search the triangles adjacent to this one to find the
        # enclosed point
        while True:
            if self.enclosed(u, v, w, pt):
                return u, v, w

            # Otherwise, find the next triangle that has u in it
            x = self.adjacent(u, w)
            if x is None:
                break

            # Update the v's and w's to march over the triangle
            v = w
            w = x

            if w == winit:
                break

        # Walk the mesh from here...
        for t in self.tris.values():
            if self.enclosed(t[0], t[1], t[2], pt):
                return t

        return None

    def dig_cavity(self, u, v, w):
        """
        u is the index of a new point to be inserted. Is the
        oriented triangle (u,v,w) Delaunay or not?
        """
        if (w, v) in self.pslg_edges:
            self.add_triangle(u, v, w)
            return

        # Get the adjacent node
        x = self.adjacent(w, v)

        if x is not None:
            if self.incircle(self.pts[u], (w, v, x)):
                self.delete_triangle(w, v, x)
                self.dig_cavity(u, v, x)
                self.dig_cavity(u, x, w)
                return

        self.add_triangle(u, v, w)

        return

    def add_vertex(self, pt):
        """
        Add the vertex to the underlying Delaunay triangularization
        """

        # Find the enclosing triangle
        t = self.find_enclosing(pt)

        # Add the point to the quadtree and point list
        u = len(self.pts)  # Pt index
        self.root.add_point(u, pt)
        self.pts.append(pt)

        if t is not None:
            v, w, x = t
            self.delete_triangle(v, w, x)
            self.dig_cavity(u, v, w)
            self.dig_cavity(u, w, x)
            self.dig_cavity(u, x, v)

        return

    def add_vertex_frontal(self, pt, t):
        """Add the vertex contained within the specified triangle"""
        u = len(self.pts)  # Pt index
        self.root.add_point(u, pt)
        self.pts.append(pt)

        if t is not None:
            v, w, x = t
            self.delete_triangle(v, w, x)
            self.dig_cavity(u, v, w)
            self.dig_cavity(u, w, x)
            self.dig_cavity(u, x, v)

        return

    def circumcircle(self, t):
        """Compute the circumcircle radii for this set of triangles"""
        p = [self.pts[v] for v in t]
        r1 = np.sqrt((p[1][0] - p[0][0]) ** 2 + (p[1][1] - p[0][1]) ** 2)
        r2 = np.sqrt((p[2][0] - p[1][0]) ** 2 + (p[2][1] - p[1][1]) ** 2)
        r3 = np.sqrt((p[2][0] - p[0][0]) ** 2 + (p[2][1] - p[0][1]) ** 2)
        return max(r1, r2, r3)

    def computeIntersection(self, m, e, u, v, w):
        """
        Compute the distance to the intersection of the line m + alpha*e with the
        two lines from u to w and v to w
        """

        # m + alpha*e = pts[u] + beta*(pts[w] - pts[u])
        # => alpha*e + beta*(pts[u] - pts[w]) = pts[u] - m
        a2 = [self.pts[u][i] - self.pts[w][i] for i in range(2)]

        # Check if the line is orthogonal to this direction
        if e[0] * a2[1] - e[1] * a2[0] != 0.0:
            b = [self.pts[u][i] - m[i] for i in range(2)]
            A = np.vstack((e, a2)).T
            ans = np.linalg.solve(A, b)

            # If beta is on the interval [0,1], alpha is the distance
            if ans[1] >= 0.0 and ans[1] <= 1.0:
                return ans[0]

        # If we get here and there's no solution, then we have a degenerate triangle
        b = [self.pts[v][i] - m[i] for i in range(2)]
        a2 = [self.pts[v][i] - self.pts[w][i] for i in range(2)]
        A = np.vstack((e, a2)).T
        ans = np.linalg.solve(A, b)

        # If beta is on the interval [0,1], alpha is the distance
        if ans[1] >= 0.0 and ans[1] <= 1.0:
            return ans[0]

        return -1.0

    def frontal(self, h, plot=False, freq=50):
        """Refine the triangles using the frontal method"""

        # The length of the edge of the triangle
        de = 0.5 * np.sqrt(3.0) * h

        # Add all of the triangles to the active set
        status = {}
        circumcircles = {}
        for key in self.tris:
            status[key] = "waiting"

            # Check whether the triangle should be labeled as active
            t = self.tris[key]
            u, v, w = t
            for e in [(u, v), (v, w), (w, u)]:
                if e in self.pslg_edges:
                    status[key] = "active"

        # Compute the circumcircles
        for key in self.tris:
            circumcircles[key] = self.circumcircle(self.tris[key])

        # Compute the circumcircle radii of all triangles
        itr = 0
        while "active" in status.values():
            if itr % freq == 0:
                print("iteration = %d" % (itr))
                if plot:
                    self.status = status
                    self.plot()
                    plt.show()
            itr += 1

            # Find the wors ttriangle
            rmax = 0.0
            tri = -1
            # If we maintained a sorted list of the worst offenders, this
            # could be faster
            for key in circumcircles:
                if status[key] == "active" and circumcircles[key] > rmax:
                    rmax = circumcircles[key]
                    tri = key

            # Now insert a new point based on the Voronoi criteria
            # Determine the triangle's accepted or boundary edge
            t = self.tris[tri]

            # Check if we are along an edge of the PSLG
            edge = None
            u, v, w = t
            for e in [(u, v), (v, w), (w, u)]:
                if e in self.pslg_edges:
                    edge = e

            # Check whether one of the adjacent edges is done
            if edge is None:
                for e in [(u, v), (v, w), (w, u)]:
                    index = self.get_triangle(e[1], e[0])
                    if index is not None and status[index] == "done":
                        edge = e
                        break

            # Compute the location of the new point
            # | i   j   k |
            # | dx  dy  0 | = i*dy - j*dx
            # | 0   0   1 |
            u, v = edge
            w = self.adjacent(u, v)
            m = 0.5 * np.array(
                [self.pts[u][0] + self.pts[v][0], self.pts[u][1] + self.pts[v][1]]
            )
            e = -np.array(
                [self.pts[v][1] - self.pts[u][1], self.pts[u][0] - self.pts[v][0]]
            )

            # Compute half the distance between the u and v points
            p = 0.5 * np.sqrt(np.sum(e**2))
            e = 0.5 * e / p

            # Try the size-optimal point first
            pt = m + de * e

            # Find the enclosing triangle for the
            if not self.enclosed(t[0], t[1], t[2], pt):
                t = self.find_enclosing(pt)

                if t is None:
                    # Pick a point that is actually in the current triangle
                    q = self.computeIntersection(m, e, u, v, w)
                    rho = 0.5 * q
                    pt = m + rho * e
                    t = self.tris[tri]

            # Use the edge for this triangle to determine where to insert
            # the new point
            old_tri = 1 * self.tri_key
            self.deleted_tris = []
            self.add_vertex_frontal(pt, t)
            new_tri = self.tri_key

            # Free the deleted triangles
            for key in self.deleted_tris:
                circumcircles.pop(key)
                status.pop(key)

            # Compute the circumcircles of the new triangles and check
            # whether they belong in the done category or not...
            for key in range(old_tri, new_tri):
                circumcircles[key] = self.circumcircle(self.tris[key])
                if circumcircles[key] <= 1.5 * h:
                    status[key] = "done"
                else:
                    status[key] = "waiting"

            # Mark the triangle with the edge that was just computed
            # as being done, regardless of whether it is or not...
            # Otherwise the same extrapolation point will be used which
            # will duplicate points within the domain leading to chaos!
            if (u, v) in self.edge_to_tris:
                tri = self.edge_to_tris[(u, v)]
                status[tri] = "done"

            # Label the new triangles with active status
            for key in range(old_tri, new_tri):
                if status[key] == "done":
                    continue

                # Extract the triangle
                t = self.tris[key]
                u, v, w = t

                # Check the adjacent triangles
                for edge in [(u, v), (v, w), (w, u)]:
                    index = self.get_triangle(edge[1], edge[0])
                    if index is not None and status[index] == "done":
                        status[key] = "active"
                        break
                    if edge in self.pslg_edges:
                        status[key] = "active"
                        break

        self.status = status

        return

    def plot(self):
        """Plot the Delaunay triangularization"""

        plt.figure()

        for edge in self.pslg_edges:
            x = [self.pts[v][0] for v in edge]
            y = [self.pts[v][1] for v in edge]
            plt.plot(x, y, "b", linewidth=2)

        for key in self.tris:
            t = self.tris[key]
            t2 = [v for v in t]
            t2.append(t[0])
            x = [self.pts[v][0] for v in t2]
            y = [self.pts[v][1] for v in t2]

            if hasattr(self, "status"):
                if key in self.status.keys():
                    if self.status[key] == "active":
                        plt.fill(x, y, "b", alpha=0.2, edgecolor="r")
                    elif self.status[key] == "done":
                        plt.fill(x, y, "r", alpha=0.2, edgecolor="r")
                    elif self.status[key] == "waiting":
                        plt.fill(x, y, "g", alpha=0.2, edgecolor="r")
                else:
                    plt.fill(x, y, "y", alpha=0.2, edgecolor="k")
            else:
                plt.fill(x, y, "g", alpha=0.2, edgecolor="k")

        for pt in self.pts[4:]:
            plt.plot(pt[0], pt[1], "ro")

        return


# Create a list of points that forms
h = 0.75
r1 = 1.5
r2 = 0.75
r3 = 0.2
hole = None

# The list of points and segments
pts = []
segs = []

if "circle" in sys.argv:
    # Create the points for the outer circle
    npts = (int)((2 * r1 * np.pi) / h)
    h = 2 * np.pi * r1 / npts
    u = np.linspace(0, 2 * np.pi, npts + 1)[:-1]
    for i in xrange(npts):
        pts.append([r1 * np.cos(u[i]), r1 * np.sin(u[i])])
        if i == npts - 1:
            segs.append([i, 0])
        else:
            segs.append([i, i + 1])

elif "annulus" in sys.argv:
    hole = [0, 0]
    # Create the points for the outer circle
    npts = (int)((2 * r1 * np.pi) / h)
    h = 2 * np.pi * r1 / npts
    u = np.linspace(0, 2 * np.pi, npts + 1)[:-1]
    for i in xrange(npts):
        pts.append([r1 * np.cos(u[i]), r1 * np.sin(u[i])])
        if i == npts - 1:
            segs.append([i, 0])
        else:
            segs.append([i, i + 1])

    offset = npts

    # Create the points for the inner circle
    npts = (int)((2 * r2 * np.pi) / h)
    u = np.linspace(0, 2 * np.pi, npts + 1)[:-1]
    for i in xrange(npts):
        pts.append([r2 * np.cos(u[-1 - i]), r2 * np.sin(u[-1 - i])])
        if i == npts - 1:
            segs.append([offset + i, offset])
        else:
            segs.append([offset + i, offset + i + 1])

elif "triangle" in sys.argv:
    npts = (int)(2 * r1 / h)
    u = np.linspace(-r1, r1, npts + 1)

    for i in xrange(npts):
        pts.append([u[i], -r1])
    for i in xrange(npts):
        pts.append([r1, u[i]])
    for i in xrange(npts):
        v = u[-1 - i]
        x = v + 0.1 * (v - r1) * (v + r1)
        y = v - 0.1 * (v - r1) * (v + r1)
        pts.append([x, y])

    for i in xrange(3 * npts):
        segs.append([i, i + 1])
    segs[-1][1] = 0

else:
    npts = (int)(2 * r1 / h)
    u = np.linspace(-r1, r1, npts + 1)
    for i in xrange(npts):
        pts.append([u[i], -r1])
    for i in xrange(npts):
        pts.append([r1, u[i]])
    for i in xrange(npts):
        pts.append([u[-1 - i], r1])
    for i in xrange(npts):
        pts.append([-r1, u[-1 - i]])

    for i in xrange(4 * npts):
        segs.append([i, i + 1])
    segs[-1][1] = 0

    hole = [0, 0]

    offset = 4 * npts

    # Create the points for the inner circle
    npts = (int)((2 * r2 * np.pi) / h)
    u = np.linspace(0, 2 * np.pi, npts + 1)[:-1]
    for i in xrange(npts):
        pts.append([r2 * np.cos(u[-1 - i]), r2 * np.sin(u[-1 - i])])
        if i == npts - 1:
            segs.append([offset + i, offset])
        else:
            segs.append([offset + i, offset + i + 1])

# Create the triangularization
tri = triangluate(pts, segs, hole=hole)
# cProfile.run('tri.frontal(h)')
tri.frontal(h)

tri.plot()
plt.savefig("bowyer_watson.pdf")
plt.show()
