# distutils: language=c++

# This file is part of the package TMR for adaptive mesh refinement.

# Copyright (C) 2015 Georgia Tech Research Corporation.
# Additional copyright (C) 2015 Graeme Kennedy.
# All rights reserved.

# TMR is licensed under the Apache License, Version 2.0 (the "License");
# you may not use this software except in compliance with the License.
# You may obtain a copy of the License at

#   http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function, division

# For the use of MPI
from mpi4py.libmpi cimport *
cimport mpi4py.MPI as MPI

# Import numpy
cimport numpy as np
import numpy as np

# Import random
import random

cdef tmr_init():
    if not TMRIsInitialized():
        TMRInitialize()

# Ensure that numpy is initialized
np.import_array()

# Initialize the MPI libraries in TMR (if not already done)
tmr_init()

# Import the definition required for const strings
from libc.string cimport const_char
from libc.stdlib cimport malloc, free

# Import C methods for python
from cpython cimport PyObject, Py_INCREF

# Import TACS and ParOpt
from tacs.TACS cimport *
from tacs.constitutive cimport *
from tacs.functions cimport *
from paropt.ParOpt cimport *
from tacs import TACS

# Import the definitions
from TMR cimport *

# Include the mpi4py header
cdef extern from "mpi-compat.h":
    pass

# Max level
MAX_LEVEL = TMR_MAX_LEVEL

# Set the different types of meshes
NO_MESH = TMR_NO_MESH
STRUCTURED = TMR_STRUCTURED
UNSTRUCTURED = TMR_UNSTRUCTURED
TRIANGLE = TMR_TRIANGLE

# Set the type of interpolation to use
UNIFORM_POINTS = TMR_UNIFORM_POINTS
GAUSS_LOBATTO_POINTS = TMR_GAUSS_LOBATTO_POINTS
BERNSTEIN_POINTS = TMR_BERNSTEIN_POINTS

cdef class Vertex:
    """
    The vertex class is used to store both the point and to
    represent the underlying geometry
    """
    cdef TMRVertex *ptr
    def __cinit__(self):
        self.ptr = NULL

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def evalPoint(self):
        """
        evalPoint(self)

        Evaluate the point on the parametric surface and returns the node
        location
        """
        cdef TMRPoint pt
        self.ptr.evalPoint(&pt)
        return np.array([pt.x, pt.y, pt.z])

    def setName(self, aname):
        """
        setName(self, aname)

        Set the name associated with the vertex

        Parameters
        ----------
        aname: str
            name associated with the vertex
        """
        cdef char *name = tmr_convert_str_to_chars(aname)
        if self.ptr:
            self.ptr.setName(name)

    def getName(self):
        """
        getName(self)

        Get the name associated with the vertex

        Returns
        ----------
        name: str
            name associated with the entity
        """
        cdef const char *name = NULL
        if self.ptr:
            name = self.ptr.getName()
            return tmr_convert_char_to_str(name)
        return None

    def getEntityId(self):
        """
        getEntityId(self)

        Get the entity id associated with the vertex

        Returns
        ----------
        id: int
            id number of the model
        """
        if self.ptr:
            return self.ptr.getEntityId()
        return -1

    def setNodeNum(self, num):
        """
        setNodeNum(self, num)
        Set the node number associated with the vertex

        Parameters
        ----------
        num: int
            node number associated with the vertex
        """
        cdef int n = num
        self.ptr.setNodeNum(&n)
        return n

cdef _init_Vertex(TMRVertex *ptr):
    vertex = Vertex()
    vertex.ptr = ptr
    if ptr != NULL:
        vertex.ptr.incref()
    return vertex

cdef class Edge:
    """
    The edge class is used to store both the points on an edge and to
    represent the underlying geometry
    """
    cdef TMREdge *ptr
    def __cinit__(self):
        self.ptr = NULL

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def getRange(self):
        cdef double tmin, tmax
        self.ptr.getRange(&tmin, &tmax)
        return tmin, tmax

    def evalPoint(self, double t):
        """
        evalPoint(self, t)

        Evaluate a node location on a curve with a single
        parametric argument *t*. Provides access to the first and second vertices
        that begin/end the edge. Given a parametric location t and a face, return
        the (*u,v*) location on the face.

        Parameters
        ----------
        t: double
          parameteric location
        """
        cdef TMRPoint pt
        self.ptr.evalPoint(t, &pt)
        return np.array([pt.x, pt.y, pt.z])

    def setName(self, aname):
        cdef char *name = tmr_convert_str_to_chars(aname)
        if self.ptr:
            self.ptr.setName(name)

    def getName(self):
        cdef const char *name = NULL
        if self.ptr:
            name = self.ptr.getName()
            return tmr_convert_char_to_str(name)
        return None

    def getEntityId(self):
        if self.ptr:
            return self.ptr.getEntityId()
        return -1

    def setVertices(self, Vertex v1, Vertex v2):
        self.ptr.setVertices(v1.ptr, v2.ptr)

    def getVertices(self):
        cdef TMRVertex *v1 = NULL
        cdef TMRVertex *v2 = NULL
        self.ptr.getVertices(&v1, &v2)
        return _init_Vertex(v1), _init_Vertex(v2)

    def isDegenerate(self):
        return self.ptr.isDegenerate()

    def setSource(self, Edge e):
        self.ptr.setSource(e.ptr)

    def getSource(self):
        cdef TMREdge *e
        self.ptr.getSource(&e)
        return _init_Edge(e)

    def setMesh(self, EdgeMesh mesh):
        self.ptr.setMesh(mesh.ptr)

    def writeToVTK(self, fname):
        cdef char *filename = tmr_convert_str_to_chars(fname)
        self.ptr.writeToVTK(filename)

cdef _init_Edge(TMREdge *ptr):
    edge = Edge()
    edge.ptr = ptr
    if ptr != NULL:
        edge.ptr.incref()
    return edge

cdef class EdgeLoop:
    """
    Contains an oriented set of edges that enclose a surface. Provides the edge
    list and their relative orientations in the loop
    """
    cdef TMREdgeLoop *ptr
    def __cinit__(self, list edges=None, list dirs=None):
        cdef int nedges = 0
        cdef TMREdge **e = NULL
        cdef int *d = NULL
        self.ptr = NULL

        # If no directions are specified, try to find them
        # based on the node locations
        if edges is not None and dirs is None:
            edge_list = edges[:]
            nedges = len(edge_list)
            dirs = [1]
            elist = [edge_list.pop()]
            v1, vnext = elist[0].getVertices()
            for k in range(nedges-1):
                for i, edge in enumerate(edge_list):
                    v1, v2 = edge.getVertices()
                    if v1.getEntityId() == vnext.getEntityId():
                        dirs.append(1)
                        elist.append(edge)
                        edge_list.remove(edge)
                        vnext = v2
                        break
                    elif v2.getEntityId() == vnext.getEntityId():
                        dirs.append(-1)
                        elist.append(edge)
                        edge_list.remove(edge)
                        vnext = v1
                        break
            if len(elist) == nedges:
                edges = elist
            else:
                errmsg = 'EdgeLoop: Could not compute directions from edges'
                raise ValueError(errmsg)

        if (edges is not None and dirs is not None and
            len(edges) == len(dirs)):
            nedges = len(edges)
            e = <TMREdge**>malloc(nedges*sizeof(TMREdge*))
            d = <int*>malloc(nedges*sizeof(int))
            for i in range(nedges):
                e[i] = (<Edge>edges[i]).ptr
                d[i] = <int>dirs[i]

            self.ptr = new TMREdgeLoop(nedges, e, d)
            self.ptr.incref()
            free(e)
            free(d)

    def __decalloc__(self):
        if self.ptr:
            self.ptr.decref()

    def getEntityId(self):
        if self.ptr:
            return self.ptr.getEntityId()
        return -1

    def getEdgeLoop(self):
        cdef int nedges = 0
        cdef TMREdge **edges = NULL
        cdef const int *dirs = NULL
        self.ptr.getEdgeLoop(&nedges, &edges, &dirs)
        e = []
        d = []
        for i in range(nedges):
            e.append(_init_Edge(edges[i]))
            d.append(dirs[i])
        return e, d

cdef _init_EdgeLoop(TMREdgeLoop *ptr):
    loop = EdgeLoop()
    loop.ptr = ptr
    if ptr != NULL:
        loop.ptr.incref()
    return loop

cdef class Face:
    """
    The face class is used to store both the points and edges and to
    represent the underlying geometry
    """
    cdef TMRFace *ptr
    def __cinit__(self):
        self.ptr = NULL

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def getRange(self):
        cdef double umin, vmin, umax, vmax
        self.ptr.getRange(&umin, &vmin, &umax, &vmax)
        return umin, vmin, umax, vmax

    def evalPoint(self, double u, double v):
        """
        evalPoint(self, u, v)

        Evaluate a node location given a parametric location
        (*u,v*). Access the edge loops that bound the face.

        Parameters
        ----------
        u: double
         Parametric location in u
        v: double
         Parametric location in v
        """
        cdef TMRPoint pt
        self.ptr.evalPoint(u, v, &pt)
        return np.array([pt.x, pt.y, pt.z])

    def setName(self, aname):
        cdef char *name = tmr_convert_str_to_chars(aname)
        if self.ptr:
            self.ptr.setName(name)

    def getName(self):
        cdef const char *name = NULL
        if self.ptr:
            name = self.ptr.getName()
            return tmr_convert_char_to_str(name)
        return None

    def getEntityId(self):
        if self.ptr:
            return self.ptr.getEntityId()
        return -1

    def getNumEdgeLoops(self):
        return self.ptr.getNumEdgeLoops()

    def addEdgeLoop(self, int _dir, EdgeLoop loop):
        self.ptr.addEdgeLoop(_dir, loop.ptr)

    def getEdgeLoop(self, k):
        cdef TMREdgeLoop *loop = NULL
        self.ptr.getEdgeLoop(k, &loop)
        if loop:
            return _init_EdgeLoop(loop)
        return None

    def setSource(self, Volume v, Face f):
        self.ptr.setSource(v.ptr, f.ptr)

    def getSource(self):
        cdef TMRFace *f
        cdef TMRVolume *v
        cdef int d
        self.ptr.getSource(&d, &v, &f)
        return d, _init_Volume(v), _init_Face(f)

    def setMesh(self, FaceMesh mesh):
        self.ptr.setMesh(mesh.ptr)

    def getMesh(self):
        cdef TMRFaceMesh *fmesh
        self.ptr.getMesh(&fmesh)
        return _init_FaceMesh(fmesh)

    def writeToVTK(self, fname):
        cdef char *filename = tmr_convert_str_to_chars(fname)
        self.ptr.writeToVTK(filename)

cdef _init_Face(TMRFace *ptr):
    face = Face()
    face.ptr = ptr
    if ptr != NULL:
        face.ptr.incref()
    return face

cdef class Volume:
    """
    Contains an oriented collection of surfaces that enclose a volume. Provides
    a list of surfaces and their relative orientations in the volume.
    """
    cdef TMRVolume *ptr
    def __cinit__(self, list faces=None, list dirs=None):
        cdef int nfaces = 0
        cdef int *d = NULL
        cdef TMRFace **f = NULL
        self.ptr = NULL
        if faces is not None and dirs is not None:
            if len(faces) == len(dirs):
                nfaces = len(faces)
                d = <int*>malloc(nfaces*sizeof(int))
                f = <TMRFace**>malloc(nfaces*sizeof(TMRFace*))
                for i in range(nfaces):
                    f[i] = (<Face>faces[i]).ptr
                    d[i] = <int>dirs[i]
                self.ptr = new TMRVolume(nfaces, f, d)
                self.ptr.incref()
                free(d)
                free(f)
            else:
                errmsg = 'Face and direction lists must be the same length'
                raise ValueError(errmsg)

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def setName(self, aname):
        cdef char *name = tmr_convert_str_to_chars(aname)
        if self.ptr:
            self.ptr.setName(name)

    def getName(self):
        cdef const char *name = NULL
        if self.ptr:
            name = self.ptr.getName()
            return tmr_convert_char_to_str(name)
        return None

    def getEntityId(self):
        if self.ptr:
            return self.ptr.getEntityId()
        return -1

    def getFaces(self):
        cdef TMRFace **f
        cdef int num_faces = 0
        if self.ptr:
            self.ptr.getFaces(&num_faces, &f, NULL)
        faces = []
        for i in range(num_faces):
            faces.append(_init_Face(f[i]))
        return faces

    def writeToVTK(self, fname):
        cdef char *filename = tmr_convert_str_to_chars(fname)
        self.ptr.writeToVTK(filename)

cdef _init_Volume(TMRVolume *ptr):
    vol = Volume()
    vol.ptr = ptr
    if ptr != NULL:
        vol.ptr.incref()
    return vol

cdef class Curve:
    cdef TMRCurve *ptr
    def __cinit__(self):
        self.ptr = NULL

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def setName(self, aname):
        cdef char *name = tmr_convert_str_to_chars(aname)
        if self.ptr:
            self.ptr.setName(name)

    def getEntityId(self):
        if self.ptr:
            return self.ptr.getEntityId()
        return -1

    def writeToVTK(self, fname):
        cdef char *filename = tmr_convert_str_to_chars(fname)
        self.ptr.writeToVTK(filename)

    def getData(self):
        cdef int nu, ku
        cdef const double *Tu = NULL
        cdef const double *wts = NULL
        cdef const TMRPoint *pts = NULL
        cdef TMRBsplineCurve *curve = NULL
        curve = _dynamicBsplineCurve(self.ptr)

        if curve is not NULL:
            curve.getData(&nu, &ku, &Tu, &wts, &pts)
            tu = np.zeros(nu+ku)
            w = np.ones((nu))
            p = np.zeros((nu, 3))
            for i in range(nu+ku):
                tu[i] = Tu[i]
            if wts is not NULL:
                for i in range(nu):
                    w[i] = wts[i]
            for i in range(nu):
                p[i,0] = pts[i].x
                p[i,1] = pts[i].y
                p[i,2] = pts[i].z
            return ku, tu, w, p
        return None

    def split(self, double t):
        cdef TMRBsplineCurve *curve = NULL
        cdef TMRBsplineCurve *c1 = NULL
        cdef TMRBsplineCurve *c2 = NULL
        curve = _dynamicBsplineCurve(self.ptr)

        if curve is not NULL:
            curve.split(t, &c1, &c2)
            return _init_Curve(c1), _init_Curve(c2)
        return None

cdef _init_Curve(TMRCurve *ptr):
    curve = Curve()
    curve.ptr = ptr
    if ptr != NULL:
        curve.ptr.incref()
    return curve

cdef class Pcurve:
    cdef TMRPcurve *ptr
    def __cinit__(self):
        self.ptr = NULL

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def setName(self, aname):
        cdef char *name = tmr_convert_str_to_chars(aname)
        if self.ptr:
            self.ptr.setName(name)

    def getEntityId(self):
        if self.ptr:
            return self.ptr.getEntityId()
        return -1

cdef class Surface:
    cdef TMRSurface *ptr
    def __cinit__(self):
        self.ptr = NULL

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def setName(self, aname):
        cdef char *name = tmr_convert_str_to_chars(aname)
        if self.ptr:
            self.ptr.setName(name)

    def writeToVTK(self, fname):
        cdef char *filename = tmr_convert_str_to_chars(fname)
        self.ptr.writeToVTK(filename)

    def getData(self):
        cdef int nu, nv, ku, kv
        cdef const double *Tu = NULL
        cdef const double *Tv = NULL
        cdef const double *wts = NULL
        cdef const TMRPoint *pts = NULL
        cdef TMRBsplineSurface *surf = NULL
        surf = _dynamicBsplineSurface(self.ptr)

        if surf is not NULL:
            surf.getData(&nu, &nv, &ku, &kv, &Tu, &Tv, &wts, &pts)
            tu = np.zeros(nu+ku)
            tv = np.zeros(nv+kv)
            w = np.ones((nu*nv))
            p = np.zeros((nu*nv, 3))
            for i in range(nu+ku):
                tu[i] = Tu[i]
            for i in range(nv+kv):
                tv[i] = Tv[i]
            if wts is not NULL:
                for i in range(nu*nv):
                    w[i] = wts[i]
            for i in range(nv*nu):
                p[i,0] = pts[i].x
                p[i,1] = pts[i].y
                p[i,2] = pts[i].z
            return ku, kv, tu, tv, w, p
        return None

cdef _init_Surface(TMRSurface *ptr):
    surface = Surface()
    surface.ptr = ptr
    if ptr != NULL:
        surface.ptr.incref()
    return surface

cdef class BsplineCurve(Curve):
    def __cinit__(self, np.ndarray[double, ndim=2, mode='c'] pts, int k=4):
        cdef int nctl = pts.shape[0]
        cdef ku = k
        if ku > nctl:
            ku = nctl
        cdef TMRPoint* p = <TMRPoint*>malloc(nctl*sizeof(TMRPoint))
        for i in range(nctl):
            p[i].x = pts[i,0]
            p[i].y = pts[i,1]
            p[i].z = pts[i,2]
        self.ptr = new TMRBsplineCurve(nctl, ku, p)
        self.ptr.incref()
        free(p)

cdef class BsplinePcurve(Pcurve):
    def __cinit__(self, np.ndarray[double, ndim=2, mode='c'] pts,
                  np.ndarray[double, ndim=1, mode='c'] tu=None,
                  np.ndarray[double, ndim=1, mode='c'] wts=None, int k=4):
        cdef int nctl = pts.shape[0]
        cdef ku = k
        cdef double *ctu = NULL
        cdef double *cwts = NULL
        if ku > nctl:
            ku = nctl
        self.ptr = NULL
        if tu is not None and wts is not None:
            if len(tu) != nctl+ku:
                errmsg = 'Incorrect BsplinePcurve knot length'
                raise ValueError(errmsg)
            self.ptr = new TMRBsplinePcurve(nctl, ku, <double*>tu.data,
                                            <double*>wts.data,
                                            <double*>pts.data)
        elif tu is not None and wts is None:
            self.ptr = new TMRBsplinePcurve(nctl, ku, <double*>tu.data,
                                            <double*>pts.data)
        elif tu is None and wts is None:
            self.ptr = new TMRBsplinePcurve(nctl, ku, <double*>pts.data)
        else:
            errmsg = 'BsplinePcurve: must supply knots and weights'
            raise ValueError(errmsg)

        self.ptr.incref()
        return

cdef class BsplineSurface(Surface):
    def __cinit__(self, np.ndarray[double, ndim=3, mode='c'] pts,
                  int ku=4, int kv=4):
        cdef int nx = pts.shape[0]
        cdef int ny = pts.shape[1]
        cdef kx = ku
        cdef ky = kv
        cdef TMRPoint* p = <TMRPoint*>malloc(nx*ny*sizeof(TMRPoint))
        if kx > nx:
            kx = nx
        if ky > ny:
            ky = ny
        for j in range(ny):
            for i in range(nx):
                p[i + j*nx].x = pts[i,j,0]
                p[i + j*nx].y = pts[i,j,1]
                p[i + j*nx].z = pts[i,j,2]
        self.ptr = new TMRBsplineSurface(nx, ny, kx, ky, p)
        self.ptr.incref()
        free(p)

cdef class VertexFromPoint(Vertex):
    def __cinit__(self, np.ndarray[double, ndim=1, mode='c'] pt):
        cdef TMRPoint point
        point.x = pt[0]
        point.y = pt[1]
        point.z = pt[2]
        self.ptr = new TMRVertexFromPoint(point)
        self.ptr.incref()

cdef class VertexFromEdge(Vertex):
    def __cinit__(self, Edge edge, double t):
        self.ptr = new TMRVertexFromEdge(edge.ptr, t)
        self.ptr.incref()

cdef class VertexFromFace(Vertex):
    def __cinit__(self, Face face, double u, double v):
        self.ptr = new TMRVertexFromFace(face.ptr, u, v)
        self.ptr.incref()

cdef class EdgeFromFace(Edge):
    def __cinit__(self, Face face, Pcurve pcurve, int degen=0):
        self.ptr = new TMREdgeFromFace(face.ptr, pcurve.ptr, degen)
        self.ptr.incref()

    def addEdgeFromFace(self, Face face, Pcurve pcurve):
        cdef TMREdgeFromFace *ef = NULL
        ef = _dynamicEdgeFromFace(self.ptr)
        if ef:
            ef.addEdgeFromFace(face.ptr, pcurve.ptr)

cdef class EdgeFromCurve(Edge):
    def __cinit__(self, Curve curve):
        self.ptr = new TMREdgeFromCurve(curve.ptr)
        self.ptr.incref()

cdef class FaceFromSurface(Face):
    def __cinit__(self, Surface surf):
        self.ptr = new TMRFaceFromSurface(surf.ptr)
        self.ptr.incref()

cdef class TFIEdge(Edge):
    def __cinit__(self, Vertex v1, Vertex v2):
        self.ptr = new TMRTFIEdge(v1.ptr, v2.ptr)
        self.ptr.incref()

cdef class TFIFace(Face):
    def __cinit__(self, list edges, list dirs=None, list verts=None):
        cdef TMREdge *e[4]
        cdef int d[4]
        cdef TMRVertex *v[4]
        cdef list edge_list

        if len(edges) != 4:
            errmsg = 'TFIFace: Number of edges must be 4'
            raise ValueError(errmsg)

        if dirs is None or verts is None:
            edge_list = edges
            edges = [edge_list.pop()]
            dirs = [1]
            v1, vnext = edges[0].getVertices()
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

            # Remove the last entry from the list of vertices
            verts.pop()

        if len(dirs) != 4:
            errmsg = 'TFIFace: Number of edge directions must be 4'
            raise ValueError(errmsg)
        if len(verts) != 4:
            errmsg = 'TFIFace: Number of vertices must be 4'
            raise ValueError(errmsg)
        for i in range(4):
            e[i] = (<Edge>edges[i]).ptr
            v[i] = (<Vertex>verts[i]).ptr
            d[i] = <int>dirs[i]
        self.ptr = new TMRTFIFace(e, d, v)
        self.ptr.incref()

cdef class CurveInterpolation:
    cdef TMRCurveInterpolation *ptr
    def __cinit__(self, np.ndarray[double, ndim=2, mode='c'] pts):
        cdef int nctl = pts.shape[0]
        cdef TMRPoint* p = <TMRPoint*>malloc(nctl*sizeof(TMRPoint))
        for i in range(nctl):
            p[i].x = pts[i,0]
            p[i].y = pts[i,1]
            p[i].z = pts[i,2]
        self.ptr = new TMRCurveInterpolation(p, nctl)
        self.ptr.incref()
        free(p)

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def setNumControlPoints(self, int nctl):
        self.ptr.setNumControlPoints(nctl)
        return

    def createCurve(self, int ku):
        cdef TMRBsplineCurve *curve = self.ptr.createCurve(ku)
        return _init_Curve(curve)

cdef class CurveLofter:
    cdef TMRCurveLofter *ptr
    def __cinit__(self, curves):
        cdef int ncurves = len(curves)
        cdef TMRBsplineCurve **crvs = NULL
        cdef TMRBsplineCurve *bspline = NULL
        crvs = <TMRBsplineCurve**>malloc(ncurves*sizeof(TMRBsplineCurve*))
        for i in range(ncurves):
            bspline =  _dynamicBsplineCurve((<Curve>curves[i]).ptr)
            if bspline != NULL:
               crvs[i] = bspline
            else:
                errstr = 'CurveLofter: Lofting curves must be BsplineCurves'
                raise ValueError(errstr)
        self.ptr = new TMRCurveLofter(crvs, ncurves)
        self.ptr.incref()
        free(crvs)

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def createSurface(self, int kv):
        cdef TMRSurface *surf = self.ptr.createSurface(kv)
        return _init_Surface(surf)

cdef class Model:
    """
    Contains an ordered collection of vertices, edges, faces, and volumes that
    define a model geometry.
    """
    cdef TMRModel *ptr
    def __cinit__(self, verts=None, edges=None, faces=None, vols=None):
        # Set the pointer to NULL
        self.ptr = NULL
        cdef int nverts = 0
        cdef TMRVertex **v = NULL
        if verts is not None:
            nverts = len(verts)
            v = <TMRVertex**>malloc(nverts*sizeof(TMRVertex*))
            for i in range(len(verts)):
                v[i] = (<Vertex>verts[i]).ptr
                if v[i] is NULL:
                    errmsg = 'Vertex %d is NULL'%(i)
                    raise ValueError(errmsg)

        cdef int nedges = 0
        cdef TMREdge **e = NULL
        if edges is not None:
            nedges = len(edges)
            e = <TMREdge**>malloc(nedges*sizeof(TMREdge*))
            for i in range(len(edges)):
                e[i] = (<Edge>edges[i]).ptr
                if e[i] is NULL:
                    errmsg = 'Edge %d is NULL'%(i)
                    raise ValueError(errmsg)

        cdef int nfaces = 0
        cdef TMRFace **f = NULL
        if faces is not None:
            nfaces = len(faces)
            f = <TMRFace**>malloc(nfaces*sizeof(TMRFace*))
            for i in range(len(faces)):
                f[i] = (<Face>faces[i]).ptr
                if f[i] is NULL:
                    errmsg = 'Face %d is NULL'%(i)
                    raise ValueError(errmsg)

        cdef int nvols = 0
        cdef TMRVolume **b = NULL
        if vols is not None:
            nvols = len(vols)
            b = <TMRVolume**>malloc(nvols*sizeof(TMRVolume*))
            for i in range(len(vols)):
                b[i] = (<Volume>vols[i]).ptr
                if b[i] is NULL:
                    errmsg = 'Volume %d is NULL'%(i)
                    raise ValueError(errmsg)

        if v and e and f and b:
            self.ptr = new TMRModel(nverts, v, nedges, e, nfaces, f, nvols, b)
        elif v and e and f:
            self.ptr = new TMRModel(nverts, v, nedges, e, nfaces, f, 0, NULL)
        elif v and e:
            self.ptr = new TMRModel(nverts, v, nedges, e, 0, NULL, 0, NULL)

        if self.ptr:
            self.ptr.incref()

        if v: free(v)
        if e: free(e)
        if f: free(f)
        if b: free(b)
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def getVolumes(self):
        cdef TMRVolume **vol
        cdef int num_vol = 0
        if self.ptr:
            self.ptr.getVolumes(&num_vol, &vol)
        volumes = []
        for i in range(num_vol):
            volumes.append(_init_Volume(vol[i]))
        return volumes

    def getFaces(self):
        cdef TMRFace **f
        cdef int num_faces = 0
        if self.ptr:
            self.ptr.getFaces(&num_faces, &f)
        faces = []
        for i in range(num_faces):
            faces.append(_init_Face(f[i]))
        return faces

    def getEdges(self):
        cdef TMREdge **e
        cdef int num_edges = 0
        if self.ptr:
            self.ptr.getEdges(&num_edges, &e)
        edges = []
        for i in range(num_edges):
            edges.append(_init_Edge(e[i]))
        return edges

    def getVertices(self):
        cdef TMRVertex **v
        cdef int num_verts = 0
        if self.ptr:
            self.ptr.getVertices(&num_verts, &v)
        verts = []
        for i in range(num_verts):
            verts.append(_init_Vertex(v[i]))
        return verts

    def writeModelToTecplot(self, fname,
                            vlabels=True, elabels=True, flabels=True,
                            off_scale=0.0):
        '''Write a representation of the edge loops to a file'''
        fp = open(fname, 'w')
        fp.write('Variables = x, y, z, tx, ty, tz\n')

        # Write out the vertices
        verts = self.getVertices()
        index = 0
        for v in verts:
            pt = v.evalPoint()
            if vlabels:
                dx = off_scale*random.uniform(-1.0, 1.0)
                fp.write('TEXT CS=GRID3D, X=%e, Y=%e, Z=%e, T=\"Vertex %d\"\n'%(
                    pt[0]+dx, pt[1]+dx, pt[2]+dx, index))
            fp.write('Zone T = \"Vertex %d\"\n'%(index))
            fp.write('%e %e %e 0 0 0\n'%(pt[0], pt[1], pt[2]))
            index += 1

        # Write out the edges
        edges = self.getEdges()
        index = 0
        for e in edges:
            v1, v2 = e.getVertices()
            pt1 = v1.evalPoint()
            pt2 = v2.evalPoint()
            if elabels:
                dx = off_scale*random.uniform(-1.0, 1.0)
                pt = 0.5*(pt1 + pt2)
                fp.write('TEXT CS=GRID3D, X=%e, Y=%e, Z=%e, T=\"Edge %d\"\n'%(
                    pt[0]+dx, pt[1]+dx, pt[2]+dx, index))
            fp.write('Zone T = \"Edge %d\"\n'%(index))
            fp.write('%e %e %e  %e %e %e\n'%(pt1[0], pt1[1], pt1[2],
                pt2[0] - pt1[0], pt2[1] - pt1[1], pt2[2] - pt1[2]))
            fp.write('%e %e %e 0 0 0\n'%(pt2[0], pt2[1], pt2[2]))
            index += 1

        # Write out the faces
        faces = self.getFaces()
        index = 0
        for f in faces:
            xav = np.zeros(3)
            count = 0
            for k in range(f.getNumEdgeLoops()):
                fp.write('Zone T = \"Face %d Loop %d\"\n'%(index, k))

                loop = f.getEdgeLoop(k)
                e, dirs = loop.getEdgeLoop()
                pts = np.zeros((len(e)+1, 3))
                tx = np.zeros((len(e)+1, 3))
                for i in range(len(e)):
                    v1, v2 = e[i].getVertices()
                    if dirs[i] > 0:
                        pt1 = v1.evalPoint()
                        pt2 = v2.evalPoint()
                    else:
                        pt1 = v2.evalPoint()
                        pt2 = v1.evalPoint()
                    if i == 0:
                        pts[0,:] = pt1[:]
                    pts[i+1,:] = pt2[:]

                for i in range(len(e)):
                    tx[i,:] = pts[i+1,:] - pts[i,:]
                    xav[:] += 0.5*(pts[i+1,:] + pts[i,:])
                    count += 1

                for i in range(len(e)+1):
                    fp.write('%e %e %e %e %e %e\n'%(
                        pts[i,0], pts[i,1], pts[i,2],
                        tx[i,0], tx[i,1], tx[i,2]))

            if count != 0 and flabels:
                xav /= count
                dx = off_scale*random.uniform(-1.0, 1.0)
                fp.write('TEXT CS=GRID3D, X=%e, Y=%e, Z=%e, T=\"Face %d\"\n'%(
                    xav[0]+dx, xav[1]+dx, xav[2]+dx, index))

            index += 1
        return

cdef _init_Model(TMRModel* ptr):
    model = Model()
    model.ptr = ptr
    if ptr != NULL:
        model.ptr.incref()
    return model

cdef class MeshOptions:
    """
    Defines a number of options that modify the meshing algorithm.
    """
    cdef TMRMeshOptions ptr
    def __cinit__(self):
        self.ptr = TMRMeshOptions()

    def __dealloc__(self):
        return

    property mesh_type_default:
        """
        Default mesh type either structured or unstructured. In the case that it
        is set to structured, the algorithm firsts checks if it is possible to
        use a mapped mesh, and reverts to an unstructured algorithm otherwise.
        """
        def __get__(self):
            return self.ptr.mesh_type_default
        def __set__(self, TMRFaceMeshType value):
            self.ptr.mesh_type_default = value

    property num_smoothing_steps:
        """
        Number of smoothing steps to apply during both the Laplacian and quad
        smoothing algorithms

        Parameters
        ----------
        value: int
           Number of smoothing steps
        """
        def __get__(self):
            return self.ptr.num_smoothing_steps
        def __set__(self, value):
            self.ptr.num_smoothing_steps=value

    property frontal_quality_factor:
        """
        Use the mesh quality indicator to determine when to accept new triangles
        in the frontal algorithm.

        Parameters
        ----------
        value: double
          Quality factor
        """
        def __get__(self):
            return self.ptr.frontal_quality_factor
        def __set__(self, value):
            self.ptr.frontal_quality_factor = value

    property triangularize_print_level:
        """
        Print level to provide more verbosity during the triangularization
        algorithm

        Parameters
        ----------
        value: int
           Print level for the operation
        """
        def __get__(self):
            return self.ptr.triangularize_print_level
        def __set__(self, value):
            self.ptr.triangularize_print_level = value

    property triangularize_print_iter:
        def __get__(self):
            return self.ptr.triangularize_print_iter
        def __set__(self, value):
            if value >= 1:
                self.ptr.triangularize_print_iter = value

    property reset_mesh_objects:
        def __get__(self):
            return self.ptr.reset_mesh_objects
        def __set__(self, value):
            self.ptr.reset_mesh_objects = value

    property write_mesh_quality_histogram:
        """
        Write out a histogram of the mesh quality in the final smoothed
        quadrilateral mesh.

        Parameters
        ----------
        value: bool
           Whether or not to write the mesh quality histogram
        """
        def __get__(self):
            return self.ptr.write_mesh_quality_histogram
        def __set__(self, value):
            self.ptr.write_mesh_quality_histogram = value

    property write_init_domain_triangle:
        def __get__(self):
            return self.ptr.write_init_domain_triangle
        def __set__(self, value):
            self.ptr.write_init_domain_triangle = value

    property write_triangularize_intermediate:
        def __get__(self):
            return self.ptr.write_triangularize_intermediate
        def __set__(self, value):
            self.ptr.write_triangularize_intermediate = value

    property write_pre_smooth_triangle:
        def __get__(self):
            return self.ptr.write_pre_smooth_triangle
        def __set__(self, value):
            self.ptr.write_pre_smooth_triangle = value

    property write_post_smooth_triangle:
        def __get__(self):
            return self.ptr.write_post_smooth_triangle
        def __set__(self, value):
            self.ptr.write_post_smooth_triangle = value

    property write_dual_recombine:
        def __get__(self):
            return self.ptr.write_dual_recombine
        def __set__(self, value):
            self.ptr.write_dual_recombine = value

    property write_pre_smooth_quad:
        def __get__(self):
            return self.ptr.write_pre_smooth_quad
        def __set__(self, value):
            self.ptr.write_pre_smooth_quad = value

    property write_post_smooth_quad:
        def __get__(self):
            return self.ptr.write_post_smooth_quad
        def __set__(self, value):
            self.ptr.write_post_smooth_quad = value

    property write_quad_dual:
        def __get__(self):
            return self.ptr.write_quad_dual
        def __set__(self, value):
            self.ptr.write_quad_dual = value

   # @property for cython 0.26 and above
   # def num_smoothing_steps(self):
   #    return self.ptr.num_smoothing_steps
   # @num_smoothing_steps.setter
   # def num_smoothing_steps(self, value):
   #    self.ptr.num_smoothing_steps = value

cdef class ElementFeatureSize:
    """
    Base class for creating feature size
    """
    cdef TMRElementFeatureSize *ptr
    def __cinit__(self, *args, **kwargs):
        self.ptr = NULL

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def getFeatureSize(self, x):
        cdef TMRPoint pt
        pt.x = x[0]
        pt.y = x[1]
        pt.z = x[2]
        return self.ptr.getFeatureSize(pt)

cdef class ConstElementSize(ElementFeatureSize):
    def __cinit__(self, double h):
        self.ptr = new TMRElementFeatureSize(h)
        self.ptr.incref()

cdef class LinearElementSize(ElementFeatureSize):
    def __cinit__(self, double hmin, double hmax,
                  double c=0.0, double ax=0.0,
                  double ay=0.0, double az=0.0):
        self.ptr = new TMRLinearElementSize(hmin, hmax, c, ax, ay, az)
        self.ptr.incref()

cdef class BoxFeatureSize(ElementFeatureSize):
    cdef TMRBoxFeatureSize *bptr
    def __cinit__(self, xlow, xhigh, double hmin, double hmax):
        cdef TMRPoint p1
        cdef TMRPoint p2
        p1.x = xlow[0]
        p1.y = xlow[1]
        p1.z = xlow[2]
        p2.x = xhigh[0]
        p2.y = xhigh[1]
        p2.z = xhigh[2]
        self.bptr = new TMRBoxFeatureSize(p1, p2, hmin, hmax)
        self.ptr = self.bptr
        self.ptr.incref()

    def addBox(self, xlow, xhigh, double h):
        cdef TMRPoint p1
        cdef TMRPoint p2
        p1.x = xlow[0]
        p1.y = xlow[1]
        p1.z = xlow[2]
        p2.x = xhigh[0]
        p2.y = xhigh[1]
        p2.z = xhigh[2]
        self.bptr.addBox(p1, p2, h)

cdef class PointFeatureSize(ElementFeatureSize):
    def __cinit__(self, np.ndarray[double, ndim=2, mode='c'] X,
                  np.ndarray[double, ndim=1, mode='c'] hvals,
                  double hmin, double hmax, int num_sample_pts=16):
        cdef int npts = 0
        cdef TMRPoint *pts
        if X.shape[0] != hvals.shape[0]:
            errmsg = 'PointFeatureSize arrays must be same size'
            raise ValueError(errmsg)
        elif X.shape[1] != 3:
            errmsg = 'PointFeatureSize expecting point (n,3) array'
            raise ValueError(errmsg)

        # Convert to input data
        npts = hvals.shape[0]
        pts = <TMRPoint*>malloc(npts*sizeof(TMRPoint))
        for i in range(npts):
            pts[i].x = X[i,0]
            pts[i].y = X[i,1]
            pts[i].z = X[i,2]
        self.ptr = new TMRPointFeatureSize(npts, pts, <double*>hvals.data,
                                           hmin, hmax, num_sample_pts)
        self.ptr.incref()
        free(pts)
        return

cdef class PointLocator:
    def __cinit__(self, np.ndarray[double, ndim=2, mode='c'] X):
        cdef int npts = 0
        cdef TMRPoint *pts
        npts = X.shape[0]
        pts = <TMRPoint*>malloc(npts*sizeof(TMRPoint))
        for i in range(npts):
            pts[i].x = X[i,0]
            pts[i].y = X[i,1]
            pts[i].z = X[i,2]
        self.ptr = new TMRPointLocator(npts, pts)
        self.ptr.incref()
        free(pts)
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def locateClosest(self, x,
                      np.ndarray[int, ndim=1] index,
                      np.ndarray[double, ndim=1] dist):
        if index.shape[0] != dist.shape[0]:
            errmsg = 'PointLocator expects equal length input/output arrays'
            raise ValueError(errmsg)

        cdef int num_found = 0
        cdef int K = index.shape[0]
        cdef TMRPoint pt
        pt.x = x[0]
        pt.y = x[1]
        pt.z = x[2]
        self.ptr.locateClosest(K, pt, &num_found,
                               <int*>index.data, <double*>dist.data)
        return num_found

cdef class Mesh:
    """
    Mesh the geometry model. This class handles the meshing for surface objects
    without any additional information. For hexahedral meshes, the model must
    have a source or target.
    """
    cdef TMRMesh *ptr
    def __cinit__(self, MPI.Comm comm, Model geo):
        cdef MPI_Comm c_comm = NULL
        self.ptr = NULL
        if comm is not None:
            c_comm = comm.ob_mpi
            self.ptr = new TMRMesh(c_comm, geo.ptr)
            self.ptr.incref()

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def mesh(self, double h=1.0, MeshOptions opts=None,
             ElementFeatureSize fs=None):
        """
        mesh(self, h, opts, fs)

        Mesh the model with the provided mesh spacing *h*, default meshing
        options :class:`~TMR.MeshOptions` opts if it is not provided or given
        :class:`~TMR.ElementFeatureSize` fs

        Parameters
        ----------
        h: double
           Mesh spacing

        opts: :class:`~TMR.MeshOptions`
           Meshing options

        fs: :class:`~TMR.ElementFeatureSize`
           Element Feature Size options
        """
        cdef TMRMeshOptions default
        if fs is not None:
            if opts is None:
                self.ptr.mesh(default, fs.ptr)
            else:
                self.ptr.mesh(opts.ptr, fs.ptr)
        else:
            if opts is None:
                self.ptr.mesh(default, h)
            else:
                self.ptr.mesh(opts.ptr, h)

    def getMeshPoints(self):
        """
        getMeshPoints(self)

        Retrieve a global array of the mesh locations
        """
        cdef TMRPoint *X
        cdef int npts = 0
        npts = self.ptr.getMeshPoints(&X)
        Xp = np.zeros((npts, 3), dtype=np.double)
        for i in range(npts):
            Xp[i,0] = X[i].x
            Xp[i,1] = X[i].y
            Xp[i,2] = X[i].z
        return Xp

    def getQuadConnectivity(self):
        """
        getQuadConnectivity(self)

        Retrieve the global connectivity from the quadrilateral mesh
        """
        cdef const int *quads = NULL
        cdef int nquads = 0

        self.ptr.getQuadConnectivity(&nquads, &quads)
        q = np.zeros((nquads, 4), dtype=np.intc)
        for i in range(nquads):
            q[i,0] = quads[4*i]
            q[i,1] = quads[4*i+1]
            q[i,2] = quads[4*i+2]
            q[i,3] = quads[4*i+3]
        return q

    def getTriConnectivity(self):
        """
        getTriConnectivity(self)

        Retrieve the global connectivity from the triangular mesh
        """
        cdef const int *tris = NULL
        cdef int ntris = 0

        self.ptr.getTriConnectivity(&ntris, &tris)
        t = np.zeros((ntris, 3), dtype=np.intc)
        for i in range(ntris):
            t[i,0] = tris[3*i]
            t[i,1] = tris[3*i+1]
            t[i,2] = tris[3*i+2]
        return t

    def getHexConnectivity(self):
        """
        getHexConnectivity(self)

        Retrieve the global connectivity from the hexahedral mesh
        """
        cdef const int *hex = NULL
        cdef int nhex = 0

        self.ptr.getHexConnectivity(&nhex, &hex)
        he = np.zeros((nhex,8), dtype=np.intc)
        for i in range(nhex):
            he[i,0] = hex[8*i]
            he[i,1] = hex[8*i+1]
            he[i,2] = hex[8*i+2]
            he[i,3] = hex[8*i+3]
            he[i,4] = hex[8*i+4]
            he[i,5] = hex[8*i+5]
            he[i,6] = hex[8*i+6]
            he[i,7] = hex[8*i+7]
        return he

    def createModelFromMesh(self):
        """
        createModelFromMesh(self)

        Create a geometry model based on the input mesh
        """
        cdef TMRModel *model = NULL
        model = self.ptr.createModelFromMesh()
        return _init_Model(model)

    def writeToBDF(self, fname, outtype=None):
        """
        writeToBDF(self, fname, outtype=None)

        Write both the quadrilateral and hexahedral mesh to a BDF file

        Parameters
        ----------
        fname: str
               filename

        outtype: str
                 Type of mesh to output to BDF file i.e. quad or hex
        """
        cdef char *filename = tmr_convert_str_to_chars(fname)
        cdef int flag = 3
        if outtype is None:
            flag = 3
        elif outtype == 'quad':
            flag = 1
        elif outtype == 'hex':
            flag = 2
        self.ptr.writeToBDF(filename, flag)

    def writeToVTK(self, fname, outtype=None):
        """
        writeToVTK(self, fname, outtype=None)

        Write both the quadrilateral and hexahedral mesh to a VTK file

        Parameters
        ----------
        fname: str
               filename

        outtype: str
                 Type of mesh to output to VTK file i.e. quad or hex
        """
        cdef char *filename = tmr_convert_str_to_chars(fname)
        cdef int flag = 3
        if outtype is None:
            flag = 3
        elif outtype == 'quad':
            flag = 1
        elif outtype == 'hex':
            flag = 2
        self.ptr.writeToVTK(filename, flag)

cdef class EdgeMesh:
    """
    This is the class that stores the node numbers along an edge
    """
    cdef TMREdgeMesh *ptr
    def __cinit__(self, MPI.Comm comm, Edge e,
                  np.ndarray[double, ndim=2, mode='c'] _X=None):
        cdef MPI_Comm c_comm = comm.ob_mpi
        cdef TMRPoint *X = NULL
        cdef int npts = 0
        if _X is not None:
            npts = _X.shape[0]
            X = <TMRPoint*>malloc(npts*sizeof(TMRPoint))
            for i in range(npts):
                X[i].x = _X[i,0]
                X[i].y = _X[i,1]
                X[i].z = _X[i,2]
        self.ptr = new TMREdgeMesh(c_comm, e.ptr, X, npts)
        if X:
            free(X)
        self.ptr.incref()

    def __dealloc__(self):
        pass

    def mesh(self, double h, MeshOptions opts=None):
        cdef TMRMeshOptions options
        cdef TMRElementFeatureSize *fs = NULL
        fs = new TMRElementFeatureSize(h)
        fs.incref()
        if opts is None:
            self.ptr.mesh(options, fs)
        else:
            self.ptr.mesh(opts.ptr, fs)
        fs.decref()

cdef class FaceMesh:
    """
    This is the class that stores the structured/unstructured surface mesh
    """
    cdef TMRFaceMesh *ptr
    def __cinit__(self, MPI.Comm comm=None, Face f=None,
                  np.ndarray[double, ndim=2, mode='c'] _X=None,
                  np.ndarray[int, ndim=2, mode='c'] _quads=None):
        cdef MPI_Comm c_comm = NULL
        cdef TMRPoint *X = NULL
        cdef int *quads = NULL
        cdef int npts = 0
        cdef int nquads = 0
        if _X is not None and _quads is not None:
            c_comm = comm.ob_mpi
            npts = _X.shape[0]
            nquads = _quads.shape[0]
            X = <TMRPoint*>malloc(npts*sizeof(TMRPoint))
            quads = <int*>malloc(4*nquads*sizeof(int))
            for i in range(npts):
                X[i].x = _X[i,0]
                X[i].y = _X[i,1]
                X[i].z = _X[i,2]
            for i in range(nquads):
                quads[4*i] = _quads[i,0]
                quads[4*i+1] = _quads[i,1]
                quads[4*i+2] = _quads[i,2]
                quads[4*i+3] = _quads[i,3]
            self.ptr = new TMRFaceMesh(c_comm, f.ptr, X, npts, quads, nquads)
            if X:
                free(X)
            if quads:
                free(quads)
            self.ptr.incref()

    def __dealloc__(self):
        pass

    def mesh(self, double h, MeshOptions opts=None):
        cdef TMRMeshOptions options
        cdef TMRElementFeatureSize *fs = NULL
        fs = new TMRElementFeatureSize(h)
        fs.incref()
        if opts is None:
            self.ptr.mesh(options, fs)
        else:
            self.ptr.mesh(opts.ptr, fs)
        fs.decref()

    def writeToVTK(self, fname):
        cdef char *filename = tmr_convert_str_to_chars(fname)
        self.ptr.writeToVTK(filename)

cdef _init_FaceMesh(TMRFaceMesh *ptr):
    fm = FaceMesh()
    fm.ptr = ptr
    if ptr != NULL:
        fm.ptr.incref()
    return fm

cdef class VolumeMesh:
    """
    This is the class that stores the hexahedral mesh
    """
    cdef TMRVolumeMesh *ptr
    def __cinit__(self, MPI.Comm comm=None, Volume v=None):
        cdef MPI_Comm c_comm = NULL
        c_comm = comm.ob_mpi
        self.ptr = new TMRVolumeMesh(c_comm, v.ptr)
        self.ptr.incref()

    def __dealloc__(self):
        pass

cdef class Topology:
    """
    The main topology class that contains the objects used to build the
    underlying mesh.

    This class takes in a general :class:`~TMR.Model`, but there are additional
    requirements that are placed on the model to create a proper
    topology object. These requirements are as follows:
        #. No edge can degenerate to a vertex.
        #. No face can degenerate to an edge or vertex.
        #. All faces must be surrounded by a single edge loop with 4
        #. All volumes must contain 6 non-degenerate faces that are
           ordered in coordinate ordering as shown below. Furthermore, all
           volumes must be of type :class:`~TMR.TFIVolume`.
    """
    cdef TMRTopology *ptr
    def __cinit__(self, MPI.Comm comm=None, Model m=None):
        cdef MPI_Comm c_comm = NULL
        cdef TMRModel *model = NULL
        self.ptr = NULL
        if comm is not None and m is not None:
            c_comm = comm.ob_mpi
            model = m.ptr
            self.ptr = new TMRTopology(c_comm, model)
            self.ptr.incref()

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def getVolume(self, int index):
        cdef TMRVolume *volume
        self.ptr.getVolume(index, &volume)
        return _init_Volume(volume)

    def getFace(self, int index):
        cdef TMRFace *face
        self.ptr.getFace(index, &face)
        return _init_Face(face)

    def getEdge(self, int index):
        cdef TMREdge *edge
        self.ptr.getEdge(index, &edge)
        return _init_Edge(edge)

    def getvertex(self, int index):
        cdef TMRVertex *vert
        self.ptr.getVertex(index, &vert)
        return _init_Vertex(vert)

cdef _init_Topology(TMRTopology *ptr):
    topo = Topology()
    topo.ptr = ptr
    if ptr != NULL:
        topo.ptr.incref()
    return topo

cdef class QuadrantArray:
    cdef TMRQuadrantArray *ptr
    cdef int self_owned
    def __cinit__(self):
        self.ptr = NULL

    def __dealloc__(self):
        if self.ptr and self.self_owned:
            del self.ptr

    def __len__(self):
        cdef int size = 0
        self.ptr.getArray(NULL, &size)
        return size

    def __getitem__(self, int k):
        cdef int size = 0
        cdef TMRQuadrant *array
        self.ptr.getArray(&array, &size)
        if k < 0 or k >= size:
            errmsg = 'Quadrant array index %d out of range [0,%d)'%(k, size)
            raise IndexError(errmsg)
        quad = Quadrant()
        quad.x = array[k].x
        quad.y = array[k].y
        quad.level = array[k].level
        quad.info = array[k].info
        quad.face = array[k].face
        quad.tag = array[k].tag
        return quad

    def __setitem__(self, int k, Quadrant quad):
        cdef int size = 0
        cdef TMRQuadrant *array
        self.ptr.getArray(&array, &size)
        if k < 0 or k >= size:
            errmsg = 'Quadrant array index %d out of range [0,%d)'%(k, size)
            raise IndexError(errmsg)
        array[k].x = quad.x
        array[k].y = quad.y
        array[k].level = quad.level
        array[k].info = quad.info
        array[k].face = quad.face
        array[k].tag = quad.tag
        return

    def findIndex(self, Quadrant quad, use_nodes=False):
        cdef int size = 0
        cdef TMRQuadrant *array
        cdef TMRQuadrant *t
        cdef int index = 0
        cdef int _use_nodes = 0
        if use_nodes:
            _use_nodes = 1
        self.ptr.getArray(&array, &size)
        t = self.ptr.contains(&quad.quad, _use_nodes)
        if t == NULL:
            return None
        index = t - array
        return index

cdef _init_QuadrantArray(TMRQuadrantArray *array, int self_owned):
    arr = QuadrantArray()
    arr.ptr = array
    arr.self_owned = self_owned
    return arr

cdef class Quadrant:
    cdef TMRQuadrant quad
    def __cinit__(self):
        self.quad.face = 0
        self.quad.x = 0
        self.quad.y = 0
        self.quad.tag = 0
        self.quad.level = 0
        self.quad.info = 0

    property face:
        def __get__(self):
            return self.quad.face
        def __set__(self, value):
            self.quad.face = value

    property x:
        def __get__(self):
            return self.quad.x
        def __set__(self, value):
            self.quad.x = value

    property y:
        def __get__(self):
            return self.quad.y
        def __set__(self, value):
            self.quad.y = value

    property tag:
        def __get__(self):
            return self.quad.tag
        def __set__(self, value):
            self.quad.tag = value

    property level:
        def __get__(self):
            return self.quad.level
        def __set__(self, value):
            self.quad.level = value

    property info:
        def __get__(self):
            return self.quad.info
        def __set__(self, value):
            self.quad.info = value

cdef class QuadForest:
    """
    A parallel forest of quadtrees

    This class defines a parallel forest of quadrtrees. The connectivity
    between quadtrees is defined at a global level. The quadrants can
    easily be redistributed across processors using the repartition()
    call.
    """
    cdef TMRQuadForest *ptr
    def __cinit__(self, MPI.Comm comm=None, int order=2,
                  TMRInterpolationType interp=GAUSS_LOBATTO_POINTS):
        cdef MPI_Comm c_comm = NULL
        self.ptr = NULL
        if comm is not None:
            c_comm = comm.ob_mpi
            self.ptr = new TMRQuadForest(c_comm, order, interp)
            self.ptr.incref()

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def setMeshOrder(self, int order,
                     TMRInterpolationType interp=GAUSS_LOBATTO_POINTS):
        self.ptr.setMeshOrder(order, interp)

    def getMeshOrder(self):
        return self.ptr.getMeshOrder()

    def setTopology(self, Topology topo):
        self.ptr.setTopology(topo.ptr)

    def getTopology(self):
        cdef TMRTopology *topo = self.ptr.getTopology()
        if topo is not NULL:
            return _init_Topology(topo)
        return None

    def repartition(self):
        self.ptr.repartition()

    def createTrees(self, int depth=0):
        self.ptr.createTrees(depth)

    def createRandomTrees(self, int nrand=10, int min_lev=0, int max_lev=8):
        self.ptr.createRandomTrees(nrand, min_lev, max_lev)

    def refine(self, np.ndarray[int, ndim=1, mode='c'] refine=None,
               int min_lev=0, int max_lev=MAX_LEVEL):
        if refine is not None:
            self.ptr.refine(<int*>refine.data, min_lev, max_lev)
        else:
            self.ptr.refine(NULL, min_lev, max_lev)
        return

    def duplicate(self):
        cdef TMRQuadForest *dup = NULL
        dup = self.ptr.duplicate()
        return _init_QuadForest(dup)

    def coarsen(self):
        cdef TMRQuadForest *dup = NULL
        dup = self.ptr.coarsen()
        return _init_QuadForest(dup)

    def balance(self, int btype):
        self.ptr.balance(btype)

    def createNodes(self):
        self.ptr.createNodes()

    def getQuadsWithName(self, aname):
        cdef char *name = tmr_convert_str_to_chars(aname)
        cdef TMRQuadrantArray *array = NULL
        array = self.ptr.getQuadsWithName(name)
        return _init_QuadrantArray(array, 1)

    def getNodesWithName(self, aname):
        cdef char *name = tmr_convert_str_to_chars(aname)
        cdef int size = 0
        cdef int *nodes = NULL
        size = self.ptr.getNodesWithName(name, &nodes)
        array = np.zeros(size, dtype=np.intc)
        for i in range(size):
            array[i] = nodes[i]
        _deleteMe(nodes)
        return array

    def getQuadrants(self):
        cdef TMRQuadrantArray *array = NULL
        self.ptr.getQuadrants(&array)
        return _init_QuadrantArray(array, 0)

    def getPoints(self):
        cdef TMRPoint *X = NULL
        cdef int npts = 0
        npts = self.ptr.getPoints(&X)
        if X != NULL:
            Xp = np.zeros((npts, 3), dtype=np.double)
            for i in range(npts):
                Xp[i,0] = X[i].x
                Xp[i,1] = X[i].y
                Xp[i,2] = X[i].z
            return Xp
        else:
            errmsg = 'TMRQuadForest: No node locations'
            raise RuntimeError(errmsg)

    def getLocalNodeNumber(self, int node):
        return self.ptr.getLocalNodeNumber(node)

    def getNodeRange(self):
        cdef int size = 0
        cdef const int *node_range = NULL
        size = self.ptr.getOwnedNodeRange(&node_range)
        if node_range != NULL:
            r = np.zeros(size+1, dtype=np.intc)
            for i in range(size+1):
                r[i] = node_range[i]
            return r
        else:
            errmsg = 'TMRQuadForest: No node range'
            raise RuntimeError(errmsg)

    def getMeshConn(self):
        """
        getMeshConn(self)
        Return the distributed mesh connectivity using a global ordering
        """
        cdef const int *conn
        cdef int nelems
        cdef int order = self.ptr.getMeshOrder()
        self.ptr.getNodeConn(&conn, &nelems)
        if conn != NULL:
            quads = np.zeros((order*order*nelems), dtype=np.intc)
            for i in range(order*order*nelems):
                quads[i] = conn[i]
            return quads.reshape((nelems, order*order))
        else:
            errmsg = 'TMRQuadForest: No mesh connectivity'
            raise RuntimeError(errmsg)

    def getDepNodeConn(self):
        cdef int ndep = 0
        cdef const int *_ptr = NULL
        cdef const int *_conn = NULL
        cdef const double *_weights = NULL
        ndep = self.ptr.getDepNodeConn(&_ptr, &_conn, &_weights)
        ptr = np.zeros(ndep+1, dtype=np.intc)
        conn = np.zeros(_ptr[ndep], dtype=np.intc)
        weights = np.zeros(_ptr[ndep], dtype=np.double)
        for i in range(ndep+1):
            ptr[i] = _ptr[i]
        for i in range(ptr[ndep]):
            conn[i] = _conn[i]
            weights[i] = _weights[i]
        return ptr, conn, weights

    def writeToVTK(self, fname):
        cdef char *filename = tmr_convert_str_to_chars(fname)
        self.ptr.writeToVTK(filename)

    def writeForestToVTK(self, fname):
        cdef char *filename = tmr_convert_str_to_chars(fname)
        self.ptr.writeForestToVTK(filename)

    def createInterpolation(self, QuadForest forest, VecInterp vec):
        if self.ptr == forest.ptr:
            errmsg = 'Cannot interpolate between the same object'
            raise ValueError(errmsg)
        self.ptr.createInterpolation(forest.ptr, vec.ptr)

cdef _init_QuadForest(TMRQuadForest* ptr):
    forest = QuadForest()
    forest.ptr = ptr
    if ptr != NULL:
        forest.ptr.incref()
    return forest

cdef class OctantArray:
    cdef TMROctantArray *ptr
    cdef int self_owned
    def __cinit__(self):
        self.self_owned = 0
        self.ptr = NULL

    def __dealloc__(self):
        if self.ptr and self.self_owned:
            del self.ptr

    def __len__(self):
        cdef int size = 0
        self.ptr.getArray(NULL, &size)
        return size

    def __getitem__(self, int k):
        cdef int size = 0
        cdef TMROctant *array
        self.ptr.getArray(&array, &size)
        if k < 0 or k >= size:
            errmsg = 'Octant array index %d out of range [0,%d)'%(k, size)
            raise IndexError(errmsg)
        oc = Octant()
        oc.x = array[k].x
        oc.y = array[k].y
        oc.z = array[k].z
        oc.level = array[k].level
        oc.info = array[k].info
        oc.block = array[k].block
        oc.tag = array[k].tag
        return oc

    def __setitem__(self, int k, Octant oc):
        cdef int size = 0
        cdef TMROctant *array
        self.ptr.getArray(&array, &size)
        if k < 0 or k >= size:
            errmsg = 'Octant array index %d out of range [0,%d)'%(k, size)
            raise IndexError(errmsg)
        array[k].x = oc.x
        array[k].y = oc.y
        array[k].z = oc.z
        array[k].level = oc.level
        array[k].info = oc.info
        array[k].block = oc.block
        array[k].tag = oc.tag
        return

    def findIndex(self, Octant oc, use_nodes=False):
        cdef int size = 0
        cdef TMROctant *array
        cdef TMROctant *t
        cdef int index = 0
        cdef int _use_nodes = 0
        if use_nodes:
            _use_nodes = 1
        self.ptr.getArray(&array, &size)
        t = self.ptr.contains(&oc.octant, _use_nodes)
        if t == NULL:
            return None
        index = t - array
        return index

cdef _init_OctantArray(TMROctantArray *array, int self_owned):
    arr = OctantArray()
    arr.ptr = array
    arr.self_owned = self_owned
    return arr

cdef class Octant:
    cdef TMROctant octant
    def __cinit__(self):
        self.octant.block = 0
        self.octant.x = 0
        self.octant.y = 0
        self.octant.z = 0
        self.octant.tag = 0
        self.octant.level = 0
        self.octant.info = 0

    property block:
        def __get__(self):
            return self.octant.block
        def __set__(self, value):
            self.octant.block = value

    property x:
        def __get__(self):
            return self.octant.x
        def __set__(self, value):
            self.octant.x = value

    property y:
        def __get__(self):
            return self.octant.y
        def __set__(self, value):
            self.octant.y = value

    property z:
        def __get__(self):
            return self.octant.z
        def __set__(self, value):
            self.octant.z = value

    property tag:
        def __get__(self):
            return self.octant.tag
        def __set__(self, value):
            self.octant.tag = value

    property level:
        def __get__(self):
            return self.octant.level
        def __set__(self, value):
            self.octant.level = value

    property info:
        def __get__(self):
            return self.octant.info
        def __set__(self, value):
            self.octant.info = value

cdef class OctForest:
    """
    This class defines a forest of octrees. The octrees within the
    forest can be distributed across processors. The connectivity
    between octrees is defined on all processors by setting a octree to
    node connectivity.

    The octrees can be redistributed across processors by using the
    repartition function. This destroys the nodes that may have been
    created (but can easily be recomputed).
    """
    cdef TMROctForest *ptr
    def __cinit__(self, MPI.Comm comm=None, int order=2,
                  TMRInterpolationType interp=GAUSS_LOBATTO_POINTS):
        cdef MPI_Comm c_comm = NULL
        self.ptr = NULL
        if comm is not None:
            c_comm = comm.ob_mpi
            self.ptr = new TMROctForest(c_comm, order, interp)
            self.ptr.incref()

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def setMeshOrder(self, int order,
                     TMRInterpolationType interp=GAUSS_LOBATTO_POINTS):
        self.ptr.setMeshOrder(order, interp)

    def getMeshOrder(self):
        return self.ptr.getMeshOrder()

    def setTopology(self, Topology topo):
        """
        setTopology(self, topo)
        Set the topology (and determine the connectivity)

        Parameters
        ----------
        topo:  :class:`~TMR.Topology`
          Set topo as the underlying topology object
        """
        self.ptr.setTopology(topo.ptr)

    def getTopology(self):
        cdef TMRTopology *topo = self.ptr.getTopology()
        if topo is not NULL:
            return _init_Topology(topo)
        return None

    def setConnectivity(self, np.ndarray[int, ndim=2, mode='c'] conn):
        '''Set non-degenerate connectivity'''
        cdef int num_nodes = 0
        cdef int num_blocks = 0
        if conn.shape[1] != 8:
            errmsg = 'Expected array of shape (n, 8)'
            raise ValueError(errmsg)
        num_blocks = conn.shape[0]
        num_nodes = np.max(conn)+1
        self.ptr.setConnectivity(num_nodes, <int*>conn.data, num_blocks)

    def repartition(self, int max_rank=-1):
        """
        repartition(self, max_rank=-1)
        Repartition the mesh across processors. This redistributes the elements
        so that there are an equal, or nearly equal, number of elements on each
        processor.

        Parameters
        ----------
        max_rank: int
          Number of processors to distribute the mesh across. If negative, the
          mesh is distribute across all processors
        """
        self.ptr.repartition(max_rank)

    def createTrees(self, int depth=0):
        """
        createTrees(self, depth_level)
        Create all the octrees in the mesh with the given level of refinement.

        Parameters
        ----------
        refine_level: int
          Octree of refinement level depth
        """
        self.ptr.createTrees(depth)

    def createRandomTrees(self, int nrand=10, int min_lev=0, int max_lev=8):
        """
        createRandomTrees(self, nrand, min_level, max_level)
        Create a set of trees with random levels of refinement. This is useful
        for testing purposes.

        Parameters
        ----------
        nrand: int
          Number of random octrees
        min_level: int
          Minimum refinement level
        max_level: int
          Maximum refinement level
        """
        self.ptr.createRandomTrees(nrand, min_lev, max_lev)

    def refine(self, np.ndarray[int, ndim=1, mode='c'] refine=None,
               int min_lev=0, int max_lev=MAX_LEVEL):
        """
        refine(self, refine, min_level, max_level)
        Refine the elements in the mesh by the specified number of levels. If a
        negative number is supplied, coarsen the element.

        Parameters
        ----------
        refine: array of int
          Array of integers
        min_level: int
          Minimum level of refinement
        max_level: int
          Maximum refinement level
        """
        if refine is not None:
            self.ptr.refine(<int*>refine.data, min_lev, max_lev)
        else:
            self.ptr.refine(NULL, min_lev, max_lev)
        return

    def duplicate(self):
        """
        duplicate(self)
        Duplicate a forest by copying connectivity, topology and elements (but
        not duplicating nodes)

        Returns
        -------
        dup: :class:`~TMR.OctForest`
         The duplicate OctForest
        """
        cdef TMROctForest *dup = NULL
        dup = self.ptr.duplicate()
        return _init_OctForest(dup)

    def coarsen(self):
        """
        coarsen(self)
        Create a new forest object by coarsening all the elements within the
        mesh by one level, if possible. Does not create new nodes.

        Returns
        -------
        dup: :class:`~TMR.OctForest`
         The coarsened OctForest
        """
        cdef TMROctForest *dup = NULL
        dup = self.ptr.coarsen()
        return _init_OctForest(dup)

    def balance(self, int btype):
        """
        balance(self, btype)
        Balance all the elements in the mesh to achieve a 2-to-1 balance

        Parameters
        ----------
        btype: bool
          Whether or not to balance corners
        """
        self.ptr.balance(btype)

    def createNodes(self):
        """
        createNodes(self)
        Create all the nodes within the mesh
        """
        self.ptr.createNodes()

    def getOctsWithName(self, aname):
        """
        getOctsWithName(self, aname)
        Get an array of octants with the desired name

        Parameters
        ----------
        aname: str
          Get the octants with aname

        Returns
        -------
        array: :class:`~TMR.OctantArray`
           Array of octants with name
        """
        cdef char *name = tmr_convert_str_to_chars(aname)
        cdef TMROctantArray *array = NULL
        array = self.ptr.getOctsWithName(name)
        return _init_OctantArray(array, 1)

    def getNodesWithName(self, aname):
        """
        getNodesWithName(self, aname)
        Get an array of nodes with the desired name

        Parameters
        ----------
        aname: str
          Get the nodes with aname

        Returns
        -------
        array: array of int
           Array of node numbers with name
        """
        cdef char *name = tmr_convert_str_to_chars(aname)
        cdef int size = 0
        cdef int *nodes = NULL
        size = self.ptr.getNodesWithName(name, &nodes)
        array = np.zeros(size, dtype=np.intc)
        for i in range(size):
            array[i] = nodes[i]
        _deleteMe(nodes)
        return array

    def getOctants(self):
        cdef TMROctantArray *array = NULL
        self.ptr.getOctants(&array)
        return _init_OctantArray(array, 0)

    def getPoints(self):
        cdef TMRPoint *X = NULL
        cdef int npts = 0
        npts = self.ptr.getPoints(&X)
        Xp = np.zeros((npts, 3), dtype=np.double)
        for i in range(npts):
            Xp[i,0] = X[i].x
            Xp[i,1] = X[i].y
            Xp[i,2] = X[i].z
        return Xp

    def getNodeRange(self):
        cdef int size = 0
        cdef const int *node_range = NULL
        size = self.ptr.getOwnedNodeRange(&node_range)
        r = np.zeros(size+1, dtype=np.intc)
        for i in range(size+1):
            r[i] = node_range[i]
        return r

    def getMeshConn(self):
        """
        getMeshConn(self)
        Return the distributed mesh connectivity using a global ordering
        """
        cdef const int *conn
        cdef int nelems
        cdef int order = self.ptr.getMeshOrder()
        self.ptr.getNodeConn(&conn, &nelems)
        octs = np.zeros((order*order*order*nelems), dtype=np.intc)
        for i in range(order*order*order*nelems):
            octs[i] = conn[i]
        return octs.reshape((nelems, order*order*order))

    def getDepNodeConn(self):
        """
        getDepNodeConn(self)
        Return the dependent node connectivity, weights and number of dependent
        nodes from the octree object

        Return
        ------
        ptr: array of int
          Array of dependent nodes
        conn: array of int
          Array of connectivity of dependent nodes
        weights: array of double
          Array of weights associated with the dependent nodes
        """
        cdef int ndep = 0
        cdef const int *_ptr = NULL
        cdef const int *_conn = NULL
        cdef const double *_weights = NULL
        ndep = self.ptr.getDepNodeConn(&_ptr, &_conn, &_weights)
        ptr = np.zeros(ndep+1, dtype=np.intc)
        conn = np.zeros(_ptr[ndep], dtype=np.intc)
        weights = np.zeros(_ptr[ndep], dtype=np.double)
        for i in range(ndep+1):
            ptr[i] = _ptr[i]
        for i in range(ptr[ndep]):
            conn[i] = _conn[i]
            weights[i] = _weights[i]
        return ptr, conn, weights

    def writeToVTK(self, fname):
        cdef char *filename = tmr_convert_str_to_chars(fname)
        self.ptr.writeToVTK(filename)

    def writeForestToVTK(self, fname):
        cdef char *filename = tmr_convert_str_to_chars(fname)
        self.ptr.writeForestToVTK(filename)

    def createInterpolation(self, OctForest forest, VecInterp vec):
        """
        createInterpolation(self, forest, vec)
        Create an interpolation object between two octrees that share a common
        topology.

        Parameters
        ----------
        forest: :class:`~TMR.OctForest`
          The octree forest that has the common topology
        vec: :class:`~TACS.VecInterp`
          The interpolation operator created
        """
        if self.ptr == forest.ptr:
            errmsg = 'Cannot interpolate between the same object'
            raise ValueError(errmsg)
        self.ptr.createInterpolation(forest.ptr, vec.ptr)

cdef _init_OctForest(TMROctForest* ptr):
    forest = OctForest()
    forest.ptr = ptr
    if ptr != NULL:
        forest.ptr.incref()
    return forest

def LoadModel(fname, int print_lev=0):
    """
    LoadModel(fname, print_lev)

    Initialization of the OpenCascade geometry from an IGES/STEP files

    Parameters
    ----------
    fname: str
      Name of the geometry file
    print_level: int
      Print level for operation
    """
    cdef char *filename = tmr_convert_str_to_chars(fname)
    cdef TMRModel *model = NULL
    if fname.lower().endswith(('step', 'stp')):
        model = TMR_LoadModelFromSTEPFile(filename, print_lev)
    elif fname.lower().endswith(('igs', 'iges')):
        model = TMR_LoadModelFromIGESFile(filename, print_lev)
    elif fname.lower().endswith(('egads')):
        model = TMR_LoadModelFromEGADSFile(filename, print_lev)
    if model is NULL:
        errmsg = 'Error loading model. File %s does not exist?'%(fname)
        raise RuntimeError(errmsg)
    return _init_Model(model)

cdef class BoundaryConditions:
    cdef TMRBoundaryConditions* ptr
    def __cinit__(self):
        self.ptr = new TMRBoundaryConditions()
        self.ptr.incref()

    def __dealloc__(self):
        self.ptr.decref()

    def getNumBoundaryConditions(self):
        return self.ptr.getNumBoundaryConditions()

    def addBoundaryCondition(self, aname,
                             list bc_nums=None, list bc_vals=None):
        cdef char *name = tmr_convert_str_to_chars(aname)
        cdef int *nums = NULL
        cdef double *vals = NULL
        cdef int num_bcs = 0
        if bc_nums is not None and bc_vals is not None:
            if len(bc_nums) != len(bc_vals):
                errstr = 'Boundary condition lists must be the same length'
                raise ValueError(errstr)
            num_bcs = len(bc_nums)
            nums = <int*>malloc(num_bcs*sizeof(int))
            vals = <double*>malloc(num_bcs*sizeof(double))
            for i in range(len(bc_nums)):
                nums[i] = <int>bc_nums[i]
                vals[i] = <double>bc_vals[i]
            self.ptr.addBoundaryCondition(name, num_bcs, nums, vals)
            free(nums)
            free(vals)
        elif bc_nums is not None:
            num_bcs = len(bc_nums)
            nums = <int*>malloc(num_bcs*sizeof(int))
            for i in range(len(bc_nums)):
                nums[i] = <int>bc_nums[i]
            self.ptr.addBoundaryCondition(name, num_bcs, nums, NULL)
            free(nums)
        else:
            self.ptr.addBoundaryCondition(name, 0, NULL, NULL)
        return

cdef TACSElement* _createQuadElement(void *_self, int order,
                                     TMRQuadrant *quad):
    cdef TACSElement *elem = NULL
    q = Quadrant()
    q.quad.x = quad.x
    q.quad.y = quad.y
    q.quad.level = quad.level
    q.quad.info = quad.info
    q.quad.face = quad.face
    q.quad.tag = quad.tag
    e = (<object>_self).createElement(order, q)
    if e is not None:
        (<Element>e).ptr.incref()
        elem = (<Element>e).ptr
        return elem
    return NULL

cdef class QuadCreator:
    cdef TMRCyQuadCreator *ptr
    def __cinit__(self, BoundaryConditions bcs, *args, **kwargs):
        self.ptr = new TMRCyQuadCreator(bcs.ptr)
        self.ptr.incref()
        self.ptr.setSelfPointer(<void*>self)
        self.ptr.setCreateQuadElement(_createQuadElement)
        return

    def __dealloc__(self):
        self.ptr.decref()

    def createTACS(self, QuadForest forest,
                   OrderingType ordering=TACS.PY_NATURAL_ORDER):
        cdef TACSAssembler *assembler = NULL
        assembler = self.ptr.createTACS(forest.ptr, ordering)
        return _init_Assembler(assembler)

cdef TACSElement* _createOctElement(void *_self, int order,
                                    TMROctant *octant):
    cdef TACSElement *elem = NULL
    q = Octant()
    q.octant.x = octant.x
    q.octant.y = octant.y
    q.octant.z = octant.z
    q.octant.level = octant.level
    q.octant.info = octant.info
    q.octant.block = octant.block
    q.octant.tag = octant.tag
    e = (<object>_self).createElement(order, q)
    if e is not None:
        (<Element>e).ptr.incref()
        elem = (<Element>e).ptr
        return elem
    return NULL

cdef class OctCreator:
    cdef TMRCyOctCreator *ptr
    def __cinit__(self, BoundaryConditions bcs, *args, **kwargs):
        self.ptr = new TMRCyOctCreator(bcs.ptr)
        self.ptr.incref()
        self.ptr.setSelfPointer(<void*>self)
        self.ptr.setCreateOctElement(_createOctElement)
        return

    def __dealloc__(self):
        self.ptr.decref()

    def createTACS(self, OctForest forest,
                   OrderingType ordering=TACS.PY_NATURAL_ORDER):
        cdef TACSAssembler *assembler = NULL
        assembler = self.ptr.createTACS(forest.ptr, ordering)
        return _init_Assembler(assembler)

cdef TACSElement* _createQuadTopoElement(void *_self, int order,
                                         TMRQuadrant *quad,
                                         TMRIndexWeight *weights,
                                         int nweights):
    cdef TACSElement *elem = NULL
    q = Quadrant()
    q.quad.x = quad.x
    q.quad.y = quad.y
    q.quad.level = quad.level
    q.quad.info = quad.info
    q.quad.face = quad.face
    q.quad.tag = quad.tag
    idx = []
    wvals = []
    for i in range(nweights):
        idx.append(weights[i].index)
        wvals.append(weights[i].weight)
    e = (<object>_self).createElement(order, q, idx, wvals)
    if e is not None:
        (<Element>e).ptr.incref()
        elem = (<Element>e).ptr
        return elem
    return NULL

cdef class QuadTopoCreator:
    cdef TMRCyTopoQuadCreator *ptr
    def __cinit__(self, BoundaryConditions bcs, QuadForest filt, *args, **kwargs):
        self.ptr = new TMRCyTopoQuadCreator(bcs.ptr, filt.ptr)
        self.ptr.incref()
        self.ptr.setSelfPointer(<void*>self)
        self.ptr.setCreateQuadTopoElement(_createQuadTopoElement)
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def createTACS(self, QuadForest forest,
                   OrderingType ordering=TACS.PY_NATURAL_ORDER,scale=1.0):
        cdef TACSAssembler *assembler = NULL
        assembler = self.ptr.createTACS(forest.ptr, ordering, scale)
        return _init_Assembler(assembler)

    def getFilter(self):
        cdef TMRQuadForest *filtr = NULL
        self.ptr.getFilter(&filtr)
        return _init_QuadForest(filtr)

    def getMap(self):
        cdef TACSVarMap *vmap = NULL
        self.ptr.getMap(&vmap)
        return _init_VarMap(vmap)

    def getIndices(self):
        cdef TACSBVecIndices *indices = NULL
        self.ptr.getIndices(&indices)
        return _init_VecIndices(indices)

cdef TACSElement* _createQuadConformTopoElement( void *_self, int order,
                                                 TMRQuadrant *quad,
                                                 int *index,
                                                 int nweights,
                                                 TMRQuadForest *filtr ):
    cdef TACSElement *elem = NULL
    q = Quadrant()
    q.quad.x = quad.x
    q.quad.y = quad.y
    q.quad.level = quad.level
    q.quad.info = quad.info
    q.quad.face = quad.face
    q.quad.tag = quad.tag
    idx = []

    qf = _init_QuadForest(filtr)
    for i in range(nweights):
        idx.append(index[i])
    e = (<object>_self).createElement(order, q, idx, qf)
    if e is not None:
        (<Element>e).ptr.incref()
        elem = (<Element>e).ptr
        return elem
    return NULL

cdef class QuadConformTopoCreator:
    cdef TMRCyTopoQuadConformCreator *ptr
    def __cinit__(self, BoundaryConditions bcs, QuadForest forest,
                  int order=-1, TMRInterpolationType interp=GAUSS_LOBATTO_POINTS,
                  *args, **kwargs):
        self.ptr = new TMRCyTopoQuadConformCreator(bcs.ptr, forest.ptr, \
                                                   order, interp)
        self.ptr.incref()
        self.ptr.setSelfPointer(<void*>self)
        self.ptr.setCreateQuadTopoElement(_createQuadConformTopoElement)
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def createTACS(self, QuadForest forest,
                   OrderingType ordering=TACS.PY_NATURAL_ORDER,scale=1.0):
        cdef TACSAssembler *assembler = NULL
        assembler = self.ptr.createTACS(forest.ptr, ordering, scale)
        return _init_Assembler(assembler)

    def getFilter(self):
        cdef TMRQuadForest *filtr = NULL
        self.ptr.getFilter(&filtr)
        return _init_QuadForest(filtr)

cdef TACSElement* _createOctTopoElement(void *_self, int order,
                                        TMROctant *octant,
                                        TMRIndexWeight *weights,
                                        int nweights):
    cdef TACSElement *elem = NULL
    oct = Octant()
    oct.octant.x = octant.x
    oct.octant.y = octant.y
    oct.octant.z = octant.z
    oct.octant.level = octant.level
    oct.octant.info = octant.info
    oct.octant.block = octant.block
    oct.octant.tag = octant.tag
    idx = []
    wvals = []
    for i in range(nweights):
        idx.append(weights[i].index)
        wvals.append(weights[i].weight)
    e = (<object>_self).createElement(order, oct, idx, wvals)
    if e is not None:
        (<Element>e).ptr.incref()
        elem = (<Element>e).ptr
        return elem
    return NULL

cdef class OctTopoCreator:
    cdef TMRCyTopoOctCreator *ptr
    def __cinit__(self, BoundaryConditions bcs, OctForest filt,
                  *args, **kwargs):
        self.ptr = NULL
        self.ptr = new TMRCyTopoOctCreator(bcs.ptr, filt.ptr)
        self.ptr.incref()
        self.ptr.setSelfPointer(<void*>self)
        self.ptr.setCreateOctTopoElement(_createOctTopoElement)
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def createTACS(self, OctForest forest,
                   OrderingType ordering=TACS.PY_NATURAL_ORDER,
                   double scale=1.0):
        cdef TACSAssembler *assembler = NULL
        assembler = self.ptr.createTACS(forest.ptr, ordering, scale)
        return _init_Assembler(assembler)

    def getFilter(self):
        cdef TMROctForest *filtr = NULL
        self.ptr.getFilter(&filtr)
        return _init_OctForest(filtr)

    def getMap(self):
        cdef TACSVarMap *vmap = NULL
        self.ptr.getMap(&vmap)
        return _init_VarMap(vmap)

    def getIndices(self):
        cdef TACSBVecIndices *indices = NULL
        self.ptr.getIndices(&indices)
        return _init_VecIndices(indices)

cdef TACSElement* _createOctConformTopoElement( void *_self, int order,
                                                TMROctant *octant,
                                                int *index,
                                                int nweights,
                                                TMROctForest *filtr):
    cdef TACSElement *elem = NULL
    Oct = Octant()
    Oct.octant.x = octant.x
    Oct.octant.y = octant.y
    Oct.octant.z = octant.z
    Oct.octant.level = octant.level
    Oct.octant.info = octant.info
    Oct.octant.block = octant.block
    Oct.octant.tag = octant.tag
    idx = []

    of = _init_OctForest(filtr)
    for i in range(nweights):
        idx.append(index[i])
    e = (<object>_self).createElement(order, Oct, idx, of)
    if e is not None:
        (<Element>e).ptr.incref()
        elem = (<Element>e).ptr
        return elem
    return NULL

cdef class OctConformTopoCreator:
    cdef TMRCyTopoOctConformCreator *ptr
    def __cinit__(self, BoundaryConditions bcs, OctForest forest,
                  int order=-1, TMRInterpolationType interp=GAUSS_LOBATTO_POINTS,
                  *args, **kwargs):
        self.ptr = new TMRCyTopoOctConformCreator(bcs.ptr, forest.ptr,
                                                  order, interp)
        self.ptr.incref()
        self.ptr.setSelfPointer(<void*>self)
        self.ptr.setCreateOctTopoElement(_createOctConformTopoElement)
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def createTACS(self, OctForest forest,
                   OrderingType ordering=TACS.PY_NATURAL_ORDER,
                   double scale=1.0):
        cdef TACSAssembler *assembler = NULL
        assembler = self.ptr.createTACS(forest.ptr, ordering, scale)
        return _init_Assembler(assembler)

    def getFilter(self):
        cdef TMROctForest *filtr = NULL
        self.ptr.getFilter(&filtr)
        return _init_OctForest(filtr)

cdef class StiffnessProperties:
    cdef TMRStiffnessProperties *ptr
    def __cinit__(self, list _rho, list _E, list _nu, list _ys=None,
                  list _aT=None, list _kcond=None,
                  double q=5.0, double qtemp=0.0, double qcond=0.0,
                  double eps=0.3, double k0=1e-6,
                  double beta=15.0, double xoffset=0.5, int use_project=0):
        cdef TacsScalar *rho = NULL
        cdef TacsScalar *E = NULL
        cdef TacsScalar *nu = NULL
        cdef TacsScalar *ys = NULL
        cdef TacsScalar *aT = NULL
        cdef TacsScalar *kcond = NULL
        cdef int nmats = 0
        if len(_rho) != len(_E) or len(_rho) != len(_nu):
            errmsg = 'Must specify the same number of properties'
            raise ValueError(errmsg)
        nmats = len(_E)
        rho = <TacsScalar*>malloc(nmats*sizeof(TacsScalar));
        E = <TacsScalar*>malloc(nmats*sizeof(TacsScalar));
        nu = <TacsScalar*>malloc(nmats*sizeof(TacsScalar));
        if (_ys):
            ys = <TacsScalar*>malloc(nmats*sizeof(TacsScalar));
        if (_aT):
            aT = <TacsScalar*>malloc(nmats*sizeof(TacsScalar));
        if (_kcond):
            kcond = <TacsScalar*>malloc(nmats*sizeof(TacsScalar));
        for i in range(nmats):
            rho[i] = <TacsScalar>_rho[i]
            E[i] = <TacsScalar>_E[i]
            nu[i] = <TacsScalar>_nu[i]
            if (_ys):
                ys[i] = <TacsScalar>_ys[i]
            if (_aT):
                aT[i]= <TacsScalar>_aT[i]
            if (_kcond):
                kcond[i] = <TacsScalar>_kcond[i]
        self.ptr = new TMRStiffnessProperties(nmats, q, eps, k0, beta, xoffset,
                                              rho, E, nu, ys, aT, kcond,
                                              qtemp, qcond, use_project)
        self.ptr.incref()
        free(rho)
        free(E)
        free(nu)
        if (ys):
            free(ys)
        if (aT):
            free(aT)
        if (kcond):
            free(kcond)
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def getNumMaterials(self):
        return self.ptr.nmats

    property q:
        def __get__(self):
            return self.ptr.q
        def __set__(self, value):
            self.ptr.q = value

    property eps:
        def __get__(self):
            return self.ptr.eps
        def __set__(self, value):
            self.ptr.eps = value

cdef class QuadStiffnessProperties:
    cdef TMRQuadStiffnessProperties *ptr
    def __cinit__(self, list _rho, list _E, list _nu, list _ys=None,
                  list _aT=None, list _kcond=None,
                  double q=5.0, double qtemp=0.0, double qcond=0.0,
                  double eps=0.3, double k0=1e-6,
                  double beta=15.0, double xoffset=0.5, int use_project=0):
        cdef TacsScalar *rho = NULL
        cdef TacsScalar *E = NULL
        cdef TacsScalar *nu = NULL
        cdef TacsScalar *ys = NULL
        cdef TacsScalar *aT = NULL
        cdef TacsScalar *kcond = NULL
        cdef int nmats = 0
        if len(_rho) != len(_E) or len(_rho) != len(_nu):
            errmsg = 'Must specify the same number of properties'
            raise ValueError(errmsg)
        nmats = len(_E)
        rho = <TacsScalar*>malloc(nmats*sizeof(TacsScalar));
        E = <TacsScalar*>malloc(nmats*sizeof(TacsScalar));
        nu = <TacsScalar*>malloc(nmats*sizeof(TacsScalar));
        if (_ys):
            ys = <TacsScalar*>malloc(nmats*sizeof(TacsScalar));
        if (_aT):
            aT = <TacsScalar*>malloc(nmats*sizeof(TacsScalar));
        if (_kcond):
            kcond = <TacsScalar*>malloc(nmats*sizeof(TacsScalar));
        for i in range(nmats):
            rho[i] = <TacsScalar>_rho[i]
            E[i] = <TacsScalar>_E[i]
            nu[i] = <TacsScalar>_nu[i]
            if (_ys):
                ys[i] = <TacsScalar>_ys[i]
            if (_aT):
                aT[i]= <TacsScalar>_aT[i]
            if (_kcond):
                kcond[i] = <TacsScalar>_kcond[i]
        self.ptr = new TMRQuadStiffnessProperties(nmats, q, eps,
                                                  k0, beta, xoffset,
                                                  rho, E, nu, ys, aT, kcond,
                                                  qtemp, qcond, use_project)
        self.ptr.incref()
        free(rho)
        free(E)
        free(nu)
        if (ys):
            free(ys)
        if (aT):
            free(aT)
        if (kcond):
            free(kcond)
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def getNumMaterials(self):
        return self.ptr.nmats

    property q:
        def __get__(self):
            return self.ptr.q
        def __set__(self, value):
            self.ptr.q = value

    property eps:
        def __get__(self):
            return self.ptr.eps
        def __set__(self, value):
            self.ptr.eps = value

cdef class AnisotropicProperties:
    cdef TMRAnisotropicProperties *ptr
    def __cinit__(self, list _rho, list _C,
                  double q=5.0, double k0=1e-6,
                  double beta=15.0, double xoffset=0.5, int use_project=0):
        cdef TacsScalar *rho = NULL
        cdef TacsScalar *C = NULL
        cdef int nmats = 0
        if len(_rho) != len(_C):
            errmsg = 'Must specify the correct number of properties'
            raise ValueError(errmsg)
        nmats = len(_rho)
        rho = <TacsScalar*>malloc(nmats*sizeof(TacsScalar));
        C = <TacsScalar*>malloc(21*nmats*sizeof(TacsScalar));

        for i in range(nmats):
            rho[i] = <TacsScalar>_rho[i]
        for i in range(nmats):
            for j in range(21):
                C[21*i+j] = <TacsScalar>_C[i][j]
        self.ptr = new TMRAnisotropicProperties(nmats, q, k0, beta, xoffset,
                                                rho, C, use_project)
        self.ptr.incref()
        free(rho)
        free(C)
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def getNumMaterials(self):
        return self.ptr.nmats

    property q:
        def __get__(self):
            return self.ptr.q
        def __set__(self, value):
            self.ptr.q = value

cdef class OctStiffness(SolidStiff):
    def __cinit__(self, StiffnessProperties props,
                  list index=None, list weights=None):
        cdef TMRIndexWeight *w = NULL
        cdef int nw = 0
        self.ptr = NULL
        if weights is None or index is None:
            errmsg = 'Must define weights and indices'
            raise ValueError(errmsg)
        if len(weights) != len(index):
            errmsg = 'Weights and index list lengths must be the same'
            raise ValueError(errmsg)

        # Extract the weights
        nw = len(weights)
        w = <TMRIndexWeight*>malloc(nw*sizeof(TMRIndexWeight));
        for i in range(nw):
            w[i].weight = <double>weights[i]
            w[i].index = <int>index[i]

        # Create the constitutive object
        self.ptr = new TMROctStiffness(w, nw, props.ptr)
        self.ptr.incref()
        free(w)
        return

cdef class ThermoOctStiffness(CoupledSolid):
    def __cinit__(self, StiffnessProperties props,
                  list index=None, list weights=None,
                  OctForest filtr=None):
        cdef TMRIndexWeight *w = NULL
        cdef int* ind = NULL
        cdef int nw = 0
        self.ptr = NULL
        cdef TMROctForest *of = NULL
        # if weights is None or index is None:
        #     errmsg = 'Must define weights and indices'
        #     raise ValueError(errmsg)
        # if len(weights) != len(index):
        #     errmsg = 'Weights and index list lengths must be the same'
        #     raise ValueError(errmsg)
        # Extract the weights
        if weights:
            nw = len(weights)
            w = <TMRIndexWeight*>malloc(nw*sizeof(TMRIndexWeight));
            for i in range(nw):
                w[i].weight = <double>weights[i]
                w[i].index = <int>index[i]
        else:
            nw = len(index)
            ind = <int*>malloc(nw*sizeof(int));
            for i in range(nw):
                ind[i] = <int>index[i]
        if filtr:
            of = filtr.ptr
        # Create the constitutive object
        self.ptr = new TMRCoupledThermoOctStiffness(w, nw, props.ptr, of, ind)
        self.ptr.incref()
        if w:
            free(w)
        if ind:
            free(ind)
        return

cdef class QuadStiffness(PlaneStress):
    def __cinit__(self, QuadStiffnessProperties props,
                  list index=None, list weights=None):
        cdef TMRIndexWeight *w = NULL
        cdef int nw = 0
        self.ptr = NULL
        if weights is None or index is None:
            errmsg = 'Must define weights and indices'
            raise ValueError(errmsg)
        if len(weights) != len(index):
            errmsg = 'Weights and index list lengths must be the same'
            raise ValueError(errmsg)

        # # Check that the lengths are less than 4
        # if len(weights) > 4:
        #     errmsg = 'Weight/index lists too long > 4'
        #     raise ValueError(errmsg)

        # Extract the weights
        nw = len(weights)
        w = <TMRIndexWeight*>malloc(nw*sizeof(TMRIndexWeight));
        for i in range(nw):
            w[i].weight = <double>weights[i]
            w[i].index = <int>index[i]

        # Create the constitutive object
        self.ptr = new TMRQuadStiffness(w, nw, props.ptr)
        self.ptr.incref()
        free(w)
        return

cdef class ThermoQuadStiffness(CoupledPlaneStress):
    def __cinit__(self, QuadStiffnessProperties props,
                  list index=None, list weights=None,
                  QuadForest filtr=None ):
        cdef TMRIndexWeight *w = NULL
        cdef int* ind = NULL
        cdef int nw = 0
        self.ptr = NULL
        cdef TMRQuadForest *qf = NULL
        # if weights is None or index is None:
        #     errmsg = 'Must define weights and indices'
        #     raise ValueError(errmsg)
        # if len(weights) != len(index):
        #     errmsg = 'Weights and index list lengths must be the same'
        #     raise ValueError(errmsg)

        # # Check that the lengths are less than 4
        # if len(weights) > 4:
        #     errmsg = 'Weight/index lists too long > 4'
        #     raise ValueError(errmsg)

        # Extract the weights
        if weights:
            nw = len(weights)
            w = <TMRIndexWeight*>malloc(nw*sizeof(TMRIndexWeight));
            for i in range(nw):
                w[i].weight = <double>weights[i]
                w[i].index = <int>index[i]
        else:
            nw = len(index)
            ind = <int*>malloc(nw*sizeof(int));
            for i in range(nw):
                ind[i] = <int>index[i]

        if filtr:
            qf = filtr.ptr
        # Create the constitutive object
        self.ptr =new TMRCoupledThermoQuadStiffness(w, nw, props.ptr, qf, ind)
        self.ptr.incref()
        if w:
            free(w)
        if ind:
            free(ind)
        return

cdef class AnisotropicStiffness(SolidStiff):
    def __cinit__(self, AnisotropicProperties props,
                  list index=None, list weights=None):
        cdef TMRIndexWeight *w = NULL
        cdef int nw = 0
        self.ptr = NULL
        if weights is None or index is None:
            errmsg = 'Must define weights and indices'
            raise ValueError(errmsg)
        if len(weights) != len(index):
            errmsg = 'Weights and index list lengths must be the same'
            raise ValueError(errmsg)

        # Check that the lengths are less than 8
        if len(weights) > 8:
            errmsg = 'Weight/index lists too long > 8'
            raise ValueError(errmsg)

        # Extract the weights
        nw = len(weights)
        w = <TMRIndexWeight*>malloc(nw*sizeof(TMRIndexWeight));
        for i in range(nw):
            w[i].weight = <double>weights[i]
            w[i].index = <int>index[i]

        # Create the constitutive object
        self.ptr = new TMRAnisotropicStiffness(w, nw, props.ptr)
        self.ptr.incref()
        free(w)

def createMg(list assemblers, list forests, double omega=1.0,
             use_coarse_direct_solve=True,
             use_chebyshev_smoother=False):
    cdef int nlevels = 0
    cdef TACSAssembler **assm = NULL
    cdef TMRQuadForest **qforest = NULL
    cdef TMROctForest **oforest = NULL
    cdef TACSMg *mg = NULL
    cdef int isqforest = 0
    cdef int coarse_direct = 0
    cdef int use_cheb = 0
    if use_coarse_direct_solve:
        coarse_direct = 1
    if use_chebyshev_smoother:
        use_cheb = 1

    if len(assemblers) != len(forests):
        errstr = 'Number of Assembler and Forest objects must be equal'
        raise ValueError(errstr)
    nlevels = len(assemblers)

    for i in range(nlevels):
        if isinstance(forests[i], QuadForest):
            isqforest = 1
        elif isinstance(forests[i], OctForest):
            isqforest = 0

    assm = <TACSAssembler**>malloc(nlevels*sizeof(TACSAssembler*))
    if isqforest:
        qforest = <TMRQuadForest**>malloc(nlevels*sizeof(TMRQuadForest*))
        for i in range(nlevels):
            assm[i] = (<Assembler>assemblers[i]).ptr
            qforest[i] = (<QuadForest>forests[i]).ptr
        TMR_CreateTACSMg(nlevels, assm, qforest, &mg, omega,
                         coarse_direct, use_cheb)
        free(qforest)
    else:
        oforest = <TMROctForest**>malloc(nlevels*sizeof(TMROctForest*))
        for i in range(nlevels):
            assm[i] = (<Assembler>assemblers[i]).ptr
            oforest[i] = (<OctForest>forests[i]).ptr
        TMR_CreateTACSMg(nlevels, assm, oforest, &mg, omega,
                         coarse_direct, use_cheb)
        free(oforest)
    free(assm)
    if mg != NULL:
        return _init_Mg(mg)
    return None

def strainEnergyError(forest, Assembler coarse,
                      forest_refined, Assembler refined):
    """
    strainEnergyError(forest, coarse_assembler, forest_refined,
                      refined_assembler)

    The following function performs a mesh refinement based on a strain
    energy criteria. It is based on the following relationship for
    linear finite-element analysis

    .. math::
         \Pi(u-u_h,u-u_h) = \Pi(u,u) - \Pi(u_h,u_h)

    where a(u,u) is the trilinear strain energy functional, u is the
    exact solution, and uh is the discretized solution at any mesh
    level. This relies on the relationship that :math:`\Pi(u_h, u - u_h) = 0`
    which is satisfied due to the method of Galerkin/Ritz.

    The following function computes a localized error indicator using
    the element-wise strain energy. The code computes a higher-order
    reconstructed solution using a cubic enrichment functions. These
    enrichment functions expand the original solution space and are
    computed based on a least-squares approximation with nodal gradient
    values. The localized error indicator is evaluated as follows:

    .. math::
         err = \sum_{i=1}^{4} \Pi_e(u_{ce}, u_{ce}) - \Pi_e(u_e, u_e)

    where :math:`u_{ce}` is the element-wise cubic element reconstruction
    projected onto a uniformly refined mesh.

    Parameters
    -----------
    forest_coarse: :class:`~TMR.OctForest` or :class:`~TMR.QuadForest`
      Forest for current mesh level
    coarse_assembler: :class:`~TACS.Assembler`
      Finite assembler class for associated with forest
    forest_refined: :class:`~TMR.OctForest` or :class:`~TMR.QuadForest`
      Forest for refined mesh level
    refined_assembler: :class:`~TACS.Assembler`
      Finite assembler class for associated with forest_refined

    Returns
    --------
    ans: double
      Total strain energy error

    err: array of double
      Elemental strain energy error
    """
    cdef double ans = 0.0
    cdef TMROctForest *oct_forest = NULL
    cdef TMROctForest *oct_forest_refined = NULL
    cdef TMRQuadForest *quad_forest = NULL
    cdef TMRQuadForest *quad_forest_refined = NULL
    cdef np.ndarray err = None
    err = np.zeros(coarse.ptr.getNumElements(), dtype=np.double)
    if isinstance(forest, OctForest):
        oct_forest = (<OctForest>forest).ptr
        oct_forest_refined = (<OctForest>forest_refined).ptr
        ans = TMR_StrainEnergyErrorEst(oct_forest, coarse.ptr,
                                       oct_forest_refined, refined.ptr,
                                       <double*>err.data)
    elif isinstance(forest, QuadForest):
        quad_forest = (<QuadForest>forest).ptr
        quad_forest_refined = (<QuadForest>forest_refined).ptr
        ans = TMR_StrainEnergyErrorEst(quad_forest, coarse.ptr,
                                       quad_forest_refined, refined.ptr,
                                       <double*>err.data)
    return ans, err

def adjointError(forest, Assembler coarse,
                 forest_refined, Assembler refined,
                 Vec solution, Vec adjoint):
    """
    adjointError(forest, coarse_assembler, forest_refined, refined_assembler,
                 solution_refined, adjoint_refined)

    Refine the mesh using the original solution and the adjoint solution

    Parameters
    -----------
    forest: :class:`~TMR.OctForest` or :class:`~TMR.QuadForest`
      Forest for current mesh level
    coarse_assembler: :class:`~TACS.Assembler`
      Finite assembler class for associated with forest
    forest_refined: :class:`~TMR.OctForest` or :class:`~TMR.QuadForest`
      Higher-order forest for refined mesh level
    refined_assembler: :class:`~TACS.Assembler`
      Higher-order finite assembler class for associated with forest_refined
    solution_refined: :class:`~TACS.Vec`
      The higher-order solution (or approximation)
    adjoint_refined: :class:`~TACS.Vec`
      The difference between the refined and coarse adjoint solutions computed
      in some manner

    Returns
    -------
    ans: double
      Total strain energy error

    err: array of double
      Elemental strain energy error

    adj_corr: TacsScalar
      Adjoint-based functional correction
    """
    cdef TacsScalar ans = 0.0
    cdef TacsScalar adj_corr = 0.0
    cdef TMROctForest *oct_forest = NULL
    cdef TMROctForest *oct_forest_refined = NULL
    cdef TMRQuadForest *quad_forest = NULL
    cdef TMRQuadForest *quad_forest_refined = NULL
    cdef np.ndarray err = None
    cdef TacsScalar err_est = 0.0
    err = np.zeros(coarse.ptr.getNumElements(), dtype=np.double)
    if isinstance(forest, OctForest):
        oct_forest = (<OctForest>forest).ptr
        oct_forest_refined = (<OctForest>forest_refined).ptr
        err_est = TMR_AdjointErrorEst(oct_forest, coarse.ptr,
                                      oct_forest_refined, refined.ptr,
                                      solution.ptr, adjoint.ptr,
                                      <double*>err.data, &adj_corr)
    elif isinstance(forest, QuadForest):
        quad_forest = (<QuadForest>forest).ptr
        quad_forest_refined = (<QuadForest>forest_refined).ptr
        err_est = TMR_AdjointErrorEst(quad_forest, coarse.ptr,
                                      quad_forest_refined, refined.ptr,
                                      solution.ptr, adjoint.ptr,
                                      <double*>err.data, &adj_corr)
    return err_est, adj_corr, err

def computeInterpSolution(forest, Assembler coarse,
                          forest_refined, Assembler refined,
                          Vec uvec=None, Vec uvec_refined=None,
                          compute_diff=False):
    cdef TMROctForest *oct_forest = NULL
    cdef TMROctForest *oct_forest_refined = NULL
    cdef TMRQuadForest *quad_forest = NULL
    cdef TMRQuadForest *quad_forest_refined = NULL
    cdef TACSBVec *uvec_ptr = NULL
    cdef TACSBVec *uvec_refined_ptr = NULL
    if uvec is not None:
        uvec_ptr = uvec.ptr
    if uvec_refined is not None:
        uvec_refined_ptr = uvec_refined.ptr
    if isinstance(forest, OctForest):
        oct_forest = (<OctForest>forest).ptr
        oct_forest_refined = (<OctForest>forest_refined).ptr
        TMR_ComputeInterpSolution(oct_forest, coarse.ptr,
                                  oct_forest_refined, refined.ptr,
                                  uvec_ptr, uvec_refined_ptr)
    elif isinstance(forest, QuadForest):
        quad_forest = (<QuadForest>forest).ptr
        quad_forest_refined = (<QuadForest>forest_refined).ptr
        TMR_ComputeInterpSolution(quad_forest, coarse.ptr,
                                  quad_forest_refined, refined.ptr,
                                  uvec_ptr, uvec_refined_ptr)
    return

def computeReconSolution(forest, Assembler coarse,
                         forest_refined, Assembler refined,
                         Vec uvec=None, Vec uvec_refined=None,
                         compute_diff=False):
    cdef TMROctForest *oct_forest = NULL
    cdef TMROctForest *oct_forest_refined = NULL
    cdef TMRQuadForest *quad_forest = NULL
    cdef TMRQuadForest *quad_forest_refined = NULL
    cdef TACSBVec *uvec_ptr = NULL
    cdef TACSBVec *uvec_refined_ptr = NULL
    cdef int diff = 0
    if compute_diff:
        diff = 1
    if uvec is not None:
        uvec_ptr = uvec.ptr
    if uvec_refined is not None:
        uvec_refined_ptr = uvec_refined.ptr
    if isinstance(forest, OctForest):
        oct_forest = (<OctForest>forest).ptr
        oct_forest_refined = (<OctForest>forest_refined).ptr
        TMR_ComputeReconSolution(oct_forest, coarse.ptr,
                                 oct_forest_refined, refined.ptr,
                                 uvec_ptr, uvec_refined_ptr, diff)
    elif isinstance(forest, QuadForest):
        quad_forest = (<QuadForest>forest).ptr
        quad_forest_refined = (<QuadForest>forest_refined).ptr
        TMR_ComputeReconSolution(quad_forest, coarse.ptr,
                                 quad_forest_refined, refined.ptr,
                                 uvec_ptr, uvec_refined_ptr, diff)
    return

def writeSTLToBin(fname, OctForest forest,
                  Vec x, int offset=0, double cutoff=0.5):
    cdef char *filename = tmr_convert_str_to_chars(fname)
    TMR_GenerateBinFile(filename, forest.ptr, x.ptr, offset, cutoff)
    return

cdef class StressConstraint:
    cdef TMRStressConstraint *ptr
    def __cinit__(self, OctForest oct, Assembler assembler,
                  TacsScalar ks_weight):
        self.ptr = NULL
        self.ptr = new TMRStressConstraint(oct.ptr,
                                           assembler.ptr,
                                           ks_weight)
        self.ptr.incref()

    def __dealloc__(self):
        if self.ptr:
           self.ptr.decref()

    def evalCon(self, Vec uvec):
        return self.ptr.evalConstraint(uvec.ptr)

    def evalConDeriv(self,
                     np.ndarray[double, ndim=1, mode='c'] dfdx,
                     Vec dfdu):
        cdef int size = len(dfdx)
        self.ptr.evalConDeriv(<TacsScalar*>dfdx.data, size,
                              dfdu.ptr)
        return

    def writeReconToTec(self, Vec uvec, fname, TacsScalar ys=1e6):
        cdef char *filename = tmr_convert_str_to_chars(fname)
        self.ptr.writeReconToTec(uvec.ptr, filename, ys)
        return

cdef class LagrangeFilter(TopoFilter):
    def __cinit__(self, list assemblers, list filters,
                  list varmaps=[], list varindices=[], int vars_per_node=1):
        cdef int nlevels = 0
        cdef int isqforest = 0
        cdef TACSAssembler **assemb = NULL
        cdef TMROctForest **ofiltr = NULL
        cdef TMRQuadForest **qfiltr = NULL
        cdef TACSVarMap **vmaps = NULL
        cdef TACSBVecIndices **vindex = NULL

        if (len(assemblers) != len(filters) or
            len(assemblers) != len(varmaps) or
            len(assemblers) != len(varindices)):
            errmsg = 'LagrangeFilter must have equal number of objects in lists'
            raise ValueError(errmsg)

        nlevels = len(assemblers)
        for i in range(nlevels):
            if isinstance(filters[i], QuadForest):
                isqforest = 1
            elif isinstance(filters[i], OctForest):
                isqforest = 0

        assemb = <TACSAssembler**>malloc(nlevels*sizeof(TACSAssembler*))
        vmaps = <TACSVarMap**>malloc(nlevels*sizeof(TACSVarMap*))
        vindex = <TACSBVecIndices**>malloc(nlevels*sizeof(TACSBVecIndices*))

        for i in range(nlevels):
            vmaps[i] = (<VarMap>varmaps[i]).ptr
            vindex[i] = (<VecIndices>varindices[i]).ptr

        if isqforest:
            qfiltr = <TMRQuadForest**>malloc(nlevels*sizeof(TMRQuadForest*))
            for i in range(nlevels):
                qfiltr[i] = (<QuadForest>filters[i]).ptr
                assemb[i] = (<Assembler>assemblers[i]).ptr
            self.ptr = new TMRLagrangeFilter(nlevels, assemb, qfiltr,
                                             vmaps, vindex, vars_per_node)
            self.ptr.incref()
            free(qfiltr)
        else:
            ofiltr = <TMROctForest**>malloc(nlevels*sizeof(TMROctForest*))
            for i in range(nlevels):
                ofiltr[i] = (<OctForest>filters[i]).ptr
                assemb[i] = (<Assembler>assemblers[i]).ptr
            self.ptr = new TMRLagrangeFilter(nlevels, assemb, ofiltr,
                                             vmaps, vindex, vars_per_node)
            self.ptr.incref()
            free(ofiltr)

        free(assemb)
        free(vmaps)
        free(vindex)
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def getMap(self):
        return _init_VarMap(self.ptr.getDesignVarMap())

cdef class ConformFilter(TopoFilter):
    def __cinit__(self, list assemblers, list filters, int vars_per_node=1):
        cdef int nlevels = 0
        cdef int isqforest = 0
        cdef TACSAssembler **assemb = NULL
        cdef TMROctForest **ofiltr = NULL
        cdef TMRQuadForest **qfiltr = NULL

        if len(assemblers) != len(filters):
            errmsg = 'ConformFilter must have equal number of objects in lists'
            raise ValueError(errmsg)

        nlevels = len(assemblers)
        for i in range(nlevels):
            if isinstance(filters[i], QuadForest):
                isqforest = 1
            elif isinstance(filters[i], OctForest):
                isqforest = 0

        assemb = <TACSAssembler**>malloc(nlevels*sizeof(TACSAssembler*))
        if isqforest:
            qfiltr = <TMRQuadForest**>malloc(nlevels*sizeof(TMRQuadForest*))
            for i in range(nlevels):
                qfiltr[i] = (<QuadForest>filters[i]).ptr
                assemb[i] = (<Assembler>assemblers[i]).ptr
            self.ptr = new TMRConformFilter(nlevels, assemb, qfiltr,
                                            vars_per_node)
            self.ptr.incref()
            free(qfiltr)
        else:
            ofiltr = <TMROctForest**>malloc(nlevels*sizeof(TMROctForest*))
            for i in range(nlevels):
                ofiltr[i] = (<OctForest>filters[i]).ptr
                assemb[i] = (<Assembler>assemblers[i]).ptr
            self.ptr = new TMRConformFilter(nlevels, assemb, ofiltr,
                                            vars_per_node)
            self.ptr.incref()
            free(ofiltr)

        free(assemb)
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def getMap(self):
        return _init_VarMap(self.ptr.getDesignVarMap())

cdef class HelmholtzFilter(TopoFilter):
    def __cinit__(self, double radius, list assemblers,
                  list filters, int vars_per_node=1):
        cdef int nlevels = 0
        cdef int isqforest = 0
        cdef TACSAssembler **assemb = NULL
        cdef TMROctForest **ofiltr = NULL
        cdef TMRQuadForest **qfiltr = NULL

        if len(assemblers) != len(filters):
            errmsg = 'HelmholtzFilter must have equal number of objects in lists'
            raise ValueError(errmsg)

        nlevels = len(assemblers)
        for i in range(nlevels):
            if isinstance(filters[i], QuadForest):
                isqforest = 1
            elif isinstance(filters[i], OctForest):
                isqforest = 0

        assemb = <TACSAssembler**>malloc(nlevels*sizeof(TACSAssembler*))
        if isqforest:
            qfiltr = <TMRQuadForest**>malloc(nlevels*sizeof(TMRQuadForest*))
            for i in range(nlevels):
                qfiltr[i] = (<QuadForest>filters[i]).ptr
                assemb[i] = (<Assembler>assemblers[i]).ptr
            self.ptr = new TMRHelmholtzFilter(radius, nlevels, assemb, qfiltr,
                                              vars_per_node)
            self.ptr.incref()
            free(qfiltr)
        else:
            ofiltr = <TMROctForest**>malloc(nlevels*sizeof(TMROctForest*))
            for i in range(nlevels):
                ofiltr[i] = (<OctForest>filters[i]).ptr
                assemb[i] = (<Assembler>assemblers[i]).ptr
            self.ptr = new TMRHelmholtzFilter(radius, nlevels, assemb, ofiltr,
                                              vars_per_node)
            self.ptr.incref()
            free(ofiltr)

        free(assemb)
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def getMap(self):
        return _init_VarMap(self.ptr.getDesignVarMap())

cdef class MatrixFilter(TopoFilter):
    def __cinit__(self, double s, int N, list assemblers,
                  list filters, int vars_per_node=1):
        cdef int nlevels = 0
        cdef int isqforest = 0
        cdef TACSAssembler **assemb = NULL
        cdef TMROctForest **ofiltr = NULL
        cdef TMRQuadForest **qfiltr = NULL

        if len(assemblers) != len(filters):
            errmsg = 'MatrixFilter must have equal number of objects in lists'
            raise ValueError(errmsg)

        nlevels = len(assemblers)
        for i in range(nlevels):
            if isinstance(filters[i], QuadForest):
                isqforest = 1
            elif isinstance(filters[i], OctForest):
                isqforest = 0

        assemb = <TACSAssembler**>malloc(nlevels*sizeof(TACSAssembler*))
        if isqforest:
            qfiltr = <TMRQuadForest**>malloc(nlevels*sizeof(TMRQuadForest*))
            for i in range(nlevels):
                qfiltr[i] = (<QuadForest>filters[i]).ptr
                assemb[i] = (<Assembler>assemblers[i]).ptr
            self.ptr = new TMRMatrixFilter(s, N, nlevels, assemb, qfiltr,
                                           vars_per_node)
            self.ptr.incref()
            free(qfiltr)
        else:
            ofiltr = <TMROctForest**>malloc(nlevels*sizeof(TMROctForest*))
            for i in range(nlevels):
                ofiltr[i] = (<OctForest>filters[i]).ptr
                assemb[i] = (<Assembler>assemblers[i]).ptr
            self.ptr = new TMRMatrixFilter(s, N, nlevels, assemb, ofiltr,
                                           vars_per_node)
            self.ptr.incref()
            free(ofiltr)

        free(assemb)
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def getMap(self):
        return _init_VarMap(self.ptr.getDesignVarMap())

cdef class TopoProblem(pyParOptProblemBase):
    def __cinit__(self, TopoFilter fltr, Pc pc,
                  int gmres_subspace=50, double rtol=1e-9):
        cdef TACSMg *mg = NULL

        # Check for a multigrid preconditioner
        mg = _dynamicTACSMg(pc.ptr)
        if mg == NULL:
            raise ValueError('TopoProblem requires a TACSMg preconditioner')

        self.ptr = new TMRTopoProblem(fltr.ptr, mg, gmres_subspace, rtol)
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def setF5OutputFlags(self, int freq, ElementType elem_type,
                         int flag):
        cdef TMRTopoProblem *prob = NULL
        prob = _dynamicTopoProblem(self.ptr)
        if prob == NULL:
            errmsg = 'Expected TMRTopoProblem got other type'
            raise ValueError(errmsg)
        prob.setF5OutputFlags(freq, elem_type, flag)
        return

    def setF5EigenOutputFlags(self, int freq, ElementType elem_type,
                             int flag):
        cdef TMRTopoProblem *prob = NULL
        prob = _dynamicTopoProblem(self.ptr)
        if prob == NULL:
            errmsg = 'Expected TMRTopoProblem got other type'
            raise ValueError(errmsg)
        prob.setF5EigenOutputFlags(freq, elem_type, flag)
        return

    def setLoadCases(self, list forces):
        cdef TACSBVec **f = NULL
        cdef int nforces = len(forces)
        cdef TMRTopoProblem *prob = NULL
        prob = _dynamicTopoProblem(self.ptr)
        if prob == NULL:
            errmsg = 'Expected TMRTopoProblem got other type'
            raise ValueError(errmsg)
        f = <TACSBVec**>malloc(nforces*sizeof(TACSBVec*))
        for i in range(nforces):
            f[i] = (<Vec>forces[i]).ptr
        prob.setLoadCases(f, nforces)
        free(f)
        return

    def getNumLoadCases(self):
        cdef TMRTopoProblem *prob = NULL
        prob = _dynamicTopoProblem(self.ptr)
        if prob == NULL:
            errmsg = 'Expected TMRTopoProblem got other type'
            raise ValueError(errmsg)
        return prob.getNumLoadCases()

    def addConstraints(self, int case, list funcs, list offset, list scale):
        cdef int nfuncs = 0
        cdef TacsScalar *_offset = NULL
        cdef TacsScalar *_scale = NULL
        cdef TACSFunction **f = NULL
        cdef TMRTopoProblem *prob = NULL
        prob = _dynamicTopoProblem(self.ptr)
        if prob == NULL:
            errmsg = 'Expected TMRTopoProblem got other type'
            raise ValueError(errmsg)
        if case < 0 or case >= prob.getNumLoadCases():
            errmsg = 'Load case out of expected range'
            raise ValueError(errmsg)
        if len(funcs) != len(offset) or len(funcs) != len(scale):
            errmsg = 'Expected equal function, offset and scale counts'
            raise ValueError(errmsg)

        nfuncs = len(funcs)
        f = <TACSFunction**>malloc(nfuncs*sizeof(TACSFunction*))
        _offset = <TacsScalar*>malloc(nfuncs*sizeof(TacsScalar))
        _scale = <TacsScalar*>malloc(nfuncs*sizeof(TacsScalar))
        for i in range(nfuncs):
            f[i] = (<Function>funcs[i]).ptr
            _offset[i] = <TacsScalar>offset[i]
            _scale[i] = <TacsScalar>scale[i]
        prob.addConstraints(case, f, _offset, _scale, nfuncs)
        free(f)
        free(_offset)
        free(_scale)
        return

    def addStressConstraint(self, int case, StressConstraint sc,
                            TacsScalar offset=1.0, TacsScalar scale=1.0):
        cdef TMRTopoProblem *prob = NULL
        prob = _dynamicTopoProblem(self.ptr)
        if prob == NULL:
            errmsg = 'Expected TMRTopoProblem got other type'
            raise ValueError(errmsg)
        prob.addStressConstraint(case, sc.ptr, offset, scale)
        return

    def addLinearConstraints(self, list vecs, list offset):
        cdef int nvecs
        cdef TacsScalar *_offset = NULL
        cdef ParOptVec **_vecs = NULL
        cdef TMRTopoProblem *prob = NULL
        prob = _dynamicTopoProblem(self.ptr)
        if prob == NULL:
            errmsg = 'Expected TMRTopoProblem got other type'
            raise ValueError(errmsg)
        if len(vecs) != len(offset):
            errmsg = 'Expected equal vector and offset counts'
            raise ValueError(errmsg)

        nvecs = len(vecs)
        _offset = <TacsScalar*>malloc(nvecs*sizeof(TacsScalar))
        _vecs = <ParOptVec**>malloc(nvecs*sizeof(ParOptVec*))
        for i in range(nvecs):
            _offset[i] = <TacsScalar>offset[i]
            _vecs[i] = (<PVec>vecs[i]).ptr
        prob.addLinearConstraints(_vecs, _offset, nvecs)
        free(_offset)
        free(_vecs)
        return

    def addFrequencyConstraint(self, double sigma, int num_eigvals,
                               TacsScalar ks_weight=30.0,
                               TacsScalar offset=0.0,
                               TacsScalar scale=0.0,
                               int max_subspace_size=100,
                               double eigtol=1e-8,
                               int use_jd=0, int fgmres_size=5,
                               double eig_rtol=1e-12,
                               double eig_atol=1e-30,
                               int num_recycle=0,
                               JDRecycleType recycle_type=JD_NUM_RECYCLE,
                               int track_eigen_iters=0):
        '''
        Add buckling/natural frequency constraints
        '''
        prob = _dynamicTopoProblem(self.ptr)
        if prob == NULL:
            errmsg = 'Expected TMRTopoProblem got other type'
            raise ValueError(errmsg)
        prob.addFrequencyConstraint(sigma, num_eigvals,
                                    ks_weight, offset,
                                    scale, max_subspace_size,
                                    eigtol, use_jd, fgmres_size,
                                    eig_rtol, eig_atol, num_recycle,
                                    recycle_type, track_eigen_iters)
        return

    def addBucklingConstraint(self, double sigma, int num_eigvals,
                              TacsScalar ks_weight=30.0,
                              TacsScalar offset=0.0, TacsScalar scale=0.0,
                              int max_lanczos=100, double eigtol=1e-8):
        '''
        Add buckling/natural frequency constraints
        '''
        prob = _dynamicTopoProblem(self.ptr)
        if prob == NULL:
            errmsg = 'Expected TMRTopoProblem got other type'
            raise ValueError(errmsg)
        prob.addBucklingConstraint(sigma, num_eigvals, ks_weight,
                                   offset, scale, max_lanczos, eigtol)
        return

    def setObjective(self, list weights, list funcs=None):
        cdef int lenw = 0
        cdef TacsScalar *w = NULL
        cdef TMRTopoProblem *prob = NULL
        prob = _dynamicTopoProblem(self.ptr)
        if prob == NULL:
            errmsg = 'Expected TMRTopoProblem got other type'
            raise ValueError(errmsg)
        lenw = len(weights)
        if lenw != prob.getNumLoadCases():
            errmsg = 'Incorrect number of weights'
            raise ValueError(errmsg)
        w = <TacsScalar*>malloc(lenw*sizeof(TacsScalar))
        for i in range(lenw):
            w[i] = weights[i]
        # Check if list of functions are provided
        cdef int nfuncs = 0
        cdef TACSFunction **f = NULL
        if funcs:
            # Get the objective function associated with each load case
            nfuncs = len(funcs)
            f = <TACSFunction**>malloc(nfuncs*sizeof(TACSFunction*))
            for i in range(nfuncs):
                f[i] = (<Function>funcs[i]).ptr
            prob.setObjective(w,f)
        else:
            prob.setObjective(w)
        free(w)
        if (f):
            free(f)
        return

    def initialize(self):
        cdef TMRTopoProblem *prob = NULL
        prob = _dynamicTopoProblem(self.ptr)
        if prob == NULL:
            errmsg = 'Expected TMRTopoProblem got other type'
            raise ValueError(errmsg)
        prob.initialize()
        return

    def setPrefix(self, _prefix):
        cdef char *prefix = tmr_convert_str_to_chars(_prefix)
        cdef TMRTopoProblem *prob = NULL
        prob = _dynamicTopoProblem(self.ptr)
        if prob == NULL:
            errmsg = 'Expected TMRTopoProblem got other type'
            raise ValueError(errmsg)
        prob.setPrefix(prefix)
        return

    def setIterationCounter(self, int count):
        cdef TMRTopoProblem *prob = NULL
        prob = _dynamicTopoProblem(self.ptr)
        if prob == NULL:
            errmsg = 'Expected TMRTopoProblem got other type'
            raise ValueError(errmsg)
        prob.setIterationCounter(count)
        return

    def convertPVecToVec(self, PVec pvec):
        cdef ParOptBVecWrap *new_vec = NULL
        new_vec = _dynamicParOptBVecWrap(pvec.ptr)
        if new_vec == NULL:
            errmsg = 'Expected ParOptBVecWrap got other type'
            raise ValueError(errmsg)
        return _init_Vec(new_vec.vec)

    def setInitDesignVars(self, PVec pvec, PVec lbvec=None,
                          PVec ubvec=None):
        cdef TMRTopoProblem *prob = NULL
        prob = _dynamicTopoProblem(self.ptr)
        if prob == NULL:
            errmsg = 'Expected TMRTopoProblem got other type'
            raise ValueError(errmsg)
        cdef ParOptVec *lb =NULL
        cdef ParOptVec *ub = NULL
        if lbvec is not None:
            lb = lbvec.ptr
        if ubvec is not None:
            ub = ubvec.ptr
        prob.setInitDesignVars(pvec.ptr,lb,ub)

    def setUseRecycledSolution(self, int truth):
        cdef TMRTopoProblem *prob = NULL
        prob = _dynamicTopoProblem(self.ptr)
        if prob == NULL:
            errmsg = 'Expected TMRTopoProblem got other type'
            raise ValueError(errmsg)
        prob.setUseRecycledSolution(truth)
        return
