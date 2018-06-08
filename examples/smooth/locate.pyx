#!/bin/python

import numpy as np
cimport numpy as np

# Ensure that numpy is initialized
np.import_array()

cdef extern from "LocatePoint.h":
   cppclass LocatePoint:
      LocatePoint(double *, int, int)
      int locateClosest(double*)
      void locateKClosest(int, int*, double*, double*)

cdef class locate:
   cdef LocatePoint *ptr
     
   def __cinit__(self, np.ndarray[double, ndim=2, mode='c'] xpts):
       cdef int nnodes = xpts.shape[0]
       cdef int max_size = 20
       self.ptr = new LocatePoint(<double*>xpts.data, nnodes, max_size)

   def __dealloc__(self):
      del self.ptr

   def locateClosest(self, np.ndarray[double, ndim=1] x):
       return self.ptr.locateClosest(<double*>x.data)

   def locateKClosest(self, np.ndarray[double, ndim=1, mode='c'] x,
                      np.ndarray[int, ndim=1] index,
                      np.ndarray[double, ndim=1] dist):
       cdef int K = index.shape[0]
       self.ptr.locateKClosest(K, <int*>index.data, 
                               <double*>dist.data,
                               <double*>x.data)
       return
