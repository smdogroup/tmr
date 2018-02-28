import numpy as np
from mpi4py import MPI
from tmr import TMR
from tacs import TACS

conn = np.array([[0, 1, 3, 4, 6, 7, 9, 10],
                 [8, 11, 2, 5, 7, 10, 1, 4]], dtype=np.intc)

comm = MPI.COMM_WORLD

fine = TMR.OctForest(comm)
fine.setConnectivity(conn)
fine.createTrees(0)

refine = np.array([2, 0], dtype=np.intc)
fine.refine(refine)
fine.balance(1)
fine.createNodes()

coarse = fine.duplicate()
coarse.refine()
coarse.balance(1)
coarse.createNodes()

coarse_range = coarse.getNodeRange()
nc = coarse_range[comm.rank+1] - coarse_range[comm.rank]
coarse_map = TACS.VarMap(comm, nc)

fine_range = fine.getNodeRange()
nf = fine_range[comm.rank+1] - fine_range[comm.rank]
fine_map = TACS.VarMap(comm, nf)

interp = TACS.VecInterp(coarse_map, fine_map, 1)
fine.createInterpolation(coarse, interp)