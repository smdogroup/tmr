from utils import cantilever_egads, lbracket_egads
import argparse
from mpi4py import MPI

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--domain', type=str, nargs='*',
        default=['cantilever', 'lbracket'],
        choices=['cantilever', 'lbracket'])
    p.add_argument('--AR', type=float, nargs='*', default=[1.0])
    p.add_argument('--len0', type=float, default=1.0)
    args = p.parse_args()

    comm = MPI.COMM_WORLD

    if 'cantilever' in args.domain:
        for AR in args.AR:
            lx = args.len0*AR
            ly = args.len0
            lz = args.len0
            cantilever_egads(comm, lx, ly, lz)

    if 'lbracket' in args.domain:
        for AR in args.AR:
            lx = args.len0*AR
            ly = args.len0
            lz = args.len0
            lbracket_egads(comm, lx, ly, lz)