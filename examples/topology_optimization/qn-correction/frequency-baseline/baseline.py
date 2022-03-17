"""
Compute the base frequencies for different domains with concentrated mass
"""
from tmr import TMR, TopOptUtils
from tacs import TACS, constitutive

import numpy as np
from mpi4py import MPI
import argparse
import sys

from utils import CreatorCallback, MFilterCreator
from utils import create_forest

class BaseFreq:

    def __init__(self, comm, domain, AR, ratio, len0, r0_frac,
                 htarget, mg_levels, qval, eig_method, 
                 max_jd_size, max_gmres_size, max_lanczos,
                 mscale, kscale, num_eigenvalues):

        self.eig_method = eig_method
        self.ratio = ratio
        self.max_jd_size = max_jd_size
        self.max_gmres_size = max_gmres_size
        self.max_lanczos = max_lanczos
        self.num_eigenvalues = num_eigenvalues
        self.mscale = mscale
        self.kscale = kscale

        # Create geometry, material, bc, mesh
        self.lx = len0*AR
        self.ly = len0
        self.lz = len0
        material_props = constitutive.MaterialProperties(rho=2600.0, E=70e3, nu=0.3, ys=100.0)
        stiffness_props = TMR.StiffnessProperties(material_props, k0=1e-3, eps=0.2,
            q=qval, qmass=qval)
        forest = create_forest(comm, self.lx, self.ly, self.lz, ratio, htarget, mg_levels-1, domain)
        bcs = TMR.BoundaryConditions()
        if domain == 'mbb':
            bcs.addBoundaryCondition('symmetry', [0], [0.0])
            bcs.addBoundaryCondition('support', [1,2], [0.0, 0.0])
        else:
            bcs.addBoundaryCondition('fixed', [0,1,2], [0.0, 0.0, 0.0])

        # Topology analysis problem object
        mfilter = MFilterCreator(r0_frac, 20, a=len0)
        obj = CreatorCallback(bcs, stiffness_props)
        problem = TopOptUtils.createTopoProblem(forest,
                                                obj.creator_callback,
                                                mfilter.filter_callback,
                                                use_galerkin=True,
                                                nlevels=mg_levels)

        # Set members and allocate space for data
        self.fltr = problem.getTopoFilter()
        self.mg = problem.getMg()
        self.problem = problem
        self.assembler = self.fltr.getAssembler()
        self.comm = self.assembler.getMPIComm()
        self.kmat = self.assembler.createMat()
        self.mmat = self.assembler.createMat()
        self.m0mat = self.assembler.createMat()
        self.k0mat = self.assembler.createMat()

        # Create a non-design mass matrix
        self.mvec = self.assembler.createDesignVec()
        mvals = self.mvec.getArray()
        Xpts = forest.getPoints()
        n_local_nodes = Xpts.shape[0]
        n_ext_pre = forest.getExtPreOffset()
        offset = n_ext_pre

        tol = 1e-6
        depth = 0.1
        if domain == 'cantilever':
            xmin = (1-depth)*self.lx - tol
            xmax = self.lx + tol
            ymin = 0.25*self.ly - tol
            ymax = 0.75*self.ly + tol
            zmin = 0.0*self.lz - tol
            zmax = 0.2*self.lz + tol

        elif domain == 'michell':
            xmin = (1-depth)*self.lx - tol
            xmax = self.lx + tol
            ymin = 0.25*self.ly - tol
            ymax = 0.75*self.ly + tol
            zmin = 0.4*self.lz - tol
            zmax = 0.6*self.lz + tol

        elif domain == 'mbb':
            xmin = 0.0*self.lx - tol
            xmax = 0.2*self.lx + tol
            ymin = 0.25*self.ly - tol
            ymax = 0.75*self.ly + tol
            zmin = (1-depth)*self.lz - tol
            zmax = self.lz + tol

        elif domain == 'lbracket':
            RATIO = self.ratio
            xmin = (1-depth)*self.lx - tol
            xmax = self.lx + tol
            ymin = 0.25*self.ly - tol
            ymax = 0.75*self.ly + tol
            zmin = 0.5*RATIO*self.lz - tol
            zmax = 1.0*RATIO*self.lz + tol

        for i in range(offset, n_local_nodes):
            x, y, z = Xpts[i]
            if xmin < x < xmax:
                if ymin < y < ymax:
                    if zmin < z < zmax:
                        mvals[i-offset] = 1.0

        # Assemble a constant non-design mass matrix
        dv = self.assembler.createDesignVec()
        self.assembler.getDesignVars(dv)
        self.assembler.setDesignVars(self.mvec)
        self.assembler.assembleMatType(TACS.MASS_MATRIX, self.m0mat)
        self.assembler.assembleMatType(TACS.STIFFNESS_MATRIX, self.k0mat)
        self.m0mat.scale(self.mscale)
        self.k0mat.scale(self.kscale)
        self.assembler.setDesignVars(dv)

        # Set up solver
        if self.eig_method == 'jd':
            jd_oper = TACS.JDFrequencyOperator(self.assembler, self.kmat,
                self.mmat, self.mg.getMat(), self.mg)
            self.jd = TACS.JacobiDavidson(jd_oper, self.num_eigenvalues,
                self.max_jd_size, self.max_gmres_size)
            self.jd.setTolerances(eig_rtol=1e-6, eig_atol=1e-8, rtol=1e-6, atol=1e-12)
            self.jd.setRecycle(self.num_eigenvalues)

        elif self.eig_method == 'lanczos':
            ksmk = TACS.KSM(self.kmat, self.mg, 15, 0, 0)
            ep_oper = TACS.EPGeneralizedShiftInvertOp(0.0, ksmk, self.mmat)
            self.sep = TACS.SEPsolver(ep_oper, self.max_lanczos,
                                      TACS.SEP_FULL, self.assembler.getBcMap())
            self.sep.setTolerances(1e-6, TACS.SEP_SMALLEST, self.num_eigenvalues)

        return

    def solve(self, dv):
        """Solve for the GEP and return the smallest eigenvalue"""

        self.assembler.zeroVariables()

        # Set design variable
        self.assembler.setDesignVars(dv)

        # Assemble kmat
        self.assembler.assembleMatType(TACS.STIFFNESS_MATRIX,
                                    self.kmat)
        self.kmat.axpy(1.0, self.k0mat)
        self.assembler.applyMatBCs(self.kmat)

        # Assemble mmat
        self.assembler.assembleMatType(TACS.MASS_MATRIX,
                                    self.mmat)
        self.mmat.axpy(1.0, self.m0mat)
        self.assembler.applyMatBCs(self.mmat)

        # Solve the GEP
        if self.eig_method == 'jd':

            # Assemble and factor the preconditioner
            mgmat = self.mg.getMat()
            mgmat.copyValues(self.kmat)
            self.mg.assembleGalerkinMat()
            self.mg.factor()

            self.jd.solve(print_flag=True, print_level=1)
            assert(self.jd.getNumConvergedEigenvalues() == self.num_eigenvalues)
            eigvec = self.assembler.createVec()
            eig, err = self.jd.extractEigenvector(0, eigvec)

        elif self.eig_method == 'lanczos':

            # Assemble the multigrid preconditioner:
            mgmat = self.mg.getMat()

            # mgmat = K
            mgmat.copyValues(self.kmat)

            self.mg.assembleGalerkinMat()
            self.mg.factor()

            self.sep.solve(self.comm, print_flag=True)
            eigvec = self.assembler.createVec()
            eig, err = self.sep.extractEigenvector(0, eigvec)

        res1 = self.assembler.createVec()
        res2 = self.assembler.createVec()
        res3 = self.assembler.createVec()

        Kv = self.assembler.createVec()
        Mv = self.assembler.createVec()
        Av = self.assembler.createVec()
        self.kmat.mult(eigvec, Kv)
        self.mmat.mult(eigvec, Mv)

        res1.copyValues(Kv)
        res1.axpy(-eig, Mv)

        res2.copyValues(Kv)
        res2.axpy(-eig, eigvec)

        res3.copyValues(Av)
        res3.axpy(-eig, eigvec)

        print("(general) |Kv - lambda*Mv| = {:20.10e}".format(res1.norm()))

        # v1 = self.assembler.createVec()
        # M1 = self.assembler.createMat()
        # v1.setRand()
        # M1.addDiag(1.0) # Why this doesn't work??
        # print("(test) norm(v1)    = {:20.10e}".format(v1.norm()))
        # M1.mult(v1, v1)
        # print("(test) norm(M1*v1) = {:20.10e}".format(v1.norm()))

        return eig, err

def run_baseline_case(domain, AR, ratio, len0, r0_frac, htarget, mg_levels,
                      qval=5.0, mscale=50.0, kscale=1.0, eig_method='jd',
                      max_jd_size=200, max_gmres_size=30,
                      max_lanczos=50, num_eigenvals=5):
    comm = MPI.COMM_WORLD
    bf = BaseFreq(comm, domain, AR, ratio, len0,
                  r0_frac, htarget, mg_levels,
                  qval, eig_method, 
                  max_jd_size, max_gmres_size,
                  max_lanczos, mscale, kscale,
                  num_eigenvals)
    dv = TMR.convertPVecToVec(bf.problem.createDesignVec())
    dv_vals = dv.getArray()
    dv_vals[:] = 0.95
    eig, err = bf.solve(dv)

    if comm.rank == 0:
        print("--------- Summary ---------")
        print("domain:   {:>20s}".format(domain))
        print("AR:       {:>20.1f}".format(AR))
        print("ratio:    {:>20.2f}".format(ratio))
        print("len0:     {:>20.2f}".format(len0))
        print("r0-frac:  {:>20.2f}".format(r0_frac))
        print("htarget:  {:>20.2f}".format(htarget))
        print("mg-levels:{:>20d}".format(mg_levels))
        print("qval:     {:>20.2f}".format(qval))
        print("mscale:   {:>20.2f}".format(mscale))
        print("kscale:   {:>20.2f}".format(kscale))
        print("min eig:  {:>20.2f}".format(eig))
        print("method:   {:>20s}".format(eig_method))
        print("min eig:  {:>20e}".format(eig))
        print("error:    {:>20e}".format(err))

    ret_dict = {}
    ret_dict["domain"] = domain
    ret_dict["AR"] = AR
    ret_dict["ratio"] = ratio
    ret_dict["len0"] = len0
    ret_dict["r0-frac"] = r0_frac
    ret_dict["htarget"] = htarget
    ret_dict["mg-levels"] = mg_levels
    ret_dict["qval"] = qval
    ret_dict["mscale"] = mscale
    ret_dict["kscale"] = kscale
    ret_dict["method"] = eig_method
    ret_dict["min eig"] = eig
    ret_dict["error"] = err

    return ret_dict

if __name__ == "__main__":

    p = argparse.ArgumentParser()

    # Geometry
    p.add_argument('--domain', type=str, default='cantilever',
        choices=['cantilever', 'michell', 'mbb', 'lbracket'])
    p.add_argument('--AR', type=float, default=1.0)
    p.add_argument('--ratio', type=float, default=0.4)
    p.add_argument('--len0', type=float, default=1.0)
    p.add_argument('--r0-frac', type=float, default=0.05)
    p.add_argument('--htarget', type=float, default=1.0)
    p.add_argument('--mg-levels', type=int, default=3)
    p.add_argument('--qval', type=float, default=5.0)
    p.add_argument('--mscale', type=float, default=50.0)
    p.add_argument('--kscale', type=float, default=1.0)

    # Solver
    p.add_argument('--eig-method', type=str, default='jd',
        choices=['jd', 'lanczos'])
    p.add_argument('--max-jd-size', type=int, default=200)
    p.add_argument('--max-gmres-size', type=int, default=60)
    p.add_argument('--max-lanczos', type=int, default=50)
    p.add_argument('--num-eigenvalues', type=int, default=5)

    args = p.parse_args()

    run_baseline_case(args.domain, args.AR, args.ratio, args.len0,
                        args.r0_frac, args.htarget, args.mg_levels,
                        args.qval, args.mscale, args.kscale, args.eig_method,
                        args.max_jd_size, args.max_gmres_size,
                        args.max_lanczos, args.num_eigenvalues)
