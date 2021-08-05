"""
Compute the base frequencies for different domains with concentrated mass
"""
from tmr import TMR, TopOptUtils
from tacs import TACS, constitutive

import numpy as np
from mpi4py import MPI
import argparse
import sys

sys.path.append('../eigenvalue')
from utils import CreatorCallback, MFilterCreator
from utils import create_forest

class BaseFreq:

    def __init__(self, comm, domain, AR, ratio, len0, r0_frac,
                 htarget, mg_levels, qval, eig_method, sep_op_type,
                 eig_problem, shift, estimate,
                 max_jd_size, max_gmres_size, max_lanczos,
                 non_design_mass, num_eigenvalues):

        self.eig_method = eig_method
        self.sep_op_type = sep_op_type
        self.eig_problem = eig_problem
        self.shift = shift
        self.estimate = estimate
        self.max_jd_size = max_jd_size
        self.max_gmres_size = max_gmres_size
        self.max_lanczos = max_lanczos
        self.num_eigenvalues = num_eigenvalues
        self.non_design_mass = non_design_mass

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
        self.Amat = self.assembler.createMat()

        # Create a non-design mass matrix
        self.mvec = self.assembler.createDesignVec()
        mvals = self.mvec.getArray()
        Xpts = forest.getPoints()
        n_local_nodes = Xpts.shape[0]
        n_ext_pre = forest.getExtPreOffset()
        offset = n_ext_pre

        tol = 1e-6
        if domain == 'cantilever':
            xmin = self.lx - tol
            xmax = self.lx + tol
            ymin = 0.25*self.ly - tol
            ymax = 0.75*self.ly + tol
            zmin = 0.0*self.lz - tol
            zmax = 0.2*self.lz + tol

        elif domain == 'michell':
            xmin = self.lx - tol
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
            zmin = self.lz - tol
            zmax = self.lz + tol

        elif domain == 'lbracket':
            RATIO = self.ratio
            xmin = self.lx - tol
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
        self.m0mat.scale(self.non_design_mass)
        self.assembler.setDesignVars(dv)

        # Set up solver
        if self.eig_method == 'jd':
            if self.eig_problem == 'simple':
                jd_oper = TACS.JDSimpleOperator(self.assembler, self.kmat, self.mg)
            elif self.eig_problem == 'Amat':
                jd_oper = TACS.JDSimpleOperator(self.assembler, self.Amat, self.mg)
            elif self.eig_problem == 'general':
                jd_oper = TACS.JDFrequencyOperator(self.assembler, self.kmat,
                    self.mmat, self.mg.getMat(), self.mg)
            else:
                raise ValueError("Invalid eig_problem:", self.eig_problem)
            self.jd = TACS.JacobiDavidson(jd_oper, self.num_eigenvalues,
                self.max_jd_size, self.max_gmres_size)
            self.jd.setTolerances(eig_rtol=1e-6, eig_atol=1e-8, rtol=1e-6, atol=1e-12)
            self.jd.setRecycle(self.num_eigenvalues)

        elif self.eig_method == 'lanczos':
            ksmk = TACS.KSM(self.kmat, self.mg, 15, 0, 0)
            ksmA = TACS.KSM(self.Amat, self.mg, 15, 0, 0)
            if self.eig_problem == 'simple':
                if self.sep_op_type == 'regular':
                    ep_oper = TACS.EPRegularOp(self.kmat)
                elif self.sep_op_type == 'shift-invert':
                    ep_oper = TACS.EPShiftInvertOp(self.shift, ksmk)
                else:
                    raise ValueError("Invalid sep_op_type")
            elif self.eig_problem == 'Amat':
                if self.sep_op_type == 'regular':
                    ep_oper = TACS.EPRegularOp(self.Amat)
                elif self.sep_op_type == 'shift-invert':
                    ep_oper = TACS.EPShiftInvertOp(self.shift, ksmA)
                else:
                    raise ValueError("Invalid sep_op_type")
            elif self.eig_problem == 'general':
                if self.sep_op_type == 'general-shift-invert':
                    ep_oper = TACS.EPGeneralizedShiftInvertOp(self.shift, ksmk, self.mmat)
                else:
                    raise ValueError("Invalid sep_op_type")
            else:
                raise ValueError("Invalid eig_problem:", self.eig_problem)
            self.sep = TACS.SEPsolver(ep_oper, self.max_lanczos,
                                      TACS.SEP_FULL, self.assembler.getBcMap())
            self.sep.setTolerances(1e-6, TACS.SEP_SMALLEST, self.num_eigenvalues)

        else:
            raise ValueError("Invalid eig_method:", self.eig_method)

        return

    def solve(self, dv):
        """Solve for the GEP and return the smallest eigenvalue"""

        self.assembler.zeroVariables()

        # Set design variable
        self.assembler.setDesignVars(dv)

        # Assemble kmat
        self.assembler.assembleMatType(TACS.STIFFNESS_MATRIX,
                                    self.kmat)

        # Assemble mmat (No matter we use or not)
        self.assembler.assembleMatType(TACS.MASS_MATRIX,
                                    self.mmat)
        self.mmat.axpy(1.0, self.m0mat)
        self.assembler.applyMatBCs(self.mmat)

        # Assembler Amat (No matter we use or not)
        # A = K - estimate*M
        self.Amat.copyValues(self.kmat)
        self.Amat.axpy(-self.estimate, self.mmat)

        # Solve the GEP
        if self.eig_method == 'jd':

            # Assemble and factor the preconditioner
            if self.eig_problem == 'general' or self.eig_problem == 'simple':
                self.mg.assembleMatType(TACS.STIFFNESS_MATRIX)
            elif self.eig_problem == 'Amat':
                self.mg.assembleMatCombo(TACS.STIFFNESS_MATRIX, 1.0,
                    TACS.MASS_MATRIX, -self.estimate)
            else:
                raise ValueError("invalid eig_problem")

            self.mg.assembleMatType(TACS.STIFFNESS_MATRIX)
            self.mg.factor()

            self.jd.solve(print_flag=True, print_level=1)
            assert(self.jd.getNumConvergedEigenvalues() == self.num_eigenvalues)
            eigvec = self.assembler.createVec()
            eig, err = self.jd.extractEigenvector(0, eigvec)

        elif self.eig_method == 'lanczos':

            # Assemble the multigrid preconditioner:
            mgmat = self.mg.getMat()

            # mgmat = K - shift*M
            if self.eig_problem == 'general':
                mgmat.copyValues(self.kmat)
                mgmat.axpy(-self.shift, self.mmat)
                # self.kmat.axpy(-self.shift, self.mmat) # Is this necessary?
            # mgmat = K - shift*I or mgmat = K
            elif self.eig_problem == 'simple':
                if self.sep_op_type == 'regular':
                    mgmat.copyValues(self.kmat)
                elif self.sep_op_type == 'shift-invert':
                    mgmat.copyValues(self.kmat)
                    mgmat.addDiag(-self.shift)
                else:
                    raise ValueError("Invalid sep-op-type")
            # mgmat = A - shift*I or mgmat = A
            elif self.eig_problem == 'Amat':
                if self.sep_op_type == 'regular':
                    mgmat.copyValues(self.Amat)
                elif self.sep_op_type == 'shift-invert':
                    mgmat.copyValues(self.Amat)
                    mgmat.addDiag(-self.shift)
                else:
                    raise ValueError("Invalid sep-op-type")
            else:
                raise ValueError("Invalid eig_problem")

            self.mg.assembleGalerkinMat()
            self.mg.factor()

            self.sep.solve(self.comm, print_flag=True)
            eigvec = self.assembler.createVec()
            eig, err = self.sep.extractEigenvector(0, eigvec)

        else:
            raise ValueError("Invalid eig_method")

        res1 = self.assembler.createVec()
        res2 = self.assembler.createVec()
        res3 = self.assembler.createVec()

        Kv = self.assembler.createVec()
        Mv = self.assembler.createVec()
        Av = self.assembler.createVec()
        self.kmat.mult(eigvec, Kv)
        self.mmat.mult(eigvec, Mv)
        self.Amat.mult(eigvec, Av)

        res1.copyValues(Kv)
        res1.axpy(-eig, Mv)

        res2.copyValues(Kv)
        res2.axpy(-eig, eigvec)

        res3.copyValues(Av)
        res3.axpy(-eig, eigvec)

        print("(general) |Kv - lambda*Mv| = {:20.10e}".format(res1.norm()))
        print("(simple)  |Kv - lambda*v | = {:20.10e}".format(res2.norm()))
        print("(A-matrix)|Av - lambda*v | = {:20.10e}".format(res3.norm()))
        print("(eigenvec) v^T*v           = {:20.10e}".format(eigvec.dot(eigvec)))

        # v1 = self.assembler.createVec()
        # M1 = self.assembler.createMat()
        # v1.setRand()
        # M1.addDiag(1.0) # Why this doesn't work??
        # print("(test) norm(v1)    = {:20.10e}".format(v1.norm()))
        # M1.mult(v1, v1)
        # print("(test) norm(M1*v1) = {:20.10e}".format(v1.norm()))

        return eig, err

def run_baseline_case(domain, AR, ratio, len0, r0_frac, htarget, mg_levels,
                      qval=5.0, non_design_mass=10.0, eig_method='lanczos',
                      sep_op_type='general-shift-invert', eig_problem='general',
                      shift=0.0, estimate=0.0, max_jd_size=100, max_gmres_size=30,
                      max_lanczos=50, num_eigenvals=5):
    comm = MPI.COMM_WORLD
    bf = BaseFreq(comm, domain, AR, ratio, len0,
                  r0_frac, htarget, mg_levels,
                  qval, eig_method, sep_op_type, eig_problem,
                  shift, estimate,
                  max_jd_size, max_gmres_size,
                  max_lanczos, non_design_mass,
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
        print("ndv mass: {:>20.2f}".format(non_design_mass))
        print("")
        print("method:   {:>20s}".format(eig_method))
        print("problem:  {:>20s}".format(eig_problem))
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
    ret_dict["ndv mass"] = non_design_mass
    ret_dict["method"] = eig_method
    ret_dict["problem"] = eig_problem
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
    p.add_argument('--mg-levels', type=int, default=4)
    p.add_argument('--qval', type=float, default=5.0)
    p.add_argument('--non-design-mass', type=float, default=10.0)

    # Solver
    p.add_argument('--eig-method', type=str, default='jd',
        choices=['jd', 'lanczos'])
    p.add_argument('--sep-op-type', type=str, default='shift-invert',
        choices=['regular', 'shift-invert', 'general-shift-invert'])
    p.add_argument('--eig-problem', type=str, default='general',
        choices=['general', 'simple', 'Amat'],
        help='Whether we solve Ku=lu or Ku=lMu or Au=lu, where A=K-e*I')
    p.add_argument('--shift', type=float, default=0.0, help='eigenvalue shift for lanczos')
    p.add_argument('--estimate', type=float, default=0.0, help='eigenvalue estimate for Amat')
    p.add_argument('--max-jd-size', type=int, default=100)
    p.add_argument('--max-gmres-size', type=int, default=30)
    p.add_argument('--max-lanczos', type=int, default=50)
    p.add_argument('--num-eigenvalues', type=int, default=5)

    args = p.parse_args()

    run_baseline_case(args.domain, args.AR, args.ratio, args.len0,
                        args.r0_frac, args.htarget, args.mg_levels,
                        args.qval, args.non_design_mass, args.eig_method,
                        args.sep_op_type, args.eig_problem, args.shift,
                        args.estimate, args.max_jd_size, args.max_gmres_size,
                        args.max_lanczos, args.num_eigenvals)