"""
Compute the base frequencies for different domains with concentrated mass
"""
from tmr import TMR, TopOptUtils
from tacs import TACS, constitutive

import numpy as np
from mpi4py import MPI
import argparse
import sys
import os

sys.path.append('../eigenvalue')
from utils import CreatorCallback, MFilterCreator
from utils import create_forest

class BaseFreq:

    def __init__(self, prefix, comm, domain, AR, ratio, len0, r0_frac,
                 htarget, mg_levels, qval, eig_method,
                 eig_problem, shift, estimate,
                 max_jd_size, max_gmres_size, max_lanczos,
                 non_design_mass, num_eigenvalues, index_eigvec_to_visual,
                 lanczos_spectrum, lanczos_mgmat_coeff, jd_Amat_shift):

        self.prefix = prefix
        self.eig_method = eig_method
        self.eig_problem = eig_problem
        self.ratio = ratio
        self.shift = shift
        self.estimate = estimate
        self.max_jd_size = max_jd_size
        self.max_gmres_size = max_gmres_size
        self.max_lanczos = max_lanczos
        self.num_eigenvalues = num_eigenvalues
        self.non_design_mass = non_design_mass
        self.index_eigvec_to_visual = index_eigvec_to_visual
        self.lanczos_spectrum = lanczos_spectrum
        self.lanczos_mgmat_coeff = lanczos_mgmat_coeff
        self.jd_Amat_shift = jd_Amat_shift

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

        self.Imat = self.assembler.createMat()
        self.Imat.addDiag(1.0)

        # Allocate space for eigenvalues and eigenvectors
        self.eigvals = np.zeros(self.num_eigenvalues)
        self.eigenvecs = []
        for i in range(self.num_eigenvalues):
            self.eigenvecs.append(self.assembler.createVec())

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
            if self.eig_problem == 'general':
                jd_oper = TACS.JDFrequencyOperator(self.assembler, self.kmat,
                    self.mmat, self.mg.getMat(), self.mg)
            elif self.eig_problem == 'Amat':
                jd_oper = TACS.JDSimpleOperator(self.assembler, self.Amat, self.mg)
            else:
                raise ValueError("Invalid eig_problem:", self.eig_problem)

            self.jd = TACS.JacobiDavidson(jd_oper, self.num_eigenvalues,
                self.max_jd_size, self.max_gmres_size)
            self.jd.setTolerances(eig_rtol=1e-6, eig_atol=1e-8, rtol=1e-6, atol=1e-12)
            self.jd.setRecycle(self.num_eigenvalues)

        elif self.eig_method == 'lanczos':
            self.ksm_kmat = self.assembler.createMat()
            self.ksm_Amat = self.assembler.createMat()
            ksmk = TACS.KSM(self.ksm_kmat, self.mg, 15, 0, 0)
            ksmA = TACS.KSM(self.ksm_Amat, self.mg, 15, 0, 0)

            if self.eig_problem == 'general':
                ep_oper = TACS.EPGeneralizedShiftInvertOp(self.shift, ksmk, self.mmat)
            elif self.eig_problem == 'Amat':
                ep_oper = TACS.EPShiftInvertOp(self.shift, ksmA)
            else:
                raise ValueError("Invalid eig_problem:", self.eig_problem)

            self.sep = TACS.SEPsolver(ep_oper, self.max_lanczos,
                                      TACS.SEP_FULL, self.assembler.getBcMap())
            if self.lanczos_spectrum == 'smallest':
                spectrum = TACS.SEP_SMALLEST
            elif self.lanczos_spectrum == 'smallest_magnitude':
                spectrum = TACS.SEP_SMALLEST_MAGNITUDE
            else:
                raise ValueError("invalid spectrum")
            self.sep.setTolerances(1e-6, spectrum   , self.num_eigenvalues)

        else:
            raise ValueError("Invalid eig_method:", self.eig_method)

        # Set up f5 writer
        flag = (TACS.OUTPUT_CONNECTIVITY |
                TACS.OUTPUT_NODES |
                TACS.OUTPUT_DISPLACEMENTS |
                TACS.OUTPUT_EXTRAS)
        self.f5 = TACS.ToFH5(self.assembler, TACS.SOLID_ELEMENT, flag)

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
        self.Amat.addDiag(-self.jd_Amat_shift)
        self.assembler.applyMatBCs(self.Amat)

        # Get the matrix associated with the preconditioner
        mgmat = self.mg.getMat()

        # Solve the GEP
        if self.eig_method == 'jd':

            # Assemble and factor the preconditioner
            if self.eig_problem == 'general':
                # mgmat = K
                mgmat.copyValues(self.kmat)
            elif self.eig_problem == 'Amat':
                # mgmat = K - 0.95*estimate*M
                mgmat.copyValues(self.kmat)
                mgmat.axpy(-self.estimate*0.95, self.mmat)
                mgmat.addDiag(-self.jd_Amat_shift)
            else:
                raise ValueError("invalid eig_problem")

            self.mg.assembleGalerkinMat()
            self.mg.factor()

            self.jd.solve(print_flag=True, print_level=1)

            nconverged = self.jd.getNumConvergedEigenvalues()
            if (nconverged < self.num_eigenvalues):
                if self.assembler.getMPIComm().rank == 0:
                    print("[Error] not enough eigenvalues converged! {:d}/{:d}".format(
                        nconverged, self.num_eigenvalues))

            for i in range(self.num_eigenvalues):
                self.eigvals[i], err = self.jd.extractEigenvector(i, self.eigenvecs[i])
            
            # Adjust eigenvalues and shift Amat back
            for i in range(self.num_eigenvalues):
                self.eigvals[i] += self.jd_Amat_shift
            self.Amat.addDiag(self.jd_Amat_shift)
            mgmat.addDiag(self.jd_Amat_shift)

        elif self.eig_method == 'lanczos':

            # mgmat = K - shift*M
            if self.eig_problem == 'general':
                mgmat.copyValues(self.kmat)
                mgmat.axpy(-self.shift, self.mmat)
                self.ksm_kmat.copyValues(mgmat)
            # mgmat = K - 0.95*estimate*M - shift*I
            elif self.eig_problem == 'Amat':
                mgmat.copyValues(self.kmat)
                mgmat.axpy(-self.estimate*self.lanczos_mgmat_coeff, self.mmat)
                mgmat.addDiag(-self.shift)

                self.ksm_Amat.copyValues(self.Amat)
                self.ksm_Amat.addDiag(-self.shift)
            else:
                raise ValueError("Invalid eig_problem")

            self.mg.assembleGalerkinMat()
            self.mg.factor()

            self.sep.solve(self.comm, print_flag=True)

            # print("Printing orthogonality....")
            # self.sep.printOrthogonality()
            # print("checkOrthogonality() = {:20.10e}".format(self.sep.checkOrthogonality()))
            for i in range(self.num_eigenvalues):
                self.eigvals[i], err = self.sep.extractEigenvector(i, self.eigenvecs[i])

        else:
            raise ValueError("Invalid eig_method")

        # Export to f5
        self.assembler.setVariables(self.eigenvecs[self.index_eigvec_to_visual])
        self.f5.writeToFile(os.path.join(self.prefix, "baseline.f5"))

        '''
        Check residuals
        '''

        res1 = self.assembler.createVec()
        res2 = self.assembler.createVec()
        res3 = self.assembler.createVec()
        one = self.assembler.createVec()
        one_arr = one.getArray()
        one_arr[:] = 1.0

        Kv = self.assembler.createVec()
        Mv = self.assembler.createVec()
        Av = self.assembler.createVec()

        general_res = np.zeros(self.num_eigenvalues)
        Amat_res = np.zeros(self.num_eigenvalues)
        eigvec_norm = np.zeros(self.num_eigenvalues)
        eigvec_Mnorm = np.zeros(self.num_eigenvalues)
        eigvec_sum = np.zeros(self.num_eigenvalues)

        for i in range(self.num_eigenvalues):
            self.kmat.mult(self.eigenvecs[i], Kv)
            self.mmat.mult(self.eigenvecs[i], Mv)
            self.Amat.mult(self.eigenvecs[i], Av)

            res1.copyValues(Kv)
            res1.axpy(-self.eigvals[i], Mv)

            res2.copyValues(Kv)
            res2.axpy(-self.eigvals[i], self.eigenvecs[i])

            res3.copyValues(Av)
            res3.axpy(-self.eigvals[i], self.eigenvecs[i])

            general_res[i] = res1.norm()
            Amat_res[i] = res3.norm()
            eigvec_norm[i] = self.eigenvecs[i].dot(self.eigenvecs[i])
            eigvec_Mnorm[i] = self.eigenvecs[i].dot(Mv)
            eigvec_sum[i] = self.eigenvecs[i].dot(one)

        if self.assembler.getMPIComm().rank == 0:
            print("(eigenvalues)")
            for i in range(self.num_eigenvalues):
                print("{:4d}{:20.10e}".format(i, self.eigvals[i]))

            print("\n(eigenvec) v^T*1")
            for i in range(self.num_eigenvalues):
                print("{:4d}{:20.10e}".format(i, eigvec_sum[i]))

            if self.eig_problem == 'general':
                print("\n(residual) |Kv - lambda*Mv|")
                for i in range(self.num_eigenvalues):
                    print("{:4d}{:20.10e}".format(i, general_res[i]))

            if self.eig_problem == 'Amat':
                print("\n(residual)|Av - lambda*v|")
                for i in range(self.num_eigenvalues):
                    print("{:4d}{:20.10e}".format(i, Amat_res[i]))

            # if self.eig_problem == 'Amat':
            #     print("\n(orthogonality) v^T*v")
            #     for i in range(self.num_eigenvalues):
            #         print("{:4d}{:20.10e}".format(i, eigvec_norm[i]))

            # if self.eig_problem == 'general':
            #     print("\n(orthogonality) v^TMv")
            #     for i in range(self.num_eigenvalues):
            #         print("{:4d}{:20.10e}".format(i, eigvec_Mnorm[i]))

        return self.eigvals[0], err

def run_baseline_case(prefix, domain, AR, ratio, len0, r0_frac, htarget, mg_levels,
                      qval=5.0, non_design_mass=10.0, eig_method='lanczos',
                      eig_problem='general',
                      shift=0.0, estimate=0.0, max_jd_size=100, max_gmres_size=30,
                      max_lanczos=50, num_eigenvals=5, random_design=False,
                      index_eigvec=0, lanczos_spectrum='smallest',
                      lanczos_mgmat_coeff=0.95,
                      jd_Amat_shift=0.0):

    if not os.path.isdir(prefix):
        os.mkdir(prefix)

    comm = MPI.COMM_WORLD
    bf = BaseFreq(prefix, comm, domain, AR, ratio, len0,
                  r0_frac, htarget, mg_levels,
                  qval, eig_method, eig_problem,
                  shift, estimate,
                  max_jd_size, max_gmres_size,
                  max_lanczos, non_design_mass,
                  num_eigenvals, index_eigvec,
                  lanczos_spectrum,
                  lanczos_mgmat_coeff,
                  jd_Amat_shift)
    dv = TMR.convertPVecToVec(bf.problem.createDesignVec())
    dv_vals = dv.getArray()
    dv_vals[:] = 0.95
    if random_design:
        np.random.seed(0)
        dv_vals[:] = np.random.rand()

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

    # Geometry, mesh, load
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
    p.add_argument('--eig-method', type=str, default='jd', choices=['jd', 'lanczos'])
    p.add_argument('--eig-problem', type=str, default='general', choices=['general', 'Amat'],
        help='Whether we solve Ku=lMu or Au=lu, where A=K-estimate*M')
    p.add_argument('--estimate', type=float, default=0.0, help='eigenvalue estimate for Amat')
    p.add_argument('--num-eigenvalues', type=int, default=5)
    p.add_argument('--random-design', action='store_true')
    p.add_argument('--index-eigenvec', type=int, default=0, help='which eigenvector to visualize')

    # JD options
    p.add_argument('--max-jd-size', type=int, default=100)
    p.add_argument('--max-gmres-size', type=int, default=30)
    p.add_argument('--jd-Amat-shift', type=float, default=0.0, 
        help='shift of all eigenvalues of A, shifted_eig(A) = eig(A) - shift')

    # Lanczos options
    p.add_argument('--max-lanczos', type=int, default=50)
    p.add_argument('--lanczos-shift', type=float, default=0.0, help='eigenvalue shift for lanczos')
    p.add_argument('--lanczos-spectrum', type=str, default='smallest', choices=['smallest', 'smallest_magnitude'])
    p.add_argument('--lanczos-mgmat-coeff', type=float, default=0.95)

    # OS
    p.add_argument('--prefix', type=str, default='results')
    args = p.parse_args()

    run_baseline_case(args.prefix, args.domain, args.AR, args.ratio, args.len0,
                      args.r0_frac, args.htarget, args.mg_levels,
                      args.qval, args.non_design_mass, args.eig_method,
                      args.eig_problem, args.lanczos_shift,
                      args.estimate, args.max_jd_size, args.max_gmres_size,
                      args.max_lanczos, args.num_eigenvalues, args.random_design,
                      args.index_eigenvec, args.lanczos_spectrum, args.lanczos_mgmat_coeff,
                      args.jd_Amat_shift)
