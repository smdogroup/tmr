/*
  This file is part of the package TMR for adaptive mesh refinement.

  Copyright (C) 2015 Georgia Tech Research Corporation.
  Additional copyright (C) 2015 Graeme Kennedy.
  All rights reserved.

  TMR is licensed under the Apache License, Version 2.0 (the "License");
  you may not use this software except in compliance with the License.
  You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0
  
  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
*/

#ifndef TMR_TOPO_PROBLEM_H
#define TMR_TOPO_PROBLEM_H

#include "ParOpt.h"
#include "TACSAssembler.h"
#include "TACSMg.h"
#include "TMROctForest.h"
#include "TMRQuadForest.h"
#include "TMR_RefinementTools.h"
#include "StructuralMass.h"
#include "Compliance.h"
#include "KSFailure.h"
#include "TACSBuckling.h"

/*
  Wrap a TACSBVec object with the ParOptVec interface
*/
class ParOptBVecWrap : public ParOptVec {
 public:
  ParOptBVecWrap( TACSBVec *_vec );
  ~ParOptBVecWrap();

  // Perform standard operations required for linear algebra
  // -------------------------------------------------------
  void set( ParOptScalar alpha );
  void zeroEntries();
  void copyValues( ParOptVec *pvec );
  double norm();
  double maxabs();
  double l1norm();
  ParOptScalar dot( ParOptVec *pvec );
  void mdot( ParOptVec **vecs, int nvecs, ParOptScalar *output );
  void scale( ParOptScalar alpha );
  void axpy( ParOptScalar alpha, ParOptVec *pvec );
  int getArray( ParOptScalar **array );

  // The underlying TACSBVec object
  TACSBVec *vec;
};

/*
  The implementation of the ParOptProblem class
*/
class TMRTopoProblem : public ParOptProblem {
 public:
  // Create the topology optimization object
  // ---------------------------------------
  TMRTopoProblem( int _nlevels, 
                  TACSAssembler *_tacs[],
                  TMROctForest *_filter[], 
                  TACSVarMap *_filter_maps[],
                  TACSBVecIndices *_filter_indices[],
                  TACSMg *_mg,
                  int _vars_per_node=1 );
  TMRTopoProblem( int _nlevels, 
                  TACSAssembler *_tacs[],
                  TMRQuadForest *_filter[], 
                  TACSVarMap *_filter_maps[],
                  TACSBVecIndices *_filter_indices[],
                  TACSMg *_mg,
                  int _vars_per_node=1 );
  ~TMRTopoProblem();

  // Set the load cases - note that this destroys internal information
  // stored in the load case data associated with the constraints.
  // -----------------------------------------------------------------
  void setLoadCases( TACSBVec **_forces, int _num_load_cases );
  int getNumLoadCases();

  // Add constraints associated with one of the load cases
  // -----------------------------------------------------
  void addConstraints( int _load_case, TACSFunction **_funcs,
                       const TacsScalar *_func_offset, 
                       const TacsScalar *_func_scale,
                       int num_funcs );
  void addStressConstraint( int _load_case,
                            TMRStressConstraint *stress_func,
                            TacsScalar _constr_scale=1.0 );
  void addLinearConstraints( ParOptVec **vecs,
                             TacsScalar *offset,
                             int _ncon );
  void addFrequencyConstraint( double sigma, int num_eigvals,
                               TacsScalar ks_weight,
                               TacsScalar offset, TacsScalar scale,
                               int max_subspace_size, double eigtol,
                               int use_jd=0, int fgmres_size=5,
                               double eig_rtol=1e-12,
                               double eig_atol=1e-30,
                               int num_recycle=0,
                               JDRecycleType recycle_type=JD_NUM_RECYCLE,
                               int _track_eigen_iters=0 );
  void addBucklingConstraint( double sigma, int num_eigvals,
                              TacsScalar ks_weight,
                              TacsScalar offset, TacsScalar scale,
                              int max_lanczos, double eigtol );
  // Set the objective - in this case either compliance or a function
  // for one of the load cases
  // ----------------------------------------------------------------
  void setObjective( const TacsScalar *_obj_weights, TACSFunction **funcs );
  void setObjective( const TacsScalar *_obj_weights );

  // Finish the initialization tasks - assign the number of
  // constraints and variables in the problem. Allocate arrays etc.
  // --------------------------------------------------------------
  void initialize();

  // Set the output prefix for files
  // -------------------------------
  void setPrefix( const char *prefix );

  // Set the initial design variable values
  // --------------------------------------
  void setInitDesignVars( ParOptVec *vars, ParOptVec *lb = NULL,
                          ParOptVec *ub = NULL );

  // Set the output iteration counter
  // --------------------------------
  void setIterationCounter( int iter );

  // Set the linearization point and penalty parameter
  // -------------------------------------------------
  void setLinearization( double q, ParOptVec *xvec );

  // Create a TACSBVec object containing the filter element volumes
  // --------------------------------------------------------------
  TACSBVec* createVolumeVec( double Xscale=1.0 );
  
  // Create a TACSBVec object containing the filter element area
  // --------------------------------------------------------------
  TACSBVec* createAreaVec( double Xscale=1.0 );
  
  // Create a design variable vector
  // -------------------------------
  ParOptVec *createDesignVec();

  // Set the inequality flags
  // ------------------------
  int isSparseInequality();
  int isDenseInequality();
  int useLowerBounds();
  int useUpperBounds();

  // Get the initial variables and bounds
  // ------------------------------------
  void getVarsAndBounds( ParOptVec *x, 
                         ParOptVec *lb, ParOptVec *ub );

  // Evaluate the objective and constraints
  // --------------------------------------
  int evalObjCon( ParOptVec *x, 
                  ParOptScalar *fobj, ParOptScalar *cons );

  // Evaluate the objective and constraint gradients
  // -----------------------------------------------
  int evalObjConGradient( ParOptVec *x, 
                          ParOptVec *g, ParOptVec **Ac );

  // Evaluate the product of the Hessian with the given vector px
  // ------------------------------------------------------------
  int evalHvecProduct( ParOptVec *xvec,
                       ParOptScalar *z, 
                       ParOptVec *zw,
                       ParOptVec *pxvec, 
                       ParOptVec *hvec );

  // Evaluate the sparse constraints
  // -------------------------------
  void evalSparseCon( ParOptVec *x, ParOptVec *out );

  // Compute the Jacobian-vector product out = J(x)*px
  // --------------------------------------------------
  void addSparseJacobian( double alpha, ParOptVec *x,
                          ParOptVec *px, ParOptVec *out );

  // Compute the transpose Jacobian-vector product out = J(x)^{T}*pzw
  // -----------------------------------------------------------------
  void addSparseJacobianTranspose( double alpha, ParOptVec *x,
                                   ParOptVec *pzw, 
                                   ParOptVec *out );

  // Add the inner product of the constraints to the matrix such 
  // that A += J(x)*cvec*J(x)^{T} where cvec is a diagonal matrix
  // ------------------------------------------------------------
  void addSparseInnerProduct( double alpha, ParOptVec *x,
                              ParOptVec *cvec, double *A );

  // Write the output file
  // ---------------------
  void writeOutput( int iter, ParOptVec *x );

 private:
  // Get/set values from the TACSBVec object
  int getLocalValuesFromBVec( TACSBVec *vec, TacsScalar *xloc );
  void setBVecFromLocalValues( const TacsScalar *xloc, TACSBVec *vec );

  // Extract and write the eigenvectors to file
  void writeEigenVector( int iter );

  // Store the prefix
  char *prefix;

  // Set the iteration count for printing to the file
  int iter_count;

  // Set the number of variables per node (defaults to 1)
  int vars_per_node;

  // Set the load case information. In this case, these are force
  // vectors for each load case
  int num_load_cases;
  TACSBVec **vars, **forces;

  // The linear constraints -- independent of load case
  int num_linear_con;
  ParOptVec **Alinear;
  TacsScalar *linear_offset;

  // The natural frequency TACS object
  TACSFrequencyAnalysis *freq;
  
  // Set parameters used to control the frequency constraint
  double freq_eig_tol;
  int num_freq_eigvals;
  TacsScalar freq_ks_sum;
  TacsScalar freq_ks_weight;
  TacsScalar freq_offset, freq_scale;
  int track_eigen_iters;
  KSMPrint *ksm_file;
  // The buckling TACS object
  TACSLinearBuckling **buck;
  // Set parameters used to control the buckling constraint
  double buck_eig_tol;
  int num_buck_eigvals;
  TacsScalar *buck_ks_sum;
  TacsScalar buck_ks_weight;
  TacsScalar buck_offset, buck_scale;

  // The derivative of f(x,u) w.r.t. u and the adjoint variables
  TACSBVec *dfdu, *adjoint;

  // The objective weights
  TacsScalar *obj_weights;
  TACSFunction **obj_funcs;

  // Set the constraint information for each load case
  class LoadCaseInfo {
   public:
    // Information for the TACSFunction constraints
    int num_funcs;
    TacsScalar *offset;
    TacsScalar *scale;
    TACSFunction **funcs;
    TMRStressConstraint *stress_func;
    TacsScalar stress_func_scale;
  } *load_case_info;

  // Store the design variable info
  int nlevels;
  TACSAssembler **tacs;

  // Set the information about the filter at each level
  TMROctForest **oct_filter;
  TMRQuadForest **quad_filter;
  TACSVarMap **filter_maps;
  TACSBVecIndices **filter_indices;
  TACSBVecDistribute **filter_dist;
  TACSBVecDistCtx **filter_ctx;
  TACSBVecInterp **filter_interp;

  // Create the design variable values at each level
  TACSBVec **x;

  // The initial design variable values
  int max_local_size;
  TacsScalar *xlocal;

  // Store the Krylov solver and the multigrid object 
  TACSKsm *ksm;
  TACSMg *mg;

  // The initial design variable values
  TACSBVec *xinit;
  TACSBVec *xlb, *xub;
};

#endif // TMR_TOPO_PROBLEM_H
