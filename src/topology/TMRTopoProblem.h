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

#include "ParOptProblem.h"
#include "TACSAssembler.h"
#include "TACSMg.h"
#include "TMRTopoFilter.h"
#include "TMROctForest.h"
#include "TMRQuadForest.h"
#include "TACSStructuralMass.h"
#include "TACSKSFailure.h"
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
  TMRTopoProblem( TMRTopoFilter *_filter, TACSMg *_mg,
                  int gmres_iters=50, double rtol=1e-9 );
  ~TMRTopoProblem();

  // Set the load cases - note that this destroys internal information
  // stored in the load case data associated with the constraints.
  // -----------------------------------------------------------------
  void setLoadCases( TACSBVec **_forces, int _num_load_cases );
  int getNumLoadCases();

  // Set the output frequency, element type and flags for f5 files
  // -------------------------------------------------------------
  void setF5OutputFlags( int freq, ElementType elem_type, int flag );
  void setF5EigenOutputFlags( int freq, ElementType elem_type, int flag );

  // Add constraints associated with one of the load cases
  // -----------------------------------------------------
  void addConstraints( int _load_case, TACSFunction **_funcs,
                       const TacsScalar *_func_offset,
                       const TacsScalar *_func_scale,
                       int num_funcs );
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
                               JDRecycleType recycle_type=JD_NUM_RECYCLE );
  void addBucklingConstraint( double sigma, int num_eigvals,
                              TacsScalar ks_weight,
                              TacsScalar offset, TacsScalar scale,
                              int max_lanczos, double eigtol );
  void addConstraintCallback( int ncon,
                              void *con_ptr,
                              void (*confunc)(void*, TMRTopoFilter*, TACSMg*,
                                              int, TacsScalar*),
                              void *con_grad_ptr,
                              void (*congradfunc)(void*, TMRTopoFilter*, TACSMg*,
                                                  int, TACSBVec**) );

  // Accessor functions to the underlying Assembler and Oct or QuadForest
  // --------------------------------------------------------------------
  TACSAssembler* getAssembler();
  TMRQuadForest* getFilterQuadForest();
  TMROctForest* getFilterOctForest();

  // Set the objective - in this case either compliance or a function
  // for one of the load cases
  // ----------------------------------------------------------------
  void setObjective( const TacsScalar *_obj_weights, TACSFunction **funcs );
  void setObjective( const TacsScalar *_obj_weights );
  void addObjectiveCallback( void *obj_ptr,
                             void (*objfunc)(void*, TMRTopoFilter*, TACSMg*,
                                             TacsScalar*),
                             void *obj_grad_ptr,
                             void (*objgradfunc)(void*, TMRTopoFilter*, TACSMg*,
                                                 TACSBVec*) );

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

  // Create a design variable vector
  // -------------------------------
  ParOptVec *createDesignVec();

  // Set the inequality flags
  // ------------------------
  int isSparseInequality();
  int isDenseInequality();
  int useLowerBounds();
  int useUpperBounds();

  // Set the option to use the previous solution to
  // Ku=f as the starting point for the current iteration
  // ----------------------------------------------------
  void setUseRecycledSolution( int truth );

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

  void setOutputCallback( void *data,
                          void (*func)( void*, const char*, int,
                                        TMROctForest*, TMRQuadForest*,
                                        TACSBVec* ) ){
    output_callback_ptr = data;
    writeOutputCallback = func;
  }

 private:
  void *output_callback_ptr;
  void (*writeOutputCallback)( void*, const char*, int,
                               TMROctForest*, TMRQuadForest*,
                               TACSBVec* );

  int num_callback_constraints;

  void *constraint_callback_ptr;
  void (*constraintCallback)( void*, TMRTopoFilter*, TACSMg*,
                              int, TacsScalar* );
  void *constraint_gradient_callback_ptr;
  void (*constraintGradientCallback)( void*, TMRTopoFilter*, TACSMg*,
                                      int, TACSBVec** );

  void *objective_callback_ptr;
  void (*objectiveCallback)( void*, TMRTopoFilter*, TACSMg*,
                             TacsScalar* );
  void *objective_gradient_callback_ptr;
  void (*objectiveGradientCallback)( void*, TMRTopoFilter*, TACSMg*,
                                     TACSBVec* );

  // Set the design variables across all multigrid levels
  void setDesignVars( ParOptVec *xvec );

  // Store the prefix
  char *prefix;

  // Solver parameters
  int use_recyc_sol;

  // Set the iteration count for printing to the file
  int iter_count;

  // Set the number of variables per node (defaults to 1)
  int design_vars_per_node;

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
    // TMRStressConstraint *stress_func;
    TacsScalar stress_func_offset;
    TacsScalar stress_func_scale;
  } *load_case_info;

  // Store the design variable info
  TACSAssembler *assembler;
  TMRTopoFilter *filter;

  // Store the Krylov solver and the multigrid object
  TACSKsm *ksm;
  TACSMg *mg;

  // The initial design variable values
  ParOptVec *xinit;
  ParOptVec *xlb, *xub;

  // Information to control the output frequency
  int f5_frequency, f5_eigen_frequency;
  ElementType f5_element_type, f5_eigen_element_type;
  int f5_write_flag, f5_eigen_write_flag;
};

#endif // TMR_TOPO_PROBLEM_H
