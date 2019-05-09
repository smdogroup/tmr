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

#include "TMRTopoProblem.h"
#include "TMROctStiffness.h"
#include "TACSFunction.h"
#include "Solid.h"
#include "TACSToFH5.h"
#include "TMR_TACSCreator.h"


/*
  Wrap a TACSBVec object with the ParOpt vector interface
*/
ParOptBVecWrap::ParOptBVecWrap( TACSBVec *_vec ){
  vec = _vec;
  vec->incref();
}

ParOptBVecWrap::~ParOptBVecWrap(){
  vec->decref();
}

/*
  Set all the values within the vector
*/
void ParOptBVecWrap::set( ParOptScalar alpha ){
  vec->set(alpha);
}

/*
  Zero all the entries in the vector
*/
void ParOptBVecWrap::zeroEntries(){
  vec->zeroEntries();
}

/*
  Copy the vector values
*/
void ParOptBVecWrap::copyValues( ParOptVec *pvec ){
  ParOptBVecWrap *avec = dynamic_cast<ParOptBVecWrap*>(pvec);
  if (avec){
    vec->copyValues(avec->vec);
  }
}

/*
  Compute the norm
*/
double ParOptBVecWrap::norm(){
  return vec->norm();
}

/*
  Compute the maximum absolute value of any entry in the vector
*/
double ParOptBVecWrap::maxabs(){
  TacsScalar *x = NULL;
  int size = vec->getArray(&x);

  double res = 0.0;
  for ( int i = 0; i < size; i++ ){
    if (fabs(TacsRealPart(x[i])) > res){
      res = fabs(TacsRealPart(x[i]));
    }
  }

  double infty_norm = 0.0;
  MPI_Allreduce(&res, &infty_norm, 1, MPI_DOUBLE, MPI_MAX,
                vec->getMPIComm());

  return infty_norm;
}

/*
  Compute the l1 norm of the vector
*/
double ParOptBVecWrap::l1norm(){
  TacsScalar *x = NULL;
  int size = vec->getArray(&x);

  double res = 0.0;
  for ( int i = 0; i < size; i++ ){
    res += fabs(TacsRealPart(x[i]));
  }

  double l1_norm = 0.0;
  MPI_Allreduce(&res, &l1_norm, 1, MPI_DOUBLE, MPI_SUM,
                vec->getMPIComm());

  return l1_norm;
}

/*
  Compute the dot product
*/
ParOptScalar ParOptBVecWrap::dot( ParOptVec *pvec ){
  ParOptBVecWrap *avec = dynamic_cast<ParOptBVecWrap*>(pvec);
  if (avec){
    return vec->dot(avec->vec);
  }
  return 0.0;
}

/*
  Compute multiple dot products simultaneously
*/
void ParOptBVecWrap::mdot( ParOptVec **vecs, int nvecs,
                           ParOptScalar *output ){
  TACSVec **tvecs = new TACSVec*[ nvecs ];
  for ( int k = 0; k < nvecs; k++ ){
    ParOptBVecWrap *wrap = dynamic_cast<ParOptBVecWrap*>(vecs[k]);
    tvecs[k] = NULL;
    if (wrap){
      tvecs[k] = wrap->vec;
    }
  }

  vec->mdot(tvecs, output, nvecs);
  delete [] tvecs;
}

/*
  Scale the vector
*/
void ParOptBVecWrap::scale( ParOptScalar alpha ){
  vec->scale(alpha);
}

/*
  Perform an axpy operation
*/
void ParOptBVecWrap::axpy( ParOptScalar alpha, ParOptVec *pvec ){
  ParOptBVecWrap *avec = dynamic_cast<ParOptBVecWrap*>(pvec);
  if (avec){
    return vec->axpy(alpha, avec->vec);
  }
}

/*
  Get the array (return the size) of the local part of the vector
*/
int ParOptBVecWrap::getArray( ParOptScalar **array ){
  TacsScalar *_array;
  int size = 0;
  size = vec->getArray(&_array);
  if (array){
    *array = _array;
  }
  return size;
}

/*
  Create the topology optimization problem
*/
TMRTopoProblem::TMRTopoProblem( TMRTopoFilter *_filter,
                                TACSMg *_mg,
                                int gmres_iters,
                                double rtol ):
  ParOptProblem(_filter->getMPIComm()){
  // Set the prefix to NULL
  prefix = NULL;

  // Get the filter object
  filter = _filter;
  filter->incref();

  // Get the TACSAssembler object
  tacs = filter->getAssembler();
  tacs->incref();

  // The multigrid object
  mg = _mg;
  mg->incref();

  // Set the number of variables per node
  vars_per_node = filter->getVarsPerNode();

  // Set the maximum number of indices
  max_local_size = filter->getMaxNumLocalVars();

  // Set the maximum local size
  xlocal = new TacsScalar[ max_local_size ];

  // The initial design variable values (may not be set)
  xinit = NULL;
  xlb = NULL;
  xub = NULL;

  // Allocate an adjoint and df/du vector
  dfdu = tacs->createVec();
  adjoint = tacs->createVec();
  dfdu->incref();
  adjoint->incref();

  int mpi_rank;
  MPI_Comm_rank(tacs->getMPIComm(), &mpi_rank);

  // Set up the solver
  int nrestart = 5;
  int is_flexible = 0;
  atol = 1e-30;
  // double rtol = _rtol;
  use_recyc_sol = 0;
  ksm = new GMRES(mg->getMat(0), mg,
                  gmres_iters, nrestart, is_flexible);
  ksm->incref();
  ksm->setMonitor(new KSMPrintStdout("GMRES", mpi_rank, 10));
  ksm->setTolerances(rtol, atol);

  // Set the iteration count
  iter_count = 0;

  // Set the load case information
  num_load_cases = 0;
  forces = NULL;
  vars = NULL;

  // Set the load case information
  load_case_info = NULL;

  // Set the linear constraint info to NULL
  num_linear_con = 0;
  Alinear = NULL;
  linear_offset = NULL;

  // Set the objective weight information
  obj_weights = NULL;
  obj_funcs = NULL;

  // Set up the frequency constraint data
  freq = NULL;
  freq_eig_tol = 1e-8;
  num_freq_eigvals = 5;
  freq_ks_sum = 0.0;
  freq_ks_weight = 30.0;
  freq_offset = 0.0;
  freq_scale = 1.0;
  track_eigen_iters = 0;
  ksm_file = NULL;

  // Set up the buckling constraint data
  buck = NULL;
  buck_eig_tol = 1e-8;
  num_buck_eigvals = 5;
  buck_ks_sum = NULL;
  buck_ks_weight = 30.0;
  buck_offset = 0.0;
  buck_scale = 1.0;

  // Set the defaults for the output types
  f5_frequency = -1;
  f5_element_type = TACS_ELEMENT_NONE;
  f5_write_flag = 0;
  f5_eigen_frequency = -1;
  f5_eigen_element_type = TACS_ELEMENT_NONE;
  f5_eigen_write_flag = 0;
}

/*
  Free the data stored in the object
*/
TMRTopoProblem::~TMRTopoProblem(){
  if (prefix){
    delete [] prefix;
  }

  // Free the initial design variable values (if allocated)
  if (xinit){ xinit->decref(); }
  if (xlb){ xlb->decref(); }
  if (xub){ xub->decref(); }

  // Free the local temp array
  delete [] xlocal;

  // Free the solver/multigrid information
  mg->decref();
  ksm->decref();

  dfdu->decref();
  adjoint->decref();

  // Free the variables/forces
  if (forces){
    for ( int i = 0; i < num_load_cases; i++ ){
      if (forces[i]){ forces[i]->decref(); }
      vars[i]->decref();
    }
    delete [] forces;
    delete [] vars;
  }

  // Free the load case data
  for ( int i = 0; i < num_load_cases; i++ ){
    for ( int j = 0; j < load_case_info[i].num_funcs; j++ ){
      load_case_info[i].funcs[j]->decref();
    }
    if (load_case_info[i].stress_func){
      load_case_info[i].stress_func->decref();
    }
    if (load_case_info[i].funcs){
      delete [] load_case_info[i].funcs;
    }
    if (load_case_info[i].offset){
      delete [] load_case_info[i].offset;
    }
    if (load_case_info[i].scale){
      delete [] load_case_info[i].scale;
    }
  }

  // Free data for the linear constraint
  if (linear_offset){
    delete [] linear_offset;
  }
  if (Alinear){
    for ( int i = 0; i < num_linear_con; i++ ){
      if (Alinear[i]){ Alinear[i]->decref(); }
    }
  }

  // If the objective weights exist, delete them
  if (obj_weights){
    delete [] obj_weights;
  }

  // Delete the buckling and frequency TACS object
  if (freq){
    freq->decref();
  }

  // Free the array of KS weights
  if (obj_funcs){
    for ( int i = 0; i < num_load_cases; i++ ){
      if (obj_funcs[i]){ obj_funcs[i]->decref(); }
    }
  }
  if (ksm_file){
    ksm_file->decref();
  }
}

/*
  Set the output flags for any f5 output files that will be created
*/
void TMRTopoProblem::setF5OutputFlags( int freq, ElementType elem_type,
                                       unsigned int flag ){
  f5_frequency = freq;
  f5_element_type = elem_type;
  f5_write_flag = flag;
}

/*
  Set the output flags specifically for the eigenvalues
*/
void TMRTopoProblem::setF5EigenOutputFlags( int freq, ElementType elem_type,
                                            unsigned int flag ){
  f5_eigen_frequency = freq;
  f5_eigen_element_type = elem_type;
  f5_eigen_write_flag = flag;
}

/*
  Set the directory prefix to use for this load case
*/
void TMRTopoProblem::setPrefix( const char *_prefix ){
  if (prefix){ delete [] prefix; }
  prefix = new char[ strlen(_prefix)+1 ];
  strcpy(prefix, _prefix);
}

/*
  Set the load cases for each problem
*/
void TMRTopoProblem::setLoadCases( TACSBVec **_forces, int _num_load_cases ){
  // Pre-incref the input forces
  for ( int i = 0; i < _num_load_cases; i++ ){
    if (_forces[i]){
      _forces[i]->incref();
    }
  }

  // Free the forces/variables if any exist
  if (forces){
    for ( int i = 0; i < num_load_cases; i++ ){
      if (forces[i]){
        forces[i]->decref();
      }
      vars[i]->decref();
    }
    delete [] forces;
    delete [] vars;
  }

  num_load_cases = _num_load_cases;
  forces = new TACSBVec*[ num_load_cases ];
  vars = new TACSBVec*[ num_load_cases ];
  for ( int i = 0; i < num_load_cases; i++ ){
    forces[i] = _forces[i];
    vars[i] = tacs->createVec();
    vars[i]->incref();
    vars[i]->zeroEntries();
  }

  // If using the previous solution as a starting
  // point, set the atol based on rtol
  // if (use_recyc_sol){
  //   double fnorm = 0.0;
  //   for ( int i = 0; i < num_load_cases; i++ ){
  //     fnorm += forces[i]->norm();
  //   }
  //   atol = fnorm*rtol;
  //   ksm->setTolerances(rtol, atol);
  // }
  
  // Allocate the load case information
  load_case_info = new LoadCaseInfo[ num_load_cases ];
  for ( int i = 0; i < num_load_cases; i++ ){
    load_case_info[i].num_funcs = 0;
    load_case_info[i].funcs = NULL;
    load_case_info[i].offset = NULL;
    load_case_info[i].scale = NULL;
    load_case_info[i].stress_func = NULL;
    load_case_info[i].stress_func_offset = 1.0;
    load_case_info[i].stress_func_scale = 1.0;
  }
}

/*
  Get the number of load cases
*/
int TMRTopoProblem::getNumLoadCases(){
  return num_load_cases;
}

/*
  Set the constraint functions for each of the specified load cases
*/
void TMRTopoProblem::addConstraints( int load_case,
                                     TACSFunction **funcs,
                                     const TacsScalar *offset,
                                     const TacsScalar *scale,
                                     int num_funcs ){
  if (!load_case_info){
    fprintf(stderr, "TMRTopoProblem error: Must call setLoadCases() \
      before adding constraints\n");
    return;
  }
  if (load_case < 0 || load_case >= num_load_cases){
    fprintf(stderr, "TMRTopoProblem error: Load case out of range\n");
    return;
  }

  for ( int i = 0; i < num_funcs; i++ ){
    funcs[i]->incref();
  }

  // Free the load case if it has been allocated before
  if (load_case_info[load_case].num_funcs > 0){
    for ( int j = 0; j < load_case_info[load_case].num_funcs; j++ ){
      load_case_info[load_case].funcs[j]->decref();
    }
    delete [] load_case_info[load_case].funcs;
    delete [] load_case_info[load_case].offset;
    delete [] load_case_info[load_case].scale;
  }

  // Allocate the data
  load_case_info[load_case].num_funcs = num_funcs;

  if (num_funcs > 0){
    load_case_info[load_case].funcs = new TACSFunction*[ num_funcs ];
    load_case_info[load_case].offset = new TacsScalar[ num_funcs ];
    load_case_info[load_case].scale = new TacsScalar[ num_funcs ];

    // Copy over the values
    for ( int i = 0; i < num_funcs; i++ ){
      load_case_info[load_case].funcs[i] = funcs[i];
      load_case_info[load_case].offset[i] = offset[i];
      load_case_info[load_case].scale[i] = scale[i];
    }
  }
  else {
    load_case_info[load_case].funcs = NULL;
    load_case_info[load_case].stress_func = NULL;
    load_case_info[load_case].offset = NULL;
    load_case_info[load_case].scale = NULL;
  }
}

/*
  Add a stress constraint to the given load case using the TMRStressConstraint
  class
*/
void TMRTopoProblem::addStressConstraint( int load_case,
                                          TMRStressConstraint *stress_func,
                                          TacsScalar constr_offset,
                                          TacsScalar constr_scale ){
  load_case_info[load_case].stress_func = stress_func;
  load_case_info[load_case].stress_func_offset = constr_offset;
  load_case_info[load_case].stress_func_scale = constr_scale;
  load_case_info[load_case].stress_func->incref();
}

/*
  Add linear constraints to the problem.
*/
void TMRTopoProblem::addLinearConstraints( ParOptVec **vecs,
                                           TacsScalar *offset,
                                           int _ncon ){
  for ( int i = 0; i < _ncon; i++ ){
    vecs[i]->incref();
  }

  if (linear_offset){
    delete [] linear_offset;
  }
  if (Alinear){
    for ( int i = 0; i < num_linear_con; i++ ){
      Alinear[i]->decref();
    }
    delete [] Alinear;
  }

  // Allocate the new space
  num_linear_con = _ncon;
  linear_offset = new TacsScalar[ num_linear_con ];
  Alinear = new ParOptVec*[ num_linear_con ];
  for ( int i = 0; i < num_linear_con; i++ ){
    linear_offset[i] = offset[i];
    Alinear[i] = vecs[i];
  }
}

/*
  Add a natural frequency constraint
*/
void TMRTopoProblem::addFrequencyConstraint( double sigma,
                                             int num_eigvals,
                                             TacsScalar ks_weight,
                                             TacsScalar offset,
                                             TacsScalar scale,
                                             int max_subspace_size,
                                             double eigtol,
                                             int use_jd,
                                             int fgmres_size,
                                             double eig_rtol,
                                             double eig_atol,
                                             int num_recycle,
                                             JDRecycleType recycle_type,
                                             int _track_eigen_iters ){
  if (!freq){
    // Create a mass matrix for the frequency constraint
    TACSMat *mmat = tacs->createMat();
    if (use_jd){
      // Create the preconditioner matrix
      TACSMat *kmat = tacs->createMat();
      TACSMat *pcmat = mg->getMat(0);

      // Get preconditioner from Mg
      TACSPc *pc = mg;

      freq = new TACSFrequencyAnalysis(tacs, sigma, mmat,
                                       kmat, pcmat, pc,
                                       max_subspace_size, fgmres_size,
                                       num_eigvals, eigtol, eig_rtol,
                                       eig_atol, num_recycle, recycle_type);
    }
    else{
      // Create the frequency analysis object
      freq = new TACSFrequencyAnalysis(tacs, sigma, mmat,
                                       mg->getMat(0), ksm,
                                       max_subspace_size, num_eigvals,
                                       eigtol);
    }
    freq->incref();
  }

  // Set a parameters that control how the natural frequency
  // constraint is implemented
  freq_eig_tol = eigtol;
  num_freq_eigvals = num_eigvals;
  freq_ks_sum = 0.0;
  freq_ks_weight = ks_weight;
  freq_offset = offset;
  freq_scale = scale;
  track_eigen_iters = _track_eigen_iters;
  if (track_eigen_iters){
    // Get the processor rank
    int mpi_rank;
    char line[256];
    MPI_Comm_rank(tacs->getMPIComm(), &mpi_rank);
    if (use_jd){
      sprintf(line, "eigen_iteration_jd_recycle%02d_res%d.dat", num_recycle,
	      track_eigen_iters);

    }
    else {
      sprintf(line, "eigen_iteration_lanczos_res%d.dat",track_eigen_iters);
    }
    ksm_file = new KSMPrintFile(line,
                                "KSM", mpi_rank, 1);
    ksm_file->incref();
  }
}

/*
  Add a buckling constraint
*/
void TMRTopoProblem::addBucklingConstraint( double sigma,
                                            int num_eigvals,
                                            TacsScalar ks_weight,
                                            TacsScalar offset,
                                            TacsScalar scale,
                                            int max_lanczos,
                                            double eigtol ){
  if (!buck){

    // Create a geometric stiffness matrix for buckling constraint
    TACSMat *gmat = tacs->createMat();
    TACSMat *kmat = tacs->createMat();
    TACSMat *aux_mat;
    ksm->getOperators(&aux_mat, NULL);

    buck = new TACSLinearBuckling*[ num_load_cases ];
    buck_ks_sum = new TacsScalar[ num_load_cases ];
    for ( int i = 0; i < num_load_cases; i++ ){
      // Create the buckling analysis object
      buck[i] = new TACSLinearBuckling(tacs, sigma, gmat, kmat,
                                       aux_mat, ksm, max_lanczos,
                                       num_eigvals, eigtol);

      buck[i]->incref();
    }
  }

  // Set a parameters that control how the natural frequency
  // constraint is implemented
  buck_eig_tol = eigtol;
  num_buck_eigvals = num_eigvals;
  memset(buck_ks_sum, 0.0, num_load_cases*sizeof(TacsScalar));
  buck_ks_weight = ks_weight;
  buck_offset = offset;
  buck_scale = scale;
}

/*
  Set the objective weight values - this indicates a compliance objective
*/
void TMRTopoProblem::setObjective( const TacsScalar *_obj_weights ){
  if (!obj_weights){
    obj_weights = new TacsScalar[ num_load_cases ];
  }
  if (obj_funcs){
    for ( int i = 0; i < num_load_cases; i++ ){
      if (obj_funcs[i]){ obj_funcs[i]->decref(); }
    }
  }
  for ( int i = 0; i < num_load_cases; i++ ){
    obj_weights[i] = _obj_weights[i];
  }
}

/*
  Set the objective weight values - this indicates a compliance objective
*/
void TMRTopoProblem::setObjective( const TacsScalar *_obj_weights,
                                   TACSFunction **_obj_funcs ){
  if (!obj_weights){
    obj_weights = new TacsScalar[ num_load_cases ];
  }
  if (!obj_funcs){
    obj_funcs = new TACSFunction*[ num_load_cases ];
  }
  else {
    for ( int i = 0; i < num_load_cases; i++ ){
      if (obj_funcs[i]){ obj_funcs[i]->decref(); }
    }
  }
  for ( int i = 0; i < num_load_cases; i++ ){
    obj_weights[i] = _obj_weights[i];
    obj_funcs[i] = _obj_funcs[i];
    obj_funcs[i]->incref();
  }
}

/*
  Initialize the problem sizes
*/
void TMRTopoProblem::initialize(){
  // Add up the total number of constraints
  int num_constraints = num_linear_con;
  for ( int i = 0; i < num_load_cases; i++ ){
    num_constraints += load_case_info[i].num_funcs;
    if (load_case_info[i].stress_func){
      num_constraints++;
    }
  }

  if (freq){
    num_constraints++;
  }
  if (buck){
    for ( int i = 0; i < num_load_cases; i++ ){
      num_constraints++;
    }
  }

  // Set the problem sizes
  int nvars = filter->getNumLocalVars();
  int nw = 0;
  int nwblock = 0;
  if (vars_per_node > 1){
    nw = nvars/vars_per_node;
    nwblock = 1;
  }

  setProblemSizes(nvars, num_constraints, nw, nwblock);
}

/*
  Set the initial design variable values
*/
void TMRTopoProblem::setInitDesignVars( ParOptVec *xvars,
                                        ParOptVec *lb,
                                        ParOptVec *ub ){
  // Create the design vector if it doesn't exist
  if (xvars){
    if (!xinit){
      xinit = createDesignVec();
      xinit->incref();
    }

    // Copy over the local values
    xinit->copyValues(xvars);
    setDesignVars(xinit);
  }

  // Copy over the lower bound
  if (lb){
    if (!xlb){
      xlb = createDesignVec();
      xlb->incref();
    }
    xlb->copyValues(lb);
  }

  // Copy over the upper bound
  if (ub){
    if (!xub){
      xub = createDesignVec();
      xub->incref();
    }
    xub->copyValues(ub);
  }
}

/*
  Set the iteration count
*/
void TMRTopoProblem::setIterationCounter( int iter ){
  iter_count = iter;
}

/*
  Create a design variable vector
*/
ParOptVec *TMRTopoProblem::createDesignVec(){
  return new ParOptBVecWrap(filter->createVec());
}

/*
  Set whether these should be considered sparse inequalities
*/
int TMRTopoProblem::isSparseInequality(){
  // These are sparse equality constraints
  return 0;
}

/*
  Use the inequality constraint - this seems to work better
*/
int TMRTopoProblem::isDenseInequality(){
  // These are sparse inequality constraints
  return 1;
}

/*
  Always use the lower bound variables
*/
int TMRTopoProblem::useLowerBounds(){
  return 1;
}

/*
  We impose an upper bound on the variables even though in the case
  when we have reciprocal variables they are not really defeind.
*/
int TMRTopoProblem::useUpperBounds(){
  if (vars_per_node > 1){
    return 0;
  }
  return 1;
}

/*
  If True, use the solution to Ku=f from the previous iteration
  as the starting point to the current iterative solve
*/
void TMRTopoProblem::setUseRecycledSolution( int truth ){
  use_recyc_sol = truth;
}

/*
  Set the initial design variables
*/
void TMRTopoProblem::getVarsAndBounds( ParOptVec *xvec,
                                       ParOptVec *lbvec,
                                       ParOptVec *ubvec ){
  int mpi_rank;
  MPI_Comm_rank(tacs->getMPIComm(), &mpi_rank);
  if (xvec){
    // If the initial design vector exists, copy it directly
    if (xinit){
      xvec->copyValues(xinit);
    }
    else {
      ParOptBVecWrap *wrap = dynamic_cast<ParOptBVecWrap*>(xvec);
      if (wrap){
        // Get the design variable values from TACS directly
        memset(xlocal, 0, max_local_size*sizeof(TacsScalar));
        tacs->getDesignVars(xlocal, max_local_size);
        filter->setValues(xlocal, wrap->vec);
      }
    }
  }

  // Handle the lower and upper bounds
  if (lbvec || ubvec){
    int has_lower = 0, has_upper = 0;
    if (lbvec && xlb){
      has_lower = 1;
      lbvec->copyValues(xlb);
    }
    if (ubvec && xub){
      has_upper = 1;
      ubvec->copyValues(xub);
    }

    if (!(has_lower && has_upper)){
      // Get the design variable values from TACS
      TacsScalar *upper = new TacsScalar[ max_local_size ];
      memset(xlocal, 0, max_local_size*sizeof(TacsScalar));
      memset(upper, 0, max_local_size*sizeof(TacsScalar));
      tacs->getDesignVarRange(xlocal, upper, max_local_size);

      if (!has_lower){
        ParOptBVecWrap *lbwrap = dynamic_cast<ParOptBVecWrap*>(lbvec);
        if (lbwrap){
          lbwrap->vec->zeroEntries();
          filter->setValues(xlocal, lbwrap->vec);
        }
      }
      if (!has_upper){
        ParOptBVecWrap *ubwrap = dynamic_cast<ParOptBVecWrap*>(ubvec);
        if (ubwrap){
          ubwrap->vec->zeroEntries();
          filter->setValues(upper, ubwrap->vec);
        }
      }
      delete [] upper;
    }
  }
}

/*
  Set the design variable values consistently across all multigrid
  levels by interpolating the design variables using the filter
  interpolation objects.
*/
void TMRTopoProblem::setDesignVars( ParOptVec *pxvec ){
  ParOptBVecWrap *wrap = dynamic_cast<ParOptBVecWrap*>(pxvec);

  if (wrap){
    // Set the design variable values
    filter->setDesignVars(wrap->vec);
  }
}

/*
  Evaluate the objective and constraints
*/
int TMRTopoProblem::evalObjCon( ParOptVec *pxvec,
                                ParOptScalar *fobj,
                                ParOptScalar *cons ){
  // Get the rank of comm
  int mpi_rank;
  MPI_Comm_rank(tacs->getMPIComm(), &mpi_rank);

  // Set the design variable values on all mesh levels
  setDesignVars(pxvec);

  // Zero the variables
  tacs->zeroVariables();

  // Assemble the Jacobian on each level
  double alpha = 1.0, beta = 0.0, gamma = 0.0;
  mg->assembleJacobian(alpha, beta, gamma, NULL);
  mg->factor();

  // Set the objective value
  *fobj = 0.0;

  // Keep track of the constraint number
  int count = 0;

  // Compute the linear constraint
  for ( int i = 0; i < num_linear_con; i++, count++ ){
    cons[count] = linear_offset[i] + Alinear[i]->dot(pxvec);
  }

  for ( int i = 0; i < num_load_cases; i++ ){
    if (forces[i]){
      // Solve the system: K(x)*u = forces
      if (use_recyc_sol){
        ksm->solve(forces[i], vars[i], 0);
      }
      else {
	ksm->solve(forces[i], vars[i]);
      }
      tacs->setBCs(vars[i]);

      // Set the variables into TACSAssembler
      tacs->setVariables(vars[i]);

      // Add the contribution to the objective
      if (obj_funcs){
        if (obj_funcs[i]){
          TacsScalar fobj_val;
          tacs->evalFunctions(&obj_funcs[i], 1, &fobj_val);
          *fobj += obj_weights[i]*fobj_val;
        }
      }
      // If we are dealing with a compliance objective
      else {
        *fobj += obj_weights[i]*vars[i]->dot(forces[i]);
      }

      // Evaluate the constraints
      int num_funcs = load_case_info[i].num_funcs;
      if (num_funcs > 0){
        tacs->evalFunctions(load_case_info[i].funcs,
                            num_funcs, &cons[count]);

        if (mpi_rank == 0){
          for ( int k = 0; k < num_funcs; k++ ){
            printf("%-30s %25.10e\n",
                   load_case_info[i].funcs[k]->functionName(),
                   cons[count + k]);
            TACSKSFailure *ks_fail =
              dynamic_cast<TACSKSFailure*>(load_case_info[i].funcs[k]);
            if (ks_fail){
              TacsScalar max_fail = ks_fail->getMaximumFailure();
              printf("%-30s %25.10e\n", "TACSKSFailure max stress",
                     max_fail);
            }
          }
        }

        // Scale and offset the constraints that we just evaluated
        for ( int j = 0; j < num_funcs; j++ ){
          TacsScalar offset = load_case_info[i].offset[j];
          TacsScalar scale = load_case_info[i].scale[j];
          cons[count + j] = scale*(cons[count+j] + offset);
        }
        count += num_funcs;
      }

      // Evaluate the stress constraint
      if (load_case_info[i].stress_func){
        TacsScalar con_offset = load_case_info[i].stress_func_offset;

        cons[count] =
          con_offset - load_case_info[i].stress_func->evalConstraint(vars[i]);
        cons[count] *= load_case_info[i].stress_func_scale;
        count++;
      }
    }
  }

  // Compute the natural frequency constraint, if any
  if (freq){
    // Keep track of the number of eigenvalues with unacceptable
    // error. If more than one exists, re-solve the eigenvalue
    // problem again.
    int err_count = 1;

    // Keep track of the smallest eigenvalue
    double shift = 0.95;
    TacsScalar smallest_eigval = 0.0;
    while (err_count > 0){
      // Set the error counter to zero
      err_count = 0;
      // Solve the eigenvalue problem
      freq->solve(new KSMPrintStdout("KSM", mpi_rank, 1),
                  ksm_file);

      // Extract the first k eigenvalues
      for ( int k = 0; k < num_freq_eigvals; k++ ){
        TacsScalar error;
        TacsScalar eigval = freq->extractEigenvalue(k, &error);
        if (eigval < 0.0){
          eigval *= -1.0;
        }

        if (k == 0){
          smallest_eigval = eigval;
        }
        if (eigval < smallest_eigval){
          smallest_eigval = eigval;
        }
        if (error > freq_eig_tol){
          err_count++;
        }
      }
      err_count = 0;
      // If there is significant error in computing the eigenvalues,
      // reset the buckling computation
      if (err_count > 0){
        double sigma = shift*smallest_eigval;
        shift += 0.5*(1.0 - shift) + shift;
        freq->setSigma(sigma);
      }
    }

    // Evaluate the KS function of the lowest eigenvalues
    freq_ks_sum = 0.0;

    // Weight on the KS function
    for (int k = 0; k < num_freq_eigvals; k++){
      TacsScalar error;
      TacsScalar eigval = freq->extractEigenvalue(k, &error);
      if (eigval < 0.0){
        eigval *= -1.0;
      }

      // Add up the contribution to the ks function
      freq_ks_sum += exp(-freq_ks_weight*(eigval - smallest_eigval));
    }

    // Evaluate the KS function of the aggregation of the eigenvalues
    cons[count] = (smallest_eigval - log(freq_ks_sum)/freq_ks_weight);
    cons[count] = freq_scale*(cons[count] + freq_offset);
    count++;
  }
  // Compute the buckling constraint, if any
  if (buck){
    for ( int i = 0; i < num_load_cases; i++ ){
      if (forces[i]){
        // Keep track of the number of eigenvalues with unacceptable
        // error. If more than one exists, re-solve the eigenvalue
        // problem again.
        int err_count = 1;

        // Keep track of the smallest eigenvalue
        double shift = 0.95;
        TacsScalar smallest_eigval = 0.0;

        while (err_count > 0){
          // Set the error counter to zero
          err_count = 0;

          // Solve the eigenvalue problem
          buck[i]->solve(forces[i],
                         new KSMPrintStdout("KSM", mpi_rank, 1));

          // Extract the first k eigenvalues
          for ( int k = 0; k < num_buck_eigvals; k++ ){
            TacsScalar error;
            TacsScalar eigval = buck[i]->extractEigenvalue(k, &error);
            if (eigval < 0.0){
              eigval *= -1.0;
            }

            if (k == 0){
              smallest_eigval = eigval;
            }
            if (error > buck_eig_tol){
              err_count++;
            }
          }

          // If there is significant error in computing the eigenvalues,
          // reset the buckling computation
          if (err_count > 0){
            double sigma = shift*smallest_eigval;
            shift += 0.5*(1.0 - shift) + shift;
            buck[i]->setSigma(sigma);
          }
        }

        // Evaluate the KS function of the lowest eigenvalues
        buck_ks_sum[i] = 0.0;

        // Weight on the KS function
        for (int k = 0; k < num_buck_eigvals; k++){
          TacsScalar error;
          TacsScalar eigval = buck[i]->extractEigenvalue(k, &error);
          if (eigval < 0.0){
            eigval *= -1.0;
          }

          // Add up the contribution to the ks function
          buck_ks_sum[i] += exp(-buck_ks_weight*(eigval - smallest_eigval));
        }

        // Evaluate the KS function of the aggregation of the eigenvalues
        cons[count] = (smallest_eigval - log(buck_ks_sum[i])/buck_ks_weight);
        cons[count] = buck_scale*(cons[count] + buck_offset);
        count++;
      }
    }
  }

  return 0;
}

/*
  Evaluate the objective and constraint gradients
*/
int TMRTopoProblem::evalObjConGradient( ParOptVec *xvec,
                                        ParOptVec *gvec,
                                        ParOptVec **Acvec ){
  int mpi_rank;
  MPI_Comm_rank(tacs->getMPIComm(), &mpi_rank);

  // Evaluate the derivative of the weighted compliance with
  // respect to the design variables
  ParOptBVecWrap *wrap = dynamic_cast<ParOptBVecWrap*>(gvec);
  if (wrap){
    TACSBVec *g = wrap->vec;
    g->zeroEntries();

    // Evaluate the gradient of the objective - weighted sum of
    // compliances
    memset(xlocal, 0, max_local_size*sizeof(TacsScalar));

    if (obj_funcs){
      for ( int i = 0; i < num_load_cases; i++ ){
        tacs->setVariables(vars[i]);
        int use_adjoint = 1;
        if (dynamic_cast<TACSStructuralMass*>(obj_funcs[i])){
          use_adjoint = 0;
        }

        if (use_adjoint){
          dfdu->zeroEntries();
          double alpha = 1.0, beta = 0.0, gamma = 0.0;
          mg->assembleJacobian(alpha, beta, gamma, NULL, TRANSPOSE);
          mg->factor();
          tacs->addSVSens(alpha, beta, gamma, &obj_funcs[i], 1, &dfdu);
          tacs->applyBCs(dfdu);

          // Solve the system of adjoint equations
          ksm->solve(dfdu, adjoint);
          tacs->addDVSens(obj_weights[i], &obj_funcs[i],
                          1, xlocal, max_local_size);
          tacs->addAdjointResProducts(-obj_weights[i], &adjoint,
                                      1, xlocal, max_local_size);
        }
        else {
          tacs->addDVSens(obj_weights[i], &obj_funcs[i], 1, xlocal,
                          max_local_size);
        }
      }

      filter->addValues(xlocal, g);
    }
    else { // For compliance objective
      for ( int i = 0; i < num_load_cases; i++ ){
        tacs->setVariables(vars[i]);
        tacs->addAdjointResProducts(-obj_weights[i], &vars[i],
                                    1, xlocal, max_local_size);
      }

      filter->addValues(xlocal, g);
    }
  }
  else {
    return 1;
  }

  // Keep track of the constraint gradient number
  int count = 0;

  // Set the linear constraint
  for ( int i = 0; i < num_linear_con; i++, count++ ){
    Acvec[count]->copyValues(Alinear[i]);
  }

  // Compute the derivative of the constraint functions
  for ( int i = 0; i < num_load_cases; i++ ){
    tacs->setVariables(vars[i]);

    // Get the number of functions for each load case
    int num_funcs = load_case_info[i].num_funcs;

    for ( int j = 0; j < num_funcs; j++ ){
      TACSFunction *func = load_case_info[i].funcs[j];
      TacsScalar scale = load_case_info[i].scale[j];

      // Try to unwrap the vector
      wrap = dynamic_cast<ParOptBVecWrap*>(Acvec[count + j]);

      if (wrap){
        TACSBVec *A = wrap->vec;
        A->zeroEntries();

        // If the function is the structural mass, then do not
        // use the adjoint, otherwise assume that we should use the
        // adjoint method to compute the gradient.
        int use_adjoint = 1;
        if (dynamic_cast<TACSStructuralMass*>(func)){
          use_adjoint = 0;
        }
        if (use_adjoint){
          // Evaluate the right-hand-side
          dfdu->zeroEntries();
          double alpha = 1.0, beta = 0.0, gamma = 0.0;
          tacs->addSVSens(alpha, beta, gamma, &func, 1, &dfdu);
          tacs->applyBCs(dfdu);

          // Solve the system of equations
          ksm->solve(dfdu, adjoint);

          // Compute the total derivative using the adjoint
          memset(xlocal, 0, max_local_size*sizeof(TacsScalar));
          tacs->addDVSens(scale, &func, 1, xlocal, max_local_size);
          tacs->addAdjointResProducts(-scale, &adjoint,
                                      1, xlocal, max_local_size);
        }
        else {
          memset(xlocal, 0, max_local_size*sizeof(TacsScalar));
          tacs->addDVSens(scale, &func, 1, xlocal, max_local_size);
        }

        filter->addValues(xlocal, A);
      }
    }
    count += num_funcs;

    // Compute the gradient with respect to the stress-reconstruction
    // ks functional
    if (load_case_info[i].stress_func){
      // Try to unwrap the vector
      wrap = dynamic_cast<ParOptBVecWrap*>(Acvec[count]);

      if (wrap){
        // Get the underlying TACS vector for the design variables
        TACSBVec *A = wrap->vec;

        // Evaluate the partial derivatives required for the adjoint
        load_case_info[i].stress_func->evalConDeriv(xlocal,
                                                    max_local_size,
                                                    dfdu);
        tacs->applyBCs(dfdu);

        // Solve the system of equations
        ksm->solve(dfdu, adjoint);

        // Compute the total derivative using the adjoint
        tacs->addAdjointResProducts(-1.0, &adjoint,
                                    1, xlocal, max_local_size);

        // Add the local values to obtain the filtered sensitivity
        filter->addValues(xlocal, A);

        // Scale the constraint by -1 since the constraint is
        // formulated as 1 - c(x, u) > 0.0
        A->scale(-load_case_info[i].stress_func_scale);
      }

      count++;
    }

    if (freq && i == 0){
      // Try to unwrap the vector
      wrap = dynamic_cast<ParOptBVecWrap*>(Acvec[count]);
      if (wrap){
        // Get the underlying TACS vector for the design variables
        TACSBVec *A = wrap->vec;

        // Zero the local components
        memset(xlocal, 0, max_local_size*sizeof(TacsScalar));

        TacsScalar *temp = new TacsScalar[ max_local_size ];
        memset(temp, 0, max_local_size*sizeof(TacsScalar));

        // Add the contribution from each eigenvalue derivative
        TacsScalar smallest_eigval = 0.0;

        // Extract the first k eigenvalues to find smallest eigval
        for ( int k = 0; k < num_freq_eigvals; k++ ){
          TacsScalar error;
          TacsScalar eigval = freq->extractEigenvalue(k, &error);
          if (eigval < 0.0){
            eigval *= -1.0;
          }

          if (k == 0){
            smallest_eigval = eigval;
          }
          if (eigval < smallest_eigval){
            smallest_eigval = eigval;
          }
        }

        for ( int k = 0; k < num_freq_eigvals; k++ ){
          // Compute the derivaive of the eigenvalue
          freq->evalEigenDVSens(k, temp, max_local_size);

          // Extract the eigenvalue itself
          TacsScalar error;
          TacsScalar ks_grad_weight = 1.0;
          TacsScalar eigval = freq->extractEigenvalue(k, &error);
          if (eigval < 0.0){
            eigval *= -1.0;
            ks_grad_weight = -1.0;
          }

          // Evaluate the weight on the gradient
          ks_grad_weight *=
            exp(-freq_ks_weight*(eigval - smallest_eigval))/freq_ks_sum;

          // Scale the constraint by the frequency scaling value
          ks_grad_weight *= freq_scale;

          // Add contribution to eigenvalue gradient
          for ( int j = 0; j < max_local_size; j++ ){
            xlocal[j] += ks_grad_weight*temp[j];
          }
        }

        // Free the data
        delete [] temp;

        filter->addValues(xlocal, A);
      }

      count++;
    }
    if (buck){
      // Compute the derivative of the buckling functions
      // Try to unwrap the vector
      wrap = dynamic_cast<ParOptBVecWrap*>(Acvec[count]);
      if (wrap){
        // Get the vector
        TACSBVec *A = wrap->vec;

        // Zero the local components
        memset(xlocal, 0, max_local_size*sizeof(TacsScalar));

        TacsScalar *temp = new TacsScalar[max_local_size];
        memset(temp, 0, max_local_size*sizeof(TacsScalar));

        // Add the contribution from each eigenvalue derivative
        TacsScalar smallest_eigval = 0.0;
        for ( int k = 0; k < num_buck_eigvals; k++ ){
          // Compute the derivaive of the eigenvalue
          buck[i]->evalEigenDVSens(k, temp, max_local_size);

          // Extract the eigenvalue itself
          TacsScalar error;
          TacsScalar ks_grad_weight = 1.0;
          TacsScalar eigval = buck[i]->extractEigenvalue(k, &error);
          if (eigval < 0.0){
            eigval *= -1.0;
            ks_grad_weight = -1.0;
          }
          if (k == 0){
            smallest_eigval = eigval;
          }

          // Evaluate the weight on the gradient
          ks_grad_weight *=
            exp(-buck_ks_weight*(eigval - smallest_eigval))/buck_ks_sum[i];

          // Scale the constraint by the buckling scaling value
          ks_grad_weight *= buck_scale;

          // Add contribution to eigenvalue gradient
          for ( int j = 0; j < max_local_size; j++ ){
            xlocal[j] += ks_grad_weight*temp[j];
          }
        }

        // Free the data
        delete [] temp;

        filter->addValues(xlocal, A);
      }
      count++;
    }
  } // end num_load_cases

  return 0;
}

/*
  Evaluate the product of the Hessian with the given vector px
*/
int TMRTopoProblem::evalHvecProduct( ParOptVec *xvec,
                                     ParOptScalar *z,
                                     ParOptVec *zw,
                                     ParOptVec *pxvec,
                                     ParOptVec *hvec ){
  return 0;
}

// Evaluate the sparse constraints
// ------------------------
void TMRTopoProblem::evalSparseCon( ParOptVec *xvec,
                                    ParOptVec *outvec ){
  if (vars_per_node > 1){
    // Dynamically cast the vectors to the ParOptBVecWrap class
    TacsScalar *x, *out;
    int size = xvec->getArray(&x);
    outvec->getArray(&out);

    // Compute the weighting constraints
    int n = size/vars_per_node;
    for ( int i = 0; i < n; i++ ){
      out[i] = -1.0;
      for ( int j = 0; j < vars_per_node; j++ ){
        out[i] += x[vars_per_node*i + j];
      }
    }
  }
}

// Compute the Jacobian-vector product out = J(x)*px
// --------------------------------------------------
void TMRTopoProblem::addSparseJacobian( double alpha,
                                        ParOptVec *xvec,
                                        ParOptVec *pxvec,
                                        ParOptVec *outvec ){
  if (vars_per_node > 1){
    TacsScalar *px, *out;
    int size = pxvec->getArray(&px);
    outvec->getArray(&out);

    // Compute the matrix-vector product
    int n = size/vars_per_node;
    for ( int i = 0; i < n; i++ ){
      for ( int j = 0; j < vars_per_node; j++ ){
        out[i] += alpha*px[vars_per_node*i + j];
      }
    }
  }
}

// Compute the transpose Jacobian-vector product out = J(x)^{T}*pzw
// -----------------------------------------------------------------
void TMRTopoProblem::addSparseJacobianTranspose( double alpha,
                                                 ParOptVec *x,
                                                 ParOptVec *pzwvec,
                                                 ParOptVec *outvec ){
  if (vars_per_node > 1){
    TacsScalar *pzw, *out;
    pzwvec->getArray(&pzw);
    int size = outvec->getArray(&out);

    int n = size/vars_per_node;
    for ( int i = 0; i < n; i++ ){
      for ( int j = 0; j < vars_per_node; j++ ){
        out[vars_per_node*i + j] += alpha*pzw[i];
      }
    }
  }
}

// Add the inner product of the constraints to the matrix such
// that A += J(x)*cvec*J(x)^{T} where cvec is a diagonal matrix
// ------------------------------------------------------------
void TMRTopoProblem::addSparseInnerProduct( double alpha,
                                            ParOptVec *x,
                                            ParOptVec *cvec,
                                            double *A ){
  if (vars_per_node > 1){
    TacsScalar *c;
    int size = cvec->getArray(&c);

    int n = size/vars_per_node;
    for ( int i = 0; i < n; i++ ){
      for ( int j = 0; j < vars_per_node; j++ ){
        A[i] += alpha*c[vars_per_node*i + j];
      }
    }
  }
}

/*
  Write the output file
*/
void TMRTopoProblem::writeOutput( int iter, ParOptVec *xvec ){
  // Print out the binary STL file for later visualization
  if (prefix){
    // Write out the file at a cut off of 0.25
    char *filename = new char[ strlen(prefix) + 100 ];

    for ( int k = 0; k < filter->getVarsPerNode(); k++ ){
      double cutoff = 0.5;
      sprintf(filename, "%s/levelset05_var%d_binary%04d.bstl",
              prefix, k, iter_count);

      // Write the STL file
      filter->writeSTLFile(k, cutoff, filename);
    }

    delete [] filename;
  }

  if (prefix){
    if (f5_frequency > 0 && (iter % f5_frequency) == 0){
      // Create the filename
      char *filename = new char[ strlen(prefix) + 100 ];
      sprintf(filename, "%s/tacs_output%04d.f5", prefix, iter_count);

      TACSToFH5 *f5 = new TACSToFH5(tacs, f5_element_type, f5_write_flag);
      f5->incref();
      f5->writeToFile(filename);
      f5->decref();
      delete [] filename;
    }

    if ((buck || freq) && f5_eigen_frequency > 0 &&
        (iter % f5_eigen_frequency) == 0){
      char *filename = new char[ strlen(prefix) + 100 ];
      TACSToFH5 *f5 = new TACSToFH5(tacs, f5_eigen_element_type,
                                    f5_eigen_write_flag);
      f5->incref();

      TACSBVec *tmp = tacs->createVec();
      tmp->incref();
      if (buck){
        for ( int i = 0; i < num_load_cases; i++ ){
          // Extract the first k eigenvectors for ith load
          for ( int k = 0; k < num_buck_eigvals; k++ ){
            TacsScalar error;
            buck[i]->extractEigenvector(k, tmp, &error);
            tacs->setVariables(tmp);
            sprintf(filename, "%s/load%d_eigenvector%02d_output%d.f5",
                    prefix, i, k, iter);
            f5->writeToFile(filename);
          }
        }
      }
      else {
        for ( int k = 0; k < num_freq_eigvals; k++ ){
          TacsScalar error;
          freq->extractEigenvector(k, tmp, &error);
          tacs->setVariables(tmp);
          sprintf(filename, "%s/freq_eigenvector%02d_output%d.f5",
                  prefix, k, iter);
          f5->writeToFile(filename);
        }
      }
      tmp->decref();

      f5->decref();
      delete [] filename;
    }
  }

  // Update the iteration count
  iter_count++;
}
