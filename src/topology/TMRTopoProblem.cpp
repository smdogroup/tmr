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

#include "TACSFunction.h"
#include "TACSToFH5.h"
#include "TMR_TACSCreator.h"

/*
  Wrap a TACSBVec object with the ParOpt vector interface
*/
ParOptBVecWrap::ParOptBVecWrap(TACSBVec *_vec) {
  vec = _vec;
  vec->incref();
}

ParOptBVecWrap::~ParOptBVecWrap() { vec->decref(); }

/*
  Set all the values within the vector
*/
void ParOptBVecWrap::set(ParOptScalar alpha) { vec->set(alpha); }

/*
  Zero all the entries in the vector
*/
void ParOptBVecWrap::zeroEntries() { vec->zeroEntries(); }

/*
  Copy the vector values
*/
void ParOptBVecWrap::copyValues(ParOptVec *pvec) {
  ParOptBVecWrap *avec = dynamic_cast<ParOptBVecWrap *>(pvec);
  if (avec) {
    vec->copyValues(avec->vec);
  }
}

/*
  Compute the norm
*/
double ParOptBVecWrap::norm() { return vec->norm(); }

/*
  Compute the maximum absolute value of any entry in the vector
*/
double ParOptBVecWrap::maxabs() {
  TacsScalar *x = NULL;
  int size = vec->getArray(&x);

  double res = 0.0;
  for (int i = 0; i < size; i++) {
    if (fabs(TacsRealPart(x[i])) > res) {
      res = fabs(TacsRealPart(x[i]));
    }
  }

  double infty_norm = 0.0;
  MPI_Allreduce(&res, &infty_norm, 1, MPI_DOUBLE, MPI_MAX, vec->getMPIComm());

  return infty_norm;
}

/*
  Compute the l1 norm of the vector
*/
double ParOptBVecWrap::l1norm() {
  TacsScalar *x = NULL;
  int size = vec->getArray(&x);

  double res = 0.0;
  for (int i = 0; i < size; i++) {
    res += fabs(TacsRealPart(x[i]));
  }

  double l1_norm = 0.0;
  MPI_Allreduce(&res, &l1_norm, 1, MPI_DOUBLE, MPI_SUM, vec->getMPIComm());

  return l1_norm;
}

/*
  Compute the dot product
*/
ParOptScalar ParOptBVecWrap::dot(ParOptVec *pvec) {
  ParOptBVecWrap *avec = dynamic_cast<ParOptBVecWrap *>(pvec);
  if (avec) {
    return vec->dot(avec->vec);
  }
  return 0.0;
}

/*
  Compute multiple dot products simultaneously
*/
void ParOptBVecWrap::mdot(ParOptVec **vecs, int nvecs, ParOptScalar *output) {
  TACSVec **tvecs = new TACSVec *[nvecs];
  for (int k = 0; k < nvecs; k++) {
    ParOptBVecWrap *wrap = dynamic_cast<ParOptBVecWrap *>(vecs[k]);
    tvecs[k] = NULL;
    if (wrap) {
      tvecs[k] = wrap->vec;
    }
  }

  vec->mdot(tvecs, output, nvecs);
  delete[] tvecs;
}

/*
  Scale the vector
*/
void ParOptBVecWrap::scale(ParOptScalar alpha) { vec->scale(alpha); }

/*
  Perform an axpy operation
*/
void ParOptBVecWrap::axpy(ParOptScalar alpha, ParOptVec *pvec) {
  ParOptBVecWrap *avec = dynamic_cast<ParOptBVecWrap *>(pvec);
  if (avec) {
    return vec->axpy(alpha, avec->vec);
  }
}

/*
  Get the array (return the size) of the local part of the vector
*/
int ParOptBVecWrap::getArray(ParOptScalar **array) {
  TacsScalar *_array;
  int size = 0;
  size = vec->getArray(&_array);
  if (array) {
    *array = _array;
  }
  return size;
}

/*
  Create the topology optimization problem
*/
TMRTopoProblem::TMRTopoProblem(TMRTopoFilter *_filter, TACSMg *_mg,
                               int gmres_iters, double rtol)
    : ParOptProblem(_filter->getAssembler()->getMPIComm()) {
  // Set the prefix to NULL
  prefix = NULL;

  // Get the filter object
  filter = _filter;
  filter->incref();

  // Get the TACSAssembler object
  assembler = filter->getAssembler();
  assembler->incref();

  // The multigrid object
  mg = _mg;
  if (mg) {
    mg->incref();
  }

  // Set the number of variables per node
  design_vars_per_node = assembler->getDesignVarsPerNode();

  // The initial design variable values (may not be set)
  xinit = NULL;
  xlb = NULL;
  xub = NULL;

  // Allocate an adjoint and df/du vector
  dfdu = assembler->createVec();
  adjoint = assembler->createVec();
  dfdu->incref();
  adjoint->incref();

  int mpi_rank;
  MPI_Comm_rank(assembler->getMPIComm(), &mpi_rank);

  // Set up the solver
  if (mg) {
    int nrestart = 5;
    int is_flexible = 0;
    double atol = 1e-30;
    ksm = new GMRES(mg->getMat(0), mg, gmres_iters, nrestart, is_flexible);
    ksm->incref();
    ksm->setMonitor(new KSMPrintStdout("GMRES", mpi_rank, 10));
    ksm->setTolerances(rtol, atol);
  } else {
    ksm = NULL;
  }

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

  // Callback function information
  output_callback_ptr = NULL;
  writeOutputCallback = NULL;

  num_callback_constraints = 0;
  num_callback_ineq_constraints = 0;
  constraint_callback_ptr = NULL;
  constraintCallback = NULL;
  constraint_gradient_callback_ptr = NULL;
  constraintGradientCallback = NULL;
  objective_callback_ptr = NULL;
  objectiveCallback = NULL;
  objective_gradient_callback_ptr = NULL;
  objectiveGradientCallback = NULL;
  qn_correction_callback_ptr = NULL;
  qn_correction_callback_zlen = 0;
  qnCorrectionCallback = NULL;

  // By default, we don't use quasi-Newton update correction
  use_qn_correction_comp_obj = 0;

  // Set default finite difference step
  dh_Kmat_2nd_deriv = 1e-6;

  proj_deriv = NULL;
  x_h = NULL;
  s_temp = NULL;
}

/*
  Free the data stored in the object
*/
TMRTopoProblem::~TMRTopoProblem() {
  if (prefix) {
    delete[] prefix;
  }

  // Free the initial design variable values (if allocated)
  if (xinit) {
    xinit->decref();
  }
  if (xlb) {
    xlb->decref();
  }
  if (xub) {
    xub->decref();
  }

  // Free the solver/multigrid information
  if (mg) {
    mg->decref();
  }
  if (ksm) {
    ksm->decref();
  }

  dfdu->decref();
  adjoint->decref();
  assembler->decref();

  // Free the variables/forces
  if (forces) {
    for (int i = 0; i < num_load_cases; i++) {
      if (forces[i]) {
        forces[i]->decref();
      }
      vars[i]->decref();
    }
    delete[] forces;
    delete[] vars;
  }

  // Free the load case data
  if (load_case_info) {
    for (int i = 0; i < num_load_cases; i++) {
      for (int j = 0; j < load_case_info[i].num_funcs; j++) {
        load_case_info[i].funcs[j]->decref();
      }
      if (load_case_info[i].funcs) {
        delete[] load_case_info[i].funcs;
      }
      if (load_case_info[i].offset) {
        delete[] load_case_info[i].offset;
      }
      if (load_case_info[i].scale) {
        delete[] load_case_info[i].scale;
      }
    }
    delete[] load_case_info;
  }

  // Free data for the linear constraint
  if (linear_offset) {
    delete[] linear_offset;
  }
  if (Alinear) {
    for (int i = 0; i < num_linear_con; i++) {
      if (Alinear[i]) {
        Alinear[i]->decref();
      }
    }
  }

  // If the objective weights exist, delete them
  if (obj_weights) {
    delete[] obj_weights;
  }

  // Delete the buckling and frequency TACS object
  if (freq) {
    freq->decref();
  }

  // Free the array of KS weights
  if (obj_funcs) {
    for (int i = 0; i < num_load_cases; i++) {
      if (obj_funcs[i]) {
        obj_funcs[i]->decref();
      }
    }
  }

  // Free extra vectors for qn update correction
  if (use_qn_correction_comp_obj) {
    proj_deriv->decref();
    x_h->decref();
    s_temp->decref();
  }
}

/*
  Get the TACSAssembler object associated with the filter
*/
TACSAssembler *TMRTopoProblem::getAssembler() { return filter->getAssembler(); }

/*
  Get the TMRQuadForest object associated with the filter (may be NULL)
*/
TMRQuadForest *TMRTopoProblem::getFilterQuadForest() {
  return filter->getFilterQuadForest();
}

/*
  Get the TMRQuadForest object associated with the filter (may be NULL)
*/
TMROctForest *TMRTopoProblem::getFilterOctForest() {
  return filter->getFilterOctForest();
}

/*
  Get the TMRTopoFilter object associated with the TopoProblem
*/
TMRTopoFilter *TMRTopoProblem::getTopoFilter() { return filter; }

/*
  Get the TACSMg object associated with the TopoProblem
*/
TACSMg *TMRTopoProblem::getMg() { return mg; }

/*
  Set the output flags for any f5 output files that will be created
*/
void TMRTopoProblem::setF5OutputFlags(int freq, ElementType elem_type,
                                      int flag) {
  f5_frequency = freq;
  f5_element_type = elem_type;
  f5_write_flag = flag;
}

/*
  Set the output flags specifically for the eigenvalues
*/
void TMRTopoProblem::setF5EigenOutputFlags(int freq, ElementType elem_type,
                                           int flag) {
  f5_eigen_frequency = freq;
  f5_eigen_element_type = elem_type;
  f5_eigen_write_flag = flag;
}

/*
  Set the directory prefix to use for this load case
*/
void TMRTopoProblem::setPrefix(const char *_prefix) {
  if (prefix) {
    delete[] prefix;
  }
  prefix = new char[strlen(_prefix) + 1];
  strcpy(prefix, _prefix);
}

/*
  Set the load cases for each problem
*/
void TMRTopoProblem::setLoadCases(TACSBVec **_forces, int _num_load_cases) {
  if (!mg) {
    fprintf(stderr,
            "TMRTopoProblem: Cannot call setLoadCases, multigrid object not "
            "defined\n");
    return;
  }

  // Pre-incref the input forces
  for (int i = 0; i < _num_load_cases; i++) {
    if (_forces[i]) {
      _forces[i]->incref();
    }
  }

  // Free the forces/variables if any exist
  if (forces) {
    for (int i = 0; i < num_load_cases; i++) {
      if (forces[i]) {
        forces[i]->decref();
      }
      vars[i]->decref();
    }
    delete[] forces;
    delete[] vars;
  }

  // Deallocate the load case data (if it exists)
  if (load_case_info) {
    for (int i = 0; i < num_load_cases; i++) {
      for (int j = 0; j < load_case_info[i].num_funcs; j++) {
        load_case_info[i].funcs[j]->decref();
      }
      if (load_case_info[i].funcs) {
        delete[] load_case_info[i].funcs;
      }
      if (load_case_info[i].offset) {
        delete[] load_case_info[i].offset;
      }
      if (load_case_info[i].scale) {
        delete[] load_case_info[i].scale;
      }
    }
    delete[] load_case_info;
  }

  num_load_cases = _num_load_cases;
  forces = new TACSBVec *[num_load_cases];
  vars = new TACSBVec *[num_load_cases];
  for (int i = 0; i < num_load_cases; i++) {
    forces[i] = _forces[i];
    assembler->setBCs(forces[i]);
    vars[i] = assembler->createVec();
    vars[i]->incref();
    vars[i]->zeroEntries();
  }

  // Allocate the load case information
  load_case_info = NULL;
  if (num_load_cases > 0) {
    load_case_info = new LoadCaseInfo[num_load_cases];
    for (int i = 0; i < num_load_cases; i++) {
      load_case_info[i].num_funcs = 0;
      load_case_info[i].funcs = NULL;
      load_case_info[i].offset = NULL;
      load_case_info[i].scale = NULL;
    }
  }
}

/*
  Get the number of load cases
*/
int TMRTopoProblem::getNumLoadCases() { return num_load_cases; }

/*
  Set the constraint functions for each of the specified load cases
*/
void TMRTopoProblem::addConstraints(int load_case, TACSFunction **funcs,
                                    const TacsScalar *offset,
                                    const TacsScalar *scale, int num_funcs) {
  if (!load_case_info) {
    fprintf(stderr,
            "TMRTopoProblem error: Must call setLoadCases() "
            "before adding constraints\n");
    return;
  }
  if (load_case < 0 || load_case >= num_load_cases) {
    fprintf(stderr, "TMRTopoProblem error: Load case out of range\n");
    return;
  }

  for (int i = 0; i < num_funcs; i++) {
    funcs[i]->incref();
  }

  // Free the load case if it has been allocated before
  if (load_case_info[load_case].num_funcs > 0) {
    for (int j = 0; j < load_case_info[load_case].num_funcs; j++) {
      load_case_info[load_case].funcs[j]->decref();
    }
    delete[] load_case_info[load_case].funcs;
    delete[] load_case_info[load_case].offset;
    delete[] load_case_info[load_case].scale;
  }

  // Allocate the data
  load_case_info[load_case].num_funcs = num_funcs;

  if (num_funcs > 0) {
    load_case_info[load_case].funcs = new TACSFunction *[num_funcs];
    load_case_info[load_case].offset = new TacsScalar[num_funcs];
    load_case_info[load_case].scale = new TacsScalar[num_funcs];

    // Copy over the values
    for (int i = 0; i < num_funcs; i++) {
      load_case_info[load_case].funcs[i] = funcs[i];
      load_case_info[load_case].offset[i] = offset[i];
      load_case_info[load_case].scale[i] = scale[i];
    }
  } else {
    load_case_info[load_case].funcs = NULL;
    load_case_info[load_case].offset = NULL;
    load_case_info[load_case].scale = NULL;
  }
}

/*
  Add linear constraints to the problem.
*/
void TMRTopoProblem::addLinearConstraints(ParOptVec **vecs, TacsScalar *offset,
                                          int _ncon) {
  for (int i = 0; i < _ncon; i++) {
    vecs[i]->incref();
  }

  if (linear_offset) {
    delete[] linear_offset;
  }
  if (Alinear) {
    for (int i = 0; i < num_linear_con; i++) {
      Alinear[i]->decref();
    }
    delete[] Alinear;
  }

  // Allocate the new space
  num_linear_con = _ncon;
  linear_offset = new TacsScalar[num_linear_con];
  Alinear = new ParOptVec *[num_linear_con];
  for (int i = 0; i < num_linear_con; i++) {
    linear_offset[i] = offset[i];
    Alinear[i] = vecs[i];
  }
}

/*
  Add a natural frequency constraint
*/
void TMRTopoProblem::addFrequencyConstraint(
    double sigma, int num_eigvals, TacsScalar ks_weight, TacsScalar offset,
    TacsScalar scale, int max_subspace_size, double eigtol, int use_jd,
    int fgmres_size, double eig_rtol, double eig_atol, int num_recycle,
    JDRecycleType recycle_type) {
  if (!mg) {
    fprintf(stderr,
            "TMRTopoProblem: Cannot call addFrequencyConstraint, multigrid "
            "object not defined\n");
    return;
  }

  if (!freq) {
    // Create a mass matrix for the frequency constraint
    TACSMat *mmat = assembler->createMat();
    if (use_jd) {
      // Create the preconditioner matrix
      TACSMat *kmat = assembler->createMat();
      TACSMat *pcmat = mg->getMat(0);

      // Get preconditioner from Mg
      TACSPc *pc = mg;

      freq = new TACSFrequencyAnalysis(assembler, sigma, mmat, kmat, pcmat, pc,
                                       max_subspace_size, fgmres_size,
                                       num_eigvals, eigtol, eig_rtol, eig_atol,
                                       num_recycle, recycle_type);
    } else {
      // Create the frequency analysis object
      freq =
          new TACSFrequencyAnalysis(assembler, sigma, mmat, mg->getMat(0), ksm,
                                    max_subspace_size, num_eigvals, eigtol);
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
}

/*
  Add a buckling constraint
*/
void TMRTopoProblem::addBucklingConstraint(double sigma, int num_eigvals,
                                           TacsScalar ks_weight,
                                           TacsScalar offset, TacsScalar scale,
                                           int max_lanczos, double eigtol) {
  if (!mg) {
    fprintf(stderr,
            "TMRTopoProblem: Cannot call addBucklingConstraint, multigrid "
            "object not defined\n");
    return;
  }

  if (!buck) {
    // Create a geometric stiffness matrix for buckling constraint
    TACSMat *gmat = assembler->createMat();
    TACSMat *kmat = assembler->createMat();
    TACSMat *aux_mat;
    ksm->getOperators(&aux_mat, NULL);

    buck = new TACSLinearBuckling *[num_load_cases];
    buck_ks_sum = new TacsScalar[num_load_cases];
    for (int i = 0; i < num_load_cases; i++) {
      // Create the buckling analysis object
      buck[i] = new TACSLinearBuckling(assembler, sigma, gmat, kmat, aux_mat,
                                       ksm, max_lanczos, num_eigvals, eigtol);

      buck[i]->incref();
    }
  }

  // Set a parameters that control how the natural frequency
  // constraint is implemented
  buck_eig_tol = eigtol;
  num_buck_eigvals = num_eigvals;
  memset(buck_ks_sum, 0.0, num_load_cases * sizeof(TacsScalar));
  buck_ks_weight = ks_weight;
  buck_offset = offset;
  buck_scale = scale;
}

/*
  Add callback constraint and constraint gradient calls
*/
void TMRTopoProblem::addConstraintCallback(
    int ncon, int nineq, void *con_ptr,
    void (*confunc)(void *, TMRTopoFilter *, TACSMg *, int, TacsScalar *),
    void *con_grad_ptr,
    void (*congradfunc)(void *, TMRTopoFilter *, TACSMg *, int, TACSBVec **)) {
  if (ncon > 0 && confunc && congradfunc) {
    num_callback_constraints = ncon;
    num_callback_ineq_constraints = nineq;
    constraint_callback_ptr = con_ptr;
    constraintCallback = confunc;
    constraint_gradient_callback_ptr = con_grad_ptr;
    constraintGradientCallback = congradfunc;
  }
}

/*
  Add quasi-Newton update correction callback calls
*/
void TMRTopoProblem::addQnCorrectionCallback(
    int ncon, void *callback_ptr,
    void (*callback_fun)(int, void *, ParOptVec *, ParOptScalar *, ParOptVec *,
                         ParOptVec *, ParOptVec *)) {
  if (callback_ptr) {
    qn_correction_callback_zlen = ncon;
    qn_correction_callback_ptr = callback_ptr;
    qnCorrectionCallback = callback_fun;
  }
}

/*
  Set the objective weight values - this indicates a compliance objective
*/
void TMRTopoProblem::setObjective(const TacsScalar *_obj_weights) {
  if (!obj_weights) {
    obj_weights = new TacsScalar[num_load_cases];
  }
  if (obj_funcs) {
    for (int i = 0; i < num_load_cases; i++) {
      if (obj_funcs[i]) {
        obj_funcs[i]->decref();
      }
    }
  }
  for (int i = 0; i < num_load_cases; i++) {
    obj_weights[i] = _obj_weights[i];
  }
}

/*
  Set the objective weight values - this indicates a compliance objective
*/
void TMRTopoProblem::setObjective(const TacsScalar *_obj_weights,
                                  TACSFunction **_obj_funcs) {
  if (!obj_weights) {
    obj_weights = new TacsScalar[num_load_cases];
  }
  if (!obj_funcs) {
    obj_funcs = new TACSFunction *[num_load_cases];
  } else {
    for (int i = 0; i < num_load_cases; i++) {
      if (obj_funcs[i]) {
        obj_funcs[i]->decref();
      }
    }
  }
  for (int i = 0; i < num_load_cases; i++) {
    obj_weights[i] = _obj_weights[i];
    obj_funcs[i] = _obj_funcs[i];
    obj_funcs[i]->incref();
  }
}

/*
  Add callback for objective and objective gradient calls
*/
void TMRTopoProblem::addObjectiveCallback(
    void *obj_ptr,
    void (*objfunc)(void *, TMRTopoFilter *, TACSMg *, TacsScalar *),
    void *obj_grad_ptr,
    void (*objgradfunc)(void *, TMRTopoFilter *, TACSMg *, TACSBVec *)) {
  if (objfunc && objgradfunc) {
    objective_callback_ptr = obj_ptr;
    objectiveCallback = objfunc;
    objective_gradient_callback_ptr = obj_grad_ptr;
    objectiveGradientCallback = objgradfunc;
  }
}

/*
  Initialize the problem sizes
*/
void TMRTopoProblem::initialize() {
  // Add up the total number of constraints
  int num_constraints = num_linear_con;

  for (int i = 0; i < num_load_cases; i++) {
    num_constraints += load_case_info[i].num_funcs;
  }

  if (freq) {
    num_constraints++;
  }
  if (buck) {
    for (int i = 0; i < num_load_cases; i++) {
      num_constraints++;
    }
  }
  num_constraints += num_callback_constraints;

  // Get the design node map
  TACSNodeMap *designMap = assembler->getDesignNodeMap();

  // Set the problem sizes
  int nvars = design_vars_per_node * designMap->getNumNodes();
  int nw = 0;
  int nwblock = 0;
  if (design_vars_per_node > 1) {
    nw = nvars / design_vars_per_node;
    nwblock = 1;
  }

  int num_eq = num_callback_constraints - num_callback_ineq_constraints;
  int num_ineq = num_constraints - num_eq;

  setProblemSizes(nvars, num_constraints, nw);
  setNumInequalities(num_ineq, nwblock);
}

/*
  Set the initial design variable values
*/
void TMRTopoProblem::setInitDesignVars(ParOptVec *xvars, ParOptVec *lb,
                                       ParOptVec *ub) {
  // Create the design vector if it doesn't exist
  if (xvars) {
    if (!xinit) {
      xinit = createDesignVec();
      xinit->incref();
    }

    // Copy over the local values
    xinit->copyValues(xvars);
    setDesignVars(xinit);
  }

  // Copy over the lower bound
  if (lb) {
    if (!xlb) {
      xlb = createDesignVec();
      xlb->incref();
    }
    xlb->copyValues(lb);
  }

  // Copy over the upper bound
  if (ub) {
    if (!xub) {
      xub = createDesignVec();
      xub->incref();
    }
    xub->copyValues(ub);
  }
}

/*
  Set the iteration count
*/
void TMRTopoProblem::setIterationCounter(int iter) { iter_count = iter; }

/*
  Create a design variable vector
*/
ParOptVec *TMRTopoProblem::createDesignVec() {
  return new ParOptBVecWrap(assembler->createDesignVec());
}

/*
  Set whether these should be considered sparse inequalities
*/
int TMRTopoProblem::isSparseInequality() {
  // These are sparse equality constraints
  return 0;
}

/*
  Always use the lower bound variables
*/
int TMRTopoProblem::useLowerBounds() { return 1; }

/*
  We impose an upper bound on the variables even though in the case
  when we have reciprocal variables they are not really defeind.
*/
int TMRTopoProblem::useUpperBounds() {
  if (design_vars_per_node > 1) {
    return 0;
  }
  return 1;
}

// TODO: Is this right?
ParOptQuasiDefMat *TMRTopoProblem::createQuasiDefMat() {
  int nwblock = 0;
  if (design_vars_per_node > 1) {
    nwblock = 1;
  }
  return new ParOptQuasiDefBlockMat(this, nwblock);
}

/*
  Set the initial design variables
*/
void TMRTopoProblem::getVarsAndBounds(ParOptVec *xvec, ParOptVec *lbvec,
                                      ParOptVec *ubvec) {
  int mpi_rank;
  MPI_Comm_rank(assembler->getMPIComm(), &mpi_rank);
  if (xvec) {
    // If the initial design vector exists, copy it directly
    if (xinit) {
      xvec->copyValues(xinit);
    } else {
      ParOptBVecWrap *wrap = dynamic_cast<ParOptBVecWrap *>(xvec);
      if (wrap) {
        assembler->getDesignVars(wrap->vec);
      }
    }
  }

  // Handle the lower and upper bounds
  if (lbvec || ubvec) {
    int has_lower = 0, has_upper = 0;
    if (lbvec && xlb) {
      has_lower = 1;
      lbvec->copyValues(xlb);
    }
    if (ubvec && xub) {
      has_upper = 1;
      ubvec->copyValues(xub);
    }

    if (!(has_lower && has_upper)) {
      // Get the design variable values from TACS
      ParOptBVecWrap *lbwrap = dynamic_cast<ParOptBVecWrap *>(lbvec);
      ParOptBVecWrap *ubwrap = dynamic_cast<ParOptBVecWrap *>(ubvec);
      if (lbwrap && ubwrap) {
        assembler->getDesignVarRange(lbwrap->vec, ubwrap->vec);
      }
    }
  }
}

/*
  Set the design variable values consistently across all multigrid
  levels by interpolating the design variables using the filter
  interpolation objects.
*/
void TMRTopoProblem::setDesignVars(ParOptVec *pxvec) {
  ParOptBVecWrap *wrap = dynamic_cast<ParOptBVecWrap *>(pxvec);

  if (wrap) {
    // Set the design variable values
    filter->setDesignVars(wrap->vec);
  }
}

/*
  Evaluate the objective and constraints
*/
int TMRTopoProblem::evalObjCon(ParOptVec *pxvec, ParOptScalar *fobj,
                               ParOptScalar *cons) {
  // Get the rank of comm
  int mpi_rank;
  MPI_Comm_rank(assembler->getMPIComm(), &mpi_rank);

  // Set the design variable values on all mesh levels
  setDesignVars(pxvec);

  // Set the objective value
  *fobj = 0.0;

  // Keep track of the constraint number
  int count = 0;

  // Compute the linear constraint
  for (int i = 0; i < num_linear_con; i++, count++) {
    cons[count] = linear_offset[i] + Alinear[i]->dot(pxvec);
  }

  // Perform the objective function callback
  if (objectiveCallback) {
    TacsScalar fcallback = 0.0;
    objectiveCallback(objective_callback_ptr, filter, mg, &fcallback);
    *fobj = fcallback;
  }

  if (num_load_cases > 0 && mg) {
    // Zero the variables
    assembler->zeroVariables();

    // Assemble the Jacobian on each level
    double alpha = 1.0, beta = 0.0, gamma = 0.0;
    mg->assembleJacobian(alpha, beta, gamma, NULL);
    mg->factor();

    for (int i = 0; i < num_load_cases; i++) {
      if (forces[i]) {
        // Solve the system: K(x)*u = forces
        ksm->solve(forces[i], vars[i]);
        assembler->setBCs(vars[i]);

        // Set the variables into TACSAssembler
        assembler->setVariables(vars[i]);

        // Add the contribution to the objective
        if (obj_funcs) {
          if (obj_funcs[i]) {
            TacsScalar fobj_val;
            assembler->evalFunctions(1, &obj_funcs[i], &fobj_val);

            if (mpi_rank == 0) {
              printf("%-30s %25.10e\n", obj_funcs[i]->getObjectName(),
                     fobj_val);
              TACSKSFailure *ks_fail =
                  dynamic_cast<TACSKSFailure *>(obj_funcs[i]);
              if (ks_fail) {
                TacsScalar max_fail = ks_fail->getMaximumFailure();
                printf("%-30s %25.10e\n", "TACSKSFailure max stress", max_fail);
              }
            }

            *fobj += obj_weights[i] * fobj_val;
          }
        }
        // If we are dealing with a compliance objective
        else {
          *fobj += obj_weights[i] * vars[i]->dot(forces[i]);
        }

        // Evaluate the constraints
        int num_funcs = load_case_info[i].num_funcs;
        if (num_funcs > 0) {
          assembler->evalFunctions(num_funcs, load_case_info[i].funcs,
                                   &cons[count]);

          if (mpi_rank == 0) {
            for (int k = 0; k < num_funcs; k++) {
              printf("%-30s %25.10e\n",
                     load_case_info[i].funcs[k]->getObjectName(),
                     cons[count + k]);
              TACSKSFailure *ks_fail =
                  dynamic_cast<TACSKSFailure *>(load_case_info[i].funcs[k]);
              if (ks_fail) {
                TacsScalar max_fail = ks_fail->getMaximumFailure();
                printf("%-30s %25.10e\n", "TACSKSFailure max stress", max_fail);
              }
            }
          }

          // Scale and offset the constraints that we just evaluated
          for (int j = 0; j < num_funcs; j++) {
            TacsScalar offset = load_case_info[i].offset[j];
            TacsScalar scale = load_case_info[i].scale[j];
            cons[count + j] = scale * (cons[count + j] + offset);
          }
          count += num_funcs;
        }
      }
    }
  }

  // Compute the natural frequency constraint, if any
  if (freq) {
    // Keep track of the number of eigenvalues with unacceptable
    // error. If more than one exists, re-solve the eigenvalue
    // problem again.
    int err_count = 1;

    // Keep track of the smallest eigenvalue
    double shift = 0.95;
    TacsScalar smallest_eigval = 0.0;
    while (err_count > 0) {
      // Set the error counter to zero
      err_count = 0;
      // Solve the eigenvalue problem
      freq->solve(new KSMPrintStdout("KSM", mpi_rank, 1));

      // Extract the first k eigenvalues
      for (int k = 0; k < num_freq_eigvals; k++) {
        TacsScalar error;
        TacsScalar eigval = freq->extractEigenvalue(k, &error);
        if (eigval < 0.0) {
          eigval *= -1.0;
        }

        if (k == 0) {
          smallest_eigval = eigval;
        }
        if (eigval < smallest_eigval) {
          smallest_eigval = eigval;
        }
        if (error > freq_eig_tol) {
          err_count++;
        }
      }
      err_count = 0;

      // If there is significant error in computing the eigenvalues,
      // reset the buckling computation
      if (err_count > 0) {
        double sigma = shift * smallest_eigval;
        shift += 0.5 * (1.0 - shift) + shift;
        freq->setSigma(sigma);
      }
    }

    // Evaluate the KS function of the lowest eigenvalues
    freq_ks_sum = 0.0;

    // Weight on the KS function
    for (int k = 0; k < num_freq_eigvals; k++) {
      TacsScalar error;
      TacsScalar eigval = freq->extractEigenvalue(k, &error);
      if (eigval < 0.0) {
        eigval *= -1.0;
      }

      // Add up the contribution to the ks function
      freq_ks_sum += exp(-freq_ks_weight * (eigval - smallest_eigval));
    }

    // Evaluate the KS function of the aggregation of the eigenvalues
    cons[count] = (smallest_eigval - log(freq_ks_sum) / freq_ks_weight);
    cons[count] = freq_scale * (cons[count] + freq_offset);
    count++;
  }

  // Compute the buckling constraint, if any
  if (buck) {
    for (int i = 0; i < num_load_cases; i++) {
      if (forces[i]) {
        // Keep track of the number of eigenvalues with unacceptable
        // error. If more than one exists, re-solve the eigenvalue
        // problem again.
        int err_count = 1;

        // Keep track of the smallest eigenvalue
        double shift = 0.95;
        TacsScalar smallest_eigval = 0.0;

        while (err_count > 0) {
          // Set the error counter to zero
          err_count = 0;

          // Solve the eigenvalue problem
          TACSBVec *u0 = NULL;  // TODO: Is this right? do we need a u0? --Aaron
          buck[i]->solve(forces[i], u0, new KSMPrintStdout("KSM", mpi_rank, 1));

          // Extract the first k eigenvalues
          for (int k = 0; k < num_buck_eigvals; k++) {
            TacsScalar error;
            TacsScalar eigval = buck[i]->extractEigenvalue(k, &error);
            if (eigval < 0.0) {
              eigval *= -1.0;
            }

            if (k == 0) {
              smallest_eigval = eigval;
            }
            if (error > buck_eig_tol) {
              err_count++;
            }
          }

          // If there is significant error in computing the eigenvalues,
          // reset the buckling computation
          if (err_count > 0) {
            double sigma = shift * smallest_eigval;
            shift += 0.5 * (1.0 - shift) + shift;
            buck[i]->setSigma(sigma);
          }
        }

        // Evaluate the KS function of the lowest eigenvalues
        buck_ks_sum[i] = 0.0;

        // Weight on the KS function
        for (int k = 0; k < num_buck_eigvals; k++) {
          TacsScalar error;
          TacsScalar eigval = buck[i]->extractEigenvalue(k, &error);
          if (eigval < 0.0) {
            eigval *= -1.0;
          }

          // Add up the contribution to the ks function
          buck_ks_sum[i] += exp(-buck_ks_weight * (eigval - smallest_eigval));
        }

        // Evaluate the KS function of the aggregation of the eigenvalues
        cons[count] = (smallest_eigval - log(buck_ks_sum[i]) / buck_ks_weight);
        cons[count] = buck_scale * (cons[count] + buck_offset);
        count++;
      }
    }
  }

  // Evaluate the callback constraints
  if (constraintCallback) {
    constraintCallback(constraint_callback_ptr, filter, mg,
                       num_callback_constraints, &cons[count]);
  }

  return 0;
}

/*
  Evaluate the objective and constraint gradients
*/
int TMRTopoProblem::evalObjConGradient(ParOptVec *xvec, ParOptVec *gvec,
                                       ParOptVec **Acvec) {
  int mpi_rank;
  MPI_Comm_rank(assembler->getMPIComm(), &mpi_rank);

  // Evaluate the derivative of the weighted compliance with
  // respect to the design variables
  ParOptBVecWrap *wrap = dynamic_cast<ParOptBVecWrap *>(gvec);
  if (wrap) {
    TACSBVec *g = wrap->vec;
    g->zeroEntries();

    // // Add non-design mass to mass matrix, if specified
    // if (m0vec){
    //   m0vec->axpy(1.0, xvec);
    //   setDesignVars(m0vec);
    // }

    // Evaluate the gradient of the objective function from the callback
    if (objectiveGradientCallback) {
      objectiveGradientCallback(objective_gradient_callback_ptr, filter, mg, g);
    }

    // // Reset design variable if having non-design mass
    // if (m0vec){
    //   setDesignVars(xvec);
    // }

    // Evaluate the gradient of the objective. If no objective functions are
    // set, the weighted sum of the compliance is used. Otherwise the weighted
    // sum of the objective functions from each load case are used
    if (obj_funcs && mg) {
      for (int i = 0; i < num_load_cases; i++) {
        assembler->setVariables(vars[i]);

        // In the case of the structural mass, no objective is needed
        int use_adjoint = 1;
        if (dynamic_cast<TACSStructuralMass *>(obj_funcs[i])) {
          use_adjoint = 0;
        }

        if (use_adjoint) {
          // Assemble the transpose of the Jacobian matrix
          double alpha = 1.0, beta = 0.0, gamma = 0.0;
          mg->assembleJacobian(alpha, beta, gamma, NULL, TACS_MAT_TRANSPOSE);
          mg->factor();

          // Compute the right-hand-side
          dfdu->zeroEntries();
          assembler->addSVSens(alpha, beta, gamma, 1, &obj_funcs[i], &dfdu);
          assembler->applyBCs(dfdu);

          // Solve the system of adjoint equations
          ksm->solve(dfdu, adjoint);
          assembler->addDVSens(obj_weights[i], 1, &obj_funcs[i], &g);
          assembler->addAdjointResProducts(-obj_weights[i], 1, &adjoint, &g);
        } else {
          assembler->addDVSens(obj_weights[i], 1, &obj_funcs[i], &g);
        }
      }
    } else {  // For compliance objective
      for (int i = 0; i < num_load_cases; i++) {
        assembler->setVariables(vars[i]);
        assembler->addAdjointResProducts(-obj_weights[i], 1, &vars[i], &g);
      }
    }

    double gnorm = g->norm();
    filter->addValues(g);  // Apply filter transpose to the gradient
    double fgnorm = g->norm();
    if (mpi_rank == 0) {
      printf("[TMRTopoProblem]unfiltered objective gradient norm: %20.10e\n",
             gnorm);
      printf("[TMRTopoProblem]filtered objective gradient norm:   %20.10e\n",
             fgnorm);
    }
  } else {
    return 1;
  }

  // Keep track of the constraint gradient number
  int count = 0;

  // Set the linear constraint
  for (int i = 0; i < num_linear_con; i++, count++) {
    Acvec[count]->copyValues(Alinear[i]);
  }

  // Compute the derivative of the constraint functions
  for (int i = 0; i < num_load_cases; i++) {
    assembler->setVariables(vars[i]);

    // Get the number of functions for each load case
    int num_funcs = load_case_info[i].num_funcs;

    for (int j = 0; j < num_funcs; j++) {
      TACSFunction *func = load_case_info[i].funcs[j];
      TacsScalar scale = load_case_info[i].scale[j];

      // Try to unwrap the vector
      wrap = dynamic_cast<ParOptBVecWrap *>(Acvec[count + j]);

      if (wrap) {
        TACSBVec *A = wrap->vec;
        A->zeroEntries();

        // If the function is the structural mass, then do not
        // use the adjoint, otherwise assume that we should use the
        // adjoint method to compute the gradient.
        int use_adjoint = 1;
        if (dynamic_cast<TACSStructuralMass *>(func)) {
          use_adjoint = 0;
        }
        if (use_adjoint) {
          // Evaluate the right-hand-side
          dfdu->zeroEntries();
          double alpha = 1.0, beta = 0.0, gamma = 0.0;
          assembler->addSVSens(alpha, beta, gamma, 1, &func, &dfdu);
          assembler->applyBCs(dfdu);

          // Solve the system of equations
          ksm->solve(dfdu, adjoint);

          // Compute the total derivative using the adjoint
          assembler->addDVSens(scale, 1, &func, &A);
          assembler->addAdjointResProducts(-scale, 1, &adjoint, &A);
        } else {
          assembler->addDVSens(scale, 1, &func, &A);
        }

        filter->addValues(A);
      }
    }
    count += num_funcs;

    if (freq && i == 0) {
      // Try to unwrap the vector
      wrap = dynamic_cast<ParOptBVecWrap *>(Acvec[count]);
      if (wrap) {
        // Get the underlying TACS vector for the design variables
        TACSBVec *A = wrap->vec;
        A->zeroEntries();

        // Create a temporary design vector
        TACSBVec *temp = assembler->createDesignVec();
        temp->incref();

        // Add the contribution from each eigenvalue derivative
        TacsScalar smallest_eigval = 0.0;

        // Extract the first k eigenvalues to find smallest eigval
        for (int k = 0; k < num_freq_eigvals; k++) {
          TacsScalar error;
          TacsScalar eigval = freq->extractEigenvalue(k, &error);
          if (eigval < 0.0) {
            eigval *= -1.0;
          }

          if (k == 0) {
            smallest_eigval = eigval;
          }
          if (eigval < smallest_eigval) {
            smallest_eigval = eigval;
          }
        }

        for (int k = 0; k < num_freq_eigvals; k++) {
          // Compute the derivaive of the eigenvalue
          freq->evalEigenDVSens(k, temp);

          // Extract the eigenvalue itself
          TacsScalar error;
          TacsScalar ks_grad_weight = 1.0;
          TacsScalar eigval = freq->extractEigenvalue(k, &error);
          if (eigval < 0.0) {
            eigval *= -1.0;
            ks_grad_weight = -1.0;
          }

          // Evaluate the weight on the gradient
          ks_grad_weight *=
              exp(-freq_ks_weight * (eigval - smallest_eigval)) / freq_ks_sum;

          // Scale the constraint by the frequency scaling value
          ks_grad_weight *= freq_scale;

          // Add contribution to eigenvalue gradient
          A->axpy(ks_grad_weight, temp);
        }

        // Free the data
        temp->decref();

        filter->addValues(A);
      }

      count++;
    }
    if (buck) {
      // Compute the derivative of the buckling functions
      // Try to unwrap the vector
      wrap = dynamic_cast<ParOptBVecWrap *>(Acvec[count]);
      if (wrap) {
        // Get the vector
        TACSBVec *A = wrap->vec;
        A->zeroEntries();

        // Create a temporary design vector
        TACSBVec *temp = assembler->createDesignVec();

        // Add the contribution from each eigenvalue derivative
        TacsScalar smallest_eigval = 0.0;
        for (int k = 0; k < num_buck_eigvals; k++) {
          // Compute the derivaive of the eigenvalue
          buck[i]->evalEigenDVSens(k, temp);

          // Extract the eigenvalue itself
          TacsScalar error;
          TacsScalar ks_grad_weight = 1.0;
          TacsScalar eigval = buck[i]->extractEigenvalue(k, &error);
          if (eigval < 0.0) {
            eigval *= -1.0;
            ks_grad_weight = -1.0;
          }
          if (k == 0) {
            smallest_eigval = eigval;
          }

          // Evaluate the weight on the gradient
          ks_grad_weight *= exp(-buck_ks_weight * (eigval - smallest_eigval)) /
                            buck_ks_sum[i];

          // Scale the constraint by the buckling scaling value
          ks_grad_weight *= buck_scale;

          // Add contribution to eigenvalue gradient
          A->axpy(ks_grad_weight, temp);
        }

        // Free the data
        temp->decref();

        filter->addValues(A);
      }
      count++;
    }
  }  // end num_load_cases

  // Evaluate the callback constraints
  if (constraintGradientCallback) {
    TACSBVec **vecs = new TACSBVec *[num_callback_constraints];
    for (int i = 0; i < num_callback_constraints; i++) {
      vecs[i] = NULL;
      wrap = dynamic_cast<ParOptBVecWrap *>(Acvec[count]);
      if (wrap) {
        vecs[i] = wrap->vec;
      }
      count++;
    }
    constraintGradientCallback(constraint_gradient_callback_ptr, filter, mg,
                               num_callback_constraints, vecs);

    for (int i = 0; i < num_callback_constraints; i++) {
      double vi_norm = vecs[i]->norm();
      filter->addValues(vecs[i]);
      double fvi_norm = vecs[i]->norm();
      if (mpi_rank == 0) {
        printf("[TMRTopoProblem]unfiltered constraint gradient norm: %20.10e\n",
               vi_norm);
        printf("[TMRTopoProblem]filtered constraint gradient norm:   %20.10e\n",
               fvi_norm);
      }
    }
    delete[] vecs;
  }

  return 0;
}

// Switch on quasi-Newton correction for compliance objective
// Note that this can be called only once, otherwise will
// encounter memory leak for extra vectors.
// ----------------------------------------------------------
void TMRTopoProblem::useQnCorrectionComplianceObj() {
  // Switch on the flag
  use_qn_correction_comp_obj = 1;

  // Create space for extra vectors
  proj_deriv = createDesignVec();
  proj_deriv->incref();
  x_h = createDesignVec();
  x_h->incref();
  s_temp = createDesignVec();
  s_temp->incref();
  return;
}

// Compute a correction to the quasi-Newton update
// -----------------------------------------------
void TMRTopoProblem::computeQuasiNewtonUpdateCorrection(
    ParOptVec *x, ParOptScalar *z, ParOptVec *zw, ParOptVec *s, ParOptVec *y) {
  // Get processor rank
  int rank = 0;
  MPI_Comm_rank(assembler->getMPIComm(), &rank);

  // Perform the correction via callback, if any
  if (qnCorrectionCallback) {
    qnCorrectionCallback(ncon, qn_correction_callback_ptr, x, z, zw, s, y);
  }

  /*
    Perform quasi-Newton update correction for compliance objective

    The exact Hessian of the compliance is composed of the difference
    between two contributions:

    H = P - N

    Here P is a positive-definite term while N is positive semi-definite.
    Since the true Hessian is a difference between the two, the quasi-Newton
    Hessian update can be written as:

    H*s = y = P*s - N*s

    This often leads to damped update steps as the optimization converges.
    Instead, we want to approximate just P, so  we modify y so that

    ymod ~ P*s = (H + N)*s ~ y + N*s
  */
  if (use_qn_correction_comp_obj) {
    // Zero out the vector that stores projected derivative
    proj_deriv->zeroEntries();
    ParOptBVecWrap *proj_deriv_wrap =
        dynamic_cast<ParOptBVecWrap *>(proj_deriv);

    //[Test] Update design variable back to x and solve the system again
    // setDesignVars(x);
    // for ( int i = 0; i < num_load_cases; i++ ){
    //   if (forces[i]){
    //     ksm->solve(forces[i], vars[i]);
    //     assembler->setBCs(vars[i]);
    //     assembler->setVariables(vars[i]);
    //   }
    // }

    // Compute first derivative at x
    if (proj_deriv_wrap) {
      for (int i = 0; i < num_load_cases; i++) {
        // assembler->setVariables(vars[i]);
        assembler->addAdjointResProducts(-obj_weights[i], 1, &vars[i],
                                         &proj_deriv_wrap->vec);
      }
    }

    double gx = proj_deriv->norm();
    if (rank == 0) printf("Qn correction: norm of g(x)    = %.5e\n", gx);

    // Update design variable: x <-- x + h*s
    x_h->copyValues(x);
    x_h->axpy(dh_Kmat_2nd_deriv, s);
    setDesignVars(x_h);

    // Compute first derivative at x + h*s
    if (proj_deriv_wrap) {
      for (int i = 0; i < num_load_cases; i++) {
        // assembler->setVariables(vars[i]);
        assembler->addAdjointResProducts(obj_weights[i], 1, &vars[i],
                                         &proj_deriv_wrap->vec);
      }
    }

    double dgx = proj_deriv->norm();
    if (rank == 0) printf("Qn correction: norm of dg(x)   = %.5e\n", dgx);

    // Divide by h so we get projected derivative by finite differencing
    proj_deriv->scale(1 / dh_Kmat_2nd_deriv);

    double dgxdh = proj_deriv->norm();
    if (rank == 0) printf("Qn correction: norm of dg(x)/dh= %.5e\n", dgxdh);

    // Apply filter transpose to projected derivative
    filter->addValues(proj_deriv_wrap->vec);
    // filter->applyTranspose(proj_deriv_wrap->vec, proj_deriv_wrap->vec);  //
    // Same as addValues()

    double xval = x->norm();
    double sval = s->norm();
    double yyal = y->norm();
    double yup = proj_deriv->norm();
    if (rank == 0) {
      printf("Qn correction: norm of x       = %.5e\n", xval);
      printf("Qn correction: norm of s       = %.5e\n", sval);
      printf("Qn correction: norm of y       = %.5e\n", yyal);
      printf("Qn correction: norm of yupdate = %.5e\n", yup);
    }

    // Test if y^Ts > 0
    s_temp->copyValues(s);
    ParOptBVecWrap *s_temp_wrap = dynamic_cast<ParOptBVecWrap *>(s_temp);
    filter->applyTranspose(s_temp_wrap->vec, s_temp_wrap->vec);

    ParOptScalar yTs = s_temp->dot(proj_deriv);
    if (rank == 0) printf("Qn correction: yTs  = %.5e\n", yTs);

    // Update y
    y->axpy(1.0, proj_deriv);
  }

  return;
}

/*
  Evaluate the product of the Hessian with the given vector px
*/
int TMRTopoProblem::evalHvecProduct(ParOptVec *xvec, ParOptScalar *z,
                                    ParOptVec *zw, ParOptVec *pxvec,
                                    ParOptVec *hvec) {
  return 0;
}

// Evaluate the sparse constraints
// ------------------------
void TMRTopoProblem::evalSparseCon(ParOptVec *xvec, ParOptVec *outvec) {
  if (design_vars_per_node > 1) {
    // Dynamically cast the vectors to the ParOptBVecWrap class
    TacsScalar *x, *out;
    int size = xvec->getArray(&x);
    outvec->getArray(&out);

    // Compute the weighting constraints
    int n = size / design_vars_per_node;
    for (int i = 0; i < n; i++) {
      out[i] = -1.0;
      for (int j = 0; j < design_vars_per_node; j++) {
        out[i] += x[design_vars_per_node * i + j];
      }
    }
  }
}

// Compute the Jacobian-vector product out = J(x)*px
// --------------------------------------------------
void TMRTopoProblem::addSparseJacobian(double alpha, ParOptVec *xvec,
                                       ParOptVec *pxvec, ParOptVec *outvec) {
  if (design_vars_per_node > 1) {
    TacsScalar *px, *out;
    int size = pxvec->getArray(&px);
    outvec->getArray(&out);

    // Compute the matrix-vector product
    int n = size / design_vars_per_node;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < design_vars_per_node; j++) {
        out[i] += alpha * px[design_vars_per_node * i + j];
      }
    }
  }
}

// Compute the transpose Jacobian-vector product out = J(x)^{T}*pzw
// -----------------------------------------------------------------
void TMRTopoProblem::addSparseJacobianTranspose(double alpha, ParOptVec *x,
                                                ParOptVec *pzwvec,
                                                ParOptVec *outvec) {
  if (design_vars_per_node > 1) {
    TacsScalar *pzw, *out;
    pzwvec->getArray(&pzw);
    int size = outvec->getArray(&out);

    int n = size / design_vars_per_node;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < design_vars_per_node; j++) {
        out[design_vars_per_node * i + j] += alpha * pzw[i];
      }
    }
  }
}

// Add the inner product of the constraints to the matrix such
// that A += J(x)*cvec*J(x)^{T} where cvec is a diagonal matrix
// ------------------------------------------------------------
void TMRTopoProblem::addSparseInnerProduct(double alpha, ParOptVec *x,
                                           ParOptVec *cvec, double *A) {
  if (design_vars_per_node > 1) {
    TacsScalar *c;
    int size = cvec->getArray(&c);

    int n = size / design_vars_per_node;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < design_vars_per_node; j++) {
        A[i] += alpha * c[design_vars_per_node * i + j];
      }
    }
  }
}

/*
  Write the output file
*/
void TMRTopoProblem::writeOutput(int iter, ParOptVec *xvec) {
  ParOptBVecWrap *wrap = dynamic_cast<ParOptBVecWrap *>(xvec);
  if (wrap && writeOutputCallback) {
    writeOutputCallback(output_callback_ptr, prefix, iter,
                        filter->getFilterOctForest(),
                        filter->getFilterQuadForest(), wrap->vec);
  }

  // Print out the binary STL file for later visualization
  if (prefix) {
    // Write out the file at a cut off of 0.25
    char filename[strlen(prefix) + 100];

    for (int k = 0; k < design_vars_per_node; k++) {
      double cutoff = 0.5;
      snprintf(filename, sizeof(filename),
               "%s/levelset05_var%d_binary%04d.bstl", prefix, k, iter_count);

      // Write the STL file
      filter->writeSTLFile(k, cutoff, filename);
    }
  }

  if (prefix) {
    if (f5_frequency > 0 && (iter % f5_frequency) == 0) {
      // Create the filename
      char filename[strlen(prefix) + 100];
      snprintf(filename, sizeof(filename), "%s/tacs_output%04d.f5", prefix,
               iter_count);

      TACSToFH5 *f5 = new TACSToFH5(assembler, f5_element_type, f5_write_flag);
      f5->incref();
      f5->writeToFile(filename);
      f5->decref();
    }

    if ((buck || freq) && f5_eigen_frequency > 0 &&
        (iter % f5_eigen_frequency) == 0) {
      char filename[strlen(prefix) + 100];
      TACSToFH5 *f5 =
          new TACSToFH5(assembler, f5_eigen_element_type, f5_eigen_write_flag);
      f5->incref();

      TACSBVec *tmp = assembler->createVec();
      tmp->incref();
      if (buck) {
        for (int i = 0; i < num_load_cases; i++) {
          // Extract the first k eigenvectors for ith load
          for (int k = 0; k < num_buck_eigvals; k++) {
            TacsScalar error;
            buck[i]->extractEigenvector(k, tmp, &error);
            assembler->setVariables(tmp);
            snprintf(filename, sizeof(filename),
                     "%s/load%d_eigenvector%02d_output%d.f5", prefix, i, k,
                     iter);
            f5->writeToFile(filename);
          }
        }
      } else {
        for (int k = 0; k < num_freq_eigvals; k++) {
          TacsScalar error;
          freq->extractEigenvector(k, tmp, &error);
          assembler->setVariables(tmp);
          snprintf(filename, sizeof(filename),
                   "%s/freq_eigenvector%02d_output%d.f5", prefix, k, iter);
          f5->writeToFile(filename);
        }
      }
      tmp->decref();
      f5->decref();
    }
  }

  // Update the iteration count
  iter_count++;
}
