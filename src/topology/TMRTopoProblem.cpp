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
#include "TMR_STLTools.h"
#include "TACSToFH5.h"
#include "TMRHelmholtz.h"
#include "TMR_TACSCreator.h"

/*
  Helmholtz filter creator classes
*/
class TMRQuadTACSHelmholtzCreator : public TMRQuadTACSCreator {
public:
  TMRQuadTACSHelmholtzCreator( double r ):
    TMRQuadTACSCreator(NULL){
    radius = r;
  }
  void createElements( int order,
                       TMRQuadForest *forest,
                       int num_elements,
                       TACSElement **elements ){
    TACSElement *elem = NULL;
    if (order == 2){
      elem = new TMRQuadHelmholtz<2>(radius);
    }
    else if (order == 3){
      elem = new TMRQuadHelmholtz<3>(radius);
    }
    else if (order == 4){
      elem = new TMRQuadHelmholtz<4>(radius);
    }
    else if (order == 5){
      elem = new TMRQuadHelmholtz<5>(radius);
    }
    else if (order == 6){
      elem = new TMRQuadHelmholtz<6>(radius);
    }

    for ( int i = 0; i < num_elements; i++ ){
      elements[i] = elem;
    }
  }

  double radius;
};

class TMROctTACSHelmholtzCreator : public TMROctTACSCreator {
public:
  TMROctTACSHelmholtzCreator( double r ):
    TMROctTACSCreator(NULL){
    radius = r;
  }
  void createElements( int order,
                       TMROctForest *forest,
                       int num_elements,
                       TACSElement **elements ){
    TACSElement *elem = NULL;
    if (order == 2){
      elem = new TMROctHelmholtz<2>(radius);
    }
    else if (order == 3){
      elem = new TMROctHelmholtz<3>(radius);
    }
    else if (order == 4){
      elem = new TMROctHelmholtz<4>(radius);
    }
    else if (order == 5){
      elem = new TMROctHelmholtz<5>(radius);
    }
    else if (order == 6){
      elem = new TMROctHelmholtz<6>(radius);
    }

    for ( int i = 0; i < num_elements; i++ ){
      elements[i] = elem;
    }
  }

  double radius;
};

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
TMRTopoProblem::TMRTopoProblem( int _nlevels,
                                TACSAssembler *_tacs[],
                                TMROctForest *_filter[],
                                TACSVarMap *_filter_maps[],
                                TACSBVecIndices *filter_indices[],
                                TACSMg *_mg,
                                double helmholtz_radius,
                                int _vars_per_node ):
 ParOptProblem(_tacs[0]->getMPIComm()){
  initialize(_nlevels, _tacs, _filter, NULL,
             _filter_maps, filter_indices, _mg,
             helmholtz_radius, _vars_per_node);
}

TMRTopoProblem::TMRTopoProblem( int _nlevels,
                                TACSAssembler *_tacs[],
                                TMRQuadForest *_filter[],
                                TACSVarMap *_filter_maps[],
                                TACSBVecIndices *filter_indices[],
                                TACSMg *_mg,
                                double helmholtz_radius,
                                int _vars_per_node ):
 ParOptProblem(_tacs[0]->getMPIComm()){
  initialize(_nlevels, _tacs, NULL, _filter,
             _filter_maps, filter_indices, _mg,
             helmholtz_radius, _vars_per_node);
}

void TMRTopoProblem::initialize( int _nlevels,
                                 TACSAssembler *_tacs[],
                                 TMROctForest *_oct_filter[],
                                 TMRQuadForest *_quad_filter[],
                                 TACSVarMap *_filter_maps[],
                                 TACSBVecIndices *filter_indices[],
                                 TACSMg *_mg,
                                 double helmholtz_radius,
                                 int _vars_per_node ){
  if (!_oct_filter && !_quad_filter){
    fprintf(stderr,
            "TMRTopoProblem: Quadtree and octree objects not provided.\n");
  }

  // Set the prefix to NULL
  prefix = NULL;

  // Get the processor rank
  int mpi_rank;
  MPI_Comm_rank(_tacs[0]->getMPIComm(), &mpi_rank);

  // Set the number of design variables per node i.e. multi-material or single
  // material
  vars_per_node = _vars_per_node;

  // Set the maximum number of indices
  max_local_size = 0;

  // The initial design variable values (may not be set)
  xinit = NULL;
  xlb = NULL;
  xub = NULL;

  // Set the number of levels
  nlevels = _nlevels;

  // Allocate arrays to store the assembler objects/forests
  tacs = new TACSAssembler*[ nlevels ];
  oct_filter = NULL;
  quad_filter = NULL;
  if (_oct_filter){
    oct_filter = new TMROctForest*[ nlevels ];
  }
  else {
    quad_filter = new TMRQuadForest*[ nlevels ];
  }

  for ( int k = 0; k < nlevels; k++ ){
    // Set the TACSAssembler objects for each level
    tacs[k] = _tacs[k];
    tacs[k]->incref();

    // Set the filter object
    if (_oct_filter){
      oct_filter[k] = _oct_filter[k];
      oct_filter[k]->incref();
    }
    else {
      quad_filter[k] = _quad_filter[k];
      quad_filter[k]->incref();
    }
  }

  // Allocate arrays to store the filter data
  filter_maps = new TACSVarMap*[ nlevels ];
  filter_dist = new TACSBVecDistribute*[ nlevels ];
  filter_dep_nodes = new TACSBVecDepNodes*[ nlevels ];

  // The design variable vector for each level
  x = new TACSBVec*[ nlevels ];

  // Set the helmholtz data to NULL
  helmholtz_tacs = NULL;
  helmholtz_ksm = NULL;
  helmholtz_mg = NULL;
  helmholtz_rhs = NULL;
  helmholtz_psi = NULL;
  helmholtz_vec = NULL;

  int use_helmholtz_filter = 0;
  if (helmholtz_radius > 0.0){
    use_helmholtz_filter = 1;
  }

  if (use_helmholtz_filter){
    // Allocate space for the assembler objects
    helmholtz_tacs = new TACSAssembler*[ nlevels ];

    // Create the Helmholtz creator objects
    TMROctTACSHelmholtzCreator *helmholtz_creator3d = NULL;
    TMRQuadTACSHelmholtzCreator *helmholtz_creator2d = NULL;

    if (oct_filter){
      helmholtz_creator3d = new TMROctTACSHelmholtzCreator(helmholtz_radius);
      helmholtz_creator3d->incref();
    }
    else {
      helmholtz_creator2d = new TMRQuadTACSHelmholtzCreator(helmholtz_radius);
      helmholtz_creator2d->incref();
    }

    for ( int k = 0; k < nlevels; k++ ){
      // Create the assembler object
      if (helmholtz_creator3d){
        helmholtz_tacs[k] =
          helmholtz_creator3d->createTACS(oct_filter[k],
                                          TACSAssembler::NATURAL_ORDER);
      }
      else {
        helmholtz_tacs[k] =
          helmholtz_creator2d->createTACS(quad_filter[k],
                                          TACSAssembler::NATURAL_ORDER);
      }
      helmholtz_tacs[k]->incref();

      // Extract the information that we need from the helmholtz
      // assembler object
      if (_filter_maps && filter_indices){
        filter_maps[k] = _filter_maps[k];
        filter_maps[k]->incref();

        // Create the distribution object for the design variables
        filter_dist[k] = new TACSBVecDistribute(filter_maps[k],
                                                filter_indices[k]);
        filter_dist[k]->incref();
      }
      else {
        filter_maps[k] = helmholtz_tacs[k]->getVarMap();
        filter_maps[k]->incref();

        filter_dist[k] = helmholtz_tacs[k]->getBVecDistribute();
        filter_dist[k]->incref();
      }

      filter_dep_nodes[k] = helmholtz_tacs[k]->getBVecDepNodes();
      if (filter_dep_nodes[k]){
        filter_dep_nodes[k]->incref();
      }

      // Create the design vector (can't use the helmholtz tacs
      // vectors because they have only one design variable per node)
      x[k] = new TACSBVec(filter_maps[k], vars_per_node, filter_dist[k],
                          filter_dep_nodes[k]);
      x[k]->incref();
    }

    // Create a temporary vector to store the design variables
    helmholtz_vec =
      new TACSBVec(filter_maps[0], vars_per_node, filter_dist[0],
                   filter_dep_nodes[0]);
    helmholtz_vec->incref();

    // Destroy the helmholtz creator object
    if (helmholtz_creator3d){
      helmholtz_creator3d->decref();
    }
    else {
      helmholtz_creator2d->decref();
    }

    // Create the vectors
    helmholtz_rhs = helmholtz_tacs[0]->createVec();
    helmholtz_psi = helmholtz_tacs[0]->createVec();
    helmholtz_rhs->incref();
    helmholtz_psi->incref();

    // Create the multigrid object
    double helmholtz_omega = 0.5;
    if (oct_filter){
      TMR_CreateTACSMg(nlevels, helmholtz_tacs, oct_filter,
                       &helmholtz_mg, helmholtz_omega);
    }
    else {
      TMR_CreateTACSMg(nlevels, helmholtz_tacs, quad_filter,
                       &helmholtz_mg, helmholtz_omega);
    }
    helmholtz_mg->incref();

    // Set up the solver
    int gmres_iters = 50;
    int nrestart = 5;
    int is_flexible = 0;

    // Create the GMRES object
    helmholtz_ksm = new GMRES(helmholtz_mg->getMat(0), helmholtz_mg,
                              gmres_iters, nrestart, is_flexible);
    helmholtz_ksm->incref();
    helmholtz_ksm->setMonitor(new KSMPrintStdout("Filter GMRES", mpi_rank, 10));
    helmholtz_ksm->setTolerances(1e-12, 1e-30);

    // Create the multigrid solver and factor it
    double alpha = 1.0, beta = 0.0, gamma = 0.0;
    helmholtz_mg->assembleJacobian(alpha, beta, gamma, NULL);
    helmholtz_mg->factor();

  }
  else {
    // Copy over the assembler objects and filters
    for ( int k = 0; k < nlevels; k++ ){
      filter_maps[k] = _filter_maps[k];
      filter_maps[k]->incref();

      // Create the distribution object for the design variables
      filter_dist[k] = new TACSBVecDistribute(filter_maps[k],
                                              filter_indices[k]);
      filter_dist[k]->incref();

      // By default the dependent nodes are not created unless needed
      filter_dep_nodes[k] = NULL;

      // Create the design variable vector for this level
      if ((oct_filter && oct_filter[k] &&
           oct_filter[k]->getInterpType() == TMR_BERNSTEIN_POINTS) ||
          (quad_filter && quad_filter[k] &&
           quad_filter[k]->getInterpType() == TMR_BERNSTEIN_POINTS)){
        // Extract the dependent node info from the oct filter
        int *dep_ptr, *dep_conn;
        const int *_dep_ptr, *_dep_conn;
        double *dep_weights;
        const double *_dep_weights;
        int num_dep_nodes = 0;
        if (oct_filter){
          num_dep_nodes =
            oct_filter[k]->getDepNodeConn(&_dep_ptr, &_dep_conn,
                                          &_dep_weights);
        }
        else {
          num_dep_nodes =
            quad_filter[k]->getDepNodeConn(&_dep_ptr, &_dep_conn,
                                           &_dep_weights);
        }

        // Copy over the dependent node data
        dep_ptr = new int[ num_dep_nodes+1 ];
        memcpy(dep_ptr, _dep_ptr, (num_dep_nodes+1)*sizeof(int));
        int dep_size = dep_ptr[num_dep_nodes];
        dep_conn = new int[ dep_size ];
        memcpy(dep_conn, _dep_conn, dep_size*sizeof(int));

        dep_weights = new double[ dep_size ];
        memcpy(dep_weights, _dep_weights, dep_size*sizeof(double));

        // Create the dependent node object
        filter_dep_nodes[k] = new TACSBVecDepNodes(num_dep_nodes,
                                                   &dep_ptr, &dep_conn,
                                                   &dep_weights);
        filter_dep_nodes[k]->incref();

        x[k] = new TACSBVec(filter_maps[k], vars_per_node, filter_dist[k],
                            filter_dep_nodes[k]);
      }
      else {
        // Create the design variable vector for this level
        x[k] = new TACSBVec(filter_maps[k], vars_per_node, filter_dist[k]);
      }
      x[k]->incref();
    }
  }

  // Now create the interpolation between filter levels
  filter_interp = new TACSBVecInterp*[ nlevels-1 ];

  for ( int k = 1; k < nlevels; k++ ){
    // Create the interpolation object
    filter_interp[k-1] = new TACSBVecInterp(filter_maps[k],
                                            filter_maps[k-1], vars_per_node);
    filter_interp[k-1]->incref();

    // Create the interpolation on the TMR side
    if (oct_filter){
      oct_filter[k-1]->createInterpolation(oct_filter[k], filter_interp[k-1]);
    }
    else {
      quad_filter[k-1]->createInterpolation(quad_filter[k], filter_interp[k-1]);
    }
    filter_interp[k-1]->initialize();
  }

  // Compute the max filter size
  for ( int k = 0; k < nlevels; k++ ){
    // Set the maximum local size
    const int *range;
    filter_maps[k]->getOwnerRange(&range);

    // Find the number of dependent nodes
    int ndep = 0;
    if (filter_dep_nodes[k]){
      ndep = filter_dep_nodes[k]->getDepNodes(NULL, NULL, NULL);
    }

    // Get the global indices and add its contribution to the
    // max. size of the local design variable vector
    TACSBVecIndices *indices = filter_dist[k]->getIndices();
    int size = vars_per_node*((range[mpi_rank+1] - range[mpi_rank]) +
                              indices->getNumIndices() + ndep);

    // Compute the total number of referenced nodes on this filter object
    int size2 = 0;
    if (oct_filter){
      size2 = vars_per_node*oct_filter[k]->getNodeNumbers(NULL);
    }
    else {
      size2 = vars_per_node*quad_filter[k]->getNodeNumbers(NULL);
    }

    // Update the maximum local size
    if (size > max_local_size){
      max_local_size = size;
    }
    if (size2 > max_local_size){
      max_local_size = size2;
    }
  }

  // Set the maximum local size
  xlocal = new TacsScalar[ max_local_size ];

  // The multigrid object
  mg = _mg;
  mg->incref();

  // Allocate an adjoint and df/du vector
  dfdu = tacs[0]->createVec();
  adjoint = tacs[0]->createVec();
  dfdu->incref();
  adjoint->incref();

  // Set up the solver
  int gmres_iters = 50;
  int nrestart = 5;
  int is_flexible = 0;
  ksm = new GMRES(mg->getMat(0), mg,
                  gmres_iters, nrestart, is_flexible);
  ksm->incref();
  ksm->setMonitor(new KSMPrintStdout("GMRES", mpi_rank, 10));
  ksm->setTolerances(1e-9, 1e-30);

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
}

/*
  Free the data stored in the object
*/
TMRTopoProblem::~TMRTopoProblem(){
  if (prefix){
    delete [] prefix;
  }

  // Decrease the reference counts
  for ( int k = 0; k < nlevels; k++ ){
    tacs[k]->decref();
    if (quad_filter && quad_filter[k]){
      quad_filter[k]->decref();
    }
    if (oct_filter && oct_filter[k]){
      oct_filter[k]->decref();
    }
    if (helmholtz_tacs && helmholtz_tacs[k]){
      helmholtz_tacs[k]->decref();
    }
    filter_maps[k]->decref();
    filter_dist[k]->decref();
    if (filter_dep_nodes[k]){
      filter_dep_nodes[k]->decref();
    }
    x[k]->decref();
  }
  delete [] tacs;
  if (quad_filter){
    delete [] quad_filter;
  }
  if (oct_filter){
    delete [] oct_filter;
  }

  delete [] filter_maps;
  delete [] filter_dist;
  delete [] filter_dep_nodes;
  delete [] x;

  if (helmholtz_tacs){
    delete [] helmholtz_tacs;
    helmholtz_mg->decref();
    helmholtz_ksm->decref();
    helmholtz_vec->decref();
    helmholtz_rhs->decref();
    helmholtz_psi->decref();
  }

  for ( int k = 0; k < nlevels-1; k++ ){
    filter_interp[k]->decref();
  }
  delete [] filter_interp;

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
    vars[i] = tacs[0]->createVec();
    vars[i]->incref();
  }

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
    TACSMat *mmat = tacs[0]->createMat();
    if (use_jd){
      // Create the preconditioner matrix
      TACSMat *kmat = tacs[0]->createMat();
      TACSMat *pcmat = mg->getMat(0);

      // Get preconditioner from Mg
      TACSPc *pc = mg;

      freq = new TACSFrequencyAnalysis(tacs[0], sigma, mmat,
                                       kmat, pcmat, pc,
                                       max_subspace_size, fgmres_size,
                                       num_eigvals, eigtol, eig_rtol,
                                       eig_atol, num_recycle, recycle_type);
    }
    else{
      // Create the frequency analysis object
      freq = new TACSFrequencyAnalysis(tacs[0], sigma, mmat,
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
    MPI_Comm_rank(tacs[0]->getMPIComm(), &mpi_rank);
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
    TACSMat *gmat = tacs[0]->createMat();
    TACSMat *kmat = tacs[0]->createMat();
    TACSMat *aux_mat;
    ksm->getOperators(&aux_mat, NULL);

    buck = new TACSLinearBuckling*[ num_load_cases ];
    buck_ks_sum = new TacsScalar[ num_load_cases ];
    for ( int i = 0; i < num_load_cases; i++ ){
      // Create the buckling analysis object
      buck[i] = new TACSLinearBuckling(tacs[0], sigma, gmat, kmat,
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
  int nvars = x[0]->getArray(NULL);
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
  if (quad_filter){
    return new ParOptBVecWrap(new TACSBVec(filter_maps[0],
                                           vars_per_node,
                                           filter_dist[0],
                                           filter_dep_nodes[0]));
  }
  else {
    return new ParOptBVecWrap(new TACSBVec(filter_maps[0],
                                           vars_per_node,
                                           filter_dist[0],
                                           filter_dep_nodes[0]));
  }
}

/*
  Compute the volume corresponding to each node within the filter
*/
TACSBVec* TMRTopoProblem::createVolumeVec( double Xscale ){
  // Get the dependent nodes and weight values
  const int *dep_ptr, *dep_conn;
  const double *dep_weights;
  int ndep = oct_filter[0]->getDepNodeConn(&dep_ptr, &dep_conn,
                                           &dep_weights);

  // Copy over the data
  int *dptr = new int[ ndep+1 ];
  int *dconn = new int[ dep_ptr[ndep] ];
  double *dweights = new double[ dep_ptr[ndep] ];
  memcpy(dptr, dep_ptr, (ndep+1)*sizeof(int));
  memcpy(dconn, dep_conn, dep_ptr[ndep]*sizeof(int));
  memcpy(dweights, dep_weights, dep_ptr[ndep]*sizeof(double));
  TACSBVecDepNodes *dep_nodes = new TACSBVecDepNodes(ndep, &dptr,
						     &dconn, &dweights);

  TACSBVec *vec = new TACSBVec(filter_maps[0], 1, filter_dist[0], dep_nodes);
  vec->zeroEntries();

  // Get the octants
  TMROctantArray *octants;
  oct_filter[0]->getOctants(&octants);

  // Get the array of octants
  int size;
  TMROctant *array;
  octants->getArray(&array, &size);

  // Get the nodes
  const int *conn;
  oct_filter[0]->getNodeConn(&conn);

  // Get the node locations from the filter
  TMRPoint *X;
  oct_filter[0]->getPoints(&X);

  // Allocate the memory required
  int order = oct_filter[0]->getMeshOrder();
  int num_nodes = order*order*order;
  double *N = new double[ num_nodes ];
  double *Na = new double[ num_nodes ];
  double *Nb = new double[ num_nodes ];
  double *Nc = new double[ num_nodes ];
  TMRPoint *Xpts = new TMRPoint[ num_nodes ];
  TacsScalar *area = new TacsScalar[ num_nodes ];

  // Get the quadrature points/weights
  const double *quadPts;
  const double *quadWts;
  int npts = FElibrary::getGaussPtsWts(order, &quadPts, &quadWts);

  // Loop over the elements
  for ( int i = 0; i < size; i++ ){
    // Retrieve the node numbers and x/y/z locations for element i
    // within the mesh
    for ( int kk = 0; kk < order; kk++ ){
      for ( int jj = 0; jj < order; jj++ ){
        for ( int ii = 0; ii < order; ii++ ){
          const int offset = ii + jj*order + kk*order*order;
          const int node = conn[num_nodes*i + offset];
          int index = oct_filter[0]->getLocalNodeNumber(node);

          // Copy the node location
          Xpts[offset] = X[index];
        }
      }
    }

    // Set the area
    memset(area, 0, num_nodes*sizeof(TacsScalar));

    for ( int kk = 0; kk < npts; kk++ ){
      for ( int jj = 0; jj < npts; jj++ ){
        for ( int ii = 0; ii < npts; ii++ ){
          double pt[3];
          pt[0] = quadPts[ii];
          pt[1] = quadPts[jj];
          pt[2] = quadPts[kk];

          // Compute the quadrature weight
          double wt = quadWts[ii]*quadWts[jj]*quadWts[kk];

          // Evaluate the derivative of the shape functions
          oct_filter[0]->evalInterp(pt, N, Na, Nb, Nc);

          // Compute the Jacobian transformation
          TacsScalar J[9];
          memset(J, 0, 9*sizeof(TacsScalar));
          for ( int j = 0; j < 8; j++ ){
            J[0] += Na[j]*Xpts[j].x*Xscale;
            J[1] += Nb[j]*Xpts[j].x*Xscale;
            J[2] += Nc[j]*Xpts[j].x*Xscale;

            J[3] += Na[j]*Xpts[j].y*Xscale;
            J[4] += Nb[j]*Xpts[j].y*Xscale;
            J[5] += Nc[j]*Xpts[j].y*Xscale;

            J[6] += Na[j]*Xpts[j].z*Xscale;
            J[7] += Nb[j]*Xpts[j].z*Xscale;
            J[8] += Nc[j]*Xpts[j].z*Xscale;
          }

          // Add the determinant to the area - the weights in this
          // case are just = 1.0
          TacsScalar det = FElibrary::jacobian3d(J);
          for ( int j = 0; j < num_nodes; j++ ){
            area[j] += wt*det*N[j];
          }
        }
      }
    }

    // Add the values to the vector
    vec->setValues(num_nodes, &conn[num_nodes*i], area, TACS_ADD_VALUES);
  }

  // Free the element-related data
  delete [] N;
  delete [] Na;
  delete [] Nb;
  delete [] Nc;
  delete [] Xpts;
  delete [] area;

  vec->beginSetValues(TACS_ADD_VALUES);
  vec->endSetValues(TACS_ADD_VALUES);

  return vec;
}

/*
  Compute the volume corresponding to each node within the filter
*/
TACSBVec* TMRTopoProblem::createAreaVec( double Xscale ){
  // Get the dependent nodes and weight values
  const int *dep_ptr, *dep_conn;
  const double *dep_weights;
  int ndep = quad_filter[0]->getDepNodeConn(&dep_ptr, &dep_conn,
                                            &dep_weights);

  // Copy over the data
  int *dptr = new int[ ndep+1 ];
  int *dconn = new int[ dep_ptr[ndep] ];
  double *dweights = new double[ dep_ptr[ndep] ];
  memcpy(dptr, dep_ptr, (ndep+1)*sizeof(int));
  memcpy(dconn, dep_conn, dep_ptr[ndep]*sizeof(int));
  memcpy(dweights, dep_weights, dep_ptr[ndep]*sizeof(double));
  TACSBVecDepNodes *dep_nodes = new TACSBVecDepNodes(ndep, &dptr,
						     &dconn, &dweights);

  TACSBVec *vec = new TACSBVec(filter_maps[0], 1, filter_dist[0], dep_nodes);
  vec->zeroEntries();

  // Get the quadrants
  TMRQuadrantArray *quadrants;
  quad_filter[0]->getQuadrants(&quadrants);

  // Get the array of quadrants
  int size;
  TMRQuadrant *array;
  quadrants->getArray(&array, &size);

  // Get the nodes
  const int *conn;
  quad_filter[0]->getNodeConn(&conn);

  // Get the node locations from the filter
  TMRPoint *X;
  quad_filter[0]->getPoints(&X);

  // Allocate the memory required
  int order = quad_filter[0]->getMeshOrder();
  int num_nodes = order*order;
  double *N = new double[ num_nodes ];
  double *Na = new double[ num_nodes ];
  double *Nb = new double[ num_nodes ];
  TMRPoint *Xpts = new TMRPoint[ num_nodes ];
  TacsScalar *area = new TacsScalar[ num_nodes ];

  // Get the quadrature points/weights
  const double *quadPts;
  const double *quadWts;
  int npts = FElibrary::getGaussPtsWts(order, &quadPts, &quadWts);

  // Loop over the elements
  for ( int i = 0; i < size; i++ ){
    // Retrieve the node numbers and x/y/z locations for element i
    // within the mesh
    for ( int jj = 0; jj < order; jj++ ){
      for ( int ii = 0; ii < order; ii++ ){
        const int offset = ii + jj*order;
        const int node = conn[num_nodes*i + offset];
        int index = quad_filter[0]->getLocalNodeNumber(node);

        // Copy the node location
        Xpts[offset] = X[index];
      }
    }

    // Set the area
    memset(area, 0, num_nodes*sizeof(TacsScalar));

    for ( int jj = 0; jj < npts; jj++ ){
      for ( int ii = 0; ii < npts; ii++ ){
        double pt[2];
        pt[0] = quadPts[ii];
        pt[1] = quadPts[jj];

        // Compute the quadrature weight
        double wt = quadWts[ii]*quadWts[jj];

        // Evaluate the derivative of the shape functions
        quad_filter[0]->evalInterp(pt, N, Na, Nb);

        // Compute the Jacobian transformation
        TacsScalar J[4];
        memset(J, 0, 4*sizeof(TacsScalar));
        for ( int j = 0; j < 4; j++ ){
          J[0] += Na[j]*Xpts[j].x*Xscale;
          J[1] += Nb[j]*Xpts[j].x*Xscale;

          J[2] += Na[j]*Xpts[j].y*Xscale;
          J[3] += Nb[j]*Xpts[j].y*Xscale;
        }

        // Add the determinant to the area - the weights in this
        // case are just = 1.0
        TacsScalar det = FElibrary::jacobian2d(J);
        for ( int j = 0; j < num_nodes; j++ ){
          area[j] += wt*det*N[j];
        }
      }
    }

    // Add the values to the vector
    vec->setValues(num_nodes, &conn[num_nodes*i], area, TACS_ADD_VALUES);
  }

  // Free the element-related data
  delete [] N;
  delete [] Na;
  delete [] Nb;
  delete [] Xpts;
  delete [] area;

  vec->beginSetValues(TACS_ADD_VALUES);
  vec->endSetValues(TACS_ADD_VALUES);

  return vec;
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
  Set the initial design variables
*/
void TMRTopoProblem::getVarsAndBounds( ParOptVec *xvec,
                                       ParOptVec *lbvec,
                                       ParOptVec *ubvec ){
  int mpi_rank;
  MPI_Comm_rank(tacs[0]->getMPIComm(), &mpi_rank);
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
        tacs[0]->getDesignVars(xlocal, max_local_size);

        // Set the local values into the vector
        setBVecFromLocalValues(0, xlocal, wrap->vec, TACS_INSERT_VALUES);

        // Insert the values across all processors
        wrap->vec->beginDistributeValues();
        wrap->vec->endDistributeValues();
        wrap->vec->beginSetValues(TACS_INSERT_VALUES);
        wrap->vec->endSetValues(TACS_INSERT_VALUES);        
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
      tacs[0]->getDesignVarRange(xlocal, upper, max_local_size);

      if (!has_lower){
        ParOptBVecWrap *lbwrap = dynamic_cast<ParOptBVecWrap*>(lbvec);
        if (lbwrap){
          lbwrap->vec->zeroEntries();
          // Set the values from the local array
          setBVecFromLocalValues(0, xlocal, lbwrap->vec, TACS_INSERT_VALUES);
          lbwrap->vec->beginDistributeValues();
          lbwrap->vec->endDistributeValues();
          lbwrap->vec->beginSetValues(TACS_INSERT_VALUES);
          lbwrap->vec->endSetValues(TACS_INSERT_VALUES);
        }
      }
      if (!has_upper){        
        ParOptBVecWrap *ubwrap = dynamic_cast<ParOptBVecWrap*>(ubvec);
        if (ubwrap){
          ubwrap->vec->zeroEntries();
          setBVecFromLocalValues(0, upper, ubwrap->vec, TACS_INSERT_VALUES);
          ubwrap->vec->beginDistributeValues();
          ubwrap->vec->endDistributeValues();
          ubwrap->vec->beginSetValues(TACS_INSERT_VALUES);
          ubwrap->vec->endSetValues(TACS_INSERT_VALUES);
        }
      }
      delete [] upper;
    }
  }
}

/*
  Get the local values of the design variables from the TACSBVec
  object and set them in xlocal.

  This function requires that the values already be distributed
  between processors so that the external values (owned by other
  processors) are correct.

  input:
  vec:     the design variable vector

  output:
  xlocal:  the local design variable array

  returns: the number of design variables on this processor
*/
int TMRTopoProblem::getLocalValuesFromBVec( int level,
                                            TACSBVec *vec,
                                            TacsScalar *xloc ){
  // If it is for a quadtree problem
  if (quad_filter &&
      quad_filter[level]->getInterpType() == TMR_BERNSTEIN_POINTS){
    const int *node_numbers;
    int num_nodes = quad_filter[level]->getNodeNumbers(&node_numbers);
    vec->getValues(num_nodes, node_numbers, &xloc[0]);
    return num_nodes;
  }
  else if (oct_filter &&
           oct_filter[level]->getInterpType() == TMR_BERNSTEIN_POINTS){
    const int *node_numbers;
    int num_nodes = oct_filter[level]->getNodeNumbers(&node_numbers);
    vec->getValues(num_nodes, node_numbers, &xloc[0]);
    return num_nodes;
  }
  else {
    TacsScalar *x_vals, *x_ext_vals;
    int size = vec->getArray(&x_vals);
    int ext_size = vec->getExtArray(&x_ext_vals);

    memcpy(xloc, x_vals, size*sizeof(TacsScalar));
    if (x_ext_vals){
      memcpy(&xloc[size], x_ext_vals, ext_size*sizeof(TacsScalar));
    }
    return size + ext_size;
  }
}

/*
  Set the values into the TACSBVec object using the local
  size
*/
void TMRTopoProblem::setBVecFromLocalValues( int level,
                                             const TacsScalar *xloc,
                                             TACSBVec *vec,
                                             TACSBVecOperation op ){
  if (quad_filter &&
      quad_filter[level]->getInterpType() == TMR_BERNSTEIN_POINTS){
    const int *node_numbers;
    int num_nodes = quad_filter[level]->getNodeNumbers(&node_numbers);
    vec->setValues(num_nodes, node_numbers, xloc, op);
  }
  else if (oct_filter &&
           oct_filter[level]->getInterpType() == TMR_BERNSTEIN_POINTS){
    const int *node_numbers;
    int num_nodes = oct_filter[level]->getNodeNumbers(&node_numbers);
    vec->setValues(num_nodes, node_numbers, xloc, op);
  }
  else {
    TacsScalar *x_vals, *x_ext_vals;
    int size = vec->getArray(&x_vals);
    int ext_size = vec->getExtArray(&x_ext_vals);
    memcpy(x_vals, xloc, size*sizeof(TacsScalar));
    if (x_ext_vals){
      memcpy(x_ext_vals, &xloc[size], ext_size*sizeof(TacsScalar));
    }
  }
}

/*
  Apply the Helmholtz filter to the design variables.

  Here the input/output vector are the same
*/
void TMRTopoProblem::applyHelmholtzFilter( TACSBVec *xvars ){
  if (helmholtz_tacs){
    // Get the number of local elements
    int num_elements = helmholtz_tacs[0]->getNumElements();

    // Get the maximum number of nodes per element (all elements in
    // this assembler object have the same number of nodes)
    int max_nodes = helmholtz_tacs[0]->getMaxElementNodes();

    // Get all of the elements in the filter
    TACSElement **elements = helmholtz_tacs[0]->getElements();

    // Allocate space for element-level data
    double *N = new double[ max_nodes ];
    TacsScalar *Xpts = new TacsScalar[ 3*max_nodes ];
    TacsScalar *x_values = new TacsScalar[ vars_per_node*max_nodes ];
    TacsScalar *rhs_values = new TacsScalar[ max_nodes ];

    for ( int k = 0; k < vars_per_node; k++ ){
      // Zero the entries in the RHS vector
      helmholtz_rhs->zeroEntries();
      helmholtz_tacs[0]->zeroVariables();
      
      for ( int i = 0; i < num_elements; i++ ){
        // Get the values for this element
        int len;
        const int *nodes;
        helmholtz_tacs[0]->getElement(i, &nodes, &len);
        helmholtz_tacs[0]->getElement(i, Xpts);

        // Get the values of the design variables at the nodes
        xvars->getValues(len, nodes, x_values);

        // Zero the values on the right-hand-side
        memset(rhs_values, 0, max_nodes*sizeof(TacsScalar));

        // Get the number of quadrature points
        const int num_quad_pts = elements[i]->getNumGaussPts();

        // Perform the integration over the element
        for ( int n = 0; n < num_quad_pts; n++ ){
          double pt[3];
          double wt = elements[i]->getGaussWtsPts(n, pt);
          TacsScalar h = elements[i]->getDetJacobian(pt, Xpts);
          h *= wt;

          elements[i]->getShapeFunctions(pt, N);

          for ( int ii = 0; ii < max_nodes; ii++ ){
            const double *ns = N;
            const TacsScalar *xs = &x_values[k];
            TacsScalar v = 0.0;
            for ( int jj = 0; jj < max_nodes; jj++ ){
              v += ns[0]*xs[0];
              ns++;
              xs += vars_per_node;
            }
            rhs_values[ii] += h*N[ii]*v;
          }
        }

        helmholtz_rhs->setValues(len, nodes, rhs_values, TACS_ADD_VALUES);
      }

      // Complete the assembly process
      helmholtz_rhs->beginSetValues(TACS_ADD_VALUES);
      helmholtz_rhs->endSetValues(TACS_ADD_VALUES);
      
      // Solve for the filtered values of the design variables
      helmholtz_ksm->solve(helmholtz_rhs, helmholtz_psi);
      helmholtz_tacs[0]->reorderVec(helmholtz_psi);
      helmholtz_tacs[0]->setVariables(helmholtz_psi);
  
      // Get the output array
      TacsScalar *xarr;
      xvars->getArray(&xarr);

      // Get the Helmholtz solution vector
      TacsScalar *hpsi;
      const int size = helmholtz_psi->getArray(&hpsi);

      for ( int i = 0; i < size; i++ ){
        xarr[vars_per_node*i + k] = hpsi[i];
      }      
    }
    iter_count++;

    delete [] N;
    delete [] Xpts;
    delete [] x_values;
    delete [] rhs_values;
  }
}

/*
  Compute the sensitivity w.r.t the Helmholtz filter
*/
void TMRTopoProblem::reverseHelmholtzFilter( TACSBVec *input,
                                             TACSBVec *output ){
  if (helmholtz_tacs){
    output->zeroEntries();
    // Get the number of local elements
    int num_elements = helmholtz_tacs[0]->getNumElements();

    // Get the maximum number of nodes per element (all elements in
    // this assembler object have the same number of nodes)
    int max_nodes = helmholtz_tacs[0]->getMaxElementNodes();

    // Get all of the elements in the filter
    TACSElement **elements = helmholtz_tacs[0]->getElements();

    // Allocate space for element-level data
    double *N = new double[ max_nodes ];
    TacsScalar *Xpts = new TacsScalar[ 3*max_nodes ];
    TacsScalar *x_values = new TacsScalar[ vars_per_node*max_nodes ];
    TacsScalar *psi_values = new TacsScalar[ max_nodes ];

    for ( int k = 0; k < vars_per_node; k++ ){
      // Get the input array
      TacsScalar *xarr;
      input->getArray(&xarr);

      // Get the Helmholtz solution vector
      TacsScalar *hrhs;
      const int size = helmholtz_rhs->getArray(&hrhs);

      for ( int i = 0; i < size; i++ ){
        hrhs[i] = xarr[vars_per_node*i + k];
      }
      
      // Solve for the filtered values of the design variables
      helmholtz_ksm->solve(helmholtz_rhs, helmholtz_psi);
      helmholtz_tacs[0]->reorderVec(helmholtz_psi);

      // Distribute the values from the solution
      helmholtz_psi->beginDistributeValues();
      helmholtz_psi->endDistributeValues();

      for ( int i = 0; i < num_elements; i++ ){
        // Get the values for this element
        int len;
        const int *nodes;
        helmholtz_tacs[0]->getElement(i, &nodes, &len);
        helmholtz_tacs[0]->getElement(i, Xpts);

        // Get the values of the design variables at the nodes
        helmholtz_psi->getValues(len, nodes, psi_values);

        // Zero the values on the right-hand-side
        memset(x_values, 0, vars_per_node*max_nodes*sizeof(TacsScalar));

        // Get the number of quadrature points
        const int num_quad_pts = elements[i]->getNumGaussPts();

        // Perform the integration over the element
        for ( int n = 0; n < num_quad_pts; n++ ){
          double pt[3];
          double wt = elements[i]->getGaussWtsPts(n, pt);
          TacsScalar h = elements[i]->getDetJacobian(pt, Xpts);
          h *= wt;

          elements[i]->getShapeFunctions(pt, N);

          for ( int ii = 0; ii < max_nodes; ii++ ){
            const double *ns = N;
            const TacsScalar *psi = &psi_values[0];
            TacsScalar v = 0.0;
            for ( int jj = 0; jj < max_nodes; jj++ ){
              v += ns[0]*psi[0];
              ns++;
              psi++;
            }

            x_values[vars_per_node*ii + k] += h*N[ii]*v;
          }
        }

        output->setValues(len, nodes, x_values, TACS_ADD_VALUES);
      }
    }

    delete [] N;
    delete [] Xpts;
    delete [] x_values;
    delete [] psi_values;

    output->beginSetValues(TACS_ADD_VALUES);
    output->endSetValues(TACS_ADD_VALUES);
  }
}

/*
  Set the design variable values consistently across all multigrid levels by
  interpolating the design variables using the filter interpolation objects.
*/
void TMRTopoProblem::setDesignVars( ParOptVec *pxvec ){
  ParOptBVecWrap *wrap = dynamic_cast<ParOptBVecWrap*>(pxvec);

  if (wrap){
    // Extract the TACSBVec class
    TACSBVec *xvec = wrap->vec;

    // Copy the values to the local design variable vector
    x[0]->copyValues(xvec);

    // Distribute the design variable values
    x[0]->beginDistributeValues();
    x[0]->endDistributeValues();

    if (helmholtz_tacs){
      applyHelmholtzFilter(x[0]);
    }

    // Copy the values to the local array
    int size = getLocalValuesFromBVec(0, x[0], xlocal);
    tacs[0]->setDesignVars(xlocal, size);

    // Set the design variable values on all processors
    for ( int k = 0; k < nlevels-1; k++ ){
      filter_interp[k]->multWeightTranspose(x[k], x[k+1]);
      // Distribute the design variable values
      x[k+1]->beginDistributeValues();
      x[k+1]->endDistributeValues();

      // Set the design variable values
      size = getLocalValuesFromBVec(k+1, x[k+1], xlocal);
      tacs[k+1]->setDesignVars(xlocal, size);
    }
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
  MPI_Comm_rank(tacs[0]->getMPIComm(), &mpi_rank);

  // Set the design variable values on all mesh levels
  setDesignVars(pxvec);

  // Zero the variables
  tacs[0]->zeroVariables();

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
      ksm->solve(forces[i], vars[i]);
      tacs[0]->setBCs(vars[i]);

      // Set the variables into TACSAssembler
      tacs[0]->setVariables(vars[i]);

      // Add the contribution to the objective
      if (obj_funcs){
        if (obj_funcs[i]){
          TacsScalar fobj_val;
          tacs[0]->evalFunctions(&obj_funcs[i], 1, &fobj_val);
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
        tacs[0]->evalFunctions(load_case_info[i].funcs,
                               num_funcs, &cons[count]);
        if (mpi_rank == 0){
          printf("Mass: %e\n", cons[count]);
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
  MPI_Comm_rank(tacs[0]->getMPIComm(), &mpi_rank);

  // Evaluate the derivative of the weighted compliance with
  // respect to the design variables
  gvec->zeroEntries();
  ParOptBVecWrap *wrap = dynamic_cast<ParOptBVecWrap*>(gvec);
  if (wrap){
    TACSBVec *g = wrap->vec, *gB = wrap->vec;
    if (helmholtz_vec){
      g = helmholtz_vec;
      gB = wrap->vec;
    }
    g->zeroEntries();

    // Evaluate the gradient of the objective - weighted sum of
    // compliances
    memset(xlocal, 0, max_local_size*sizeof(TacsScalar));
    if (obj_funcs){
      for ( int i = 0; i < num_load_cases; i++ ){
        tacs[0]->setVariables(vars[i]);
        int use_adjoint = 1;
        if (dynamic_cast<TACSStructuralMass*>(obj_funcs[i])){
          use_adjoint = 0;
        }

        if (use_adjoint){
          dfdu->zeroEntries();
          double alpha = 1.0, beta = 0.0, gamma = 0.0;
          mg->assembleJacobian(alpha, beta, gamma, NULL, TRANSPOSE);
          mg->factor();
          tacs[0]->addSVSens(alpha, beta, gamma, &obj_funcs[i], 1, &dfdu);
          tacs[0]->applyBCs(dfdu);

          // Solve the system of adjoint equations
          ksm->solve(dfdu, adjoint);
          tacs[0]->addDVSens(obj_weights[i], &obj_funcs[i],
                             1, xlocal, max_local_size);
          tacs[0]->addAdjointResProducts(-obj_weights[i], &adjoint,
                                         1, xlocal, max_local_size);
        }
        else {
          tacs[0]->addDVSens(obj_weights[i], &obj_funcs[i], 1, xlocal,
                             max_local_size);
        }
        setBVecFromLocalValues(0, xlocal, g, TACS_ADD_VALUES);
        g->beginSetValues(TACS_ADD_VALUES);
        g->endSetValues(TACS_ADD_VALUES);
      }

      // Apply the filter sensitivity
      if (helmholtz_vec){
        reverseHelmholtzFilter(g, gB);
      }
    }
    else { // For compliance objective
      for ( int i = 0; i < num_load_cases; i++ ){
        tacs[0]->setVariables(vars[i]);
        tacs[0]->addAdjointResProducts(-obj_weights[i], &vars[i],
                                       1, xlocal, max_local_size);
      }

      setBVecFromLocalValues(0, xlocal, g, TACS_ADD_VALUES);
      g->beginSetValues(TACS_ADD_VALUES);
      g->endSetValues(TACS_ADD_VALUES);

      // Apply the filter sensitivity
      if (helmholtz_vec){
        reverseHelmholtzFilter(g, gB);
      }
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
    tacs[0]->setVariables(vars[i]);

    // Get the number of functions for each load case
    int num_funcs = load_case_info[i].num_funcs;

    for ( int j = 0; j < num_funcs; j++ ){
      TACSFunction *func = load_case_info[i].funcs[j];
      TacsScalar scale = load_case_info[i].scale[j];

      // Try to unwrap the vector
      wrap = dynamic_cast<ParOptBVecWrap*>(Acvec[count + j]);

      if (wrap){
        // Set up the vectors so that it works for the cases with and
        // without the Helmholtz filter
        TACSBVec *A = wrap->vec, *B = NULL;
        if (helmholtz_vec){
          A = helmholtz_vec;
          B = wrap->vec;
        }
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
          tacs[0]->addSVSens(alpha, beta, gamma, &func, 1, &dfdu);
          tacs[0]->applyBCs(dfdu);

          // Solve the system of equations
          ksm->solve(dfdu, adjoint);

          // Compute the total derivative using the adjoint
          memset(xlocal, 0, max_local_size*sizeof(TacsScalar));
          tacs[0]->addDVSens(scale, &func, 1, xlocal, max_local_size);
          tacs[0]->addAdjointResProducts(-scale, &adjoint,
                                         1, xlocal, max_local_size);
        }
        else {
          memset(xlocal, 0, max_local_size*sizeof(TacsScalar));
          tacs[0]->addDVSens(scale, &func, 1, xlocal, max_local_size);
        }

        // Assemble the gradient across all procs
        setBVecFromLocalValues(0, xlocal, A, TACS_ADD_VALUES);
        A->beginSetValues(TACS_ADD_VALUES);
        A->endSetValues(TACS_ADD_VALUES);

        // Apply the filter sensitivity
        if (helmholtz_vec){
          reverseHelmholtzFilter(A, B);
        }
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
        TACSBVec *A = wrap->vec, *B = wrap->vec;
        if (helmholtz_vec){
          A = helmholtz_vec;
          B = wrap->vec;
        }

        // Evaluate the partial derivatives required for the adjoint
        load_case_info[i].stress_func->evalConDeriv(xlocal,
                                                    max_local_size,
                                                    dfdu);
        tacs[0]->applyBCs(dfdu);

        // Solve the system of equations
        ksm->solve(dfdu, adjoint);

        // Compute the total derivative using the adjoint
        tacs[0]->addAdjointResProducts(-1.0, &adjoint,
                                       1, xlocal, max_local_size);

        // Wrap the vector class
        setBVecFromLocalValues(0, xlocal, A, TACS_ADD_VALUES);
        A->beginSetValues(TACS_ADD_VALUES);
        A->endSetValues(TACS_ADD_VALUES);

        // Apply the filter sensitivity
        if (helmholtz_vec){
          reverseHelmholtzFilter(A, B);
        }

        // Scale the constraint by -1 since the constraint is
        // formulated as 1 - c(x, u) > 0.0
        B->scale(-load_case_info[i].stress_func_scale);
      }

      count++;
    }

    if (freq && i == 0){
      // Try to unwrap the vector
      wrap = dynamic_cast<ParOptBVecWrap*>(Acvec[count]);
      if (wrap){
        // Get the underlying TACS vector for the design variables
        TACSBVec *A = wrap->vec, *B = wrap->vec;
        if (helmholtz_vec){
          A = helmholtz_vec;
          B = wrap->vec;
        }

        // Zero the local components
        memset(xlocal, 0, max_local_size*sizeof(TacsScalar));

        TacsScalar *temp = new TacsScalar[max_local_size];
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

        // Add the values for each constraint gradient
        setBVecFromLocalValues(0, xlocal, A, TACS_ADD_VALUES);
        A->beginSetValues(TACS_ADD_VALUES);
        A->endSetValues(TACS_ADD_VALUES);

        // Apply the filter sensitivity
        if (helmholtz_vec){
          reverseHelmholtzFilter(A, B);
        }
      }

      count++;
    }
    if (buck){
      // Compute the derivative of the buckling functions
      // Try to unwrap the vector
      wrap = dynamic_cast<ParOptBVecWrap*>(Acvec[count]);
      if (wrap){
        // Get the vector
        TACSBVec *A = wrap->vec, *B = wrap->vec;
        if (helmholtz_vec){
          A = helmholtz_vec;
          B = wrap->vec;
        }

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

        // Add the values for each constraint gradient
        setBVecFromLocalValues(0, xlocal, A, TACS_ADD_VALUES);
        A->beginSetValues(TACS_ADD_VALUES);
        A->endSetValues(TACS_ADD_VALUES);

        // Apply the filter sensitivity
        if (helmholtz_vec){
          reverseHelmholtzFilter(A, B);
        }
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
  if (prefix && oct_filter){
    // Write out the file at a cut off of 0.25
    char *filename = new char[ strlen(prefix) + 100 ];

    for ( int k = 0; k < vars_per_node; k++ ){
      double cutoff = 0.5;
      sprintf(filename, "%s/levelset05_var%d_binary%04d.bstl",
              prefix, k, iter_count);
      TMR_GenerateBinFile(filename, oct_filter[0], x[0], k, cutoff);
    }

    delete [] filename;
  }  
  else if (prefix && quad_filter){
    // Write out the file at a cut off of 0.25
    char *filename = new char[ strlen(prefix) + 100 ];

    sprintf(filename, "%s/tacs_output%04d.f5",
            prefix, iter_count);
    // Create the visualization for the object
    unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                               TACSElement::OUTPUT_DISPLACEMENTS |
                               TACSElement::OUTPUT_EXTRAS);
    TACSToFH5 *f5 = new TACSToFH5(tacs[0], TACS_PLANE_STRESS,
                                  write_flag);
    f5->incref();
    f5->writeToFile(filename);
    f5->decref();
    delete [] filename;
  }

  if ((buck || freq) && iter_count % 250 == 0){
    writeEigenVector(iter);
  }

  if (iter_count % 50 == 0){
    // Get the processor rank
    int mpi_rank;
    MPI_Comm_rank(tacs[0]->getMPIComm(), &mpi_rank);

    // Allocate the filename array
    char *filename = new char[ strlen(prefix) + 100 ];

    for ( int i = 0; i < num_load_cases; i++ ){
      if (load_case_info[i].stress_func){
        sprintf(filename, "%s/load%d_stress_proc%d.dat",
                prefix, i, mpi_rank);
        load_case_info[i].stress_func->writeReconToTec(vars[i],
                                                       filename, 1.0);
      }
    }

    // Free the filename array
    delete [] filename;
  }

  // Update the iteration count
  iter_count++;
}

/*
  Write the eigenvector file
*/
void TMRTopoProblem::writeEigenVector( int iter ){
  // Only valid if buckling is used
  if (buck){
    char outfile[256];
    TACSBVec *tmp = tacs[0]->createVec();
    tmp->incref();
    tmp->zeroEntries();
    // Create the visualization for the object
    unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                               TACSElement::OUTPUT_DISPLACEMENTS);
    TACSToFH5 *f5 = new TACSToFH5(tacs[0], TACS_SOLID,
                                  write_flag);
    f5->incref();
    for ( int i = 0; i < num_load_cases; i++ ){
      // Extract the first k eigenvectors for ith load
      for ( int k = 0; k < num_buck_eigvals; k++ ){
        TacsScalar error;
        buck[i]->extractEigenvector(k, tmp, &error);
        tacs[0]->setVariables(tmp);
        sprintf(outfile, "%s/load%d_eigenvector%02d_output%d.f5",
                prefix, i, k, iter);
        f5->writeToFile(outfile);
      }
    }
    f5->decref();
    tmp->decref();
  }
  else {
    char outfile[256];
    TACSBVec *tmp = tacs[0]->createVec();
    tmp->incref();
    tmp->zeroEntries();
    // Create the visualization for the object
    unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                               TACSElement::OUTPUT_DISPLACEMENTS);
    TACSToFH5 *f5 = new TACSToFH5(tacs[0], TACS_SOLID,
                                  write_flag);
    f5->incref();
    // Extract the first k eigenvectors for ith load
    for ( int k = 0; k < num_freq_eigvals; k++ ){
      TacsScalar error;
      freq->extractEigenvector(k, tmp, &error);
      tacs[0]->setVariables(tmp);
      sprintf(outfile, "%s/freq_eigenvector%02d_output%d.f5",
	      prefix, k, iter);
      f5->writeToFile(outfile);
    }

    f5->decref();
    tmp->decref();
  }
}
