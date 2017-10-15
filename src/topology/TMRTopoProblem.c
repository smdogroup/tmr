#include "TMRTopoProblem.h"
#include "TMROctStiffness.h"
#include "TACSFunction.h"
#include "Solid.h"
#include "TMR_STLTools.h"
#include "TACSToFH5.h"

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
                                TACSBVecIndices *_filter_indices[],
                                TACSMg *_mg,
                                int _vars_per_node ):
 ParOptProblem(_tacs[0]->getMPIComm()){
  // Set the prefix to NULL
  prefix = NULL;

  // Get the processor rank
  int mpi_rank;
  MPI_Comm_rank(_tacs[0]->getMPIComm(), &mpi_rank);
  
  // Set the number of variables per node
  vars_per_node = _vars_per_node;

  // Set the maximum number of indices
  max_local_size = 0;

  // The initial design variable values (may not be set)
  xinit = NULL;

  // Set the number of levels
  nlevels = _nlevels;

  // Allocate the arrays
  tacs = new TACSAssembler*[ nlevels ];
  oct_filter = new TMROctForest*[ nlevels ];
  quad_filter = NULL;
  filter_maps = new TACSVarMap*[ nlevels ];
  filter_indices = new TACSBVecIndices*[ nlevels ];
  filter_dist = new TACSBVecDistribute*[ nlevels ];
  filter_ctx = new TACSBVecDistCtx*[ nlevels ];
  filter_interp = new TACSBVecInterp*[ nlevels-1 ];

  // The design variable vector for each level
  x = new TACSBVec*[ nlevels ];

  // Copy over the assembler objects and filters
  for ( int k = 0; k < nlevels; k++ ){
    // Set the TACSAssembler objects for each level
    tacs[k] = _tacs[k];
    tacs[k]->incref();

    // Set the filter object
    oct_filter[k] = _filter[k];
    oct_filter[k]->incref();
 
    // Copy over the filter information
    filter_maps[k] = _filter_maps[k];
    filter_maps[k]->incref();

    filter_indices[k] = _filter_indices[k];
    filter_indices[k]->incref();

    // Set the maximum local size
    const int *range;
    filter_maps[k]->getOwnerRange(&range);
    int size = vars_per_node*((range[mpi_rank+1] - range[mpi_rank]) +
                              filter_indices[k]->getNumIndices());

    // Update the maximum local size
    if (size > max_local_size){
      max_local_size = size;
    }

    // Create the distribution object for the design variables
    filter_dist[k] = new TACSBVecDistribute(filter_maps[k], 
                                            filter_indices[k]);
    filter_dist[k]->incref();

    // Create the transfer context
    filter_ctx[k] = filter_dist[k]->createCtx(vars_per_node);
    filter_ctx[k]->incref();

    // Create the design variable vector for this level
    x[k] = new TACSBVec(filter_maps[k], vars_per_node, filter_dist[k]);
    x[k]->incref();
  }

  // Now create the interpolation between filter levels
  for ( int k = 1; k < nlevels; k++ ){
    // Create the interpolation object
    filter_interp[k-1] = new TACSBVecInterp(filter_maps[k],
                                            filter_maps[k-1], vars_per_node);
    filter_interp[k-1]->incref();

    // Create the interpolation on the TMR side
    oct_filter[k-1]->createInterpolation(oct_filter[k], filter_interp[k-1]);
    filter_interp[k-1]->initialize();
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
  int gmres_iters = 100; 
  int nrestart = 2;
  int is_flexible = 1;
  ksm = new GMRES(mg->getMat(0), mg, 
                  gmres_iters, nrestart, is_flexible);
  ksm->incref();
  ksm->setMonitor(new KSMPrintStdout("GMRES", mpi_rank, 10));
  ksm->setTolerances(1e-10, 1e-30);

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
}

/*
  Create the topology optimization problem for 2D quad meshes
*/
TMRTopoProblem::TMRTopoProblem( int _nlevels, 
                                TACSAssembler *_tacs[],
                                TMRQuadForest *_filter[], 
                                TACSVarMap *_filter_maps[],
                                TACSBVecIndices *_filter_indices[],
                                TACSMg *_mg,
                                int _vars_per_node ):
 ParOptProblem(_tacs[0]->getMPIComm()){
  // Set the prefix to NULL
  prefix = NULL;
  
  // Get the processor rank
  int mpi_rank;
  MPI_Comm_rank(_tacs[0]->getMPIComm(), &mpi_rank);
  
  // Set the number of variables per node
  vars_per_node = _vars_per_node;

  // Set the maximum number of indices
  max_local_size = 0;

  // The initial design variable values (may not be set)
  xinit = NULL;

  // Set the number of levels
  nlevels = _nlevels;

  // Allocate the arrays
  tacs = new TACSAssembler*[ nlevels ];
  quad_filter = new TMRQuadForest*[ nlevels ];
  oct_filter = NULL;
  filter_maps = new TACSVarMap*[ nlevels ];
  filter_indices = new TACSBVecIndices*[ nlevels ];
  filter_dist = new TACSBVecDistribute*[ nlevels ];
  filter_ctx = new TACSBVecDistCtx*[ nlevels ];
  filter_interp = new TACSBVecInterp*[ nlevels-1 ];

  // The design variable vector for each level
  x = new TACSBVec*[ nlevels ];

  // Copy over the assembler objects and filters
  for ( int k = 0; k < nlevels; k++ ){
    // Set the TACSAssembler objects for each level
    tacs[k] = _tacs[k];
    tacs[k]->incref();

    // Set the filter object
    quad_filter[k] = _filter[k];
    quad_filter[k]->incref();

    // Copy over the filter information
    filter_maps[k] = _filter_maps[k];
    filter_maps[k]->incref();

    filter_indices[k] = _filter_indices[k];
    filter_indices[k]->incref();

    // Set the maximum local size
    const int *range;
    filter_maps[k]->getOwnerRange(&range);
    int size = vars_per_node*((range[mpi_rank+1] - range[mpi_rank]) +
                              filter_indices[k]->getNumIndices());

    // Update the maximum local size
    if (size > max_local_size){
      max_local_size = size;
    }

    // Create the distribution object for the design variables
    filter_dist[k] = new TACSBVecDistribute(filter_maps[k], 
                                            filter_indices[k]);
    filter_dist[k]->incref();

    // Create the transfer context
    filter_ctx[k] = filter_dist[k]->createCtx(vars_per_node);
    filter_ctx[k]->incref();

    // Create the design variable vector for this level
    x[k] = new TACSBVec(filter_maps[k], vars_per_node, filter_dist[k]);
    x[k]->incref();
  }

  // Now create the interpolation between filter levels
  for ( int k = 1; k < nlevels; k++ ){
    // Create the interpolation object
    filter_interp[k-1] = new TACSBVecInterp(filter_maps[k],
                                            filter_maps[k-1], vars_per_node);
    filter_interp[k-1]->incref();

    // Create the interpolation on the TMR side
    quad_filter[k-1]->createInterpolation(quad_filter[k], filter_interp[k-1]);
    filter_interp[k-1]->initialize();
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
  int gmres_iters = 100; 
  int nrestart = 2;
  int is_flexible = 1;
  ksm = new GMRES(mg->getMat(0), mg, 
                  gmres_iters, nrestart, is_flexible);
  ksm->incref();
  ksm->setMonitor(new KSMPrintStdout("GMRES", mpi_rank, 10));
  ksm->setTolerances(1e-10, 1e-30);

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
    filter_maps[k]->decref();
    filter_indices[k]->decref();
    filter_dist[k]->decref();
    filter_ctx[k]->decref();
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
  delete [] filter_indices;
  delete [] filter_dist;
  delete [] filter_ctx;
  delete [] x;

  for ( int k = 0; k < nlevels-1; k++ ){
    filter_interp[k]->decref();
  }
  delete [] filter_interp;

  // Free the initial design variable values (if allocated)
  if (xinit){ xinit->decref(); }

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
    load_case_info[load_case].offset = NULL;
    load_case_info[load_case].scale = NULL;
  }
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
                                             int max_lanczos,
                                             double eigtol ){
  if (!freq){
    // Create a mass matrix for the frequency constraint
    TACSMat *mmat = tacs[0]->createMat();

    // Create the frequency analysis object
    freq = new TACSFrequencyAnalysis(tacs[0], sigma, mmat, 
                                     mg->getMat(0), ksm,
                                     max_lanczos, num_eigvals, eigtol);
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
  for ( int i = 0; i < num_load_cases; i++ ){
    if (_obj_funcs[i]){ obj_funcs[i]->incref(); }
  }
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
  }

  if (freq){
    num_constraints++;
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
void TMRTopoProblem::setInitDesignVars( ParOptVec *xvars ){
  if (xvars){
    ParOptBVecWrap *wrap = dynamic_cast<ParOptBVecWrap*>(xvars);
    if (wrap){
      if (xinit){ xinit->decref(); }
      xinit = new TACSBVec(filter_maps[0], 
                           vars_per_node, filter_dist[0]);
      xinit->incref();
      xinit->copyValues(wrap->vec);
    }
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
  return new ParOptBVecWrap(new TACSBVec(filter_maps[0],
                                         vars_per_node,
                                         filter_dist[0]));
}

/*
  Compute the volume corresponding to each node within the filter
*/
TACSBVec* TMRTopoProblem::createVolumeVec(){
  // Get the dependent nodes and weight values
  const int *dep_ptr, *dep_conn;
  const double *dep_weights;
  int ndep = oct_filter[0]->getDepNodeConn(&dep_ptr, &dep_conn, &dep_weights);

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
  TMROctantArray *nodes;
  oct_filter[0]->getNodes(&nodes);
  int node_size;
  TMROctant *node_array;
  nodes->getArray(&node_array, &node_size);  

  // Get the node locations from the filter
  TMRPoint *X;
  oct_filter[0]->getPoints(&X);
  
  // Loop over the elements 
  for ( int i = 0; i < size; i++ ){
    int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level);

    // Retrieve the node numbers and x/y/z locations for element i
    // within the mesh
    int node_nums[8];
    TMRPoint Xpts[8];
    for ( int kk = 0; kk < 2; kk++ ){
      for ( int jj = 0; jj < 2; jj++ ){
        for ( int ii = 0; ii < 2; ii++ ){
          // Find the octant node (without level/tag)
          TMROctant node;
          node.block = array[i].block;
          node.x = array[i].x + ii*h;
          node.y = array[i].y + jj*h;
          node.z = array[i].z + kk*h;
          oct_filter[0]->transformNode(&node);
          
          // Find the corresponding node
          const int use_nodes = 1;
          TMROctant *t = nodes->contains(&node, use_nodes);
          
          if (t){
            // Copy the node location
            Xpts[ii + 2*jj + 4*kk] = X[t - node_array];
            
            // Copy the node number
            node_nums[ii + 2*jj + 4*kk] = t->tag;
          }
        }
      }
    }
    
    // Compute the element area
    TacsScalar Area = 0.0;
    
    // The quadrature points
    const double gaussPts[] = {-0.577350269189626, 0.577350269189626};

    for ( int kk = 0; kk < 2; kk++ ){
      for ( int jj = 0; jj < 2; jj++ ){
        for ( int ii = 0; ii < 2; ii++ ){
          double pt[3];
          pt[0] = gaussPts[ii];
          pt[1] = gaussPts[jj];
          pt[2] = gaussPts[kk];

          // Evaluate the derivative of the shape functions
          double N[8], Na[8], Nb[8], Nc[8];
          FElibrary::triLagrangeSF(N, Na, Nb, Nc, pt, 2);

          // Compute the Jacobian transformation
          TacsScalar J[9];
          memset(J, 0, 9*sizeof(TacsScalar));
          for ( int j = 0; j < 8; j++ ){
            J[0] += Na[j]*Xpts[j].x;
            J[1] += Nb[j]*Xpts[j].x;
            J[2] += Nc[j]*Xpts[j].x;
            
            J[3] += Na[j]*Xpts[j].y;
            J[4] += Nb[j]*Xpts[j].y;
            J[5] += Nc[j]*Xpts[j].y;

            J[6] += Na[j]*Xpts[j].z;
            J[7] += Nb[j]*Xpts[j].z;
            J[8] += Nc[j]*Xpts[j].z;
          }

          // Add the determinant to the area - the weights in this
          // case are just = 1.0
          TacsScalar det = FElibrary::jacobian3d(J);
          Area += det;
        }
      }
    }
  
    // Distribute the element area equally to all nodes in the
    // element
    TacsScalar a[8];
    for ( int k = 0; k < 8; k++ ){
      a[k] = 0.125*Area;
    }
    
    // Add the values to the vector
    vec->setValues(8, node_nums, a, TACS_ADD_VALUES);
  }

  vec->beginSetValues(TACS_ADD_VALUES);
  vec->endSetValues(TACS_ADD_VALUES);

  return vec;
}

/*
  Compute the volume corresponding to each node within the filter
*/
TACSBVec* TMRTopoProblem::createAreaVec(){
  // Get the dependent nodes and weight values
  const int *dep_ptr, *dep_conn;
  const double *dep_weights;
  int ndep = quad_filter[0]->getDepNodeConn(&dep_ptr, &dep_conn, &dep_weights);

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
  TMRQuadrantArray *nodes;
  quad_filter[0]->getNodes(&nodes);
  int node_size;
  TMRQuadrant *node_array;
  nodes->getArray(&node_array, &node_size);  

  // Get the node locations from the filter
  TMRPoint *X;
  quad_filter[0]->getPoints(&X);
  
  // Loop over the elements 
  for ( int i = 0; i < size; i++ ){
    int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level);

    // Retrieve the node numbers and x/y/z locations for element i
    // within the mesh
    int node_nums[4];
    TMRPoint Xpts[4];
    
    for ( int jj = 0; jj < 2; jj++ ){
      for ( int ii = 0; ii < 2; ii++ ){
        // Find the quadrant node (without level/tag)
        TMRQuadrant node;
        node.face = array[i].face;
        node.x = array[i].x + ii*h;
        node.y = array[i].y + jj*h;
        quad_filter[0]->transformNode(&node);
          
        // Find the corresponding node
        const int use_nodes = 1;
        TMRQuadrant *t = nodes->contains(&node, use_nodes);
          
        if (t){
          // Copy the node location
          Xpts[ii + 2*jj] = X[t - node_array];
            
          // Copy the node number
          node_nums[ii + 2*jj] = t->tag;
        }
      }
    }

    // Compute the element area
    TacsScalar Area = 0.0;
    
    // The quadrature points
    const double gaussPts[] = {-0.577350269189626, 0.577350269189626};
    for ( int jj = 0; jj < 2; jj++ ){
      for ( int ii = 0; ii < 2; ii++ ){
        double pt[2];
        pt[0] = gaussPts[ii];
        pt[1] = gaussPts[jj];
        
        // Evaluate the derivative of the shape functions
        double N[4], Na[4], Nb[4];
        FElibrary::biLagrangeSF(N, Na, Nb, pt, 2);

        // Compute the Jacobian transformation
        TacsScalar J[4];
        memset(J, 0, 4*sizeof(TacsScalar));
        for ( int j = 0; j < 4; j++ ){
          J[0] += Na[j]*Xpts[j].x;
          J[1] += Nb[j]*Xpts[j].x;
          
          J[2] += Na[j]*Xpts[j].y;
          J[3] += Nb[j]*Xpts[j].y;
        }
        // Add the determinant to the area - the weights in this
        // case are just = 1.0
        TacsScalar det = FElibrary::jacobian2d(J);
        Area += det;
      }
    }

    // Distribute the element area equally to all nodes in the
    // element
    TacsScalar a[4];
    for ( int k = 0; k < 4; k++ ){
      a[k] = 0.25*Area;
    }
    
    // Add the values to the vector
    vec->setValues(4, node_nums, a, TACS_ADD_VALUES);
  }

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
  // Get the values of the design variables from the inner-most
  // version of TACS
  if (xvec){ 
    ParOptBVecWrap *wrap = dynamic_cast<ParOptBVecWrap*>(xvec);
    if (wrap){
      if (xinit){
        wrap->vec->copyValues(xinit);
      }
      else {
        memset(xlocal, 0, max_local_size*sizeof(TacsScalar));
        tacs[0]->getDesignVars(xlocal, max_local_size);

        // Set the local values into the vector
        setBVecFromLocalValues(xlocal, wrap->vec);
        wrap->vec->beginSetValues(TACS_INSERT_NONZERO_VALUES);
        wrap->vec->endSetValues(TACS_INSERT_NONZERO_VALUES);
      }
    }
  }

  if (lbvec || ubvec){
    TacsScalar *upper = new TacsScalar[ max_local_size ];
    memset(xlocal, 0, max_local_size*sizeof(TacsScalar));
    memset(upper, 0, max_local_size*sizeof(TacsScalar));
    tacs[0]->getDesignVarRange(xlocal, upper, max_local_size);

    if (lbvec){
      ParOptBVecWrap *lbwrap = dynamic_cast<ParOptBVecWrap*>(lbvec);
      if (lbwrap){
        setBVecFromLocalValues(xlocal, lbwrap->vec);
        lbwrap->vec->beginSetValues(TACS_INSERT_NONZERO_VALUES);
        lbwrap->vec->endSetValues(TACS_INSERT_NONZERO_VALUES);
      }
    }
    if (ubvec){
      ParOptBVecWrap *ubwrap = dynamic_cast<ParOptBVecWrap*>(ubvec);
      if (ubwrap){
        setBVecFromLocalValues(upper, ubwrap->vec);
        ubwrap->vec->beginSetValues(TACS_INSERT_NONZERO_VALUES);
        ubwrap->vec->endSetValues(TACS_INSERT_NONZERO_VALUES);
      }
    }
    delete [] upper;
  }
}

/*
  Set the linearization point when the analysis/design problem uses
  the TMRLinearOctStiffness class
*/
void TMRTopoProblem::setLinearization( double q, ParOptVec *xvec ){
  ParOptBVecWrap *wrap = dynamic_cast<ParOptBVecWrap*>(xvec);

  // Set the design variable values
  if (wrap){
    // Copy the values to the local design variable vector
    x[0]->copyValues(wrap->vec);
  
    // Distribute the design variable values
    x[0]->beginDistributeValues();
    x[0]->endDistributeValues();

    // Copy the values to the local array
    int size = getLocalValuesFromBVec(x[0], xlocal);

    // Get the res vector
    int num_elements = tacs[0]->getNumElements();
    for ( int k = 0; k < num_elements; k++ ){
      TACSElement *element = tacs[0]->getElement(k, NULL, NULL, 
                                                 NULL, NULL);
      
      // Get the constitutive object
      TACSConstitutive *constitutive = element->getConstitutive();
      
      if (constitutive){
        TMRLinearOctStiffness *con =
          dynamic_cast<TMRLinearOctStiffness*>(constitutive);

        if (con){
          con->setLinearization(q, xlocal, size);
        }
      }
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
int TMRTopoProblem::getLocalValuesFromBVec( TACSBVec *vec,
                                            TacsScalar *xloc ){
  TacsScalar *x_vals, *x_ext_vals;
  int size = vec->getArray(&x_vals);
  int ext_size = vec->getExtArray(&x_ext_vals);
  memcpy(xloc, x_vals, size*sizeof(TacsScalar));
  if (x_ext_vals){
    memcpy(&xloc[size], x_ext_vals, ext_size*sizeof(TacsScalar));
  }
  return size + ext_size;
}

/*
  Set the values into the TACSBVec object using the local
  size
*/
void TMRTopoProblem::setBVecFromLocalValues( const TacsScalar *xloc,
                                             TACSBVec *vec ){
  TacsScalar *x_vals, *x_ext_vals;
  int size = vec->getArray(&x_vals);
  int ext_size = vec->getExtArray(&x_ext_vals);
  memcpy(x_vals, xloc, size*sizeof(TacsScalar));
  if (x_ext_vals){
    memcpy(x_ext_vals, &xloc[size], ext_size*sizeof(TacsScalar));
  }
}

/*
  Evaluate the objective and constraints
*/
int TMRTopoProblem::evalObjCon( ParOptVec *pxvec,
                                ParOptScalar *fobj, 
                                ParOptScalar *cons ){
  ParOptBVecWrap *wrap = dynamic_cast<ParOptBVecWrap*>(pxvec);

  if (wrap){
    TACSBVec *xvec = wrap->vec;

    // Copy the values to the local design variable vector
    x[0]->copyValues(xvec);

    // Distribute the design variable values
    x[0]->beginDistributeValues();
    x[0]->endDistributeValues();

    // Copy the values to the local array
    int size = getLocalValuesFromBVec(x[0], xlocal);
    tacs[0]->setDesignVars(xlocal, size);

    // Set the design variable values on all processors
    for ( int k = 0; k < nlevels-1; k++ ){
      filter_interp[k]->multWeightTranspose(x[k], x[k+1]);

      // Distribute the design variable values
      x[k+1]->beginDistributeValues();
      x[k+1]->endDistributeValues();

      // Set the design variable values
      size = getLocalValuesFromBVec(x[k+1], xlocal);
      tacs[k+1]->setDesignVars(xlocal, size);
    }

    // Zero the variables
    tacs[0]->zeroVariables();

    // Assemble the Jacobian on each level
    double alpha = 1.0, beta = 0.0, gamma = 0.0;
    mg->assembleJacobian(alpha, beta, gamma, NULL);
    mg->factor();
    
    // Set the objective value
    *fobj = 0.0;

    // Get the rank of comm
    int mpi_rank;
    MPI_Comm_rank(tacs[0]->getMPIComm(), &mpi_rank);

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
        tacs[0]->applyBCs(vars[i]);
      
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
        tacs[0]->evalFunctions(load_case_info[i].funcs, 
                               num_funcs, &cons[count]);

        // Scale and offset the constraints that we just evaluated
        for ( int j = 0; j < num_funcs; j++ ){
          TacsScalar offset = load_case_info[i].offset[j];
          TacsScalar scale = load_case_info[i].scale[j];
          cons[count + j] = scale*(cons[count+j] + offset);
        }
        count += num_funcs;
      }
    }
  
    // Compute the natural frequency constraint, if any
    if (freq){
      // Keep track of the number of eigenvalues with unacceptable
      // error. If more than one exists, re-solve the eigenvalue
      // problem again.
      int err_count = 1;

      // Keep track of the smallest eigenvalue
      TacsScalar smallest_eigval = 0.0;
      while (err_count > 0){
        // Set the error counter to zero
        err_count = 0;
          
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
          if (error > freq_eig_tol){
            err_count++;
          }
        }

        // If there is significant error in computing the eigenvalues,
        // reset the buckling computation
        if (err_count > 0){
          double sigma = 0.95*smallest_eigval;
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
  }
  else {
    return 1;
  }

  return 0;
}

/*
  Evaluate the objective and constraint gradients
*/
int TMRTopoProblem::evalObjConGradient( ParOptVec *xvec, 
                                        ParOptVec *gvec, 
                                        ParOptVec **Acvec ){
  // Evaluate the derivative of the weighted compliance with
  // respect to the design variables
  gvec->zeroEntries();
  ParOptBVecWrap *wrap = dynamic_cast<ParOptBVecWrap*>(gvec);
  if (wrap){
    TACSBVec *g = wrap->vec;

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

        setBVecFromLocalValues(xlocal, g);
        g->beginSetValues(TACS_ADD_VALUES);
        g->endSetValues(TACS_ADD_VALUES);
      }      
    }  
    else { // For compliance objective
      for ( int i = 0; i < num_load_cases; i++ ){
        tacs[0]->setVariables(vars[i]);
        tacs[0]->addAdjointResProducts(-obj_weights[i], &vars[i], 
                                       1, xlocal, max_local_size);
      }
      setBVecFromLocalValues(xlocal, g);
      g->beginSetValues(TACS_ADD_VALUES);
      g->endSetValues(TACS_ADD_VALUES);
    }   
  } // end if wrap
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
        // Get the vector
        TACSBVec *A = wrap->vec;

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

        // Wrap the vector class
        setBVecFromLocalValues(xlocal, A);
        A->beginSetValues(TACS_ADD_VALUES);
        A->endSetValues(TACS_ADD_VALUES);
      }
    } //  num_funcs
    count += num_funcs;

    if (freq){
      // Try to unwrap the vector
      wrap = dynamic_cast<ParOptBVecWrap*>(Acvec[count]);
      if (wrap){
        // Get the vector
        TACSBVec *A = wrap->vec;
        memset(xlocal, 0, max_local_size*sizeof(TacsScalar));

        TacsScalar *temp = new TacsScalar[max_local_size];
        memset(temp, 0, max_local_size*sizeof(TacsScalar));

        // Add the contribution from each eigenvalue derivative
        TacsScalar smallest_eigval = 0.0;
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
          if (k == 0){
            smallest_eigval = eigval;
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
        setBVecFromLocalValues(xlocal, A);
        A->beginSetValues(TACS_ADD_VALUES);
        A->endSetValues(TACS_ADD_VALUES);
      }
      
      count++;
    }
  }

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
  /*
  ParOptBVecWrap *wrap = dynamic_cast<ParOptBVecWrap*>(pxvec);
  ParOptBVecWrap *hwrap = dynamic_cast<ParOptBVecWrap*>(hvec);

  if (wrap && hwrap){
    TACSBVec *px = wrap->vec;

    // Distribute the design variable values
    px->beginDistributeValues();
    px->endDistributeValues();

    // Copy the values to the local array
    int size = getLocalValuesFromBVec(px, xlocal);

    // Zero the entries in the RHS
    res->zeroEntries();

    // Set the maximum number of nodes/variables/stresses
    static const int NUM_NODES = 8;
    static const int NUM_STRESSES = 6;
    static const int NUM_VARIABLES = 3*NUM_NODES;
  
    // Get the res vector
    int num_elements = tacs[0]->getNumElements();
    for ( int k = 0; k < num_elements; k++ ){
      TACSElement *element = tacs[0]->getElement(k, NULL, NULL, 
                                                 NULL, NULL);
      
      // Dynamically cast the element to the 2D element type
      TACS3DElement<NUM_NODES> *elem = 
        dynamic_cast<TACS3DElement<NUM_NODES>*>(element);

      if (elem){
        // Get the element data
        TacsScalar Xpts[3*NUM_NODES];
        TacsScalar vars[NUM_VARIABLES], dvars[NUM_VARIABLES];
        TacsScalar ddvars[NUM_VARIABLES];
        tacs[0]->getElement(k, Xpts, vars, dvars, ddvars);
        
        // Get the constitutive object
        TACSConstitutive *constitutive = elem->getConstitutive();
      
        TMRLinearOctStiffness *con =
          dynamic_cast<TMRLinearOctStiffness*>(constitutive);

        // Get the number of variables
        int nvars = elem->numVariables();

        if (con){
          TacsScalar elemRes[NUM_VARIABLES];
          memset(elemRes, 0, nvars*sizeof(TacsScalar));

          // The shape functions associated with the element
          double N[NUM_NODES];
          double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
          
          // The derivative of the stress with respect to the strain
          TacsScalar B[NUM_STRESSES*NUM_VARIABLES];
        
          // Get the number of quadrature points
          int numGauss = elem->getNumGaussPts();

          for ( int n = 0; n < numGauss; n++ ){
            // Retrieve the quadrature points and weight
            double pt[3];
            double weight = elem->getGaussWtsPts(n, pt);

            // Compute the element shape functions
            elem->getShapeFunctions(pt, N, Na, Nb, Nc);

            // Compute the derivative of X with respect to the
            // coordinate directions
            TacsScalar X[3], Xa[9];
            elem->solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);

            // Compute the determinant of Xa and the transformation
            TacsScalar J[9];
            TacsScalar h = FElibrary::jacobian3d(Xa, J);
            h = h*weight;

            // Compute the strain
            TacsScalar strain[NUM_STRESSES];
            elem->evalStrain(strain, J, Na, Nb, Nc, vars);
 
            // Compute the corresponding stress
            TacsScalar stress[NUM_STRESSES];
            con->calcStressDVProject(pt, strain, xlocal, size, stress);
       
            // Get the derivative of the strain with respect to the
            // nodal displacements
            elem->getBmat(B, J, Na, Nb, Nc, vars);

            TacsScalar *b = B;
            for ( int i = 0; i < nvars; i++ ){
              elemRes[i] += h*(b[0]*stress[0] + b[1]*stress[1] + 
                               b[2]*stress[2] + b[3]*stress[3] + 
                               b[4]*stress[4] + b[5]*stress[5]);
              b += NUM_STRESSES;
            }
          }

          // Get the local element ordering
          int len;
          const int *nodes;
          tacs[0]->getElement(k, &nodes, &len);

          // Add the residual values
          res->setValues(len, nodes, elemRes, TACS_ADD_VALUES);
        }
      }
    }

    // Add the residual values
    res->beginSetValues(TACS_ADD_VALUES);
    res->endSetValues(TACS_ADD_VALUES);

    // Set the boundary conditions
    tacs[0]->applyBCs(res);

    // Solve the system: K(x)*u = res
    ksm->solve(res, vars);

    // Set the design variable values
    memset(xlocal, 0, size*sizeof(TacsScalar));
    tacs[0]->addAdjointResProducts(2.0*obj_scale, &vars, 1, xlocal, size);

    // Set the variables from the product
    hvec->zeroEntries();
    setBVecFromLocalValues(xlocal, hwrap->vec);

    // Begin setting the values
    hwrap->vec->beginSetValues(TACS_ADD_VALUES);
    hwrap->vec->endSetValues(TACS_ADD_VALUES);
  }
  */
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

// Write the output file
void TMRTopoProblem::writeOutput( int iter, ParOptVec *xvec ){
  // Print out the binary STL file for later visualization
  if (prefix && oct_filter){
    // Write out the file at a cut off of 0.25
    char *filename = new char[ strlen(prefix) + 100 ];

    for ( int k = 0; k < vars_per_node; k++ ){
      double cutoff = 0.5;
      sprintf(filename, "%s/levelset05_var%d_binary%04d.bstl", prefix, k, iter_count);
      TMR_GenerateBinFile(filename, oct_filter[0], x[0], k, cutoff);
    }
    
    delete [] filename;
  }

  // Update the iteration count
  iter_count++;
}

