#include "TMRTopoProblem.h"
#include "TACSFunction.h"
#include "TMR_STLTools.h"

ParOptBVecWrap::ParOptBVecWrap( TACSBVec *_vec ){
  vec = _vec;
  vec->incref();
}

ParOptBVecWrap::~ParOptBVecWrap(){
  vec->decref();
}

// Perform standard operations required for linear algebra
// -------------------------------------------------------
void ParOptBVecWrap::set( ParOptScalar alpha ){ 
  vec->set(alpha); 
}

void ParOptBVecWrap::zeroEntries(){ 
  vec->zeroEntries(); 
}

void ParOptBVecWrap::copyValues( ParOptVec *pvec ){
  ParOptBVecWrap *avec = dynamic_cast<ParOptBVecWrap*>(pvec);
  if (avec){
    vec->copyValues(avec->vec);
  }
}
 
double ParOptBVecWrap::norm(){ 
  return vec->norm(); 
}

double ParOptBVecWrap::maxabs(){
  TacsScalar *x = NULL;
  int size = vec->getArray(&x);
  
  double res = 0.0;
  for ( int i = 0; i < size; i++ ){
    if (fabs(RealPart(x[i])) > res){
      res = fabs(RealPart(x[i]));
    }
  }
  
  double infty_norm = 0.0;
  MPI_Allreduce(&res, &infty_norm, 1, MPI_DOUBLE, MPI_MAX, vec->getMPIComm());
  
  return infty_norm;
}

ParOptScalar ParOptBVecWrap::dot( ParOptVec *pvec ){
  ParOptBVecWrap *avec = dynamic_cast<ParOptBVecWrap*>(pvec);
  if (avec){
    return vec->dot(avec->vec);
  }
  return 0.0;
}

void ParOptBVecWrap::mdot( ParOptVec **vecs, int nvecs, ParOptScalar *output ){
  TACSVec **tvecs = new TACSVec*[ nvecs ];
  for ( int k = 0; k < nvecs; k++ ){
    ParOptBVecWrap *wrap = dynamic_cast<ParOptBVecWrap*>(vecs[k]);
    tvecs[k] = NULL;
    if (wrap){
      tvecs[k] = wrap->vec;
    }
  }
  
  vec->mdot(tvecs, output, nvecs);
  delete tvecs;
}

void ParOptBVecWrap::scale( ParOptScalar alpha ){ 
  vec->scale(alpha); 
}

void ParOptBVecWrap::axpy( ParOptScalar alpha, ParOptVec *pvec ){
  ParOptBVecWrap *avec = dynamic_cast<ParOptBVecWrap*>(pvec);
  if (avec){
    return vec->axpy(alpha, avec->vec);
  }
}

int ParOptBVecWrap::getArray( ParOptScalar **array ){
  TacsScalar *_array;
  vec->getArray(&_array);
  *array = _array;    
}



TMRTopoProblem::TMRTopoProblem( int _nlevels, 
                                TACSAssembler *_tacs[],
                                TACSBVec *_force,
                                TMROctForest *_filter[], 
                                TACSVarMap *_filter_maps[],
                                TACSBVecIndices *_filter_indices[],
                                TACSMg *_mg,
                                double _target_mass,
                                const char *_prefix ):
ParOptProblem(_tacs[0]->getMPIComm()){
  // Copy over the prefix - max of 128 characters
  strncpy(prefix, _prefix, sizeof(prefix));
  
  // Set the number of levels
  nlevels = _nlevels;
  if (nlevels > MAX_NUM_LEVELS){
    nlevels = MAX_NUM_LEVELS;
  }
  
  // Get the processor rank
  int mpi_rank;
  MPI_Comm_rank(_tacs[0]->getMPIComm(), &mpi_rank);
  
  // Set the maximum number of indices
  int max_local_size = 0;
  
  // Copy over the assembler objects and filters
  for ( int k = 0; k < nlevels; k++ ){
    // Set the TACSAssembler objects for each level
    tacs[k] = _tacs[k];
    tacs[k]->incref();

    // Set the filter object
    filter[k] = _filter[k];

    // Copy over the filter information
    filter_maps[k] = _filter_maps[k];
    filter_maps[k]->incref();

    filter_indices[k] = _filter_indices[k];
    filter_indices[k]->incref();

    // Set the maximum local size
    const int *range;
    filter_maps[k]->getOwnerRange(&range);
    int size = (range[mpi_rank+1] - range[mpi_rank]) +
      filter_indices[k]->getNumIndices();

    // Update the maximum local size
    if (size > max_local_size){
      max_local_size = size;
    }

    // Create the distribution object for the design variables
    filter_dist[k] = new TACSBVecDistribute(filter_maps[k], filter_indices[k]);
    filter_dist[k]->incref();

    // Create the transfer context
    filter_ctx[k] = filter_dist[k]->createCtx(1);
    filter_ctx[k]->incref();

    // Create the design variable vector for this level
    x[k] = new TACSBVec(filter_maps[k], 1, NULL, filter_dist[k]);
  }

  // Now create the interpolation between filter levels
  for ( int k = 1; k < nlevels; k++ ){
    // Create the interpolation on the TMR side
    int *ptr, *conn;
    double *weights;
    filter[k-1]->createInterpolation(filter[k],
                                     &ptr, &conn, &weights);

    // Create the interpolation object
    filter_interp[k-1] = new TACSBVecInterp(filter_maps[k],
                                            filter_maps[k-1], 1);
    filter_interp[k-1]->incref();

    // Get the range of nodes
    const int *node_range;
    filter[k-1]->getOwnedNodeRange(&node_range);

    // Add all the values in the interpolation
    for ( int node = node_range[mpi_rank], i = 0;
          node < node_range[mpi_rank+1]; node++, i++ ){
      int len = ptr[i+1] - ptr[i];
      filter_interp[k-1]->addInterp(node, &weights[ptr[i]], 
                                    &conn[ptr[i]], len);
    }

    filter_interp[k-1]->initialize();
    
    delete [] ptr;
    delete [] conn;
    delete [] weights;
  }
  
  // Set the maximum local size
  xlocal = new TacsScalar[ max_local_size ];

  // The multigrid object
  mg = _mg;
  mg->incref();

  // Set up the solver
  int gmres_iters = 100; 
  int nrestart = 2;
  int is_flexible = 1;
  ksm = new GMRES(mg->getMat(0), mg, 
                  gmres_iters, nrestart, is_flexible);
  ksm->incref();
  ksm->setMonitor(new KSMPrintStdout("GMRES", mpi_rank, 10));
  ksm->setTolerances(1e-10, 1e-30);

  // Allocate variables/residual
  force = _force;
  force->incref();

  vars = tacs[0]->createVec();
  vars->incref();

  // Create the compliance and mass functions
  compliance = new TACSCompliance(tacs[0]);
  mass = new TACSStructuralMass(tacs[0]);
  compliance->incref();
  mass->incref();

  // Set the target mass
  target_mass = _target_mass;
  obj_scale = -1.0;

  // Set the problem sizes
  int nvars = x[0]->getArray(NULL);
  setProblemSizes(nvars, 1, 0, 0);
}

/*
  Free the data stored in the object
*/
TMRTopoProblem::~TMRTopoProblem(){
  for ( int k = 0; k < nlevels; k++ ){
    tacs[k]->decref();
   
    // Decrease the reference counts
    filter_maps[k]->decref();
    filter_indices[k]->decref();
    filter_dist[k]->decref();
    filter_ctx[k]->decref();
  }

  for ( int k = 0; k < nlevels-1; k++ ){
    filter_interp[k]->decref();
  }

  delete [] xlocal;
  mg->decref();
  ksm->decref();
  vars->decref();
  force->decref();

  compliance->decref();
  mass->decref();
}

/*
  Get the objective scaling factor
*/
ParOptScalar TMRTopoProblem::getObjectiveScaling(){
  return obj_scale;
}

/*
  Set the objective scaling factor
*/
void TMRTopoProblem::setObjectiveScaling( ParOptScalar scale ){
  obj_scale = scale;
}

/*
  Create a design variable vector
*/
ParOptVec *TMRTopoProblem::createDesignVec(){
  return new ParOptBVecWrap(new TACSBVec(filter_maps[0], 1, NULL, filter_dist[0]));
}

/*
  Set whether these should be considered sparse inequalities
*/
int TMRTopoProblem::isSparseInequality(){ return 0; }
int TMRTopoProblem::isDenseInequality(){ return 0; }
int TMRTopoProblem::useLowerBounds(){ return 1; }
int TMRTopoProblem::useUpperBounds(){ return 1; }

/*
  Set the initial design variables
*/ 
void TMRTopoProblem::getVarsAndBounds( ParOptVec *x, 
                                       ParOptVec *lb, 
                                       ParOptVec *ub ){
  // Get the values of the design variables from the inner-most
  // version of TACS
  if (x){ x->set(0.95); }

  // Set the lower and upper bounds on the design variables
  if (lb){ lb->set(0.001); }
  if (ub){ ub->set(1.0); }
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
  memcpy(&xloc[size], x_ext_vals, ext_size*sizeof(TacsScalar));
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
  memcpy(x_ext_vals, &xloc[size], ext_size*sizeof(TacsScalar));
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

    // Assemble the Jacobian on each level
    double alpha = 1.0, beta = 0.0, gamma = 0.0;
    mg->assembleJacobian(alpha, beta, gamma, NULL);
    mg->factor();
    
    // Solve the system: K(x)*u = force
    ksm->solve(force, vars);
    tacs[0]->setVariables(vars);
  
    // Evaluate the functions in parallel
    TacsScalar fvals[2];
    TACSFunction *funcs[2];
    funcs[0] = compliance;
    funcs[1] = mass;
    tacs[0]->evalFunctions(funcs, 2, fvals);
  
    // Set the scaling for the objective function
    if (obj_scale < 0.0){
      obj_scale = 1.0/fvals[0];
    }

    // Set the compliance objective and the mass constraint
    *fobj = obj_scale*fvals[0];
    cons[0] = 1.0 - fvals[1]/target_mass;
  }
  else {
    return 1;
  }

  return 0;
}

/*
  Evaluate the objective and constraint gradients
*/
int TMRTopoProblem::evalObjConGradient( ParOptVec *x, 
                                        ParOptVec *g, 
                                        ParOptVec **Ac ){
  // Evaluate the derivative of the mass and compliance with
  // respect to the design variables
  g->zeroEntries();
  Ac[0]->zeroEntries();

  // Compute the derivative of the mass and compliance w.r.t the design variables
  ParOptBVecWrap *wrap1 = NULL, *wrap2 = NULL;
  wrap1 = dynamic_cast<ParOptBVecWrap*>(g);
  wrap2 = dynamic_cast<ParOptBVecWrap*>(Ac[0]);
  if (wrap1 && wrap2){
    TACSBVec *g_vec = wrap1->vec;
    TACSBVec *m_vec = wrap2->vec;
    int size = g_vec->getArray(NULL) + g_vec->getExtArray(NULL);

    memset(xlocal, 0, size*sizeof(TacsScalar));
    tacs[0]->addAdjointResProducts(-obj_scale, &vars, 1, xlocal, size);
    setBVecFromLocalValues(xlocal, g_vec);
    g_vec->beginSetValues(ADD_VALUES);

    memset(xlocal, 0, size*sizeof(TacsScalar));
    tacs[0]->addDVSens(-1.0/target_mass, &mass, 1, xlocal, size);
    setBVecFromLocalValues(xlocal, m_vec);
    m_vec->beginSetValues(ADD_VALUES);

    // Finsh adding the values
    g_vec->endSetValues(ADD_VALUES);
    m_vec->endSetValues(ADD_VALUES);
  }
  else {
    return 1;
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
                                     ParOptVec *hvec ){}

// Evaluate the sparse constraints
// ------------------------
void TMRTopoProblem::evalSparseCon( ParOptVec *x, 
                                    ParOptVec *out ){}

// Compute the Jacobian-vector product out = J(x)*px
// --------------------------------------------------
void TMRTopoProblem::addSparseJacobian( double alpha, 
                                        ParOptVec *x,
                                        ParOptVec *px, 
                                        ParOptVec *out ){}

// Compute the transpose Jacobian-vector product out = J(x)^{T}*pzw
// -----------------------------------------------------------------
void TMRTopoProblem::addSparseJacobianTranspose( double alpha, 
                                                 ParOptVec *x,
                                                 ParOptVec *pzw, 
                                                 ParOptVec *out ){}

// Add the inner product of the constraints to the matrix such 
// that A += J(x)*cvec*J(x)^{T} where cvec is a diagonal matrix
// ------------------------------------------------------------
void TMRTopoProblem::addSparseInnerProduct( double alpha, 
                                            ParOptVec *x,
                                            ParOptVec *cvec, 
                                            double *A ){}

// Write the output file
void TMRTopoProblem::writeOutput( int iter, ParOptVec *x ){
  ParOptBVecWrap *wrap = dynamic_cast<ParOptBVecWrap*>(x);

  if (wrap){
    // Print out the binary STL file for later visualization
    int var_offset = 0;
    double cutoff = 0.25;
    char filename[256];
    sprintf(filename, "%s/levelset%.2f_binary%04d.bstl", prefix, cutoff, iter);
    int fail = TMR_GenerateBinFile(filename, filter[0], wrap->vec, var_offset, cutoff);
  }
}
