#ifdef TMR_HAS_PAROPT

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

// Perform standard operations required for linear algebra
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
    if (fabs(TacsRealPart(x[i])) > res){
      res = fabs(TacsRealPart(x[i]));
    }
  }
  
  double infty_norm = 0.0;
  MPI_Allreduce(&res, &infty_norm, 1, MPI_DOUBLE, MPI_MAX, 
                vec->getMPIComm());
  
  return infty_norm;
}

ParOptScalar ParOptBVecWrap::dot( ParOptVec *pvec ){
  ParOptBVecWrap *avec = dynamic_cast<ParOptBVecWrap*>(pvec);
  if (avec){
    return vec->dot(avec->vec);
  }
  return 0.0;
}

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
  int size = 0;
  if (array){
    size = vec->getArray(&_array);
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
  max_local_size = 0;

  // The initial design variable values (may not be set)
  xinit = NULL;
  
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
    filter_dist[k] = new TACSBVecDistribute(filter_maps[k], 
                                            filter_indices[k]);
    filter_dist[k]->incref();

    // Create the transfer context
    filter_ctx[k] = filter_dist[k]->createCtx(1);
    filter_ctx[k]->incref();

    // Create the design variable vector for this level
    x[k] = new TACSBVec(filter_maps[k], 1, filter_dist[k]);
    x[k]->incref();
  }

  // Now create the interpolation between filter levels
  for ( int k = 1; k < nlevels; k++ ){
    // Create the interpolation object
    filter_interp[k-1] = new TACSBVecInterp(filter_maps[k],
                                            filter_maps[k-1], 1);
    filter_interp[k-1]->incref();

    // Create the interpolation on the TMR side
    filter[k-1]->createInterpolation(filter[k], filter_interp[k-1]);
    filter_interp[k-1]->initialize();
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

  // Do not use inverse variables by default
  use_inverse_vars = 0;

  // Allocate the variables/residual
  vars = tacs[0]->createVec();
  res = tacs[0]->createVec();
  vars->incref();
  res->incref();

  // Create the compliance and mass functions
  compliance = new TACSCompliance(tacs[0]);
  mass = new TACSStructuralMass(tacs[0]);
  compliance->incref();
  mass->incref();

  // Set the target mass
  target_mass = _target_mass;
  obj_scale = -1.0;
  mass_scale = 1.0/target_mass;

  // Set the problem sizes
  int nvars = x[0]->getArray(NULL);
  setProblemSizes(nvars, 1, 0, 0);

  // Set the iteration count
  iter_count = 0;
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
    x[k]->decref();
  }

  for ( int k = 0; k < nlevels-1; k++ ){
    filter_interp[k]->decref();
  }

  // Free the initial design variable values (if allocated)
  if (xinit){ xinit->decref(); }

  // Free the local temp array
  delete [] xlocal;

  // Free the solver/multigrid information
  mg->decref();
  ksm->decref();
  vars->decref();
  res->decref();

  // Free the function information
  compliance->decref();
  mass->decref();
}

/*
  Set the initial design variable values
*/
void TMRTopoProblem::setInitDesignVars( ParOptVec *vars ){
  if (vars){
    ParOptBVecWrap *wrap = dynamic_cast<ParOptBVecWrap*>(vars);
    if (wrap){
      xinit = new TACSBVec(filter_maps[0], 1, filter_dist[0]);
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
  Get the mass scaling factor
*/
ParOptScalar TMRTopoProblem::getMassScaling(){
  return mass_scale;
}

/*
  Set the mass constraint scaling factor
*/
void TMRTopoProblem::setMassScaling( ParOptScalar scale ){
  mass_scale = scale;
}

/*
  Create a design variable vector
*/
ParOptVec *TMRTopoProblem::createDesignVec(){
  return new ParOptBVecWrap(new TACSBVec(filter_maps[0], 1, 
                                         filter_dist[0]));
}

/*
  Set whether these should be considered sparse inequalities
*/
int TMRTopoProblem::isSparseInequality(){ 
  return 0; 
}

/*
  Use the inequality constraint - this seems to work better
*/
int TMRTopoProblem::isDenseInequality(){ 
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
  return 1; 
}

/*
  Set the initial design variables
*/ 
void TMRTopoProblem::getVarsAndBounds( ParOptVec *x, 
                                       ParOptVec *lb, 
                                       ParOptVec *ub ){
  if (use_inverse_vars){
    if (x){ x->set(1.05); }
    if (lb){ lb->set(1.0); }
    if (ub){ ub->set(1000.0); }
  }
  else {
    // Get the values of the design variables from the inner-most
    // version of TACS
    if (x){ 
      ParOptBVecWrap *wrap = dynamic_cast<ParOptBVecWrap*>(x);
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

    if (lb || ub){
      TacsScalar *upper = new TacsScalar[ max_local_size ];
      memset(xlocal, 0, max_local_size*sizeof(TacsScalar));
      memset(upper, 0, max_local_size*sizeof(TacsScalar));
      tacs[0]->getDesignVarRange(xlocal, upper, max_local_size);

      if (lb){
        ParOptBVecWrap *lbwrap = dynamic_cast<ParOptBVecWrap*>(lb);
        if (lbwrap){
          setBVecFromLocalValues(xlocal, lbwrap->vec);
          lbwrap->vec->beginSetValues(TACS_INSERT_NONZERO_VALUES);
          lbwrap->vec->endSetValues(TACS_INSERT_NONZERO_VALUES);
        }
      }
      if (ub){
        ParOptBVecWrap *ubwrap = dynamic_cast<ParOptBVecWrap*>(ub);
        if (ubwrap){
          setBVecFromLocalValues(upper, ubwrap->vec);
          ubwrap->vec->beginSetValues(TACS_INSERT_NONZERO_VALUES);
          ubwrap->vec->endSetValues(TACS_INSERT_NONZERO_VALUES);
        }
      }
      delete [] upper;
    }
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

    // If we're using inverse variables, convert from the inverse
    // variables to the usual variables
    if (use_inverse_vars){
      TacsScalar *xvals;
      int size = x[0]->getArray(&xvals);
      for ( int i = 0; i < size; i++ ){
        xvals[i] = 1.0/xvals[i];
      }
    }
  
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

    tacs[0]->zeroVariables();

    // Assemble the Jacobian on each level
    double alpha = 1.0, beta = 0.0, gamma = 0.0;
    mg->assembleJacobian(alpha, beta, gamma, res);
    mg->factor();
    
    // Solve the system: K(x)*u = -res
    ksm->solve(res, vars);
    vars->scale(-1.0);

    tacs[0]->applyBCs(vars);
    tacs[0]->setVariables(vars);
    
    // Evaluate the compliance
    TacsScalar compliance_value = -vars->dot(res);
  
    // Evaluate the mass
    TacsScalar mass_value;
    tacs[0]->evalFunctions(&mass, 1, &mass_value);
    
    // Set the scaling for the objective function
    if (obj_scale < 0.0){
      obj_scale = 1.0/compliance_value;
    }

    // Set the compliance objective and the mass constraint
    *fobj = obj_scale*compliance_value;
    cons[0] = (target_mass - mass_value)*mass_scale;
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
                                        ParOptVec *g, 
                                        ParOptVec **Ac ){
  // Evaluate the derivative of the mass and compliance with
  // respect to the design variables
  g->zeroEntries();
  Ac[0]->zeroEntries();

  // Compute the derivative of the mass and compliance w.r.t the
  // design variables
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
    g_vec->beginSetValues(TACS_ADD_VALUES);

    memset(xlocal, 0, size*sizeof(TacsScalar));
    tacs[0]->addDVSens(-mass_scale, &mass, 1, xlocal, size);
    setBVecFromLocalValues(xlocal, m_vec);
    m_vec->beginSetValues(TACS_ADD_VALUES);

    // Finsh adding the values
    g_vec->endSetValues(TACS_ADD_VALUES);
    m_vec->endSetValues(TACS_ADD_VALUES);

    // Check if we're using reciprocal variables
    if (use_inverse_vars){
      ParOptScalar *xvals;
      TacsScalar *gvals, *mvals;
      int xsize = xvec->getArray(&xvals);
      g_vec->getArray(&gvals);
      m_vec->getArray(&mvals);
      for ( int i = 0; i < xsize; i++ ){
        gvals[i] = -gvals[i]/(xvals[i]*xvals[i]);
        mvals[i] = -mvals[i]/(xvals[i]*xvals[i]);
      }
    }
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
                                     ParOptVec *hvec ){
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
}

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
void TMRTopoProblem::writeOutput( int iter, ParOptVec *xvec ){
  // Print out the binary STL file for later visualization
  int var_offset = 0;
  double cutoff = 0.25;
  char filename[256];
  sprintf(filename, "%s/levelset025_binary%04d.bstl", prefix, iter_count);
  TMR_GenerateBinFile(filename, filter[0], x[0], var_offset, cutoff);

  cutoff = 0.5;
  sprintf(filename, "%s/levelset05_binary%04d.bstl", prefix, iter_count);
  TMR_GenerateBinFile(filename, filter[0], x[0], var_offset, cutoff);

  // Update the iteration count
  iter_count++;
}

#endif // TMR_HAS_PAROPT
