#include "TMRTopoOpt.h"
#include "Compliance.h"
#include "StructuralMass.h"

TMRTopoOpt::TMRTopoOpt( int _nlevels, 
                        TACSAssembler *_tacs[],
                        TMROctree *_filter[], 
                        TACSMg *_mg, BVec *_force,
                        double _target_mass,
                        TACSToFH5 *_f5,
                        const char *_prefix ):
ParOptProblem(_tacs[0]->getMPIComm()){
  // Copy over the prefix - max of 128 characters
  strncpy(prefix, _prefix, sizeof(prefix));
  
  // Set the number of levels
  nlevels = _nlevels;
  if (nlevels > MAX_NUM_LEVELS){
    nlevels = MAX_NUM_LEVELS;
  }

  // Copy over the assembler objects and filters
  for ( int k = 0; k < nlevels; k++ ){
    // Set the TACSAssembler objects for each level
    tacs[k] = _tacs[k];
    tacs[k]->incref();

    // Set the filter object
    filter[k] = _filter[k];
  
    // Get the number of design variables for each level
    filter[k]->getMesh(&num_design_vars[k], NULL, NULL, NULL);
    xglobal[k] = new TacsScalar[ num_design_vars[k] ];
  }

  // Set the temporary vector
  global = new TacsScalar[ num_design_vars[0] ];

  // The multigrid object
  mg = _mg;
  mg->incref();

  // The force vector
  force = _force;
  force->incref();

  // Set up the solver
  int gmres_iters = 100; 
  int nrestart = 2;
  int is_flexible = 1;
  ksm = new GMRES(mg->getMat(0), mg, 
                  gmres_iters, nrestart, is_flexible);
  ksm->incref();
  ksm->setTolerances(1e-8, 1e-30);

  // Allocate variables
  vars = tacs[0]->createVec();
  vars->incref();

  // Create the compliance and mass functions
  compliance = new Compliance(tacs[0]);
  mass = new StructuralMass(tacs[0]);
  compliance->incref();
  mass->incref();

  // Set the target mass
  target_mass = _target_mass;

  f5 = _f5;
  f5->incref();

  // Set the restriction operators
  for ( int k = 0; k < nlevels-1; k++ ){
    filter[k]->createRestriction(filter[k+1],
                                 &restrict_ptr[k], &restrict_conn[k],
                                 &restrict_weights[k]);
  }

  // Collect the design variable information
  int mpi_size, mpi_rank;
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);

  // Allocate the recv counts and recv displacements
  recvcount = new int[ mpi_size ];
  recvdisp = new int[ mpi_size ];

  // Set the number of variables for each processor
  int vars_per_proc = num_design_vars[0]/mpi_size;

  recvdisp[0] = 0;
  for ( int i = 0; i < mpi_size-1; i++ ){
    recvcount[i] = vars_per_proc;
    recvdisp[i+1] = recvdisp[i] + recvcount[i];
  }
  recvcount[mpi_size-1] = 
    num_design_vars[0] - vars_per_proc*(mpi_size-1);

  // Set the start offset
  start_index = recvdisp[mpi_rank];

  // Set the problem sizes
  const int ncon = 1;
  const int nwcon = 0, nwblock = 0;
  setProblemSizes(recvcount[mpi_rank], ncon, nwcon, nwblock);

  // Create a vector
  xinit = new ParOptVec(this->comm, this->nvars);
  xinit->set(0.95);
}

/*
  Free the data stored in the object
*/
TMRTopoOpt::~TMRTopoOpt(){
  for ( int k = 0; k < nlevels; k++ ){
    tacs[k]->decref();
    delete [] xglobal[k];
  }

  for ( int k = 0; k < nlevels-1; k++ ){
    delete [] restrict_ptr[k];
    delete [] restrict_conn[k];
    delete [] restrict_weights[k];
  }

  delete [] global;
  mg->decref();
  force->decref();
  ksm->decref();
  vars->decref();
  f5->decref();
  compliance->decref();
  mass->decref();
  delete xinit;
  delete [] recvcount;
  delete [] recvdisp;
}

/*
  Set the initial design variable components
*/
void TMRTopoOpt::setInitDesignVars( ParOptScalar x[] ){
  setLocalComponents(x, xinit);
}

/*
  Set whether these should be considered sparse inequalities
*/
int TMRTopoOpt::isSparseInequality(){ return 0; }
int TMRTopoOpt::isDenseInequality(){ return 1; }
int TMRTopoOpt::useLowerBounds(){ return 1; }
int TMRTopoOpt::useUpperBounds(){ return 1; }

/*
  Set the initial design variables
*/ 
void TMRTopoOpt::getVarsAndBounds( ParOptVec *x, 
                                   ParOptVec *lb, 
                                   ParOptVec *ub ){
  // Get the values of the design variables from the inner-most
  // version of TACS
  x->copyValues(xinit);

  // Set the lower and upper bounds on the design variables
  lb->set(0.0);
  ub->set(1.0);
}

/*
  Evaluate the objective and constraints
*/
int TMRTopoOpt::evalObjCon( ParOptVec *x, 
                            ParOptScalar *fobj, 
                            ParOptScalar *cons ){
  // Set the values of the design variables
  setGlobalComponents(x, xglobal[0]);
  
  // Set the deseign variables in TACS
  tacs[0]->setDesignVars(xglobal[0], num_design_vars[0]);
  
  // Transfer the design variables 
  for ( int k = 0; k < nlevels-1; k++ ){
    // Restrict desgin variables
    restrictDesignVars(k, xglobal[k], xglobal[k+1]);

    // Set the design variables for each level
    tacs[k+1]->setDesignVars(xglobal[k+1], num_design_vars[k+1]);
  }

  // Assemble the Jacobian on each level
  double alpha = 1.0, beta = 0.0, gamma = 0.0;
  mg->assembleJacobian(NULL, alpha, beta, gamma);
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
  
  // Set the compliance objective and the mass constraint
  *fobj = fvals[0];
  cons[0] = 1.0 - fvals[1]/target_mass;
  
  return 0;
}

/*
  Evaluate the objective and constraint gradients
*/
int TMRTopoOpt::evalObjConGradient( ParOptVec *x, 
                                    ParOptVec *g, ParOptVec **Ac ){
  // Evaluate the derivative of the mass and compliance with
  // respect to the design variables
  tacs[0]->evalDVSens(&mass, 1, global, num_design_vars[0]);
  
  // Set the local components of the compliance gradient
  setLocalComponents(global, Ac[0]);
  Ac[0]->scale(-1.0/target_mass);

  tacs[0]->evalAdjointResProducts(&vars, 1,
                                  global, num_design_vars[0]);

  // Set the local components of the objective gradient
  setLocalComponents(global, g);
  g->scale(-1.0);

  return 0;
}

/*
  Evaluate the product of the Hessian with the given vector px
*/
int TMRTopoOpt::evalHvecProduct( ParOptVec *xvec,
                                 ParOptScalar *z, 
                                 ParOptVec *zw,
                                 ParOptVec *pxvec, 
                                 ParOptVec *hvec ){}

/*
  Given the input vector that is the same on all processors,
  set the local vector components in the ParOpt object
*/
void TMRTopoOpt::setLocalComponents( ParOptScalar *xin, 
                                     ParOptVec *xout ){
  double *x;
  xout->getArray(&x);
  memcpy(x, &xin[start_index], nvars*sizeof(double)); 
}

/*
  Get the distributed vector, set the values on all other processor
*/
void TMRTopoOpt::setGlobalComponents( ParOptVec *xin, 
                                      ParOptScalar *xout ){
  double *x;
  xin->getArray(&x);
  MPI_Allgatherv(x, nvars, MPI_DOUBLE, 
                 xout, recvcount, recvdisp, MPI_DOUBLE, comm);
}

// Evaluate the sparse constraints
// ------------------------
void TMRTopoOpt::evalSparseCon( ParOptVec *x, 
                                ParOptVec *out ){}

// Compute the Jacobian-vector product out = J(x)*px
// --------------------------------------------------
void TMRTopoOpt::addSparseJacobian( double alpha, ParOptVec *x,
                                    ParOptVec *px, ParOptVec *out ){}

// Compute the transpose Jacobian-vector product out = J(x)^{T}*pzw
// -----------------------------------------------------------------
void TMRTopoOpt::addSparseJacobianTranspose( double alpha, 
                                             ParOptVec *x,
                                             ParOptVec *pzw, 
                                             ParOptVec *out ){}

// Add the inner product of the constraints to the matrix such 
// that A += J(x)*cvec*J(x)^{T} where cvec is a diagonal matrix
// ------------------------------------------------------------
void TMRTopoOpt::addSparseInnerProduct( double alpha, 
                                        ParOptVec *x,
                                        ParOptVec *cvec, 
                                        double *A ){}

// Restrict the design variables to the next lowest level
// ------------------------------------------------------
void TMRTopoOpt::restrictDesignVars( int lev, 
                                     ParOptScalar *xfine,
                                     ParOptScalar *xcoarse ){
  // Pull out the restriction data 
  const int len = num_design_vars[lev+1];
  const int *ptr = restrict_ptr[lev];
  const int *conn = restrict_conn[lev];
  const double *weights = restrict_weights[lev];
  
  // Compute the new design variables on the coarser mesh
  for ( int i = 0; i < len; i++ ){
    xcoarse[0] = 0.0;
    for ( int jp = ptr[i]; jp < ptr[i+1]; jp++ ){
      xcoarse[0] += weights[jp]*xfine[conn[0]];
      conn++;
    }
    xcoarse++;
  }
}

// Write the output file
void TMRTopoOpt::writeOutput( int iter, ParOptVec *x ){
  char filename[256];
  sprintf(filename, "%soctree_iter%d.f5", prefix, iter);
  f5->writeToFile(filename);
}
