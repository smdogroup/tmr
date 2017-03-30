#ifndef TMR_TOPO_PROBLEM_H
#define TMR_TOPO_PROBLEM_H

#ifdef TMR_HAS_PAROPT

#include "ParOpt.h"
#include "TACSAssembler.h"
#include "TACSMg.h"
#include "TMROctForest.h"
#include "StructuralMass.h"
#include "Compliance.h"

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
  static const int MAX_NUM_LEVELS = 5;

  // Create the topology optimization object
  // ---------------------------------------
  TMRTopoProblem( int _nlevels, 
                  TACSAssembler *_tacs[],
                  TMROctForest *_filter[], 
                  TACSVarMap *_filter_maps[],
                  TACSBVecIndices *_filter_indices[],
                  TACSMg *_mg,
                  double _target_mass,
                  const char *_prefix );
  ~TMRTopoProblem();

  // Set/get the objective and mass constraint scaling factors
  // ---------------------------------------------------------
  ParOptScalar getObjectiveScaling();
  void setObjectiveScaling( ParOptScalar scale );
  ParOptScalar getMassScaling();
  void setMassScaling( ParOptScalar scale );

  // Use the reciprocal variable values
  // ----------------------------------
  void setUseReciprocalVariables(){
    use_inverse_vars = 1;
  }

  // Set the linearization point and penalty parameter
  // -------------------------------------------------
  void setLinearization( double q, ParOptVec *xvec );

  // Set the iteration count
  // -----------------------
  void setIterationCount( int iter ){
    iter_count = iter;
  }

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

  // Store the prefix
  char prefix[256];

  // Set the iteration count for printing to the file
  int iter_count;

  // Store the design variable info
  int nlevels;
  TACSAssembler *tacs[MAX_NUM_LEVELS];

  // Set the information about the filter at each level
  TMROctForest *filter[MAX_NUM_LEVELS];
  TACSVarMap *filter_maps[MAX_NUM_LEVELS];
  TACSBVecIndices *filter_indices[MAX_NUM_LEVELS];
  TACSBVecDistribute *filter_dist[MAX_NUM_LEVELS];
  TACSBVecDistCtx *filter_ctx[MAX_NUM_LEVELS];
  TACSBVecInterp *filter_interp[MAX_NUM_LEVELS-1];

  // Create the design variable values at each level
  TACSBVec *x[MAX_NUM_LEVELS];

  // The initial design variable values
  int max_local_size;
  TacsScalar *xlocal;

  // Store the Krylov solver and the multigrid object 
  TACSKsm *ksm;
  TACSMg *mg;

  // The compliance and mass functions
  TACSFunction *compliance, *mass;

  // The target mass
  double target_mass;

  // The scaling for the objective and constraint
  double obj_scale;
  double mass_scale;

  // Flag to indicate whether to use inverse variables or not
  int use_inverse_vars;

  // The load case information
  TACSBVec *res, *vars, *force;
};

#endif // TMR_HAS_PAROPT
#endif // TMR_TOPO_PROBLEM_H
