#ifndef TMR_TOPO_PROBLEM_H
#define TMR_TOPO_PROBLEM_H

#include "ParOpt.h"
#include "TACSAssembler.h"
#include "TACSMg.h"
#include "TMROctForest.h"
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
                  TACSMg *_mg );
  ~TMRTopoProblem();

  // Set the load cases - note that this destroys internal information stored
  // in the load case data associated with the constraints.
  // ------------------------------------------------------------------------
  void setLoadCases( TACSBVec **_forces, int _num_load_cases );
  int getNumLoadCases();

  // Add constraints associated with one of the load cases
  // -----------------------------------------------------
  void addConstraints( int _load_case, TACSFunction **_funcs,
                       const TacsScalar *_func_offset, 
                       const TacsScalar *_func_scale,
                       int num_funcs );
  
  void addConstraints( int _load_case, int _buckling,
                       int _freq, double sigma, int _num_eigvals,
                       TacsScalar _func_offset,
                       TacsScalar _func_scale );

  // Set the objective - in this case either compliance or a function
  // for one of the load cases
  // ----------------------------------------------------------------------
  void setObjective( const TacsScalar *_obj_weights, TACSFunction **funcs );
  void setObjective( const TacsScalar *_obj_weights );

  // Finish the initialization tasks - assign the number of constraints
  // and variables in the problem. Allocate arrays etc.
  // ------------------------------------------------------------------
  void initialize();

  // Set the output prefix for files
  // -------------------------------
  void setPrefix( const char *prefix );

  // Set the initial design variable values
  // --------------------------------------
  void setInitDesignVars( ParOptVec *vars );

  // Set the output iteration counter
  // --------------------------------
  void setIterationCounter( int iter );

  // Set the linearization point and penalty parameter
  // -------------------------------------------------
  void setLinearization( double q, ParOptVec *xvec );

  // Create a TACSBVec object containing the filter element volumes
  // --------------------------------------------------------------
  TACSBVec* createVolumeVec();

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
  char *prefix;

  // Set the iteration count for printing to the file
  int iter_count;

  // Set the load case information. In this case, these are force
  // vectors for each load case
  int num_load_cases;
  TACSBVec **vars, **forces;

  // The buckling and natural frequency TACS object
  TACSLinearBuckling **buckling;
  TACSFrequencyAnalysis **freq;
  int reset_count;
  TacsScalar *ks_a;
  // The derivative of f(x,u) w.r.t. u and the adjoint variables
  TACSBVec *dfdu, *adjoint;

  // The objective weights
  TacsScalar *obj_weights;
  TACSFunction **obj_funcs;

  // Set the constraint information for each load case
  class LoadCaseInfo {
   public:
    int num_funcs;
    TacsScalar *offset;
    TacsScalar *scale;
    TACSFunction **funcs;
    int num_eigvals;
    TacsScalar bf_offset, bf_scale;
  } *load_case_info;

  // Store the design variable info
  int nlevels;
  TACSAssembler **tacs;

  // Set the information about the filter at each level
  TMROctForest **filter;
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
};

#endif // TMR_TOPO_PROBLEM_H
