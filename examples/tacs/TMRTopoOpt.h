#ifndef TMR_TOPO_OPT_H
#define TMR_TOPO_OPT_H

#include "ParOpt.h"
#include "TACSAssembler.h"
#include "TACSToFH5.h"
#include "TACSMg.h"
#include "TMROctree.h"

class TMRTopoOpt : public ParOptProblem {
 public:
  static const int MAX_NUM_LEVELS = 5;

  // Create the topology optimization object
  // ---------------------------------------
  TMRTopoOpt( int nlevels, 
              TACSAssembler *_tacs[],
              TMROctree *_filter[], 
              TACSMg *_mg, BVec *_force );
  ~TMRTopoOpt();

  // Set the inequality flags
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

  // Set the local components 
  // ------------------------
  void setLocalComponents( ParOptScalar *xin, 
                           ParOptVec *xout );

  // Set the global components
  // -------------------------
  void setGlobalComponents( ParOptVec *xin, 
                            ParOptScalar *xout );

  // Write the output file
  // ---------------------
  void writeOutput( int iter, ParOptVec *x );

 private:  
  // Store the design variable info
  int nlevels;
  TACSAssembler *tacs[MAX_NUM_LEVELS];
  int num_design_vars[MAX_NUM_LEVELS];
  ParOptScalar *xglobal[MAX_NUM_LEVELS];
  ParOptScalar *global;

  // Set the filter object
  TMROctree *filter[MAX_NUM_LEVELS];

  // Set the restriction operators
  int *restrict_ptr[MAX_NUM_LEVELS];
  int *restrict_conn[MAX_NUM_LEVELS];
  double *restrict_weights[MAX_NUM_LEVELS];

  // Restrict the design variables to the next lowest level
  // ------------------------------------------------------
  void restrictDesignVars( int lev, ParOptScalar *xfine,
                           ParOptScalar *xcoarse );

  // Store the Krylov solver and the multigrid object 
  TACSKsm *ksm;
  TACSMg *mg;

  // The compliance and mass functions
  TACSFunction *compliance, *mass;

  // The target mass
  double target_mass;

  // The load case information
  BVec *res, *vars, *force;

  // The start index for the number of local vars in the global
  // design variable vector
  int start_index;
  int *recvcount, *recvdisp;

  // The output variables
  TACSToFH5 *f5;
};

#endif // TMR_TOPO_OPT_H
