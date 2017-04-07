#ifndef TMR_TACS_TOPO_CREATOR_H
#define TMR_TACS_TOPO_CREATOR_H

#include "TMR_TACSCreator.h"
#include "TACSAssembler.h"
#include "TMROctStiffness.h"
#include "TACSElement.h"
#include "MITCShell.h"

class SolidShellWrapper : public TACSElement {
 public:
  SolidShellWrapper( MITCShell<2> *_shell ){
    shell = _shell;
    shell->incref();
  }
  ~SolidShellWrapper(){
    shell->decref();
  }
  int numDisplacements(){ return 3; }
  int numNodes(){ return 8; }
  void addResidual( double time, TacsScalar res[],
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] ){
    TacsScalar X[3*4];
    X[0] = Xpts[0];   X[1] = Xpts[1];   X[2] = Xpts[2];
    X[3] = Xpts[6];   X[4] = Xpts[7];   X[5] = Xpts[8];
    X[6] = Xpts[12];  X[7] = Xpts[13];  X[8] = Xpts[14];
    X[9] = Xpts[18];  X[10] = Xpts[19]; X[11] = Xpts[20];
    shell->addResidual(time, res, X, vars, dvars, ddvars);
  }
  virtual void addJacobian( double time, TacsScalar J[],
                            double alpha, double beta, double gamma,
                            const TacsScalar Xpts[],
                            const TacsScalar vars[],
                            const TacsScalar dvars[],
                            const TacsScalar ddvars[] ){
    TacsScalar X[3*4];
    X[0] = Xpts[0];   X[1] = Xpts[1];   X[2] = Xpts[2];
    X[3] = Xpts[6];   X[4] = Xpts[7];   X[5] = Xpts[8];
    X[6] = Xpts[12];  X[7] = Xpts[13];  X[8] = Xpts[14];
    X[9] = Xpts[18];  X[10] = Xpts[19]; X[11] = Xpts[20];
    shell->addJacobian(time, J, alpha, beta, gamma, 
                       X, vars, dvars, ddvars);
  }
  void addAdjResProduct( double time, double scale,
                         TacsScalar dvSens[], int dvLen,
                         const TacsScalar psi[],
                         const TacsScalar Xpts[],
                         const TacsScalar vars[],
                         const TacsScalar dvars[],
                         const TacsScalar ddvars[] ){
    TacsScalar X[3*4];
    X[0] = Xpts[0];   X[1] = Xpts[1];   X[2] = Xpts[2];
    X[3] = Xpts[6];   X[4] = Xpts[7];   X[5] = Xpts[8];
    X[6] = Xpts[12];  X[7] = Xpts[13];  X[8] = Xpts[14];
    X[9] = Xpts[18];  X[10] = Xpts[19]; X[11] = Xpts[20];
    shell->addAdjResProduct(time, scale, dvSens, dvLen,
                            psi, X, vars, dvars, ddvars);
  }
  
 private:
  MITCShell<2> *shell;
};


class TMROctTACSTopoCreator : public TMROctTACSCreator {
 public:
  TMROctTACSTopoCreator( TMRBoundaryConditions *_bcs,
                         TMRStiffnessProperties _properties,
                         TMROctForest *_filter,
                         const char *shell_attr=NULL, 
                         SolidShellWrapper *_shell=NULL );
  ~TMROctTACSTopoCreator();

  // Create the connectivity
  void createConnectivity( int order,
                           TMROctForest *forest,
                           int **_conn, int **_ptr,
                           int *_num_elements );

  // Create the elements
  void createElements( int order,
                       TMROctForest *forest,
                       int num_elements,
                       TACSElement **elements );

  // Get the underlying objects that define the filter
  void getForest( TMROctForest **filter );
  void getMap( TACSVarMap **_map );
  void getIndices( TACSBVecIndices **_indices );  

 private:
  // Compute the weights for a given point
  void computeWeights( TMROctant *oct, TMROctant *node,
                       TMRIndexWeight *welem );

  // The stiffness properties
  TMRStiffnessProperties properties;

  // The forest that defines the filter
  TMROctForest *filter;

  // The filter map for this object. This defines how the design
  // variables are distributed across all of the processors.
  TACSVarMap *filter_map;

  // The filter indices. This defines the relationship between the
  // local design variable numbers and the global design variable
  // numbers.
  TACSBVecIndices *filter_indices;

  // Set the top/bottom attributes
  char *shell_attr;
  SolidShellWrapper *shell;
};

#endif // TMR_TACS_TOPO_CREATOR_H

