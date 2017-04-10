#ifndef TACS_SOLID_SHELL_WRAPPER_H
#define TACS_SOLID_SHELL_WRAPPER_H

#include "TACSElement.h"
#include "MITCShell.h"

class SolidShellWrapper : public TACSElement {
 public:
  SolidShellWrapper( MITCShell<2> *_shell );
  ~SolidShellWrapper();
  
  const char * displacementName( int i );
  ElementType getElementType();
  int numDisplacements();
  int numNodes();
  int numStresses();

  void addResidual( double time, TacsScalar res[],
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] );
  void addJacobian( double time, TacsScalar J[],
                    double alpha, double beta, double gamma,
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] );
  void addAdjResProduct( double time, double scale,
                         TacsScalar dvSens[], int dvLen,
                         const TacsScalar psi[],
                         const TacsScalar Xpts[],
                         const TacsScalar vars[],
                         const TacsScalar dvars[],
                         const TacsScalar ddvars[] );

  TACSConstitutive * getConstitutive();
  int getNumGaussPts();
  double getGaussWtsPts( const int num, double * pt );
  TacsScalar getDetJacobian( const double * pt, 
                             const TacsScalar Xpts[] );
  void getStrain( TacsScalar strain[], 
                  const double pt[], 
                  const TacsScalar Xpts[],
                  const TacsScalar vars[] );
  void addStrainSVSens( TacsScalar strainSVSens[], 
                        const double pt[], 
                        const TacsScalar scale,
                        const TacsScalar strainSens[], 
                        const TacsScalar Xpts[],
                        const TacsScalar vars[] );


  void addOutputCount( int *nelems, int *nnodes, int *ncsr );
  void getOutputData( unsigned int out_type, 
                      double *data, int ld_data,
                      const TacsScalar Xpts[],
                      const TacsScalar vars[] );
  void getOutputConnectivity( int *con, int node );

 private:
  static const char *dispNames[6];

  MITCShell<2> *shell;
};

#endif // TACS_SOLID_SHELL_WRAPPER_H
