#ifndef TACS_POISSON_GRAD_FUNCTION_H
#define TACS_POISSON_GRAD_FUNCTION_H

#include "TACSFunction.h"

class TACSPoissonGradFunction : public TACSFunction {
 public:
  enum PoissonGradFunctionType { DISCRETE, CONTINUOUS,
                                 PNORM_DISCRETE, PNORM_CONTINUOUS };
  
  TACSPoissonGradFunction( TACSAssembler *_tacs,
                           TacsScalar _ksWeight,
                           const TacsScalar _dir[],
                           PoissonGradFunctionType ksType=CONTINUOUS );
  ~TACSPoissonGradFunction();
  
  const char *functionName();
  void setGradFunctionType( PoissonGradFunctionType _ksType );
  TACSFunctionCtx *createFunctionCtx();

  void initEvaluation( EvaluationType ftype );
  void finalEvaluation( EvaluationType ftype );
  void initThread( double tcoef, EvaluationType ftype,
                   TACSFunctionCtx *ctx );
  void elementWiseEval( EvaluationType ftype,
                        TACSElement *element, int elemNum,
                        const TacsScalar Xpts[], const TacsScalar vars[],
                        const TacsScalar dvars[], const TacsScalar ddvars[],
                        TACSFunctionCtx *ctx );
  void finalThread( double tcoef, EvaluationType ftype,
                    TACSFunctionCtx *ctx );

  TacsScalar getFunctionValue();

  void getElementSVSens( double alpha, double beta, double gamma, 
                         TacsScalar *elemSVSens, 
                         TACSElement *element, int elemNum,
                         const TacsScalar Xpts[], const TacsScalar vars[],
                         const TacsScalar dvars[], const TacsScalar ddvars[],
                         TACSFunctionCtx *ctx );

 private:  
  // The name of the function
  static const char *funcName;

  // The direction
  TacsScalar dir[3];

  // The value of the KS weight
  double ksWeight;

  // Intermediate values in the functional evaluation
  TacsScalar ksSum;
  TacsScalar maxValue;
  TacsScalar invPnorm;

  // Set the type of constraint aggregate
  PoissonGradFunctionType ksType;

  // The max number of nodes
  int maxNumNodes;
};

#endif // TACS_POISSON_GRAD_FUNCTION_H
