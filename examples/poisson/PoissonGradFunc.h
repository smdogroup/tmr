#ifndef TACS_POISSON_GRAD_FUNCTION_H
#define TACS_POISSON_GRAD_FUNCTION_H

#include "TACSFunction.h"

class TACSPoissonGradFunction : public TACSFunction {
 public:
  enum PoissonGradFunctionType {
    DISCRETE,
    CONTINUOUS,
    PNORM_DISCRETE,
    PNORM_CONTINUOUS
  };

  TACSPoissonGradFunction(TACSAssembler *_tacs, TacsScalar _ksWeight,
                          const TacsScalar _dir[],
                          PoissonGradFunctionType ksType = CONTINUOUS);
  ~TACSPoissonGradFunction();

  const char *getObjectName();
  void setGradFunctionType(PoissonGradFunctionType _ksType);

  void initEvaluation(EvaluationType ftype);
  void finalEvaluation(EvaluationType ftype);
  void elementWiseEval(EvaluationType ftype, int elemIndex,
                       TACSElement *element, double time, TacsScalar scale,
                       const TacsScalar Xpts[], const TacsScalar vars[],
                       const TacsScalar dvars[], const TacsScalar ddvars[]);

  TacsScalar getFunctionValue();

  oid getElementSVSens(int elemIndex, TACSElement *element, double time,
                       TacsScalar alpha, TacsScalar beta, TacsScalar gamma,
                       const TacsScalar Xpts[], const TacsScalar vars[],
                       const TacsScalar dvars[], const TacsScalar ddvars[],
                       TacsScalar dfdu[]);

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

#endif  // TACS_POISSON_GRAD_FUNCTION_H
