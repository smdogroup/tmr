#include "PoissonGradFunc.h"

#include "FElibrary.h"
#include "PoissonElement.h"
#include "TACSAssembler.h"

/*
  The context for the PoissonGradFuncCtx class
*/
class PoissonGradFuncCtx : public TACSFunctionCtx {
 public:
  PoissonGradFuncCtx(int maxNodes) {
    // Allocate the working array
    N = new double[maxNodes];
    Na = new double[maxNodes];
    Nb = new double[maxNodes];
    ksSum = 0.0;
    maxValue = -1e20;
  }
  ~PoissonGradFuncCtx() {
    delete[] N;
    delete[] Na;
    delete[] Nb;
  }

  // Data to be used for the function computation
  TacsScalar ksSum;
  TacsScalar maxValue;
  double *N, *Na, *Nb;
};

TACSPoissonGradFunction::TACSPoissonGradFunction(
    TACSAssembler *_tacs, TacsScalar _ksWeight, const TacsScalar _dir[],
    PoissonGradFunctionType _ksType)
    : TACSFunction(_tacs, TACSFunction::ENTIRE_DOMAIN,
                   TACSFunction::TWO_STAGE) {
  maxNumNodes = _tacs->getMaxElementNodes();
  dir[0] = _dir[0];
  dir[1] = _dir[1];
  ksWeight = _ksWeight;
  ksType = _ksType;
  maxValue = -1e20;
  ksSum = 0.0;
  invPnorm = 0.0;
}

TACSFunctionCtx *TACSPoissonGradFunction::createFunctionCtx() {
  return new PoissonGradFuncCtx(maxNumNodes);
}

TACSPoissonGradFunction::~TACSPoissonGradFunction() {}

/*
  TACSPoissonGradFunction function name
*/
const char *TACSPoissonGradFunction::funcName = "TACSPoissonGradFunction";

/*
  Return the function name
*/
const char *TACSPoissonGradFunction::functionName() { return funcName; }

/*
  Set the displacement aggregate type
*/
void TACSPoissonGradFunction::setGradFunctionType(
    PoissonGradFunctionType _ksType) {
  ksType = _ksType;
}

TacsScalar TACSPoissonGradFunction::getFunctionValue() {
  // Compute the final value of the KS function on all processors
  if (ksType == CONTINUOUS || ksType == DISCRETE) {
    return maxValue + log(ksSum) / ksWeight;
  } else {
    return maxValue * pow(ksSum, 1.0 / ksWeight);
  }
}

void TACSPoissonGradFunction::initEvaluation(EvaluationType ftype) {
  if (ftype == TACSFunction::INITIALIZE) {
    maxValue = -1e20;
  } else if (ftype == TACSFunction::INTEGRATE) {
    ksSum = 0.0;
  }
}

void TACSPoissonGradFunction::finalEvaluation(EvaluationType ftype) {
  if (ftype == TACSFunction::INITIALIZE) {
    // Distribute the values of the KS function computed on this domain
    TacsScalar temp = maxValue;
    MPI_Allreduce(&temp, &maxValue, 1, TACS_MPI_TYPE, TACS_MPI_MAX,
                  tacs->getMPIComm());
  } else {
    // Find the sum of the ks contributions from all processes
    TacsScalar temp = ksSum;
    MPI_Allreduce(&temp, &ksSum, 1, TACS_MPI_TYPE, MPI_SUM, tacs->getMPIComm());

    // Compute the P-norm quantity if needed
    invPnorm = 0.0;
    if (ksType == PNORM_DISCRETE || ksType == PNORM_CONTINUOUS) {
      if (ksSum != 0.0) {
        invPnorm = pow(ksSum, (1.0 - ksWeight) / ksWeight);
      }
    }
  }
}

void TACSPoissonGradFunction::initThread(double tcoef, EvaluationType ftype,
                                         TACSFunctionCtx *fctx) {
  PoissonGradFuncCtx *ctx = dynamic_cast<PoissonGradFuncCtx *>(fctx);
  if (ctx) {
    if (ftype == TACSFunction::INITIALIZE) {
      ctx->maxValue = -1e20;
      ctx->ksSum = 0.0;
    }
  }
}

void TACSPoissonGradFunction::elementWiseEval(
    EvaluationType ftype, TACSElement *element, int elemNum,
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TACSFunctionCtx *fctx) {
  PoissonGradFuncCtx *ctx = dynamic_cast<PoissonGradFuncCtx *>(fctx);
  PoissonQuad<2> *quad2 = dynamic_cast<PoissonQuad<2> *>(element);
  PoissonQuad<3> *quad3 = dynamic_cast<PoissonQuad<3> *>(element);
  PoissonQuad<4> *quad4 = dynamic_cast<PoissonQuad<4> *>(element);
  PoissonQuad<5> *quad5 = dynamic_cast<PoissonQuad<5> *>(element);

  if (ctx && (quad2 || quad3 || quad4 || quad5)) {
    // Get the number of quadrature points for this element
    const int numGauss = element->getNumGaussPts();
    const int numNodes = element->numNodes();

    // With the first iteration, find the maximum over the domain
    for (int i = 0; i < numGauss; i++) {
      // Get the Gauss points one at a time
      double pt[3];
      double weight = element->getGaussWtsPts(i, pt);

      // Compute the Jacobian transformation
      TacsScalar Xd[4];

      // Get the shape functions
      if (quad2) {
        quad2->getShapeFunctions(pt, ctx->N, ctx->Na, ctx->Nb);
        quad2->getJacobianTransform(ctx->Na, ctx->Nb, Xpts, Xd);
      } else if (quad3) {
        quad3->getShapeFunctions(pt, ctx->N, ctx->Na, ctx->Nb);
        quad3->getJacobianTransform(ctx->Na, ctx->Nb, Xpts, Xd);
      } else if (quad4) {
        quad4->getShapeFunctions(pt, ctx->N, ctx->Na, ctx->Nb);
        quad4->getJacobianTransform(ctx->Na, ctx->Nb, Xpts, Xd);
      } else {
        quad5->getShapeFunctions(pt, ctx->N, ctx->Na, ctx->Nb);
        quad5->getJacobianTransform(ctx->Na, ctx->Nb, Xpts, Xd);
      }

      // Compute the inverse of the Jacobian
      TacsScalar J[4];
      TacsScalar h = FElibrary::jacobian2d(Xd, J);

      TacsScalar px = 0.0, py = 0.0;
      const double *na = ctx->Na;
      const double *nb = ctx->Nb;
      for (int j = 0; j < numNodes; j++) {
        // Compute the derivative of the shape functions
        TacsScalar Nx = na[0] * J[0] + nb[0] * J[2];
        TacsScalar Ny = na[0] * J[1] + nb[0] * J[3];
        na++;
        nb++;

        px += Nx * vars[j];
        py += Ny * vars[j];
      }

      // Compute the scalar value
      TacsScalar value = px * dir[0] + py * dir[1];

      if (ftype == TACSFunction::INITIALIZE) {
        if (ksType == CONTINUOUS || ksType == DISCRETE) {
          if (TacsRealPart(value) > TacsRealPart(ctx->maxValue)) {
            ctx->maxValue = value;
          }
        } else {
          if (fabs(TacsRealPart(value)) > TacsRealPart(ctx->maxValue)) {
            ctx->maxValue = fabs(TacsRealPart(value));
          }
        }
      } else {
        if (ksType == CONTINUOUS) {
          ctx->ksSum += h * weight * exp(ksWeight * (value - maxValue));
        } else if (ksType == DISCRETE) {
          ctx->ksSum += exp(ksWeight * (value - maxValue));
        } else if (ksType == PNORM_CONTINUOUS) {
          ctx->ksSum +=
              h * weight * pow(fabs(TacsRealPart(value / maxValue)), ksWeight);
        } else if (ksType == PNORM_DISCRETE) {
          ctx->ksSum += pow(fabs(TacsRealPart(value / maxValue)), ksWeight);
        }
      }
    }
  }
}

void TACSPoissonGradFunction::finalThread(double tcoef, EvaluationType ftype,
                                          TACSFunctionCtx *fctx) {
  PoissonGradFuncCtx *ctx = dynamic_cast<PoissonGradFuncCtx *>(fctx);
  if (ctx) {
    if (ftype == TACSFunction::INITIALIZE) {
      if (ksType == CONTINUOUS || ksType == DISCRETE) {
        if (TacsRealPart(ctx->maxValue) > TacsRealPart(maxValue)) {
          maxValue = ctx->maxValue;
        }
      } else {
        if (fabs(TacsRealPart(ctx->maxValue)) > TacsRealPart(maxValue)) {
          maxValue = fabs(TacsRealPart(ctx->maxValue));
        }
      }
    } else {
      ksSum += tcoef * ctx->ksSum;
    }
  }
}

void TACSPoissonGradFunction::getElementSVSens(
    double alpha, double beta, double gamma, TacsScalar *elemSVSens,
    TACSElement *element, int elemNum, const TacsScalar Xpts[],
    const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TACSFunctionCtx *fctx) {
  PoissonGradFuncCtx *ctx = dynamic_cast<PoissonGradFuncCtx *>(fctx);
  PoissonQuad<2> *quad2 = dynamic_cast<PoissonQuad<2> *>(element);
  PoissonQuad<3> *quad3 = dynamic_cast<PoissonQuad<3> *>(element);
  PoissonQuad<4> *quad4 = dynamic_cast<PoissonQuad<4> *>(element);
  PoissonQuad<5> *quad5 = dynamic_cast<PoissonQuad<5> *>(element);

  // Zero the derivative of the function w.r.t. the element state
  // variables
  int numVars = element->numVariables();
  memset(elemSVSens, 0, numVars * sizeof(TacsScalar));

  if (ctx && (quad2 || quad3 || quad4 || quad5)) {
    // Get the number of quadrature points for this element
    const int numGauss = element->getNumGaussPts();
    const int numNodes = element->numNodes();

    // With the first iteration, find the maximum over the domain
    for (int i = 0; i < numGauss; i++) {
      // Get the Gauss points one at a time
      double pt[3];
      double weight = element->getGaussWtsPts(i, pt);

      // Compute the Jacobian transformation
      TacsScalar Xd[4];

      // Get the shape functions
      if (quad2) {
        quad2->getShapeFunctions(pt, ctx->N, ctx->Na, ctx->Nb);
        quad2->getJacobianTransform(ctx->Na, ctx->Nb, Xpts, Xd);
      } else if (quad3) {
        quad3->getShapeFunctions(pt, ctx->N, ctx->Na, ctx->Nb);
        quad3->getJacobianTransform(ctx->Na, ctx->Nb, Xpts, Xd);
      } else if (quad4) {
        quad4->getShapeFunctions(pt, ctx->N, ctx->Na, ctx->Nb);
        quad4->getJacobianTransform(ctx->Na, ctx->Nb, Xpts, Xd);
      } else {
        quad5->getShapeFunctions(pt, ctx->N, ctx->Na, ctx->Nb);
        quad5->getJacobianTransform(ctx->Na, ctx->Nb, Xpts, Xd);
      }

      // Compute the inverse of the Jacobian
      TacsScalar J[4];
      TacsScalar h = FElibrary::jacobian2d(Xd, J);

      TacsScalar px = 0.0, py = 0.0;
      const double *na = ctx->Na;
      const double *nb = ctx->Nb;
      for (int j = 0; j < numNodes; j++) {
        // Compute the derivative of the shape functions
        TacsScalar Nx = na[0] * J[0] + nb[0] * J[2];
        TacsScalar Ny = na[0] * J[1] + nb[0] * J[3];
        na++;
        nb++;

        px += Nx * vars[j];
        py += Ny * vars[j];
      }

      // Compute the scalar value
      TacsScalar value = px * dir[0] + py * dir[1];
      TacsScalar ptWeight = 0.0;

      if (ksType == CONTINUOUS) {
        ptWeight =
            alpha * h * weight * exp(ksWeight * (value - maxValue)) / ksSum;
      } else if (ksType == DISCRETE) {
        ptWeight = alpha * exp(ksWeight * (value - maxValue)) / ksSum;
      } else if (ksType == PNORM_CONTINUOUS) {
        ptWeight =
            value * pow(fabs(TacsRealPart(value / maxValue)), ksWeight - 2.0);
        ptWeight *= alpha * h * weight * invPnorm;
      } else if (ksType == PNORM_DISCRETE) {
        ptWeight =
            value * pow(fabs(TacsRealPart(value / maxValue)), ksWeight - 2.0);
        ptWeight *= alpha * ksWeight * invPnorm;
      }

      // Add the contribution to the sensitivity vector
      TacsScalar *s = elemSVSens;
      na = ctx->Na;
      nb = ctx->Nb;
      for (int j = 0; j < numNodes; j++) {
        TacsScalar Nx = na[0] * J[0] + nb[0] * J[2];
        TacsScalar Ny = na[0] * J[1] + nb[0] * J[3];
        s[0] += ptWeight * (dir[0] * Nx + dir[1] * Ny);
        s++;
        na++;
        nb++;
      }
    }
  }
}
