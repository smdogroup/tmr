/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2019 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#include "TMRHeatFlux.h"
#include "TACSAssembler.h"
#include "CoupledThermoSolid.h"
#include "PlaneStressCoupledThermoQuad.h"
#include "CoupledThermoSolidStiffness.h"
#include "CoupledThermoPlaneStressStiffness.h"

/*
  The context for the TMRHeatFlux function
*/
class HeatFluxIntCtx : public TACSFunctionCtx {
 public:
  HeatFluxIntCtx( TACSFunction *func,
                  int maxNodes ){
    // Allocate the working array
    value = 0.0;
    // Allocate the working array
    work = new TacsScalar[2*maxStrains + 3*maxNodes];

    // Set the pointers into the work array
    strain = &work[0];
    failSens = &work[maxStrains];
    hXptSens = &work[2*maxStrains];
    
  }
  ~HeatFluxIntCtx(){
    delete [] N;
  }

  // Data to be used for the function computation
  TacsScalar value;
  TacsScalar maxFail;
  TacsScalar ksFailSum;
  TacsScalar *strain;
  TacsScalar *failSens;
  TacsScalar *hXptSens;
  TacsScalar *work;
};

TMRHeatFluxIntegral::TMRHeatFluxIntegral( TACSAssembler *_tacs,
                                          const char *_name,
                                          TMROctForest *_oforest,
                                          TMRQuadForest *_qforest ):
TACSFunction(_tacs, TACSFunction::ENTIRE_DOMAIN,
             TACSFunction::SINGLE_STAGE, 0){
  maxNumNodes = _tacs->getMaxElementNodes();
  name = _name;
  value = 0.0;
  oforest = _oforest;
  qforest = _qforest;
  
}
TMRHeatFluxIntegral::~TMRHeatFluxIntegral(){}

/*
  TMRHeatFluxIntegral function name
*/
const char * TMRHeatFluxIntegral::funcName = "TMRHeatFluxIntegral";

/*
  Return the function name
*/
const char *TMRHeatFluxIntegral::functionName(){
  return funcName;
}

/*
  Retrieve the function value
*/
TacsScalar TMRHeatFluxIntegral::getFunctionValue(){
  return value;
}

/*
  Allocate and return the function-specific context
*/
TACSFunctionCtx *TMRHeatFluxIntegral::createFunctionCtx(){
  return new HeatFluxIntCtx(this, maxNumNodes);
}

/*
  Initialize the internal values stored within the KS function
*/
void TMRHeatFluxIntegral::initEvaluation( EvaluationType ftype ){
  value = 0.0;
}

/*
  Reduce the function values across all MPI processes
*/
void TMRHeatFluxIntegral::finalEvaluation( EvaluationType ftype ){
  TacsScalar temp = value;
  MPI_Allreduce(&temp, &value, 1, TACS_MPI_TYPE,
                MPI_SUM, tacs->getMPIComm());
}

/*
  Initialize the context for either integration or initialization
*/
void TMRHeatFluxIntegral::initThread( const double tcoef,
                                      EvaluationType ftype,
                                      TACSFunctionCtx *fctx ){
  HeatFluxCtx *ctx = dynamic_cast<HeatFluxIntCtx*>(fctx);
  if (ctx){
    ctx->value = 0.0;
  }
}

/*
  Perform the element-wise evaluation of the TACSDisplacementIntegral function.
*/
void TMRHeatFluxIntegral::elementWiseEval( EvaluationType ftype,
                                           TACSElement *element,
                                           int elemNum,
                                           const TacsScalar Xpts[],
                                           const TacsScalar vars[],
                                           const TacsScalar dvars[],
                                           const TacsScalar ddvars[],
                                           TACSFunctionCtx *fctx ){
  HeatFluxIntCtx *ctx = dynamic_cast<HeatFluxIntCtx*>(fctx);
  if (ctx){
    // Get the number of quadrature points for this element
    const int numGauss = element->getNumGaussPts();
    const int numDisps = element->numDisplacements();
    const int numNodes = element->numNodes();
    
    int numGauss = element->getNumGaussPts();
    double N[numNodes];
    int order = 0;
    // Get the constitutive object for this element
    TACSConstitutive *con = NULL;
    if (numDisps == 4){
      con =
        dynamic_cast<CoupledThermoSolidStiffness*>(element->getConstitutive());
      order = cbrt(numNodes);
    }
    else {
      con =
        dynamic_cast<CoupledThermoQuadStiffness*>(element->getConstitutive());
      order = sqrt(numNodes);
    }
    if (con){
      // With the first iteration, find the maximum over the domain
      for ( int i = 0; i < numGauss; i++ ){
        // Get the Gauss points one at a time
        double pt[3];
        double weight = element->getGaussWtsPts(i, pt);
        element->getShapeFunctions(pt, N);
        // Get the strain B*u and temperature dT
        // If 3D structure
        if (numDisps == 4){
          
            CoupledThermoSolid<2>* elem =
              dynamic_cast<CoupledThermoSolid<2>*>(element);
          }
          else {

          }
          
        elem->getBT(strain, pt, Xpts, vars);
        

      }
    } // end if constitutive
    
  }
  
}
/*
  For each thread used to evaluate the function, call the
  post-evaluation code once.
*/
void TMRHeatFluxIntegral::finalThread( const double tcoef,
                                       EvaluationType ftype,
                                       TACSFunctionCtx *fctx ){
  HeatFluxIntCtx *ctx = dynamic_cast<HeatFluxIntCtx*>(fctx);
  if (ctx){
    value += ctx->value;
  }
}
