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
                                          TMRQuadForest *_qforest,
                                          int _surface ):
TACSFunction(_tacs, TACSFunction::ENTIRE_DOMAIN,
             TACSFunction::SINGLE_STAGE, 0){
  maxNumNodes = _tacs->getMaxElementNodes();
  name = _name;
  value = 0.0;
  oforest = _oforest;
  qforest = _qforest;
  surface = _surface;
  
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
    double N[numNodes], Na[numNodes], Nb[numNodes];
    TacsScalar q[numDisps-1];
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
    // Direction vector of the surface/edge
    TacsScalar dir[3];
    if (con){
      // With the first iteration, find the maximum over the domain
      for ( int i = 0; i < numGauss; i++ ){
        // Get the Gauss points one at a time
        double pt[3];
        double weight = element->getGaussWtsPts(i, pt);
        element->getShapeFunctions(pt, N, Na, Nb);

        TacsScalar Xa[3], Xb[3];
        Xa[0] = Xa[1] = Xa[2] = 0.0;
        Xb[0] = Xb[1] = Xb[2] = 0.0;

        dir[0] = dir[1] = dir[2] = 0.0;
        
        // Get the strain B*u and temperature dT
        // If 3D structure
        if (numDisps == 4){
          if (surface < 2){
            const int ii = (order-1)*(surface % 2);
            for ( int kk = 0; kk < order; kk++ ){
              for ( int jj = 0; jj < order; jj++ ){
                const int node = ii + jj*order + kk*order*order;
              
                Xa[0] += Na[jj + kk*order]*Xpts[3*node];
                Xa[1] += Na[jj + kk*order]*Xpts[3*node+1];
                Xa[2] += Na[jj + kk*order]*Xpts[3*node+2];

                Xb[0] += Nb[jj + kk*order]*Xpts[3*node];
                Xb[1] += Nb[jj + kk*order]*Xpts[3*node+1];
                Xb[2] += Nb[jj + kk*order]*Xpts[3*node+2];
              }
            }
          }
          else if (surface < 4){
            const int jj = (order-1)*(surface % 2);
            for ( int kk = 0; kk < order; kk++ ){
              for ( int ii = 0; ii < order; ii++ ){
                const int node = ii + jj*order + kk*order*order;
              
                Xa[0] += Na[ii + kk*order]*Xpts[3*node];
                Xa[1] += Na[ii + kk*order]*Xpts[3*node+1];
                Xa[2] += Na[ii + kk*order]*Xpts[3*node+2];

                Xb[0] += Nb[ii + kk*order]*Xpts[3*node];
                Xb[1] += Nb[ii + kk*order]*Xpts[3*node+1];
                Xb[2] += Nb[ii + kk*order]*Xpts[3*node+2];
              }
            }
          }
          else {
            const int kk = (order-1)*(surface % 2);
            for ( int jj = 0; jj < order; jj++ ){
              for ( int ii = 0; ii < order; ii++ ){
                const int node = ii + jj*order + kk*order*order;
              
                Xa[0] += Na[ii + jj*order]*Xpts[3*node];
                Xa[1] += Na[ii + jj*order]*Xpts[3*node+1];
                Xa[2] += Na[ii + jj*order]*Xpts[3*node+2];

                Xb[0] += Nb[ii + jj*order]*Xpts[3*node];
                Xb[1] += Nb[ii + jj*order]*Xpts[3*node+1];
                Xb[2] += Nb[ii + jj*order]*Xpts[3*node+2];
              }
            }
          }
          // Compute the normal to the element
          Tensor::crossProduct3D(dir, Xa, Xb);

          if (order == 2){
            CoupledThermoSolid<2>* elem =
              dynamic_cast<CoupledThermoSolid<2>*>(element);
            if (elem){
              elem->getBT(strain, pt, Xpts, vars);
            }
          }
          else if (order == 3){
            CoupledThermoSolid<3>* elem =
              dynamic_cast<CoupledThermoSolid<3>*>(element);
            if (elem){
              elem->getBT(strain, pt, Xpts, vars);
            }
          }
          else if (order == 4){
            CoupledThermoSolid<4>* elem =
              dynamic_cast<CoupledThermoSolid<4>*>(element);
            if (elem){
              elem->getBT(strain, pt, Xpts, vars);
            }
          }
          else if (order == 5){
            CoupledThermoSolid<5>* elem =
              dynamic_cast<CoupledThermoSolid<5>*>(element);
            if (elem){
              elem->getBT(strain, pt, Xpts, vars);
            }
          }
          else if (order == 6){
            CoupledThermoSolid<6>* elem =
              dynamic_cast<CoupledThermoSolid<6>*>(element);
            if (elem){
              elem->getBT(strain, pt, Xpts, vars);
            }
          }
          else {
            printf("Heat flux element not implemented\n");
          }
        }
        else {
          // Compute direction
          if (surface == 0 || surface == 1){
            dir[0] = 0.0;
            dir[1] = 1.0;
          }
          else {
            dir[0] = 1.0;
            dir[1] = 0.0;
          }         

          if (order == 2){
            PlaneStressCoupledThermoQuad<2>* elem =
              dynamic_cast<PlaneStressCoupledThermoQuad<2>*>(element);
            if (elem){
              elem->getBT(strain, pt, Xpts, vars);
            }
          }
          else if (order == 3){
            PlaneStressCoupledThermoQuad<3>* elem =
              dynamic_cast<PlaneStressCoupledThermoQuad<3>*>(element);
            if (elem){
              elem->getBT(strain, pt, Xpts, vars);
            }
          }
          else if (order == 4){
            PlaneStressCoupledThermoQuad<4>* elem =
              dynamic_cast<PlaneStressCoupledThermoQuad<4>*>(element);
            if (elem){
              elem->getBT(strain, pt, Xpts, vars);
            }
          }
          else if (order == 5){
            PlaneStressCoupledThermoQuad<5>* elem =
              dynamic_cast<PlaneStressCoupledThermoQuad<5>*>(element);
            if (elem){
              elem->getBT(strain, pt, Xpts, vars);
            }
          }
          else if (order == 6){
            PlaneStressCoupledThermoQuad<6>* elem =
              dynamic_cast<PlaneStressCoupledThermoQuad<6>*>(element);
            if (elem){
              elem->getBT(strain, pt, Xpts, vars);
            }
          }
          else {
            printf("Heat flux element not implemented\n");
          }
        }          
        
        con->calculateConduction(pt, strain, q);
        if (numDisps == 3){
          value += dir[0]*q[0] + dir[1]*q[1];
        }
        else {
          value += dir[0]*q[0] + dir[1]*q[1] + dir[2]*q[2];
        }
      }
      // Add up the contribution from the quadrature
      TacsScalar h = element->getDetJacobian(pt, Xpts);
      h *= weight;

      ctx->value += h*value;
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
