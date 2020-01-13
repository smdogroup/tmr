/*
  This file is part of the package TMR for adaptive mesh refinement.

  Copyright (C) 2015 Georgia Tech Research Corporation.
  Additional copyright (C) 2015 Graeme Kennedy.
  All rights reserved.

  TMR is licensed under the Apache License, Version 2.0 (the "License");
  you may not use this software except in compliance with the License.
  You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
*/

#include "TMRMatrixFilterModel.h"
#include "TMRMatrixFilter.h"
#include "TMRMatrixCreator.h"
#include "TMR_TACSCreator.h"

/*
  Create the filter matrix
*/
TMRMatrixFilter::TMRMatrixFilter( double _r, int _N,
                                  int _nlevels,
                                  TACSAssembler *_assembler[],
                                  TMROctForest *_filter[] ):
  TMRConformFilter(_nlevels, _assembler, _filter){
  initialize_matrix(_r, _N, _filter[0], NULL);
}

TMRMatrixFilter::TMRMatrixFilter( double _r, int _N,
                                  int _nlevels,
                                  TACSAssembler *_assembler[],
                                  TMRQuadForest *_filter[] ):
  TMRConformFilter(_nlevels, _assembler, _filter){
  initialize_matrix(_r, _N, NULL, _filter[0]);
}

/*
  Initialize the matrix filter.

  This code creates a TACSAssembler object (and frees it), assembles a
  mass matrix, creates the internal variables required for the filter.
*/
void TMRMatrixFilter::initialize_matrix( double _r, int _N,
                                         TMROctForest *oct_forest,
                                         TMRQuadForest *quad_forest ){
  // Create the Assembler object
  TACSAssembler *matrix_assembler = NULL;

  // Keep track of the dimension of the problem
  int d = 2;
  if (oct_filter){
    d = 3;
    TACSElementModel *model = new TMRHexaMatrixModel();
    TMROctTACSMatrixCreator *matrix_creator3d =
      new TMROctTACSMatrixCreator(model);
    matrix_creator3d->incref();

    matrix_assembler = matrix_creator3d->createTACS(oct_forest,
                                                    TACSAssembler::NATURAL_ORDER);
    matrix_assembler->incref();
    matrix_creator3d->decref();
  }
  else {
    d = 2;
    TACSElementModel *model = new TMRQuadMatrixModel();
    TMRQuadTACSMatrixCreator *matrix_creator2d =
      new TMRQuadTACSMatrixCreator(model);
    matrix_creator2d->incref();

    matrix_assembler = matrix_creator2d->createTACS(quad_forest,
                                                    TACSAssembler::NATURAL_ORDER);
    matrix_assembler->incref();
    matrix_creator2d->decref();
  }

  // Create the matrix
  M = matrix_assembler->createMat();
  M->incref();

  // Allocate the vectors needed for the application of the filter
  Ainv = matrix_assembler->createVec();
  B = matrix_assembler->createVec();
  Tinv = matrix_assembler->createVec();
  t1 = matrix_assembler->createVec();
  t2 = matrix_assembler->createVec();
  y1 = matrix_assembler->createVec();
  y2 = matrix_assembler->createVec();
  Ainv->incref();
  B->incref();
  Tinv->incref();
  t1->incref();
  t2->incref();
  y1->incref();
  y2->incref();

  // Create a temporary design vector
  temp = assembler[0]->createDesignVec();
  temp->incref();

  // Assemble the mass matrix
  matrix_assembler->assembleJacobian(1.0, 0.0, 0.0, t2, M);

  // Free this version of TACS - it's not required anymore!
  matrix_assembler->decref();

  // Set the number of terms in the M filter
  N = _N;

  // Set the scalar value associated with the filter
  r = _r;

  // Set all vector values to 1.0
  y2->set(1.0);

  // Compute D_{i} = sum_{j=1}^{n} M_{ij} (stored in B temporarily)
  M->mult(y2, B);

  // Create the inverse of the diagonal matrix
  TacsScalar *Bi; // The components of the D matrix that will be over-written
  TacsScalar *Ai;
  int size = B->getArray(&Bi);
  Ainv->getArray(&Ai);
  for ( int i = 0; i < size; i++ ){
    TacsScalar D = Bi[0];
    // Check if we have a 2D or 3D problem
    if (D == 0.0){
      Ai[0] = 1.0;
      Bi[0] = 0.0;
    }
    else {
      if (d == 2){
        Ai[0] = 1.0/(1.0 + r*r/D);
        Bi[0] = Ai[0]*r*r/(D*D);
      }
      else {
        TacsScalar temp = pow(D, 2.0/3.0);
        Ai[0] = 1.0/(1.0 + r*r/temp);
        Bi[0] = Ai[0]*r*r/(temp*D);
      }
    }
    Bi++;
    Ai++;
  }

  // Apply the filter to create the normalization
  Tinv->set(1.0);
  applyFilter(y2, y1);

  // Create the inverse of the T matrix
  TacsScalar *T, *ty;
  size = Tinv->getArray(&T);
  y1->getArray(&ty);
  for ( int i = 0; i < size; i++ ){
    if (ty[0] != 0.0){
      T[0] = 1.0/ty[0];
    }
    else {
      T[0] = 0.0;
    }
    T++;
    ty++;
  }
}

/*
  Destroy the filter matrix
*/
TMRMatrixFilter::~TMRMatrixFilter(){
  t1->decref();
  t2->decref();
  M->decref();
  Ainv->decref();
  B->decref();
  y1->decref();
  y2->decref();
  temp->decref();
}

/*
  Compute the action of the filter on the input vector using Horner's
  method

  t1 = Ainv*in
  out = t1
  for n in range(N):
  .   out += t1 + B*M*out
*/
void TMRMatrixFilter::applyFilter( TACSBVec *in, TACSBVec *out ){
  // Compute t1 = Ainv*in
  t1->copyValues(in);
  kronecker(Ainv, t1);

  // Set out = Ainv*in
  out->copyValues(t1);

  // Apply Horner's method
  for ( int n = 0; n < N; n++ ){
    // Compute t2 = B*M*out
    M->mult(out, t2);
    kronecker(B, t2, out);

    // Compute out = t1 + B*M*out
    out->axpy(1.0, t1);
  }

  // Multiply by Tinv
  kronecker(Tinv, out);
}

/*
  Compute the transpose of the filter operation
*/
void TMRMatrixFilter::applyTranspose( TACSBVec *in, TACSBVec *out ){
  t1->copyValues(in);
  kronecker(Tinv, t1);

  // Copy the values from t1 to the out vector
  out->copyValues(t1);

  // Apply Horner's method
  for ( int n = 0; n < N; n++ ){
    // Compute M*B*out
    kronecker(B, out, t2);
    M->mult(t2, out);

    // Compute out = t1 + M*B*out
    out->axpy(1.0, t1);
  }

  // Multiply by Ainv
  kronecker(Ainv, out);
}

/*
  Compute the Kronecker product of the input vectors c and x and store
  the result in either y (if it is provided) or otherwise x

  y = (c o x)
*/
void TMRMatrixFilter::kronecker( TACSBVec *c, TACSBVec *x, TACSBVec *y ){
  if (c && x && y){
    TacsScalar *cvals, *xvals, *yvals;
    int size = c->getArray(&cvals);
    x->getArray(&xvals);
    y->getArray(&yvals);

    for ( int i = 0; i < size; i++ ){
      yvals[0] = cvals[0]*xvals[0];
      yvals++;
      xvals++;
      cvals++;
    }
  }
  else if (c && x){
    TacsScalar *cvals, *xvals;
    int size = c->getArray(&cvals);
    x->getArray(&xvals);

    for ( int i = 0; i < size; i++ ){
      xvals[0] *= cvals[0];
      xvals++;
      cvals++;
    }
  }
}

/*
  Set the design variables for each level
*/
void TMRMatrixFilter::setDesignVars( TACSBVec *xvec ){
  const int vpn = assembler[0]->getDesignVarsPerNode();

  if (vpn == 1){
    applyFilter(xvec, x[0]);
  }
  else {
    for ( int k = 0; k < vpn; k++ ){
      TacsScalar *xin, *xout;
      xvec->getArray(&xin);
      x[0]->getArray(&xout);

      // Set the pointers offset by the vars per node
      xin = &xin[k];
      xout = &xout[k];

      // Copy the values to the input vector (y1)
      TacsScalar *yin;
      int size = y1->getArray(&yin);
      for ( int i = 0; i < size; i++ ){
        yin[0] = xin[0];
        yin++;
        xin += vpn;
      }

      // Apply the matrix filter from the input y1 to the
      // output y2
      applyFilter(y1, y2);

      // Copy the values from y2 to the output vector x[0]
      TacsScalar *yout;
      y2->getArray(&yout);
      for ( int i = 0; i < size; i++ ){
        xout[0] = yout[0];
        yout++;
        xout += vpn;
      }
    }
  }

  // Temporarily allocate an array to store the varia
  assembler[0]->setDesignVars(x[0]);

  // Set the design variable values on all processors
  for ( int k = 0; k < nlevels-1; k++ ){
    filter_interp[k]->multWeightTranspose(x[k], x[k+1]);

    // Set the design variable values
    assembler[k+1]->setDesignVars(x[k+1]);
  }
}

/*
  Add values to the output TACSBVec
*/
void TMRMatrixFilter::addValues( TACSBVec *vec ){
  vec->beginSetValues(TACS_ADD_VALUES);
  vec->endSetValues(TACS_ADD_VALUES);

  temp->copyValues(vec);

  const int vpn = assembler[0]->getDesignVarsPerNode();
  if (vpn == 1){
    applyTranspose(temp, vec);
  }
  else {
    for ( int k = 0; k < vpn; k++ ){
      TacsScalar *xin, *xout;

      // Get the pointer to the array and offset by vars per node
      temp->getArray(&xin);
      xin = &xin[k];

      // Copy the values to the input vector (y1)
      TacsScalar *yin;
      int size = y1->getArray(&yin);
      for ( int i = 0; i < size; i++ ){
        yin[0] = xin[0];
        yin++;
        xin += vpn;
      }

      // Apply the matrix filter from the input y1 to the
      // output y2
      applyTranspose(y1, y2);

      // Add the values from y2 back to the output vector vec
      vec->getArray(&xout);
      xout = &xout[k];

      TacsScalar *yout;
      y2->getArray(&yout);
      for ( int i = 0; i < size; i++ ){
        xout[0] = yout[0];
        yout++;
        xout += vpn;
      }
    }
  }
}
