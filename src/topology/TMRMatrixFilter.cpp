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

#include "TMRMatrixFilter.h"

/*
  Create the filter matrix
*/
TMRMatrixFilter::TMRMatrixFilter( double _s, int _N, 
                                  int _nlevels,
                                  TACSAssembler *_tacs[],
                                  TMROctForest *_filter[],
                                  int _vars_per_node ):
  TMRConformFilter(_nlevels, _tacs, _filter, _vars_per_node){
  initialize_matrix(_s, _N, _filter[0], NULL);
}

TMRMatrixFilter::TMRMatrixFilter( double _s, int _N, 
                                  int _nlevels,
                                  TACSAssembler *_tacs[],
                                  TMRQuadForest *_filter[],
                                  int _vars_per_node ):
  TMRConformFilter(_nlevels, _tacs, _filter, _vars_per_node){
  initialize_matrix(_s, _N, NULL, _filter[0]);
}

TMRMatrixFilter::initialize_matrix( double _s, int _N,
                                    TMROctForest *oct_forest,
                                    TMRQuadForest *quad_forest ){
  // Store the pointer to the matrix
  M = _M;
  M->incref();

  // Allocate the vectors needed for the application of the filter
  Dinv = M->createVec();
  Tinv = M->createVec();
  t1 = M->createVec();
  t2 = M->createVec();
  Dinv->incref();
  Tinv->incref();
  t1->incref();
  t2->incref();
  
  // Set the number of terms in the M filter
  N = _N;

  // Set the scalar value associated with the filter
  s = _s;
  
  // Set all vector values to 1.0
  TACSBVec *_t2 = dynamic_cast<TACSBVec*>(t1);
  if (_t2){
    _t2->set(1.0);
  }

  // Compute D_{i} = sum_{j=1}^{n} M_{ij}
  M->mult(t2, Dinv);
  Dinv->scale(1.0/(s - 1.0));  

  // Create the inverse of the diagonal matrix
  TACSBVec *Dvec = dynamic_cast<TACSBVec*>(Dinv);
  if (Dvec){
    TacsScalar *D;
    int size = Dvec->getArray(&D);
    for ( int i = 0; i < size; i++ ){
      if (D[0] != 0.0){
        D[0] = 1.0/D[0];
      }
      else {
        D[0] = 0.0;
      }
    }
  }

  // Apply the filter to create the normalization
  applyFilter(t2, Tinv);
  
  // Create the inverse of the T matrix
  TACSBVec *Tvec = dynamic_cast<TACSBVec*>(Tinv);
  if (Tvec){
    TacsScalar *T;
    int size = Tvec->getArray(&T);
    for ( int i = 0; i < size; i++ ){
      if (T[0] != 0.0){
        T[0] = 1.0/T[0];
      }
      else {
        T[0] = 0.0;
      }
    }
  }
}

/*
  Destroy the filter matrix
*/
TMRMatrixFilter::~TMRMatrixFilter(){
  t1->decref();
  t2->decref();
  M->decref();
  Dinv->decref();
  Tinv->decref();
}

/*
  Compute the action of the filter on the input vector using Horner's
  method

  out = 1/s*in

  for n in range(N):
  .   out += (in + D^{-1}*M*out)/s

*/
void TMRMatrixFilter::applyFilter( TACSVec *in, TACSVec *out ){
  // Set out = in/s
  out->copyValues(in);
  out->scale(1.0/s);

  // Apply Horner's method
  for ( int n = 0; n < N; n++ ){
    // Compute t1 = D^{-1}*M*out
    M->mult(out, t1);
    kronecker(Dinv, t1);

    // Compute out = (in + D^{-1}*M*out)/s
    t1->axpy(1.0, in);
    out->axpy(1.0/s, t1);
  }

  // Multiply by Tinv
  kronecker(Tinv, out);
}

void TMRMatrixFilter::applyTranspose( TACSVec *in, TACSVec *out ){
  // Set out = in/s
  out->copyValues(in);
  out->scale(1.0/s);

  // Apply Horner's method
  for ( int n = 0; n < N; n++ ){
    // Compute M*D^{-1}*out
    kronecker(Dinv, out, t1);
    M->mult(t1, t2);

    // Compute out = (in + M*D^{-1}*out)/s
    t2->axpy(1.0, in);
    out->axpy(1.0/s, t2);
  }

  // Multiply by Tinv
  kronecker(Tinv, out);
}

/*
  Compute the Kronecker product of the input vectors c and x and store
  the result in either y (if it is provided) or otherwise x

  y = (c o x)
*/
void TMRMatrixFilter::kronecker( TACSVec *tc, TACSVec *tx, TACSVec *ty ){
  if (tc && tx && ty){
    TACSBVec *c = dynamic_cast<TACSBVec*>(tc);
    TACSBVec *x = dynamic_cast<TACSBVec*>(tx);
    TACSBVec *y = dynamic_cast<TACSBVec*>(ty);
    
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
  }
  else if (tc && tx){
    TACSBVec *c = dynamic_cast<TACSBVec*>(tc);
    TACSBVec *x = dynamic_cast<TACSBVec*>(tx);

    if (c && x){
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
}
