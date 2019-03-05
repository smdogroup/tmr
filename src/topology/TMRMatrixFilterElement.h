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

#ifndef TMR_MATRIX_FILTER_ELEMENT_H
#define TMR_MATRIX_FILTER_ELEMENT_H

#include "TMRInterpolation.h"
#include "TACSElement.h"
#include "FElibrary.h"

template <int order>
class TMRQuadMatrixElement : public TACSElement {
 public:
  static const int NUM_NODES = order*order;
  TMRQuadMatrixElement(){}
  ~TMRQuadMatrixElement(){}

  const char* elementName(){
    return "TMRQuadMatrixElement";
  }
  const char* displacementName( int i ){
    return "u";
  }
  const char* stressName( int i ){
    return NULL;
  }
  const char* strainName( int i ){
    return NULL;
  }
  const char* extraName( int i ){
    return NULL;
  }
  int numDisplacements(){ return 1; }
  int numNodes(){ return NUM_NODES; }
  int numStresses(){ return 0; }
  int numExtras(){ return 0; }
  ElementType getElementType(){
    return TACS_POISSON_2D_ELEMENT;
  }
  int getNumGaussPts(){
    return NUM_NODES;
  }
  double getGaussWtsPts( const int num, double *pt ){
    int n = num % order;
    int m = num/order;

    const double *pts, *wts;
    FElibrary::getGaussPtsWts(order, &pts, &wts);
    pt[0] = pts[n];
    pt[1] = pts[m];
    return wts[n]*wts[m];
  }
  TacsScalar getDetJacobian( const double *pt,
                             const TacsScalar Xpts[] ){
    // Get the shape functions
    double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
    getShapeFunctions(pt, N, Na, Nb);

    // Compute the Jacobian transformation
    TacsScalar Xd[4];
    getJacobianTransform(Na, Nb, Xpts, Xd);

    return Xd[0]*Xd[3] - Xd[1]*Xd[2];
  }
  void getShapeFunctions( const double pt[], double N[] ){
    double na[order], nb[order];
    bernstein_shape_functions(order, pt[0], na);
    bernstein_shape_functions(order, pt[1], nb);
    for ( int j = 0; j < order; j++ ){
      for ( int i = 0; i < order; i++ ){
        N[0] = na[i]*nb[j];
        N++;
      }
    }
  }
  void getShapeFunctions( const double pt[], double N[],
                          double Na[], double Nb[] ){
    double na[order], nb[order];
    double dna[order], dnb[order];
    bernstein_shape_func_derivative(order, pt[0], na, dna);
    bernstein_shape_func_derivative(order, pt[1], nb, dnb);
    for ( int j = 0; j < order; j++ ){
      for ( int i = 0; i < order; i++ ){
        N[0] = na[i]*nb[j];
        Na[0] = dna[i]*nb[j];
        Nb[0] = na[i]*dnb[j];
        N++;  Na++;  Nb++;
      }
    }
  }
  void getJacobianTransform( const double Na[], const double Nb[],
                             const TacsScalar Xpts[],
                             TacsScalar Xd[] ){
    Xd[0] = Xd[1] = Xd[2] = Xd[3] = 0.0;
    for ( int i = 0; i < NUM_NODES; i++ ){
      Xd[0] += Na[i]*Xpts[3*i];
      Xd[1] += Nb[i]*Xpts[3*i];
      Xd[2] += Na[i]*Xpts[3*i+1];
      Xd[3] += Nb[i]*Xpts[3*i+1];
    }
  }
  void addResidual( double time, TacsScalar res[],
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] ){}
  void addJacobian( double time, TacsScalar mat[],
                    double alpha, double beta, double gamma,
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] ){
    const double *pts, *wts;
    FElibrary::getGaussPtsWts(order, &pts, &wts);

    for ( int m = 0; m < order; m++ ){
      for ( int n = 0; n < order; n++ ){
        // Set the quadrature points
        double pt[2];
        pt[0] = pts[n];
        pt[1] = pts[m];

        // Get the shape functions
        double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
        getShapeFunctions(pt, N, Na, Nb);

        // Compute the Jacobian transformation
        TacsScalar Xd[4];
        getJacobianTransform(Na, Nb, Xpts, Xd);

        // Compute the inverse of the Jacobian
        TacsScalar J[4];
        TacsScalar h = FElibrary::jacobian2d(Xd, J);
        h *= alpha*wts[n]*wts[m];

        for ( int j = 0; j < NUM_NODES; j++ ){
          for ( int i = 0; i < NUM_NODES; i++ ){
            // Add to the matrix N^{T}*N
            mat[i + j*NUM_NODES] += h*N[i]*N[j];
          }
        }
      }
    }
  }
};

template <int order>
class TMROctMatrixElement : public TACSElement {
 public:
  static const int NUM_NODES = order*order*order;

  TMROctMatrixElement(){}
  ~TMROctMatrixElement(){}

  const char* elementName(){
    return "TMROctMatrixElement";
  }
  const char* displacementName( int i ){
    return "u";
  }
  const char* stressName( int i ){
    return NULL;
  }
  const char* strainName( int i ){
    return NULL;
  }
  const char* extraName( int i ){
    return NULL;
  }
  int numDisplacements(){ return 1; }
  int numNodes(){ return NUM_NODES; }
  int numStresses(){ return 0; }
  int numExtras(){ return 0; }
  ElementType getElementType(){
    return TACS_POISSON_3D_ELEMENT;
  }
  int getNumGaussPts(){
    return NUM_NODES;
  }
  double getGaussWtsPts( const int num, double *pt ){
    int p = (int)((num)/(order*order));
    int m = (int)((num - order*order*p)/order);
    int n = num - order*m - order*order*p;

    const double *pts, *wts;
    FElibrary::getGaussPtsWts(order, &pts, &wts);
    pt[0] = pts[n];
    pt[1] = pts[m];
    pt[2] = pts[p];

    return wts[n]*wts[m]*wts[p];
  }
  TacsScalar getDetJacobian( const double *pt,
                             const TacsScalar Xpts[] ){
    // Get the shape functions
    double N[NUM_NODES];
    double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
    getShapeFunctions(pt, N, Na, Nb, Nc);

    // Compute the Jacobian transformation
    TacsScalar Xd[9];
    getJacobianTransform(Na, Nb, Nc, Xpts, Xd);

    return FElibrary::jacobian3d(Xd);
  }
  void getShapeFunctions( const double pt[], double N[] ){
    double na[order], nb[order], nc[order];
    bernstein_shape_functions(order, pt[0], na);
    bernstein_shape_functions(order, pt[1], nb);
    bernstein_shape_functions(order, pt[2], nc);

    for ( int k = 0; k < order; k++ ){
      for ( int j = 0; j < order; j++ ){
        for ( int i = 0; i < order; i++ ){
          N[0] = na[i]*nb[j]*nc[k];
          N++;
        }
      }
    }
  }
  void getShapeFunctions( const double pt[], double N[],
                          double Na[], double Nb[], double Nc[] ){
    double na[order], nb[order], nc[order];
    double dna[order], dnb[order], dnc[order];
    bernstein_shape_func_derivative(order, pt[0], na, dna);
    bernstein_shape_func_derivative(order, pt[1], nb, dnb);
    bernstein_shape_func_derivative(order, pt[2], nc, dnc);
    for ( int k = 0; k < order; k++ ){
      for ( int j = 0; j < order; j++ ){
        for ( int i = 0; i < order; i++ ){
          N[0] = na[i]*nb[j]*nc[k];
          Na[0] = dna[i]*nb[j]*nc[k];
          Nb[0] = na[i]*dnb[j]*nc[k];
          Nc[0] = na[i]*nb[j]*dnc[k];
          N++;  Na++;  Nb++;  Nc++;
        }
      }
    }
  }
  void getJacobianTransform( const double Na[],
                             const double Nb[],
                             const double Nc[],
                             const TacsScalar Xpts[],
                             TacsScalar Xd[] ){
    Xd[0] = Xd[1] = Xd[2] = 0.0;
    Xd[3] = Xd[4] = Xd[5] = 0.0;
    Xd[6] = Xd[7] = Xd[8] = 0.0;

    for ( int i = 0; i < NUM_NODES; i++ ){
      Xd[0] += Xpts[0]*Na[0];
      Xd[1] += Xpts[0]*Nb[0];
      Xd[2] += Xpts[0]*Nc[0];

      Xd[3] += Xpts[1]*Na[0];
      Xd[4] += Xpts[1]*Nb[0];
      Xd[5] += Xpts[1]*Nc[0];

      Xd[6] += Xpts[2]*Na[0];
      Xd[7] += Xpts[2]*Nb[0];
      Xd[8] += Xpts[2]*Nc[0];

      Na++; Nb++; Nc++;
      Xpts += 3;
    }
  }
  void addResidual( double time, TacsScalar res[],
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] ){}
  void addJacobian( double time, TacsScalar mat[],
                    double alpha, double beta, double gamma,
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] ){
    const double *pts, *wts;
    FElibrary::getGaussPtsWts(order, &pts, &wts);

    for ( int p = 0; p < order; p++ ){
      for ( int m = 0; m < order; m++ ){
        for ( int n = 0; n < order; n++ ){
          // Set the quadrature points
          double pt[3];
          pt[0] = pts[n];
          pt[1] = pts[m];
          pt[2] = pts[p];

          // Get the shape functions
          double N[NUM_NODES];
          double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
          getShapeFunctions(pt, N, Na, Nb, Nc);

          // Compute the Jacobian transformation
          TacsScalar Xd[9];
          getJacobianTransform(Na, Nb, Nc, Xpts, Xd);

          // Compute the inverse of the Jacobian
          TacsScalar J[9];
          TacsScalar h = FElibrary::jacobian3d(Xd, J);
          h *= alpha*wts[n]*wts[m]*wts[p];

          for ( int j = 0; j < NUM_NODES; j++ ){
            for ( int i = 0; i < NUM_NODES; i++ ){
              // Add to the matrix N^{T}*N
              mat[i + j*NUM_NODES] += h*N[i]*N[j];
            }
          }
        }
      }
    }
  }
};

#endif // TMR_MATRIX_FILETER_ELEMENT_H
