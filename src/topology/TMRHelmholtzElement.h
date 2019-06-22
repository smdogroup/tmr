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

#ifndef TMR_HELMHOLTZ_ELEMENT_H
#define TMR_HELMHOLTZ_ELEMENT_H

/*
  The following file implements 2D and 3D Helmholtz filter elements
  for a spatial filter. Note that this will not neccessarily be a
  discrete partition-of-unity filter.
*/

#include "TMRInterpolation.h"
#include "TACSElement.h"
#include "FElibrary.h"

template <int order>
class TMRQuadHelmholtz : public TACSElement {
 public:
  static const int NUM_NODES = order*order;

  TMRQuadHelmholtz( TacsScalar r=1.0 ){
    r2 = r*r;
  }
  ~TMRQuadHelmholtz(){}

  const char* elementName(){
    return "TMRQuadHelmholtz";
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
  void getHelmholtzRes( const TacsScalar Xpts[],
                        const TacsScalar xdvs[],
                        TacsScalar res[] ){
    memset(res, 0, NUM_NODES*sizeof(TacsScalar));

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
        h *= wts[n]*wts[m];

        TacsScalar q = 0.0;
        for ( int i = 0; i < NUM_NODES; i++ ){
          q += N[i]*xdvs[i];
        }

        // Add the term to the residual
        for ( int i = 0; i < NUM_NODES; i++ ){
          res[i] += h*N[i]*q;
        }
      }
    }
  }
  void addResidual( double time, TacsScalar res[],
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
        h *= wts[n]*wts[m];

        // Compute the derivatives of fval
        TacsScalar px = 0.0, py = 0.0, q = 0.0;
        for ( int i = 0; i < NUM_NODES; i++ ){
          px += (Na[i]*J[0] + Nb[i]*J[2])*vars[i];
          py += (Na[i]*J[1] + Nb[i]*J[3])*vars[i];
          q += N[i]*vars[i];
        }

        // Add the term to the residual
        for ( int i = 0; i < NUM_NODES; i++ ){
          TacsScalar Nxi = Na[i]*J[0] + Nb[i]*J[2];
          TacsScalar Nyi = Na[i]*J[1] + Nb[i]*J[3];
          res[i] += h*(r2*(Nxi*px + Nyi*py) + N[i]*q);
        }
      }
    }
  }
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
          TacsScalar Nxj = Na[j]*J[0] + Nb[j]*J[2];
          TacsScalar Nyj = Na[j]*J[1] + Nb[j]*J[3];

          for ( int i = 0; i < NUM_NODES; i++ ){
            TacsScalar Nxi = Na[i]*J[0] + Nb[i]*J[2];
            TacsScalar Nyi = Na[i]*J[1] + Nb[i]*J[3];
            // Add to the matrix r^2*Nu^{T}*Nv
            mat[i + j*NUM_NODES] += h*r2*(Nxi*Nxj + Nyi*Nyj);

            // Add to the matrix N^{T}*N
            mat[i + j*NUM_NODES] += h*N[i]*N[j];
          }
        }
      }
    }
  }
  void addOutputCount( int *nelems, int *nnodes, int *ncsr ){
    *nelems += (order-1)*(order-1);
    *nnodes += NUM_NODES;
    *ncsr += 4*(order-1)*(order-1);
  }
  void getOutputData( unsigned int out_type,
                      double *data, int ld_data,
                      const TacsScalar Xpts[],
                      const TacsScalar vars[] ){
    for ( int m = 0; m < order; m++ ){
      for ( int n = 0; n < order; n++ ){
        int p = n + order*m;
        int index = 0;
        // Set the parametric point to extract the data
        double pt[2];
        pt[0] = -1.0 + 2.0*n/(order-1);
        pt[1] = -1.0 + 2.0*m/(order-1);

        // Get the shape functions
        double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
        getShapeFunctions(pt, N, Na, Nb);

        if (out_type & TACSElement::OUTPUT_NODES){
          for ( int k = 0; k < 3; k++ ){
            TacsScalar X = 0.0;
            for ( int i = 0; i < NUM_NODES; i++ ){
              X += N[i]*Xpts[3*p+k];
            }
            data[index+k] = TacsRealPart(X);
          }
          index += 3;
        }
        if (out_type & TACSElement::OUTPUT_DISPLACEMENTS){
          TacsScalar u = 0.0;
          for ( int i = 0; i < NUM_NODES; i++ ){
            u += N[i]*vars[i];
          }
          data[index] = TacsRealPart(u);
          index++;
        }
        data += ld_data;
      }
    }
  }
  void getOutputConnectivity( int *con, int node ){
    int p = 0;
    for ( int m = 0; m < order-1; m++ ){
      for ( int n = 0; n < order-1; n++ ){
        con[4*p]   = node + n   + m*order;
        con[4*p+1] = node + n+1 + m*order;
        con[4*p+2] = node + n+1 + (m+1)*order;
        con[4*p+3] = node + n   + (m+1)*order;
        p++;
      }
    }
  }

 private:
  // Filter length scale argument
  TacsScalar r2;
};

template <int order>
class TMROctHelmholtz : public TACSElement {
 public:
  static const int NUM_NODES = order*order*order;

  TMROctHelmholtz( TacsScalar r=1.0 ){
    // Set the length scale
    r2 = r*r;
  }
  ~TMROctHelmholtz(){}

  const char* elementName(){
    return "TMROctHelmholtz";
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
  void getHelmholtzRes( const TacsScalar Xpts[],
                        const TacsScalar xdvs[],
                        TacsScalar res[] ){
    memset(res, 0, NUM_NODES*sizeof(TacsScalar));

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
          h *= wts[n]*wts[m]*wts[p];

          TacsScalar q = 0.0;
          for ( int i = 0; i < NUM_NODES; i++ ){
            q += N[i]*xdvs[i];
          }

          // Add the term to the residual
          for ( int i = 0; i < NUM_NODES; i++ ){
            res[i] += h*N[i]*q;
          }
        }
      }
    }
  }
  void addResidual( double time, TacsScalar res[],
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
          h *= wts[n]*wts[m]*wts[p];

          // Compute the derivatives of fval
          TacsScalar px = 0.0, py = 0.0, pz = 0.0, q = 0.0;
          for ( int i = 0; i < NUM_NODES; i++ ){
            px += (Na[i]*J[0] + Nb[i]*J[3] + Nc[i]*J[6])*vars[i];
            py += (Na[i]*J[1] + Nb[i]*J[4] + Nc[i]*J[7])*vars[i];
            pz += (Na[i]*J[2] + Nb[i]*J[5] + Nc[i]*J[8])*vars[i];
            q += N[i]*vars[i];
          }

          // Add the term to the residual
          for ( int i = 0; i < NUM_NODES; i++ ){
            TacsScalar Nxi = Na[i]*J[0] + Nb[i]*J[3] + Nc[i]*J[6];
            TacsScalar Nyi = Na[i]*J[1] + Nb[i]*J[4] + Nc[i]*J[7];
            TacsScalar Nzi = Na[i]*J[2] + Nb[i]*J[5] + Nc[i]*J[8];

            res[i] += h*(r2*(Nxi*px + Nyi*py + Nzi*pz) + N[i]*q);
          }
        }
      }
    }
  }
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
            TacsScalar Nxj = Na[j]*J[0] + Nb[j]*J[3] + Nc[j]*J[6];
            TacsScalar Nyj = Na[j]*J[1] + Nb[j]*J[4] + Nc[j]*J[7];
            TacsScalar Nzj = Na[j]*J[2] + Nb[j]*J[5] + Nc[j]*J[8];

            for ( int i = 0; i < NUM_NODES; i++ ){
              TacsScalar Nxi = Na[i]*J[0] + Nb[i]*J[3] + Nc[i]*J[6];
              TacsScalar Nyi = Na[i]*J[1] + Nb[i]*J[4] + Nc[i]*J[7];
              TacsScalar Nzi = Na[i]*J[2] + Nb[i]*J[5] + Nc[i]*J[8];

              // Add to the matrix r^2*Nu^{T}*Nv + N^{T}*N
              mat[i + j*NUM_NODES] +=
                h*(r2*(Nxi*Nxj + Nyi*Nyj + Nzi*Nzj) + N[i]*N[j]);
            }
          }
        }
      }
    }
  }
  void addOutputCount( int *nelems, int *nnodes, int *ncsr ){
    *nelems += (order-1)*(order-1)*(order-1);
    *nnodes += NUM_NODES;
    *ncsr += 8*(order-1)*(order-1)*(order-1);
  }
  void getOutputData( unsigned int out_type,
                      double *data, int ld_data,
                      const TacsScalar Xpts[],
                      const TacsScalar vars[] ){
    for ( int p = 0; p < order; p++ ){
      for ( int m = 0; m < order; m++ ){
        for ( int n = 0; n < order; n++ ){
          int index = 0;
          // Set the parametric point to extract the data
          double pt[3];
          pt[0] = -1.0 + 2.0*n/(order-1);
          pt[1] = -1.0 + 2.0*m/(order-1);
          pt[2] = -1.0 + 2.0*p/(order-1);

          // Get the shape functions
          double N[NUM_NODES];
          double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
          getShapeFunctions(pt, N, Na, Nb, Nc);

          if (out_type & TACSElement::OUTPUT_NODES){
            for ( int k = 0; k < 3; k++ ){
              TacsScalar X = 0.0;
              for ( int i = 0; i < NUM_NODES; i++ ){
                X += N[i]*Xpts[3*i + k];
              }
              data[index+k] = TacsRealPart(X);
            }
            index += 3;
          }
          if (out_type & TACSElement::OUTPUT_DISPLACEMENTS){
            TacsScalar u = 0.0;
            for (int i = 0; i < NUM_NODES; i++ ){
              u += N[i]*vars[i];
            }
            data[index] = TacsRealPart(u);
            index++;
          }
          data += ld_data;
        }
      }
    }
  }
  void getOutputConnectivity( int *con, int node ){
    int j = 0;
    for ( int p = 0; p < order-1; p++ ){
      for ( int m = 0; m < order-1; m++ ){
        for ( int n = 0; n < order-1; n++ ){
          con[8*j]   = node + n   +     m*order +     p*order*order;
          con[8*j+1] = node + n+1 +     m*order +     p*order*order;
          con[8*j+2] = node + n+1 + (m+1)*order +     p*order*order;
          con[8*j+3] = node + n   + (m+1)*order +     p*order*order;
          con[8*j+4] = node + n   +     m*order + (p+1)*order*order;
          con[8*j+5] = node + n+1 +     m*order + (p+1)*order*order;
          con[8*j+6] = node + n+1 + (m+1)*order + (p+1)*order*order;
          con[8*j+7] = node + n   + (m+1)*order + (p+1)*order*order;
          j++;
        }
      }
    }
  }

 private:
  // Filter length scale argument
  TacsScalar r2;
};

#endif // TMR_HELMHOLTZ_H
