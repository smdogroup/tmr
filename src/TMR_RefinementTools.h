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

#ifndef TMR_REFINEMENT_TOOLS_H
#define TMR_REFINEMENT_TOOLS_H

#include "TMRQuadForest.h"
#include "TMROctForest.h"
#include "TACSAssembler.h"
#include "TACSMg.h"

/*
  Create a TACS multigrid object
*/
void TMR_CreateTACSMg( int nlevels, TACSAssembler *tacs[],
                       TMROctForest *forest[], TACSMg **_mg,
                       double omega=1.0,
                       int use_coarse_direct_solve=1,
                       int use_chebyshev_smoother=0 );
void TMR_CreateTACSMg( int nlevels, TACSAssembler *tacs[],
                       TMRQuadForest *forest[], TACSMg **_mg,
                       double omega=1.0,
                       int use_coarse_direct_solve=1,
                       int use_chebyshev_smoother=0 );

/*
  Compute a direct interpolation from a lower-order mesh to a
  higher-order one
*/
void TMR_ComputeInterpSolution( TMRQuadForest *forest,
                                TACSAssembler *tacs,
                                TMRQuadForest *forest_refined,
                                TACSAssembler *tacs_refined,
                                TACSBVec *_uvec,
                                TACSBVec *_uvec_refined );
void TMR_ComputeInterpSolution( TMROctForest *forest,
                                TACSAssembler *tacs,
                                TMROctForest *forest_refined,
                                TACSAssembler *tacs_refined,
                                TACSBVec *_uvec,
                                TACSBVec *_uvec_refined );

/*
  Compute the reconstructed solution on a uniformly refined mesh
*/
void TMR_ComputeReconSolution( TMRQuadForest *forest,
                               TACSAssembler *tacs,
                               TMRQuadForest *forest_refined,
                               TACSAssembler *tacs_refined,
                               TACSBVec *_uvec=NULL,
                               TACSBVec *_uvec_refined=NULL,
                               const int compute_difference=0 );
void TMR_ComputeReconSolution( TMROctForest *forest,
                               TACSAssembler *tacs,
                               TMROctForest *forest_refined,
                               TACSAssembler *tacs_refined,
                               TACSBVec *_uvec=NULL,
                               TACSBVec *_uvec_refined=NULL,
                               const int compute_difference=0 );

/*
  Print the error as a series of bins
*/
void TMR_PrintErrorBins( MPI_Comm comm, const double *error,
                         const int nelems, double *mean=NULL,
                         double *stddev=NULL );

/*
  Perform a mesh refinement based on the strain engery refinement
  criteria.
*/
double TMR_StrainEnergyErrorEst( TMRQuadForest *forest,
                                 TACSAssembler *tacs,
                                 TMRQuadForest *forest_refined,
                                 TACSAssembler *tacs_refined,
                                 double *error );
double TMR_StrainEnergyErrorEst( TMROctForest *forest,
                                 TACSAssembler *tacs,
                                 TMROctForest *forest_refined,
                                 TACSAssembler *tacs_refined,
                                 double *error );

/*
  Perform adjoint-based mesh refinement on the forest of quadtrees
*/
double TMR_AdjointErrorEst( TMRQuadForest *forest,
                            TACSAssembler *tacs,
                            TMRQuadForest *forest_refined,
                            TACSAssembler *tacs_refined,
                            TACSBVec *solution_refined,
                            TACSBVec *adjoint_refined,
                            double *error,
                            double *adj_corr );
double TMR_AdjointErrorEst( TMROctForest *forest,
                            TACSAssembler *tacs,
                            TMROctForest *forest_refined,
                            TACSAssembler *tacs_refined,
                            TACSBVec *solution_refined,
                            TACSBVec *adjoint_refined,
                            double *error,
                            double *adj_corr );

/*
  Evaluate a stress constraint based on a higher-order interpolation
  of the stresses in the problem.

  This makes strong assumptions about the element type and constitutive
  matrix. Be careful when using this method
*/
class TMRStressConstraint : public TMREntity {
 public:
  TMRStressConstraint( TMROctForest *_forest,
                       TACSAssembler *tacs,
                       TacsScalar _ks_weight=30.0 );
  ~TMRStressConstraint();

  // Evaluate the aggregated stress constraint across all processors
  TacsScalar evalConstraint( TACSBVec *_uvec );

  // Evaluate the terms required for the total derivative
  void evalConDeriv( TacsScalar *dfdx, int size, TACSBVec *dfdu );

  // Write the von Misses stress from the reconstruction to tecplot
  void writeReconToTec( TACSBVec *_uvec, const char *fname,
                        TacsScalar ys=1e6 );

 private:
  // Evaluate the element strain at the given point
  TacsScalar evalStrain( const double pt[],
                         const TacsScalar *Xpts,
                         const TacsScalar *vars,
                         const TacsScalar *ubar,
                         TacsScalar J[], TacsScalar e[] );

  // Evaluate the derivatives of the element strain
  TacsScalar addStrainDeriv( const double pt[], const TacsScalar J[],
                             const TacsScalar alpha, const TacsScalar dfde[],
                             TacsScalar dfdu[], TacsScalar dfdubar[] );

  // Compute the derivative terms related to the reconstructed solution
  void addEnrichDeriv( TacsScalar A[], TacsScalar dbdu[],
                       TacsScalar dubardu[],
                       TacsScalar dubar_duderiv[] );

  // The mesh order
  int order;
  TMROctForest *forest;
  TMROctForest *interp_forest;

  // The values used to compute the KS function
  TacsScalar ks_weight;
  TacsScalar ks_max_fail, ks_fail_sum;

  // The TACSAssembler
  TACSAssembler *tacs;

  // The solution vector and the derivatives at the nodes
  TACSBVec *uvec, *uderiv;

  // The temporary derivative vector
  TACSBVec *dfduderiv;

  // The weights on the local
  TACSBVec *weights;

  // Arrays to store the element-wise variables
  TacsScalar *Xpts;
  TacsScalar *vars, *dvars, *ddvars;

  // The values of the element-wise reconstructed derivatives
  // at each node
  TacsScalar *varderiv;

  // Values of the variables for each enrichment function
  TacsScalar *ubar;

  // Temporary vector
  TacsScalar *tmp;
};

/*
  Evaluate a curvature constraint
*/
class TMRCurvatureConstraint : public TMREntity {
 public:
  TMRCurvatureConstraint( TMROctForest *_forest,
                          TacsScalar _aggregate_weight );
  TMRCurvatureConstraint( TMROctForest *_forest,
                          TACSVarMap *_varmap,
                          TacsScalar _aggregate_weight );
  ~TMRCurvatureConstraint();

  // Evaluate the curvature constraint
  TacsScalar evalConstraint( TACSBVec *_xvec );

  // Write the curvature to file
  void writeCurvatureToFile( TACSBVec *_xvec, const char *filename );

 private:
  void computeNodeDeriv();
  void estimateHessian( const TacsScalar elem_Xpts[], 
                        const TacsScalar elem_vals[],
                        const TacsScalar elem_derivs[], 
                        TacsScalar *val, TacsScalar g[], TacsScalar H[] );
  TacsScalar evalCurvature( const TacsScalar val,
                            const TacsScalar g[],
                            const TacsScalar H[] );
  TacsScalar evalCurvDeriv( const TacsScalar val, const TacsScalar g[],
                            const TacsScalar H[], TacsScalar *dval,
                            TacsScalar dg[], TacsScalar dH[] );
  TacsScalar evalCurvature( const int elem_size,
                            const double N[], const double Na[],
                            const double Nb[], const double Nc[],
                            const TacsScalar J[], const TacsScalar Xpts[],
                            const TacsScalar elem_vals[],
                            const TacsScalar elem_deriv[] );
  TacsScalar addCurvDeriv( const TacsScalar alpha, const int elem_size,
                           const double N[], const double Na[],
                           const double Nb[], const double Nc[],
                           const TacsScalar J[], const TacsScalar Xpts[],
                           const TacsScalar elem_vals[],
                           const TacsScalar elem_deriv[],
                           TacsScalar dvals[], TacsScalar dderiv[] );

  void init( TMROctForest *_forest, TACSVarMap *_varmap,
             TacsScalar _aggregate_weight );

  // The variable map
  TACSVarMap *varmap;

  // The TMROctForest
  TMROctForest *forest;

  // The vector of weights for the reconstruction
  TACSBVec *weights;
  
  // The variables and their derivatives
  TACSBVec *xvec;
  TACSBVec *xderiv, *dfderiv;

  // Values for the curvature aggregation
  TacsScalar aggregate_weight;
  TacsScalar max_curvature;
  TacsScalar aggregate_numer, aggregate_denom;
};

#endif // TMR_REFINEMENT_TOOLS_H
