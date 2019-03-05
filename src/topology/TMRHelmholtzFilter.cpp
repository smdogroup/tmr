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

#include "TMRHelmholtzFilter.h"
#include "TMRHelmholtzElement.h"
#include "TMR_TACSCreator.h"
#include "TMR_RefinementTools.h"

/*
  Helmholtz filter creator classes
*/
class TMRQuadTACSHelmholtzCreator : public TMRQuadTACSCreator {
public:
  TMRQuadTACSHelmholtzCreator( double r ):
    TMRQuadTACSCreator(NULL){
    radius = r;
  }
  void createElements( int order,
                       TMRQuadForest *forest,
                       int num_elements,
                       TACSElement **elements ){
    TACSElement *elem = NULL;
    if (order == 2){
      elem = new TMRQuadHelmholtz<2>(radius);
    }
    else if (order == 3){
      elem = new TMRQuadHelmholtz<3>(radius);
    }
    else if (order == 4){
      elem = new TMRQuadHelmholtz<4>(radius);
    }
    else if (order == 5){
      elem = new TMRQuadHelmholtz<5>(radius);
    }
    else if (order == 6){
      elem = new TMRQuadHelmholtz<6>(radius);
    }

    for ( int i = 0; i < num_elements; i++ ){
      elements[i] = elem;
    }
  }

  double radius;
};

class TMROctTACSHelmholtzCreator : public TMROctTACSCreator {
public:
  TMROctTACSHelmholtzCreator( double r ):
    TMROctTACSCreator(NULL){
    radius = r;
  }
  void createElements( int order,
                       TMROctForest *forest,
                       int num_elements,
                       TACSElement **elements ){
    TACSElement *elem = NULL;
    if (order == 2){
      elem = new TMROctHelmholtz<2>(radius);
    }
    else if (order == 3){
      elem = new TMROctHelmholtz<3>(radius);
    }
    else if (order == 4){
      elem = new TMROctHelmholtz<4>(radius);
    }
    else if (order == 5){
      elem = new TMROctHelmholtz<5>(radius);
    }
    else if (order == 6){
      elem = new TMROctHelmholtz<6>(radius);
    }

    for ( int i = 0; i < num_elements; i++ ){
      elements[i] = elem;
    }
  }

  double radius;
};

/*
  Create the Helmholtz filter object
*/
TMRHelmholtzFilter::TMRHelmholtzFilter( double helmholtz_radius,
                                        int _nlevels,
                                        TACSAssembler *_tacs[],
                                        TMROctForest *_filter[],
                                        int _vars_per_node ):
  TMRConformFilter(_nlevels, _tacs, _filter, _vars_per_node){
  initialize_helmholtz(helmholtz_radius);
}

TMRHelmholtzFilter::TMRHelmholtzFilter( double helmholtz_radius,
                                        int _nlevels,
                                        TACSAssembler *_tacs[],
                                        TMRQuadForest *_filter[],
                                        int _vars_per_node ):
  TMRConformFilter(_nlevels, _tacs, _filter, _vars_per_node){
  initialize_helmholtz(helmholtz_radius);
}

/*
  Initialize the Helmholtz filter data
*/
void TMRHelmholtzFilter::initialize_helmholtz( double helmholtz_radius ){

  // Allocate space for the assembler objects
  helmholtz_tacs = new TACSAssembler*[ nlevels ];

  // Create the Helmholtz creator objects
  TMROctTACSHelmholtzCreator *helmholtz_creator3d = NULL;
  TMRQuadTACSHelmholtzCreator *helmholtz_creator2d = NULL;

  if (oct_filter){
    helmholtz_creator3d = new TMROctTACSHelmholtzCreator(helmholtz_radius);
    helmholtz_creator3d->incref();
  }
  else {
    helmholtz_creator2d = new TMRQuadTACSHelmholtzCreator(helmholtz_radius);
    helmholtz_creator2d->incref();
  }

  for ( int k = 0; k < nlevels; k++ ){
    // Create the assembler object
    if (helmholtz_creator3d){
      helmholtz_tacs[k] =
        helmholtz_creator3d->createTACS(oct_filter[k],
                                        TACSAssembler::NATURAL_ORDER);
    }
    else {
      helmholtz_tacs[k] =
        helmholtz_creator2d->createTACS(quad_filter[k],
                                        TACSAssembler::NATURAL_ORDER);
    }
    helmholtz_tacs[k]->incref();
  }

  // Create a temporary vector to store the design variables
  helmholtz_vec =
    new TACSBVec(filter_maps[0], vars_per_node, filter_dist[0],
                 filter_dep_nodes[0]);
  helmholtz_vec->incref();

  // Destroy the helmholtz creator object
  if (helmholtz_creator3d){
    helmholtz_creator3d->decref();
  }
  else {
    helmholtz_creator2d->decref();
  }

  // Create the vectors
  helmholtz_rhs = helmholtz_tacs[0]->createVec();
  helmholtz_psi = helmholtz_tacs[0]->createVec();
  helmholtz_rhs->incref();
  helmholtz_psi->incref();

  // Create the multigrid object
  double helmholtz_omega = 0.5;
  if (oct_filter){
    TMR_CreateTACSMg(nlevels, helmholtz_tacs, oct_filter,
                     &helmholtz_mg, helmholtz_omega);
  }
  else {
    TMR_CreateTACSMg(nlevels, helmholtz_tacs, quad_filter,
                     &helmholtz_mg, helmholtz_omega);
  }
  helmholtz_mg->incref();

  // Get the rank
  int mpi_rank;
  MPI_Comm_rank(getMPIComm(), &mpi_rank);

  // Set up the solver
  int gmres_iters = 50;
  int nrestart = 5;
  int is_flexible = 0;

  // Create the GMRES object
  helmholtz_ksm = new GMRES(helmholtz_mg->getMat(0), helmholtz_mg,
                            gmres_iters, nrestart, is_flexible);
  helmholtz_ksm->incref();
  helmholtz_ksm->setMonitor(new KSMPrintStdout("Filter GMRES", mpi_rank, 10));
  helmholtz_ksm->setTolerances(1e-12, 1e-30);

  // Create the multigrid solver and factor it
  double alpha = 1.0, beta = 0.0, gamma = 0.0;
  helmholtz_mg->assembleJacobian(alpha, beta, gamma, NULL);
  helmholtz_mg->factor();

  // Create a temporary vector
  temp = createVec();
  temp->incref();
}

TMRHelmholtzFilter::~TMRHelmholtzFilter(){
  helmholtz_mg->decref();
  helmholtz_ksm->decref();
  helmholtz_vec->decref();
  helmholtz_rhs->decref();
  helmholtz_psi->decref();

  for ( int k = 0; k < nlevels; k++ ){
    helmholtz_tacs[k]->decref();
  }
  delete [] helmholtz_tacs;
  temp->decref();
}

/*
  Apply the Helmholtz filter to the design variables.

  Here the input/output vector are the same
*/
void TMRHelmholtzFilter::applyFilter( TACSBVec *xvars ){
  // Get the number of local elements
  int num_elements = helmholtz_tacs[0]->getNumElements();

  // Get the maximum number of nodes per element (all elements in
  // this assembler object have the same number of nodes)
  int max_nodes = helmholtz_tacs[0]->getMaxElementNodes();

  // Get all of the elements in the filter
  TACSElement **elements = helmholtz_tacs[0]->getElements();

  // Allocate space for element-level data
  double *N = new double[ max_nodes ];
  TacsScalar *Xpts = new TacsScalar[ 3*max_nodes ];
  TacsScalar *x_values = new TacsScalar[ vars_per_node*max_nodes ];
  TacsScalar *rhs_values = new TacsScalar[ max_nodes ];

  for ( int k = 0; k < vars_per_node; k++ ){
    // Zero the entries in the RHS vector
    helmholtz_rhs->zeroEntries();
    helmholtz_tacs[0]->zeroVariables();

    for ( int i = 0; i < num_elements; i++ ){
      // Get the values for this element
      int len;
      const int *nodes;
      helmholtz_tacs[0]->getElement(i, &nodes, &len);
      helmholtz_tacs[0]->getElement(i, Xpts);

      // Get the values of the design variables at the nodes
      xvars->getValues(len, nodes, x_values);

      // Zero the values on the right-hand-side
      memset(rhs_values, 0, max_nodes*sizeof(TacsScalar));

      // Get the number of quadrature points
      const int num_quad_pts = elements[i]->getNumGaussPts();

      // Perform the integration over the element
      for ( int n = 0; n < num_quad_pts; n++ ){
        double pt[3];
        double wt = elements[i]->getGaussWtsPts(n, pt);
        TacsScalar h = elements[i]->getDetJacobian(pt, Xpts);
        h *= wt;

        elements[i]->getShapeFunctions(pt, N);

        for ( int ii = 0; ii < max_nodes; ii++ ){
          const double *ns = N;
          const TacsScalar *xs = &x_values[k];
          TacsScalar v = 0.0;
          for ( int jj = 0; jj < max_nodes; jj++ ){
            v += ns[0]*xs[0];
            ns++;
            xs += vars_per_node;
          }
          rhs_values[ii] += h*N[ii]*v;
        }
      }

      helmholtz_rhs->setValues(len, nodes, rhs_values, TACS_ADD_VALUES);
    }

    // Complete the assembly process
    helmholtz_rhs->beginSetValues(TACS_ADD_VALUES);
    helmholtz_rhs->endSetValues(TACS_ADD_VALUES);

    // Solve for the filtered values of the design variables
    helmholtz_ksm->solve(helmholtz_rhs, helmholtz_psi);
    helmholtz_tacs[0]->reorderVec(helmholtz_psi);
    helmholtz_tacs[0]->setVariables(helmholtz_psi);

    // Get the output array
    TacsScalar *xarr;
    xvars->getArray(&xarr);

    // Get the Helmholtz solution vector
    TacsScalar *hpsi;
    const int size = helmholtz_psi->getArray(&hpsi);

    for ( int i = 0; i < size; i++ ){
      xarr[vars_per_node*i + k] = hpsi[i];
    }
  }

  delete [] N;
  delete [] Xpts;
  delete [] x_values;
  delete [] rhs_values;
}

/*
  Compute the sensitivity w.r.t the Helmholtz filter
*/
void TMRHelmholtzFilter::applyTranspose( TACSBVec *input,
                                         TACSBVec *output ){
  output->zeroEntries();
  // Get the number of local elements
  int num_elements = helmholtz_tacs[0]->getNumElements();

  // Get the maximum number of nodes per element (all elements in
  // this assembler object have the same number of nodes)
  int max_nodes = helmholtz_tacs[0]->getMaxElementNodes();

  // Get all of the elements in the filter
  TACSElement **elements = helmholtz_tacs[0]->getElements();

  // Allocate space for element-level data
  double *N = new double[ max_nodes ];
  TacsScalar *Xpts = new TacsScalar[ 3*max_nodes ];
  TacsScalar *x_values = new TacsScalar[ vars_per_node*max_nodes ];
  TacsScalar *psi_values = new TacsScalar[ max_nodes ];

  for ( int k = 0; k < vars_per_node; k++ ){
    // Get the input array
    TacsScalar *xarr;
    input->getArray(&xarr);

    // Get the Helmholtz solution vector
    TacsScalar *hrhs;
    const int size = helmholtz_rhs->getArray(&hrhs);

    for ( int i = 0; i < size; i++ ){
      hrhs[i] = xarr[vars_per_node*i + k];
    }

    // Solve for the filtered values of the design variables
    helmholtz_ksm->solve(helmholtz_rhs, helmholtz_psi);
    helmholtz_tacs[0]->reorderVec(helmholtz_psi);

    // Distribute the values from the solution
    helmholtz_psi->beginDistributeValues();
    helmholtz_psi->endDistributeValues();

    for ( int i = 0; i < num_elements; i++ ){
      // Get the values for this element
      int len;
      const int *nodes;
      helmholtz_tacs[0]->getElement(i, &nodes, &len);
      helmholtz_tacs[0]->getElement(i, Xpts);

      // Get the values of the design variables at the nodes
      helmholtz_psi->getValues(len, nodes, psi_values);

      // Zero the values on the right-hand-side
      memset(x_values, 0, vars_per_node*max_nodes*sizeof(TacsScalar));

      // Get the number of quadrature points
      const int num_quad_pts = elements[i]->getNumGaussPts();

      // Perform the integration over the element
      for ( int n = 0; n < num_quad_pts; n++ ){
        double pt[3];
        double wt = elements[i]->getGaussWtsPts(n, pt);
        TacsScalar h = elements[i]->getDetJacobian(pt, Xpts);
        h *= wt;

        elements[i]->getShapeFunctions(pt, N);

        for ( int ii = 0; ii < max_nodes; ii++ ){
          const double *ns = N;
          const TacsScalar *psi = &psi_values[0];
          TacsScalar v = 0.0;
          for ( int jj = 0; jj < max_nodes; jj++ ){
            v += ns[0]*psi[0];
            ns++;
            psi++;
          }

          x_values[vars_per_node*ii + k] += h*N[ii]*v;
        }
      }

      output->setValues(len, nodes, x_values, TACS_ADD_VALUES);
    }
  }

  delete [] N;
  delete [] Xpts;
  delete [] x_values;
  delete [] psi_values;

  output->beginSetValues(TACS_ADD_VALUES);
  output->endSetValues(TACS_ADD_VALUES);
}

/*
  Set the design variables for each level
*/
void TMRHelmholtzFilter::setDesignVars( TACSBVec *xvec ){
  // Copy the values to the local design variable vector
  x[0]->copyValues(xvec);

  applyFilter(x[0]);

  // Distribute the design variable values
  x[0]->beginDistributeValues();
  x[0]->endDistributeValues();

  // Temporarily allocate an array to store the variables
  TacsScalar *xlocal = new TacsScalar[ getMaxNumLocalVars() ];

  // Copy the values to the local array
  int size = getLocalValuesFromBVec(0, x[0], xlocal);
  tacs[0]->setDesignVars(xlocal, size);

  // Set the design variable values on all processors
  for ( int k = 0; k < nlevels-1; k++ ){
    filter_interp[k]->multWeightTranspose(x[k], x[k+1]);
    // Distribute the design variable values
    x[k+1]->beginDistributeValues();
    x[k+1]->endDistributeValues();

    // Set the design variable values
    size = getLocalValuesFromBVec(k+1, x[k+1], xlocal);
    tacs[k+1]->setDesignVars(xlocal, size);
  }

  delete [] xlocal;
}

/*
  Add values to the output TACSBVec
*/
void TMRHelmholtzFilter::addValues( TacsScalar *xlocal, TACSBVec *vec ){
  setBVecFromLocalValues(0, xlocal, temp, TACS_ADD_VALUES);
  temp->beginSetValues(TACS_ADD_VALUES);
  temp->endSetValues(TACS_ADD_VALUES);

  applyTranspose(temp, vec);
}
