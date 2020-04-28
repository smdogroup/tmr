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
#include "TMRHelmholtzModel.h"
#include "TMRMatrixCreator.h"
#include "TMR_TACSCreator.h"
#include "TMR_RefinementTools.h"

/*
  Create the Helmholtz filter object
*/
TMRHelmholtzFilter::TMRHelmholtzFilter( double helmholtz_radius,
                                        int _nlevels,
                                        TACSAssembler *_assembler[],
                                        TMROctForest *_filter[] ):
  TMRConformFilter(_nlevels, _assembler, _filter){
  initialize_helmholtz(helmholtz_radius);
}

TMRHelmholtzFilter::TMRHelmholtzFilter( double helmholtz_radius,
                                        int _nlevels,
                                        TACSAssembler *_assembler[],
                                        TMRQuadForest *_filter[] ):
  TMRConformFilter(_nlevels, _assembler, _filter){
  initialize_helmholtz(helmholtz_radius);
}

/*
  Initialize the Helmholtz filter data
*/
void TMRHelmholtzFilter::initialize_helmholtz( double helmholtz_radius ){
  // Allocate space for the assembler objects
  helmholtz_assembler = new TACSAssembler*[ nlevels ];

  // Create the Helmholtz creator objects
  TMROctTACSMatrixCreator *helmholtz_creator3d = NULL;
  TMRQuadTACSMatrixCreator *helmholtz_creator2d = NULL;

  if (oct_filter){
    TACSElementModel *model = new TMRHexaHelmholtzModel(helmholtz_radius);
    helmholtz_creator3d = new TMROctTACSMatrixCreator(model);
    helmholtz_creator3d->incref();
  }
  else {
    TACSElementModel *model = new TMRQuadHelmholtzModel(helmholtz_radius);
    helmholtz_creator2d = new TMRQuadTACSMatrixCreator(model);
    helmholtz_creator2d->incref();
  }

  for ( int k = 0; k < nlevels; k++ ){
    // Create the assembler object
    if (helmholtz_creator3d){
      helmholtz_assembler[k] =
        helmholtz_creator3d->createTACS(oct_filter[k],
                                        TACSAssembler::NATURAL_ORDER);
    }
    else {
      helmholtz_assembler[k] =
        helmholtz_creator2d->createTACS(quad_filter[k],
                                        TACSAssembler::NATURAL_ORDER);
    }
    helmholtz_assembler[k]->incref();
  }

  // Create a temporary vector to store the design variables from the
  // full TACSAssembler matrix
  helmholtz_vec = assembler[0]->createDesignVec();
  helmholtz_vec->incref();

  // Destroy the helmholtz creator object
  if (helmholtz_creator3d){
    helmholtz_creator3d->decref();
  }
  else {
    helmholtz_creator2d->decref();
  }

  // Create the vectors
  helmholtz_rhs = helmholtz_assembler[0]->createVec();
  helmholtz_psi = helmholtz_assembler[0]->createVec();
  helmholtz_rhs->incref();
  helmholtz_psi->incref();

  // Create the multigrid object
  double helmholtz_omega = 0.5;
  if (oct_filter){
    TMR_CreateTACSMg(nlevels, helmholtz_assembler, oct_filter,
                     &helmholtz_mg, helmholtz_omega);
  }
  else {
    TMR_CreateTACSMg(nlevels, helmholtz_assembler, quad_filter,
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
  temp = assembler[0]->createDesignVec();
  temp->incref();
}

TMRHelmholtzFilter::~TMRHelmholtzFilter(){
  helmholtz_mg->decref();
  helmholtz_ksm->decref();
  helmholtz_vec->decref();
  helmholtz_rhs->decref();
  helmholtz_psi->decref();

  for ( int k = 0; k < nlevels; k++ ){
    helmholtz_assembler[k]->decref();
  }
  delete [] helmholtz_assembler;
  temp->decref();
}

/*
  Apply the Helmholtz filter to the design variables.

  Here the input/output vector are the same
*/
void TMRHelmholtzFilter::applyFilter( TACSBVec *xvars ){
  // Get the number of design variables per node
  const int vars_per_node = assembler[0]->getDesignVarsPerNode();

  // Get the number of local elements
  const int num_elements = helmholtz_assembler[0]->getNumElements();

  // Get the maximum number of nodes per element (all elements in
  // this assembler object have the same number of nodes)
  const int max_nodes = helmholtz_assembler[0]->getMaxElementNodes();

  // Get all of the elements in the filter
  TACSElement **elements = helmholtz_assembler[0]->getElements();

  // Allocate space for element-level data
  TacsScalar *Xpts = new TacsScalar[ 3*max_nodes ];
  TacsScalar *x_values = new TacsScalar[ vars_per_node*max_nodes ];
  TacsScalar *rhs_values = new TacsScalar[ max_nodes ];

  for ( int k = 0; k < vars_per_node; k++ ){
    // Zero the entries in the RHS vector
    helmholtz_rhs->zeroEntries();
    helmholtz_assembler[0]->zeroVariables();

    for ( int i = 0; i < num_elements; i++ ){
      // Get the values for this element
      int len;
      const int *nodes;
      helmholtz_assembler[0]->getElement(i, &len, &nodes);
      helmholtz_assembler[0]->getElement(i, Xpts);

      // Get the values of the design variables at the nodes
      xvars->getValues(len, nodes, x_values);

     // Condense the design variables so that the x_values array
      // is of length (len) and not len*vars_per_node
      if (vars_per_node > 1){
        for ( int j = 0; j < len; j++ ){
          x_values[j] = x_values[vars_per_node*j + k];
        }
      }

      // Get the element basis (it should exist in this case)
      TACSElementBasis *basis = elements[i]->getElementBasis();

      // Zero the values on the right-hand-side
      memset(rhs_values, 0, max_nodes*sizeof(TacsScalar));

      if (basis){
        // Get the number of quadrature points
        const int nquad = basis->getNumQuadraturePoints();

        for ( int n = 0; n < nquad; n++ ){
          // Get the quadrature weight/point pair
          double pt[3];
          double wt = basis->getQuadraturePoint(n, pt);

          // Compute the Jacobian transformation
          TacsScalar Xd[9], J[9];
          TacsScalar detJ = basis->getJacobianTransform(n, pt, Xpts, Xd, J);

          // Get the field gradient for the design variable value
          TacsScalar X[3], U;
          basis->getFieldValues(n, pt, Xpts, 1, x_values, X, &U);

          TacsScalar DUt[3] = {0.0, 0.0, 0.0};
          TacsScalar DUx[3] = {0.0, 0.0, 0.0};
          DUt[0] = U;
          basis->addWeakFormResidual(n, pt, wt*detJ, J, 1, DUt, DUx, rhs_values);
        }
      }

      helmholtz_rhs->setValues(len, nodes, rhs_values, TACS_ADD_VALUES);
    }

    // Complete the assembly process
    helmholtz_rhs->beginSetValues(TACS_ADD_VALUES);
    helmholtz_rhs->endSetValues(TACS_ADD_VALUES);

    // Solve for the filtered values of the design variables
    helmholtz_ksm->solve(helmholtz_rhs, helmholtz_psi);
    helmholtz_assembler[0]->reorderVec(helmholtz_psi);
    helmholtz_assembler[0]->setVariables(helmholtz_psi);

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

  // Get the number of design variables per node
  const int vars_per_node = assembler[0]->getDesignVarsPerNode();

  // Get the number of local elements
  const int num_elements = helmholtz_assembler[0]->getNumElements();

  // Get the maximum number of nodes per element (all elements in
  // this assembler object have the same number of nodes)
  const int max_nodes = helmholtz_assembler[0]->getMaxElementNodes();

  // Get all of the elements in the filter
  TACSElement **elements = helmholtz_assembler[0]->getElements();

  // Allocate space for element-level data
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
    helmholtz_assembler[0]->reorderVec(helmholtz_psi);

    // Distribute the values from the solution
    helmholtz_psi->beginDistributeValues();
    helmholtz_psi->endDistributeValues();

    for ( int i = 0; i < num_elements; i++ ){
      // Get the values for this element
      int len;
      const int *nodes;
      helmholtz_assembler[0]->getElement(i, &len, &nodes);
      helmholtz_assembler[0]->getElement(i, Xpts);

      // Get the values of the design variables at the nodes
      helmholtz_psi->getValues(len, nodes, psi_values);

      // Get the element basis (it should exist in this case)
      TACSElementBasis *basis = elements[i]->getElementBasis();

      // Zero the values on the right-hand-side
      memset(x_values, 0, max_nodes*sizeof(TacsScalar));

      if (basis){
        // Get the number of quadrature points
        const int nquad = basis->getNumQuadraturePoints();

        for ( int n = 0; n < nquad; n++ ){
          // Get the quadrature weight/point pair
          double pt[3];
          double wt = basis->getQuadraturePoint(n, pt);

          // Compute the Jacobian transformation
          TacsScalar Xd[9], J[9];
          TacsScalar detJ = basis->getJacobianTransform(n, pt, Xpts, Xd, J);

          // Get the field gradient for the design variable value
          TacsScalar X[3], U;
          basis->getFieldValues(n, pt, Xpts, 1, psi_values, X, &U);

          TacsScalar DUt[3] = {0.0, 0.0, 0.0};
          TacsScalar DUx[3] = {0.0, 0.0, 0.0};
          DUt[0] = U;
          basis->addWeakFormResidual(n, pt, wt*detJ, J, 1, DUt, DUx, x_values);
        }
      }

      if (vars_per_node > 1){
        for ( int j = len-1; j >= 0; j-- ){
          x_values[vars_per_node*j + k] = x_values[j];
          if (!(j == 0 && k == 0)){
            x_values[j] = 0.0;
          }
        }
      }

      output->setValues(len, nodes, x_values, TACS_ADD_VALUES);
    }
  }

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
  assembler[0]->setDesignVars(x[0]);

  // Set the design variable values on all processors
  for ( int k = 0; k < nlevels-1; k++ ){
    filter_interp[k]->multWeightTranspose(x[k], x[k+1]);

    assembler[k]->setDesignVars(x[k+1]);
  }
}

/*
  Add values to the output TACSBVec
*/
void TMRHelmholtzFilter::addValues( TACSBVec *vec ){
  vec->beginSetValues(TACS_ADD_VALUES);
  vec->endSetValues(TACS_ADD_VALUES);

  temp->copyValues(vec);
  applyTranspose(temp, vec);
}
