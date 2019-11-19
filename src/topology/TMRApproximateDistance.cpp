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

#include <math.h>
#include "TMRApproximateDistance.h"
#include "KSM.h"
#include "TMRHelmholtzModel.h"
#include "TACSElement2D.h"
#include "TACSElement3D.h"
#include "TACSQuadBasis.h"
#include "TACSHexaBasis.h"
#include "TACSQuadBernsteinBasis.h"
#include "TACSHexaBernsteinBasis.h"
#include "TACSToFH5.h"

void TMRApproximateDistance( TMRQuadForest *filter, int index,
                             double cutoff, double t,
                             TACSBVec *rho, const char *filename,
                             TACSBVec **_dist ){
  int mpi_rank;
  MPI_Comm comm = filter->getMPIComm();
  MPI_Comm_rank(comm, &mpi_rank);

  // Create the nodes for the filter
  filter->createNodes();

  // Get the mesh order
  int order = filter->getMeshOrder();

  // Get the connectivity
  const int *conn;
  int num_elements = 0, num_owned_nodes = 0;
  filter->getNodeConn(&conn, &num_elements, &num_owned_nodes);

  // Allocate the pointer array
  int *ptr = new int[ num_elements+1 ];
  for ( int i = 0; i <= num_elements; i++ ){
    ptr[i] = order*order*i;
  }

  // Create/retrieve the dependent node information
  const int *dep_ptr, *dep_conn;
  const double *dep_weights;
  int num_dep_nodes = filter->getDepNodeConn(&dep_ptr, &dep_conn,
                                             &dep_weights);

  // Create the elements using the virtual call
  TACSElement **elements = new TACSElement*[ num_elements ];

  // Allocate the model class
  TACSElementModel *model = new TMRQuadHelmholtzModel(t);

  // Set the basis class
  TACSElementBasis *basis = NULL;
  if (order == 2){
    basis = new TACSLinearQuadBasis();
  }
  else if (order == 3){
    basis = new TACSQuadraticQuadBernsteinBasis();
  }
  else if (order == 4){
    basis = new TACSCubicQuadBernsteinBasis();
  }

  TACSElement *elem = new TACSElement2D(model, basis);
  for ( int i = 0; i < num_elements; i++ ){
    elements[i] = elem;
  }

  // Create the associated TACSAssembler object
  TACSAssembler *assembler =
    new TACSAssembler(comm, 1,
                      num_owned_nodes, num_elements, num_dep_nodes);
  assembler->incref();

  // Set the element connectivity into TACSAssembler
  assembler->setElementConnectivity(ptr, conn);
  delete [] ptr;

  // Set the dependent node information
  assembler->setDependentNodes(dep_ptr, dep_conn, dep_weights);

  // Set the element array
  assembler->setElements(elements);
  delete [] elements;

  // Determine where to place the boundary conditions
  TacsScalar *rho_array;
  int size = rho->getArray(&rho_array);
  int bsize = rho->getBlockSize();

  // Get the ownership range of the filter
  const int *range;
  filter->getOwnedNodeRange(&range);

  // Check if this is an intermediate value of the design variables
  for ( int i = 0; i < size; i++ ){
    if (rho_array[bsize*i + index] >= cutoff &&
        rho_array[bsize*i + index] < 1.0 - cutoff){
      // Set the global number
      int node = range[mpi_rank] + i;

      int nbcs = 1;
      int vars = 0;
      TacsScalar value = 1.0;
      assembler->addBCs(1, &node, nbcs, &vars, &value);
    }
  }

  assembler->initialize();

  // Set the node locations
  const int *nodes;
  int num_local_nodes = filter->getNodeNumbers(&nodes);

  // Get the points
  TMRPoint *Xp;
  filter->getPoints(&Xp);

  TACSBVec *X = assembler->createNodeVec();
  X->incref();

  // Get the node array from the TACSBVec object
  TacsScalar *Xn;
  X->getArray(&Xn);

   // Loop over all the nodes
  for ( int i = 0; i < num_local_nodes; i++ ){
    if (nodes[i] >= range[mpi_rank] &&
        nodes[i] < range[mpi_rank+1]){
      int loc = nodes[i] - range[mpi_rank];
      Xn[3*loc] = Xp[i].x;
      Xn[3*loc+1] = Xp[i].y;
      Xn[3*loc+2] = Xp[i].z;
    }
  }

  // Reorder the vector if needed
  assembler->reorderVec(X);
  assembler->setNodes(X);
  X->decref();

  // Approximately solve the equations....
  TACSParallelMat *mat = assembler->createMat();
  int zero_guess = 1;
  TacsScalar omega = 1.0;
  int niters = 5;
  int symm = 1;
  TACSGaussSeidel *pc = new TACSGaussSeidel(mat, zero_guess, omega, niters, symm);

  // Allocate the GMRES object
  int m = 20;
  int nrestart = 4;
  int isflexible = 1;
  GMRES *gmres = new GMRES(mat, pc, m, nrestart, isflexible);
  gmres->incref();

  TACSBVec *rhs = assembler->createVec();
  rhs->incref();

  TACSBVec *dist = assembler->createVec();

  // Assemble the matrix
  assembler->assembleJacobian(1.0, 0.0, 0.0, NULL, mat);
  pc->factor();

  // Set boundary conditions for the right-hand-side (w(x) = 1 on boundary)
  assembler->setBCs(rhs);

  // Set a monitor for the Helmholtz equation
  gmres->setMonitor(new KSMPrintStdout("Helmholtz", mpi_rank, 10));

  // Solve the discrete Helmholtz equation to obtain the approximate distance
  // function
  gmres->solve(rhs, dist);

  // Deallocate the variable values
  gmres->decref();
  rhs->decref();

  // Extract the actual distance using the transformation
  TacsScalar *dist_array;
  dist->getArray(&dist_array);
  for ( int i = 0; i < size; i++ ){
    if (dist_array[i] <= 0.0){
      dist_array[i] = 1e20;
    }
    else if (dist_array[i] >= 1.0){
      dist_array[i] = 0.0;
    }
    else {
      dist_array[i] = -t*log(dist_array[i]);
    }
  }

  // Distribute the distance values
  dist->beginDistributeValues();
  dist->endDistributeValues();

  // Write out the result to an f5 file
  if (filename){
    assembler->setVariables(dist);

    int write_flag = (TACS_OUTPUT_CONNECTIVITY |
                      TACS_OUTPUT_NODES |
                      TACS_OUTPUT_DISPLACEMENTS);
    TACSToFH5 *f5 = new TACSToFH5(assembler, TACS_SCALAR_2D_ELEMENT,
                                  write_flag);
    f5->incref();
    f5->writeToFile(filename);
    f5->decref();
  }

  // Deallocate the TACSAssembler object
  assembler->decref();

  *_dist = dist;
}

void TMRApproximateDistance( TMROctForest *filter, int index,
                             double cutoff, double t,
                             TACSBVec *rho, const char *filename,
                             TACSBVec **_dist ){
  int mpi_rank;
  MPI_Comm comm = filter->getMPIComm();
  MPI_Comm_rank(comm, &mpi_rank);

  // Create the nodes for the filter
  filter->createNodes();

  // Get the mesh order
  int order = filter->getMeshOrder();

  // Get the connectivity
  const int *conn;
  int num_elements = 0, num_owned_nodes = 0;
  filter->getNodeConn(&conn, &num_elements, &num_owned_nodes);

  // Allocate the pointer array
  int *ptr = new int[ num_elements+1 ];
  for ( int i = 0; i <= num_elements; i++ ){
    ptr[i] = order*order*order*i;
  }

  // Create/retrieve the dependent node information
  const int *dep_ptr, *dep_conn;
  const double *dep_weights;
  int num_dep_nodes = filter->getDepNodeConn(&dep_ptr, &dep_conn,
                                             &dep_weights);

  // Create the elements using the virtual call
  TACSElement **elements = new TACSElement*[ num_elements ];

  // Allocate the model class
  TACSElementModel *model = new TMRHexaHelmholtzModel(t);

  // Set the basis class
  TACSElementBasis *basis = NULL;
  if (order == 2){
    basis = new TACSLinearHexaBasis();
  }
  else if (order == 3){
    basis = new TACSQuadraticHexaBernsteinBasis();
  }
  else if (order == 4){
    basis = new TACSCubicHexaBernsteinBasis();
  }

  TACSElement *elem = new TACSElement3D(model, basis);
  for ( int i = 0; i < num_elements; i++ ){
    elements[i] = elem;
  }

  // Create the associated TACSAssembler object
  TACSAssembler *assembler =
    new TACSAssembler(comm, 1,
                      num_owned_nodes, num_elements, num_dep_nodes);
  assembler->incref();

  // Set the element connectivity into TACSAssembler
  assembler->setElementConnectivity(ptr, conn);
  delete [] ptr;

  // Set the dependent node information
  assembler->setDependentNodes(dep_ptr, dep_conn, dep_weights);

  // Set the element array
  assembler->setElements(elements);
  delete [] elements;

  // Determine where to place the boundary conditions
  TacsScalar *rho_array;
  int size = rho->getArray(&rho_array);
  int bsize = rho->getBlockSize();

  // Get the ownership range of the filter
  const int *range;
  filter->getOwnedNodeRange(&range);

  // Check if this is an intermediate value of the design variables
  for ( int i = 0; i < size; i++ ){
    if (rho_array[bsize*i + index] >= cutoff &&
        rho_array[bsize*i + index] < 1.0 - cutoff){
      // Set the global number
      int node = range[mpi_rank] + i;

      int nbcs = 1;
      int vars = 0;
      TacsScalar value = 1.0;
      assembler->addBCs(1, &node, nbcs, &vars, &value);
    }
  }

  assembler->initialize();

  // Set the node locations
  const int *nodes;
  int num_local_nodes = filter->getNodeNumbers(&nodes);

  // Get the points
  TMRPoint *Xp;
  filter->getPoints(&Xp);

  TACSBVec *X = assembler->createNodeVec();
  X->incref();

  // Get the node array from the TACSBVec object
  TacsScalar *Xn;
  X->getArray(&Xn);

   // Loop over all the nodes
  for ( int i = 0; i < num_local_nodes; i++ ){
    if (nodes[i] >= range[mpi_rank] &&
        nodes[i] < range[mpi_rank+1]){
      int loc = nodes[i] - range[mpi_rank];
      Xn[3*loc] = Xp[i].x;
      Xn[3*loc+1] = Xp[i].y;
      Xn[3*loc+2] = Xp[i].z;
    }
  }

  // Reorder the vector if needed
  assembler->reorderVec(X);
  assembler->setNodes(X);
  X->decref();

  // Approximately solve the equations....
  TACSParallelMat *mat = assembler->createMat();
  int zero_guess = 1;
  TacsScalar omega = 1.0;
  int niters = 5;
  int symm = 1;
  TACSGaussSeidel *pc = new TACSGaussSeidel(mat, zero_guess, omega, niters, symm);

  // Allocate the GMRES object
  int m = 20;
  int nrestart = 4;
  int isflexible = 1;
  GMRES *gmres = new GMRES(mat, pc, m, nrestart, isflexible);
  gmres->incref();

  TACSBVec *rhs = assembler->createVec();
  rhs->incref();

  TACSBVec *dist = assembler->createVec();

  // Assemble the matrix
  assembler->assembleJacobian(1.0, 0.0, 0.0, NULL, mat);
  pc->factor();

  // Set boundary conditions for the right-hand-side (w(x) = 1 on boundary)
  assembler->setBCs(rhs);

  // Set a monitor for the Helmholtz equation
  gmres->setMonitor(new KSMPrintStdout("Helmholtz", mpi_rank, 10));

  // Solve the discrete Helmholtz equation to obtain the approximate distance
  // function
  gmres->solve(rhs, dist);

  // Deallocate the variable values
  gmres->decref();
  rhs->decref();

  // Extract the actual distance using the transformation
  TacsScalar *dist_array;
  dist->getArray(&dist_array);
  for ( int i = 0; i < size; i++ ){
    if (dist_array[i] <= 0.0){
      dist_array[i] = 1e20;
    }
    else if (dist_array[i] >= 1.0){
      dist_array[i] = 0.0;
    }
    else {
      dist_array[i] = -t*log(dist_array[i]);
    }
  }

  // Distribute the distance values
  dist->beginDistributeValues();
  dist->endDistributeValues();

  // Write out the result to an f5 file
  if (filename){
    assembler->setVariables(dist);

    int write_flag = (TACS_OUTPUT_CONNECTIVITY |
                      TACS_OUTPUT_NODES |
                      TACS_OUTPUT_DISPLACEMENTS);
    TACSToFH5 *f5 = new TACSToFH5(assembler, TACS_SCALAR_3D_ELEMENT,
                                  write_flag);
    f5->incref();
    f5->writeToFile(filename);
    f5->decref();
  }

  // Deallocate the TACSAssembler object
  assembler->decref();

  *_dist = dist;
}
