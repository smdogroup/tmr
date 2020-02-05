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
                             double *min_dist ){
  int mpi_rank;
  MPI_Comm comm = filter->getMPIComm();
  MPI_Comm_rank(comm, &mpi_rank);

  // Ensure that the nodes exist and then get the mesh order
  filter->createNodes();
  int order = filter->getMeshOrder();
  const int *filter_conn;
  int num_filter_elements = 0;
  filter->getNodeConn(&filter_conn, &num_filter_elements);

  // Get the block size
  int bsize = rho->getBlockSize();

  double *N = new double[ order*order ];
  TacsScalar *rho_values = new TacsScalar[ bsize*order*order ];
  TacsScalar *rho_local = new TacsScalar[ bsize ];

  // Distribute the design values
  rho->beginDistributeValues();
  rho->endDistributeValues();

  TMRQuadForest *forest = filter->duplicate();
  forest->incref();

  // Set the mesh order
  forest->setMeshOrder(2);

  // Set the number of refinements
  int num_refine = 0;
  if (order <= 3){
    num_refine = 1;
  }
  else if (order <= 5){
    num_refine = 2;
  }
  else {
    num_refine = 3;
  }

  // Search for those filter elements that need to be refined...
  int *refine = new int[ num_filter_elements ];

  // Flag those nodes which have low, high or intermediate values
  int *low = new int[ bsize ];
  int *high = new int[ bsize ];
  int *intermediate = new int[ bsize ];

  for ( int i = 0; i < num_filter_elements; i++ ){
    refine[i] = 0;
    rho->getValues(order*order, &filter_conn[order*order*i], rho_values);

    // Reset the flags
    for ( int kk = 0; kk < bsize; kk++ ){
      low[kk] = 0;
      high[kk] = 0;
      intermediate[kk] = 0;
    }

    // Scan through all the local node values for the filter
    if (index >= 0 && index < bsize){
      for ( int j = 0; j < order*order; j++ ){
        if (rho_values[bsize*j + index] < cutoff){
          low[0] = 1;
        }
        else if (rho_values[bsize*j + index] > 1.0 - cutoff){
          high[0] = 1;
        }
        else {
          intermediate[0] = 1;
        }
      }
    }
    else {
      for ( int j = 0; j < order*order; j++ ){
        for ( int kk = 0; kk < bsize; kk++ ){
          if (rho_values[bsize*j + kk] < cutoff){
            low[kk] = 1;
          }
          else if (rho_values[bsize*j + kk] > 1.0 - cutoff){
            high[kk] = 1;
          }
          else {
            intermediate[kk] = 1;
          }
        }
      }
    }

    if (index >= 0 && index < bsize){
      if ((low[0] && high[0]) || intermediate[0]){
        refine[i] = num_refine;
      }
    }
    else {
      for ( int kk = 0; kk < bsize; kk++ ){
        if ((low[kk] && high[kk]) || intermediate[kk]){
          refine[i] = num_refine;
        }
      }
    }
  }

  delete [] low;
  delete [] high;
  delete [] intermediate;

  // Refine the forest and balance it
  forest->refine(refine);
  delete [] refine;

  forest->balance(0);

  // Balance the forest and create the nodes
  forest->createNodes();

  // Get the connectivity
  const int *conn;
  int num_elements = 0, num_owned_nodes = 0;
  forest->getNodeConn(&conn, &num_elements, &num_owned_nodes);

  // Allocate the pointer array
  int *ptr = new int[ num_elements+1 ];
  for ( int i = 0; i <= num_elements; i++ ){
    ptr[i] = 4*i;
  }

  // Create/retrieve the dependent node information
  const int *dep_ptr, *dep_conn;
  const double *dep_weights;
  int num_dep_nodes = forest->getDepNodeConn(&dep_ptr, &dep_conn,
                                             &dep_weights);

  // Create the elements using the virtual call
  TACSElement **elements = new TACSElement*[ num_elements ];

  // Allocate the model class
  TACSElementModel *model = new TMRQuadHelmholtzModel(t);

  // Set the basis class
  TACSElementBasis *basis = new TACSLinearQuadBasis();

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

  // Filter quadrants
  TMRQuadrant *filter_quads;
  TMRQuadrantArray *filter_quad_array;
  filter->getQuadrants(&filter_quad_array);
  filter_quad_array->getArray(&filter_quads, NULL);

  // Get the array of the quadrants from the mesh
  TMRQuadrant *quads;
  TMRQuadrantArray *quad_array;
  forest->getQuadrants(&quad_array);
  quad_array->getArray(&quads, NULL);

  for ( int i = 0; i < num_elements; i++ ){
    // Set the node location
    TMRQuadrant n = quads[i];
    const int32_t hf = 1 << (TMR_MAX_LEVEL - quads[i].level);

    for ( int jj = 0; jj < 2; jj++ ){
      for ( int ii = 0; ii < 2; ii++ ){
        // Only execute the loop if the node is positive
        int node = conn[4*i + ii + 2*jj];
        if (node >= 0){
          // Set the node location
          n.info = ii + 2*jj;

          // Find the enclosing element on the filter mesh
          const double knots[] = {-1.0, 1.0};
          TMRQuadrant *filter_quad = filter->findEnclosing(2, knots, &n);

          if (filter_quad){
            // Get the values of the density at the nodes
            int filter_elem = filter_quad->tag;
            if (filter_elem >= 0 && filter_elem < num_filter_elements){
              rho->getValues(order*order,
                             &filter_conn[order*order*filter_elem],
                             rho_values);

              // Get the size of the filter element
              const int32_t h = 1 << (TMR_MAX_LEVEL - filter_quad->level);

              double pt[2];
              pt[0] = -1.0 + 2.0*(quads[i].x + hf*ii - filter_quad->x)/h;
              pt[1] = -1.0 + 2.0*(quads[i].y + hf*jj - filter_quad->y)/h;

              // Compute the interpolation
              filter->evalInterp(pt, N);

              // Evaluate the local values of rho based on the interpolation
              memset(rho_local, 0, bsize*sizeof(TacsScalar));
              for ( int k = 0; k < order*order; k++ ){
                for ( int kk = 0; kk < bsize; kk++ ){
                  rho_local[kk] += rho_values[k*bsize + kk]*N[k];
                }
              }

              // Check if this is an intermediate value of the design variables
              if (index >= 0 && index < bsize){
                if (rho_local[index] >= cutoff &&
                    rho_local[index] <= 1.0 - cutoff){
                  int nbcs = 1;
                  int vars = 0;
                  TacsScalar value = 1.0;
                  assembler->addBCs(1, &node, nbcs, &vars, &value);
                }
              }
              else {
                // Loop over all the dimension
                for ( int kk = 0; kk < bsize; kk++ ){
                  if (rho_local[kk] >= cutoff &&
                      rho_local[kk] <= 1.0 - cutoff){
                    int nbcs = 1;
                    int vars = 0;
                    TacsScalar value = 1.0;
                    assembler->addBCs(1, &node, nbcs, &vars, &value);
                    break;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Free the data
  delete [] N;
  delete [] rho_values;
  delete [] rho_local;

  assembler->initialize();

  // Get the ownership range of the forest
  const int *range;
  forest->getOwnedNodeRange(&range);

  // Set the node locations
  const int *nodes;
  int num_local_nodes = forest->getNodeNumbers(&nodes);

  // Get the points
  TMRPoint *Xp;
  forest->getPoints(&Xp);

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
  gmres->setTolerances(1e-12, 1e-30);

  TACSBVec *rhs = assembler->createVec();
  TACSBVec *dist = assembler->createVec();
  rhs->incref();
  dist->incref();

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

  dist->beginDistributeValues();
  dist->endDistributeValues();

  // Extract the actual distance using the transformation
  for ( int i = 0, index = 0; i < num_filter_elements; i++ ){
    double max_val = -1e20;

    while (index < num_elements &&
           filter_quads[i].contains(&quads[index])){
      for ( int ii = 0; ii < 4; ii++ ){
        // Only execute the loop if the node is positive
        int node = conn[4*index + ii];
        if (node >= 0){
          TacsScalar value;
          dist->getValues(1, &node, &value);
          if (value > max_val){
            max_val = value;
          }
        }
      }
      index++;
    }

    min_dist[i] = 1e20;
    if (max_val > 0.0){
      min_dist[i] = - t*log(max_val);
    }
  }

  // Write out the result to an f5 file
  if (filename){
    // Convert to a distance
    TacsScalar *dist_array;
    int size = dist->getArray(&dist_array);
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

  // Deallocate objects
  dist->decref();
  assembler->decref();
  forest->decref();
}

void TMRApproximateDistance( TMROctForest *filter, int index,
                             double cutoff, double t,
                             TACSBVec *rho, const char *filename,
                             double *min_dist ){
  int mpi_rank;
  MPI_Comm comm = filter->getMPIComm();
  MPI_Comm_rank(comm, &mpi_rank);

  // Ensure that the nodes exist and then get the mesh order
  filter->createNodes();
  int order = filter->getMeshOrder();
  const int *filter_conn;
  int num_filter_elements = 0;
  filter->getNodeConn(&filter_conn, &num_filter_elements);

  // Get the block size
  int bsize = rho->getBlockSize();

  double *N = new double[ order*order*order ];
  TacsScalar *rho_values = new TacsScalar[ bsize*order*order*order ];
  TacsScalar *rho_local = new TacsScalar[ bsize ];

  // Distribute the design values
  rho->beginDistributeValues();
  rho->endDistributeValues();

  TMROctForest *forest = filter->duplicate();
  forest->incref();

  // Set the mesh order
  forest->setMeshOrder(2);

  // Set the number of refinements
  int num_refine = 0;
  if (order <= 3){
    num_refine = 1;
  }
  else if (order <= 5){
    num_refine = 2;
  }
  else {
    num_refine = 3;
  }

  // Search for those filter elements that need to be refined...
  int *refine = new int[ num_filter_elements ];

  // Flag those nodes which have low, high or intermediate values
  int *low = new int[ bsize ];
  int *high = new int[ bsize ];
  int *intermediate = new int[ bsize ];

  for ( int i = 0; i < num_filter_elements; i++ ){
    refine[i] = 0;
    rho->getValues(order*order*order,
                   &filter_conn[order*order*order*i], rho_values);

    // Reset the flags
    for ( int kk = 0; kk < bsize; kk++ ){
      low[kk] = 0;
      high[kk] = 0;
      intermediate[kk] = 0;
    }

    // Scan through all the local node values for the filter
    if (index >= 0 && index < bsize){
      for ( int j = 0; j < order*order*order; j++ ){
        if (rho_values[bsize*j + index] < cutoff){
          low[0] = 1;
        }
        else if (rho_values[bsize*j + index] > 1.0 - cutoff){
          high[0] = 1;
        }
        else {
          intermediate[0] = 1;
        }
      }
    }
    else {
      for ( int j = 0; j < order*order*order; j++ ){
        for ( int kk = 0; kk < bsize; kk++ ){
          if (rho_values[bsize*j + kk] < cutoff){
            low[kk] = 1;
          }
          else if (rho_values[bsize*j + kk] > 1.0 - cutoff){
            high[kk] = 1;
          }
          else {
            intermediate[kk] = 1;
          }
        }
      }
    }

    if (index >= 0 && index < bsize){
      if ((low[0] && high[0]) || intermediate[0]){
        refine[i] = num_refine;
      }
    }
    else {
      for ( int kk = 0; kk < bsize; kk++ ){
        if ((low[kk] && high[kk]) || intermediate[kk]){
          refine[i] = num_refine;
        }
      }
    }
  }

  delete [] low;
  delete [] high;
  delete [] intermediate;

  // Refine the forest and balance it
  forest->refine(refine);
  delete [] refine;

  forest->balance(0);

  // Balance the forest and create the nodes
  forest->createNodes();

  // Get the connectivity
  const int *conn;
  int num_elements = 0, num_owned_nodes = 0;
  forest->getNodeConn(&conn, &num_elements, &num_owned_nodes);

  // Allocate the pointer array
  int *ptr = new int[ num_elements+1 ];
  for ( int i = 0; i <= num_elements; i++ ){
    ptr[i] = 8*i;
  }

  // Create/retrieve the dependent node information
  const int *dep_ptr, *dep_conn;
  const double *dep_weights;
  int num_dep_nodes = forest->getDepNodeConn(&dep_ptr, &dep_conn,
                                             &dep_weights);

  // Create the elements using the virtual call
  TACSElement **elements = new TACSElement*[ num_elements ];

  // Allocate the model class
  TACSElementModel *model = new TMRHexaHelmholtzModel(t);

  // Set the basis class
  TACSElementBasis *basis = new TACSLinearHexaBasis();

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

  // Filter octants
  TMROctant *filter_octs;
  TMROctantArray *filter_oct_array;
  filter->getOctants(&filter_oct_array);
  filter_oct_array->getArray(&filter_octs, NULL);

  // Get the array of the octants from the mesh
  TMROctant *octs;
  TMROctantArray *oct_array;
  forest->getOctants(&oct_array);
  oct_array->getArray(&octs, NULL);

  for ( int i = 0; i < num_elements; i++ ){
    // Set the node location
    TMROctant n = octs[i];
    const int32_t hf = 1 << (TMR_MAX_LEVEL - octs[i].level);

    for ( int kk = 0; kk < 2; kk++ ){
      for ( int jj = 0; jj < 2; jj++ ){
        for ( int ii = 0; ii < 2; ii++ ){
          // Only execute the loop if the node is positive
          int node = conn[8*i + ii + 2*jj + 4*kk];
          if (node >= 0){
            // Set the node location
            n.info = ii + 2*jj + 4*kk;

            // Find the enclosing element on the filter mesh
            const double knots[] = {-1.0, 1.0};
            TMROctant *filter_oct = filter->findEnclosing(2, knots, &n);

            if (filter_oct){
              // Get the values of the density at the nodes
              int filter_elem = filter_oct->tag;
              if (filter_elem >= 0 && filter_elem < num_filter_elements){
                rho->getValues(order*order*order,
                               &filter_conn[order*order*order*filter_elem],
                               rho_values);

                // Get the size of the filter element
                const int32_t h = 1 << (TMR_MAX_LEVEL - filter_oct->level);

                double pt[3];
                pt[0] = -1.0 + 2.0*(octs[i].x + hf*ii - filter_oct->x)/h;
                pt[1] = -1.0 + 2.0*(octs[i].y + hf*jj - filter_oct->y)/h;
                pt[2] = -1.0 + 2.0*(octs[i].z + hf*kk - filter_oct->z)/h;

                // Compute the interpolation
                filter->evalInterp(pt, N);

                // Evaluate the local values of rho based on the interpolation
                memset(rho_local, 0, bsize*sizeof(TacsScalar));
                for ( int k = 0; k < order*order*order; k++ ){
                  for ( int nk = 0; nk < bsize; nk++ ){
                    rho_local[nk] += rho_values[k*bsize + nk]*N[k];
                  }
                }

                // Check if this is an intermediate value of the design variables
                if (index >= 0 && index < bsize){
                  if (rho_local[index] >= cutoff &&
                      rho_local[index] <= 1.0 - cutoff){
                    int nbcs = 1;
                    int vars = 0;
                    TacsScalar value = 1.0;
                    assembler->addBCs(1, &node, nbcs, &vars, &value);
                  }
                }
                else {
                  // Loop over all the dimension
                  for ( int nk = 0; nk < bsize; nk++ ){
                    if (rho_local[nk] >= cutoff &&
                        rho_local[nk] <= 1.0 - cutoff){
                      int nbcs = 1;
                      int vars = 0;
                      TacsScalar value = 1.0;
                      assembler->addBCs(1, &node, nbcs, &vars, &value);
                      break;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Free the data
  delete [] N;
  delete [] rho_values;
  delete [] rho_local;

  assembler->initialize();

  // Get the ownership range of the forest
  const int *range;
  forest->getOwnedNodeRange(&range);

  // Set the node locations
  const int *nodes;
  int num_local_nodes = forest->getNodeNumbers(&nodes);

  // Get the points
  TMRPoint *Xp;
  forest->getPoints(&Xp);

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
  gmres->setTolerances(1e-12, 1e-30);

  TACSBVec *rhs = assembler->createVec();
  TACSBVec *dist = assembler->createVec();
  rhs->incref();
  dist->incref();

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

  dist->beginDistributeValues();
  dist->endDistributeValues();

  // Extract the actual distance using the transformation
  for ( int i = 0, index = 0; i < num_filter_elements; i++ ){
    double max_val = -1e20;

    while (index < num_elements &&
           filter_octs[i].contains(&octs[index])){
      for ( int ii = 0; ii < 8; ii++ ){
        // Only execute the loop if the node is positive
        int node = conn[8*index + ii];
        if (node >= 0){
          TacsScalar value;
          dist->getValues(1, &node, &value);
          if (value > max_val){
            max_val = value;
          }
        }
      }
      index++;
    }

    min_dist[i] = 1e20;
    if (max_val > 0.0){
      min_dist[i] = - t*log(max_val);
    }
  }

  // Write out the result to an f5 file
  if (filename){
    // Convert to a distance
    TacsScalar *dist_array;
    int size = dist->getArray(&dist_array);
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

  // Deallocate objects
  dist->decref();
  assembler->decref();
  forest->decref();
}
