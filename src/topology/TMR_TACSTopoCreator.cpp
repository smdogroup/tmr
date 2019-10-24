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

#include "TMR_TACSTopoCreator.h"

/*
  Compare integers for sorting
*/
static int compare_integers( const void *a, const void *b ){
  return (*(int*)a - *(int*)b);
}

/*
  Set up a creator class for the given filter problem
*/
TMROctTACSTopoCreator::TMROctTACSTopoCreator( TMRBoundaryConditions *_bcs,
                                              int _design_vars_per_node,
                                              TMROctForest *_filter ){
  int mpi_rank;
  MPI_Comm comm = filter->getMPIComm();
  MPI_Comm_rank(comm, &mpi_rank);

  // Create the nodes within the filter
  _filter->createNodes();

  initialize(_bcs, _design_vars_per_node, _filter);
}

/*
  Free the creator object
*/
TMROctTACSTopoCreator::~TMROctTACSTopoCreator(){
  if (filter_indices){ filter_indices->decref(); }
}

void TMROctTACSTopoCreator::computeWeights( const int mesh_order,
                                            const double *knots,
                                            TMROctant *node,
                                            TMROctant *oct,
                                            TMRIndexWeight *weights,
                                            double *tmp ){
  // Find the side length of the octant in the filter that contains
  // the element octant
  const int32_t h = 1 << (TMR_MAX_LEVEL - node->level);
  const int32_t hoct = 1 << (TMR_MAX_LEVEL - oct->level);

  // Compute the i, j, k location of the nod
  const int i = node->info % mesh_order;
  const int j = (node->info % (mesh_order*mesh_order))/mesh_order;
  const int k = node->info/(mesh_order*mesh_order);

  // Get the u/v/w values within the filter octant
  double pt[3];
  pt[0] = -1.0 + 2.0*((node->x % hoct) + 0.5*h*(1.0 + knots[i]))/hoct;
  pt[1] = -1.0 + 2.0*((node->y % hoct) + 0.5*h*(1.0 + knots[j]))/hoct;
  pt[2] = -1.0 + 2.0*((node->z % hoct) + 0.5*h*(1.0 + knots[k]))/hoct;

  // Get the Lagrange shape functions
  double *N = tmp;
  filter->evalInterp(pt, N);

  // Get the dependent node information for this mesh
  const int *dep_ptr, *dep_conn;
  const double *dep_weights;
  filter->getDepNodeConn(&dep_ptr, &dep_conn, &dep_weights);

  // Get the mesh order
  const int order = filter->getMeshOrder();

  // Get the connectivity
  const int *conn;
  filter->getNodeConn(&conn);
  const int *c = &conn[oct->tag*order*order*order];

  // Loop over the adjacent nodes within the filter
  int nweights = 0;
  for ( int kk = 0; kk < order; kk++ ){
    for ( int jj = 0; jj < order; jj++ ){
      for ( int ii = 0; ii < order; ii++ ){
        // Set the weights
        int offset = ii + jj*order + kk*order*order;
        double weight = N[offset];

        // Get the tag number
        if (c[offset] >= 0){
          weights[nweights].index = c[offset];
          weights[nweights].weight = weight;
          nweights++;
        }
        else {
          int node = -c[offset]-1;
          for ( int jp = dep_ptr[node]; jp < dep_ptr[node+1]; jp++ ){
            weights[nweights].index = dep_conn[jp];
            weights[nweights].weight = weight*dep_weights[jp];
            nweights++;
          }
        }
      }
    }
  }

  // Sort and sum the array of weights - there are only 8 nodes
  // per filter point at most
  nweights = TMRIndexWeight::uniqueSort(weights, nweights);
}

/*
  Create all of the elements for the topology optimization problem
*/
void TMROctTACSTopoCreator::createElements( int order,
                                            TMROctForest *forest,
                                            int num_elements,
                                            TACSElement **elements ){
  // Get the MPI communicator
  int mpi_rank, mpi_size;
  MPI_Comm comm = forest->getMPIComm();
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // Get the octants associated with the forest
  TMROctantArray *octants;
  forest->getOctants(&octants);

  // Get the array of octants from the forest
  int num_octs;
  TMROctant *octs;
  octants->getArray(&octs, &num_octs);

  // Create a queue for the external octants
  TMROctantQueue *queue = new TMROctantQueue();

  // The number of weights/element
  const int filter_order = filter->getMeshOrder();
  const int nweights = filter_order*filter_order*filter_order;

  // Allocate temp space
  double *tmp = new double[ nweights ];
  TMRIndexWeight *wtmp = new TMRIndexWeight[ nweights*filter_order*filter_order ];

  // Allocate the weights for all of the local elements
  TMRIndexWeight *weights = new TMRIndexWeight[ nweights*num_octs ];

  // Fake the information as if we have a third-order and we are
  // searching for the centeral node
  const int node_info = 13;
  const double node_knots[] = {-1.0, 0.0, 1.0};
  const int node_order = 3;

  for ( int i = 0; i < num_octs; i++ ){
    // Get the original octant from the forest
    TMROctant node = octs[i];
    node.info = node_info;

    // Find the enclosing central node
    int mpi_owner = 0;
    TMROctant *oct = filter->findEnclosing(node_order, node_knots,
                                           &node, &mpi_owner);

    if (!oct){
      // Push the octant to the external queue. We will handle these
      // cases seperately after a collective communication.
      node.tag = mpi_owner;
      queue->push(&node);
      weights[nweights*i].index = -1;
    }
    else {
      computeWeights(node_order, node_knots, &node,
                     oct, wtmp, tmp);
      memcpy(&weights[nweights*i], wtmp, nweights*sizeof(TMRIndexWeight));
    }
  }

  // Create a list of octants that are external
  TMROctantArray *nodes = queue->toArray();
  delete queue;

  // Distribute the nodes to the processors that own them
  int use_tags = 1;
  int *send_ptr, *recv_ptr;
  TMROctantArray *dist_nodes =
    filter->distributeOctants(nodes, use_tags, &send_ptr, &recv_ptr);
  delete nodes;

  // Get the external nodes that are local to this processor and
  // compute their weights
  int dist_size;
  TMROctant *dist_array;
  dist_nodes->getArray(&dist_array, &dist_size);

  // Create the distributed weights
  TMRIndexWeight *dist_weights = new TMRIndexWeight[ nweights*dist_size ];
  for ( int i = 0; i < dist_size; i++ ){
    TMROctant *oct = filter->findEnclosing(node_order, node_knots,
                                           &dist_array[i]);
    if (oct){
      computeWeights(node_order, node_knots, &dist_array[i],
                     oct, wtmp, tmp);
      memcpy(&dist_weights[nweights*i], wtmp, nweights*sizeof(TMRIndexWeight));
    }
    else {
      fprintf(stderr, "[%d] TMROctTACSTopoCreator: Node not found\n", mpi_rank);
    }
  }

  // Free the temporary space
  delete [] wtmp;
  delete [] tmp;

  // The distributed nodes are no longer required
  delete dist_nodes;

  // Compute the number of sends and recvs that were performed.
  int nsends = 0, nrecvs = 0;
  for ( int i = 0; i < mpi_size; i++ ){
    if (i != mpi_rank){
      if (send_ptr[i+1] - send_ptr[i] > 0){
        nsends++;
      }
      if (recv_ptr[i+1] - recv_ptr[i] > 0){
        nrecvs++;
      }
    }
  }

  // Now prepare to reverse the communication to distribute the
  // weights back to the processors that need them. First allocate
  // space for the requests
  MPI_Request *send_request = new MPI_Request[ nrecvs ];

  // Allocate space for the new weights
  TMRIndexWeight *new_weights =
    new TMRIndexWeight[ nweights*send_ptr[mpi_size] ];

  // Loop over all the ranks and send
  for ( int i = 0, j = 0; i < mpi_size; i++ ){
    if (i != mpi_rank &&
        recv_ptr[i+1] - recv_ptr[i] > 0){
      // Post the send to the destination
      int count = nweights*(recv_ptr[i+1] - recv_ptr[i]);
      MPI_Isend(&dist_weights[nweights*recv_ptr[i]], count,
                TMRIndexWeight_MPI_type,
                i, 0, comm, &send_request[j]);
      j++;
    }
  }

  // Loop over the recieve calls
  for ( int i = 0; i < mpi_size; i++ ){
    if (i != mpi_rank &&
        send_ptr[i+1] > send_ptr[i]){
      int count = nweights*(send_ptr[i+1] - send_ptr[i]);
      MPI_Recv(&new_weights[nweights*send_ptr[i]], count,
               TMRIndexWeight_MPI_type,
               i, 0, comm, MPI_STATUS_IGNORE);
    }
  }

  // Wait for any remaining sends to complete
  MPI_Waitall(nrecvs, send_request, MPI_STATUSES_IGNORE);

  // Now place the weights back into their original locations
  for ( int i = 0, j = 0; i < num_octs && j < send_ptr[mpi_size]; i++ ){
    if (weights[nweights*i].index == -1){
      memcpy(&weights[nweights*i], &new_weights[nweights*j],
             nweights*sizeof(TMRIndexWeight));
      j++;
    }
  }
  delete [] new_weights;
  delete [] send_request;
  delete [] send_ptr;
  delete [] recv_ptr;
  delete [] dist_weights;

  // The node numbers within the weights are global. Convert them to a
  // local node ordering and create a list of the external node
  // numbers referenced by the weights.

  // Get the node range for the filter design variables
  const int *filter_range;
  filter->getOwnedNodeRange(&filter_range);

  // The number of local nodes
  int num_filter_local = filter_range[mpi_rank+1] - filter_range[mpi_rank];

  // Get the external numbers from the filter itself
  const int *filter_ext;
  int num_filter_ext = filter->getNodeNumbers(&filter_ext);

  // Count up all the external nodes
  int num_ext = 0;
  int max_ext_nodes = nweights*num_octs + num_filter_ext;
  int *ext_nodes = new int[ max_ext_nodes ];

  // Add the external nodes from the filter
  for ( int i = 0; i < num_filter_ext; i++ ){
    int node = filter_ext[i];
    if (node >= 0 &&
        (node < filter_range[mpi_rank] ||
         node >= filter_range[mpi_rank+1])){
      ext_nodes[num_ext] = node;
      num_ext++;
    }
  }

  // Add the external nodes from the element-level connectivity
  for ( int i = 0; i < nweights*num_octs; i++ ){
    int node = weights[i].index;
    if (node < filter_range[mpi_rank] ||
        node >= filter_range[mpi_rank+1]){
      ext_nodes[num_ext] = node;
      num_ext++;
    }
  }

  // Sort the external array of nodes
  qsort(ext_nodes, num_ext, sizeof(int), compare_integers);

  // Remove duplicates from the array
  int len = 0;
  for ( int i = 0; i < num_ext; i++, len++ ){
    while ((i < num_ext-1) && (ext_nodes[i] == ext_nodes[i+1])){
      i++;
    }
    if (i != len){
      ext_nodes[len] = ext_nodes[i];
    }
  }

  // Truncate the array and delete the old array
  int num_ext_nodes = len;
  int *ext_node_nums = new int[ len ];
  memcpy(ext_node_nums, ext_nodes, len*sizeof(int));
  delete [] ext_nodes;

  // Set up the external filter indices for this filter.  The indices
  // objects steals the array for the external nodes.
  filter_indices = new TACSBVecIndices(&ext_node_nums, num_ext_nodes);
  filter_indices->incref();
  filter_indices->setUpInverse();

  // Scan through all of the weights and convert them to the local
  // ordering
  for ( int i = 0; i < nweights*num_octs; i++ ){
    int node = weights[i].index;
    if (node >= filter_range[mpi_rank] && node < filter_range[mpi_rank+1]){
      node = node - filter_range[mpi_rank];
    }
    else {
      node = num_filter_local + filter_indices->findIndex(node);
    }
    weights[i].index = node;
  }

  // Loop over the octants
  octants->getArray(&octs, &num_octs);
  for ( int i = 0; i < num_octs; i++ ){
    // Allocate the stiffness object
    elements[i] = createElement(order, &octs[i],
                                nweights, &weights[nweights*i]);
  }

  delete [] weights;
}

/*
  Set up a creator class for the given filter problem
*/
TMRQuadTACSTopoCreator::TMRQuadTACSTopoCreator( TMRBoundaryConditions *_bcs,
                                                int _design_vars_per_node,
                                                TMRQuadForest *_filter ):
TMRQuadTACSCreator(_bcs, _design_vars_per_node, _filter){
  int mpi_rank;
  MPI_Comm comm = filter->getMPIComm();
  MPI_Comm_rank(comm, &mpi_rank);

  // Create the nodes within the filter
  filter->createNodes();

  // Set the filter indices to NULL
  filter_indices = NULL;
}

/*
  Free the creator object
*/
TMRQuadTACSTopoCreator::~TMRQuadTACSTopoCreator(){
  if (filter_indices){ filter_indices->decref(); }
}

/*
  Compute the weights associated with the given quadrant
*/
void TMRQuadTACSTopoCreator::computeWeights( const int mesh_order,
                                             const double *knots,
                                             TMRQuadrant *node,
                                             TMRQuadrant *quad,
                                             TMRIndexWeight *weights,
                                             double *tmp, int sort ){
  // Find the side length of the octant in the filter that contains
  // the element octant
  const int32_t h = 1 << (TMR_MAX_LEVEL - node->level);
  const int32_t hquad = 1 << (TMR_MAX_LEVEL - quad->level);

  // Compute the i, j, k location of the node
  const int i = node->info % mesh_order;
  const int j = node->info/mesh_order;

  // Get the u/v values within the filter octant
  double pt[3];
  pt[0] = -1.0 + 2.0*((node->x % hquad) + 0.5*h*(1.0 + knots[i]))/hquad;
  pt[1] = -1.0 + 2.0*((node->y % hquad) + 0.5*h*(1.0 + knots[j]))/hquad;

  // Get the shape functions
  double *N = tmp;
  filter->evalInterp(pt, N);

  // Get the dependent node information for this mesh
  const int *dep_ptr, *dep_conn;
  const double *dep_weights;
  filter->getDepNodeConn(&dep_ptr, &dep_conn, &dep_weights);

  // Get the mesh order
  const int order = filter->getMeshOrder();

  // Get the connectivity
  const int *conn;
  filter->getNodeConn(&conn);
  const int *c = &conn[quad->tag*order*order];

  // Loop over the adjacent nodes within the filter
  int nweights = 0;
  for ( int jj = 0; jj < order; jj++ ){
    for ( int ii = 0; ii < order; ii++ ){
      // Set the weights
      int offset = ii + jj*order;
      double weight = N[offset];

      // Get the tag number
      if (c[offset] >= 0){
        weights[nweights].index = c[offset];
        weights[nweights].weight = weight;
        nweights++;
      }
      else {
        int node = -c[offset]-1;
        for ( int jp = dep_ptr[node]; jp < dep_ptr[node+1]; jp++ ){
          weights[nweights].index = dep_conn[jp];
          weights[nweights].weight = weight*dep_weights[jp];
          nweights++;
        }
      }
    }
  }
  if (sort){
    // Sort and sum the array of weights - there are only 8 nodes
    // per filter point at most
    nweights = TMRIndexWeight::uniqueSort(weights, nweights);
  }
}

/*
  Create all of the elements for the topology optimization problem
*/
void TMRQuadTACSTopoCreator::createElements( int order,
                                             TMRQuadForest *forest,
                                             int num_elements,
                                             TACSElement **elements ){
  // Get the MPI communicator
  int mpi_rank, mpi_size;
  MPI_Comm comm = forest->getMPIComm();
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // Get the quadrants associated with the forest
  TMRQuadrantArray *quadrants;
  forest->getQuadrants(&quadrants);

  // Get the array of quadrants from the forest
  int num_quads;
  TMRQuadrant *quads;
  quadrants->getArray(&quads, &num_quads);

  // Create a queue for the external octants
  TMRQuadrantQueue *queue = new TMRQuadrantQueue();

  // The number of weights/element
  const int filter_order = filter->getMeshOrder();
  const int nweights = filter_order*filter_order;

  // Allocate temp space
  double *tmp = new double[ nweights ];
  TMRIndexWeight *wtmp = new TMRIndexWeight[ nweights*filter_order*filter_order ];

  // Allocate the weights for all of the local elements
  TMRIndexWeight *weights = new TMRIndexWeight[ nweights*num_quads ];

  // Fake the information as if we have a third-order and we are
  // searching for the centeral node
  const int node_info = 4;
  const double node_knots[] = {-1.0, 0.0, 1.0};
  const int node_order = 3;

  for ( int i = 0; i < num_quads; i++ ){
    // Get the original quadrant from the forest
    TMRQuadrant node = quads[i];
    node.info = node_info;

    // Find the central node
    int mpi_owner = 0;
    TMRQuadrant *quad = filter->findEnclosing(node_order, node_knots,
                                              &node, &mpi_owner);

    if (!quad){
      // Push the quadrant to the external queue. We will handle these
      // cases seperately after a collective communication.
      node.tag = mpi_owner;
      queue->push(&node);
      weights[nweights*i].index = -1;
    }
    else {
      computeWeights(node_order, node_knots, &node,
                     quad, wtmp, tmp);
      memcpy(&weights[nweights*i], wtmp, nweights*sizeof(TMRIndexWeight));
    }
  }

  // Create a list of quadrants that are external
  TMRQuadrantArray *nodes = queue->toArray();
  delete queue;

  // Distribute the nodes to the processors that own them
  int use_tags = 1;
  int *send_ptr, *recv_ptr;
  TMRQuadrantArray *dist_nodes =
    filter->distributeQuadrants(nodes, use_tags, &send_ptr, &recv_ptr);
  delete nodes;

  // Get the external nodes that are local to this processor and
  // compute their weights
  int dist_size;
  TMRQuadrant *dist_array;
  dist_nodes->getArray(&dist_array, &dist_size);
  // Create the distributed weights
  TMRIndexWeight *dist_weights = new TMRIndexWeight[ nweights*dist_size ];
  for ( int i = 0; i < dist_size; i++ ){
    TMRQuadrant *quad = filter->findEnclosing(node_order, node_knots,
                                              &dist_array[i]);
    if (quad){
      computeWeights(node_order, node_knots, &dist_array[i],
                     quad, wtmp, tmp);
      memcpy(&dist_weights[nweights*i], wtmp, nweights*sizeof(TMRIndexWeight));
    }
  }

  // Free the tmporary space
  delete [] wtmp;
  delete [] tmp;

  // The distributed nodes are no longer required
  delete dist_nodes;

  // Compute the number of sends and recvs that were performed.
  int nsends = 0, nrecvs = 0;
  for ( int i = 0; i < mpi_size; i++ ){
    if (i != mpi_rank){
      if (send_ptr[i+1] - send_ptr[i] > 0){
        nsends++;
      }
      if (recv_ptr[i+1] - recv_ptr[i] > 0){
        nrecvs++;
      }
    }
  }
  // Now prepare to reverse the communication to distribute the
  // weights back to the processors that need them. First allocate
  // space for the requests
  MPI_Request *send_request = new MPI_Request[ nrecvs ];

  // Allocate space for the new weights
  TMRIndexWeight *new_weights =
    new TMRIndexWeight[ nweights*send_ptr[mpi_size] ];

  // Loop over all the ranks and send
  for ( int i = 0, j = 0; i < mpi_size; i++ ){
    if (i != mpi_rank &&
        recv_ptr[i+1] - recv_ptr[i] > 0){
      // Post the send to the destination
      int count = nweights*(recv_ptr[i+1] - recv_ptr[i]);
      MPI_Isend(&dist_weights[nweights*recv_ptr[i]], count,
                TMRIndexWeight_MPI_type,
                i, 0, comm, &send_request[j]);
      j++;
    }
  }

  // Loop over the recieve calls
  for ( int i = 0; i < mpi_size; i++ ){
    if (i != mpi_rank &&
        send_ptr[i+1] > send_ptr[i]){
      int count = nweights*(send_ptr[i+1] - send_ptr[i]);
      MPI_Recv(&new_weights[nweights*send_ptr[i]], count,
               TMRIndexWeight_MPI_type,
               i, 0, comm, MPI_STATUS_IGNORE);
    }
  }

  // Wait for any remaining sends to complete
  MPI_Waitall(nrecvs, send_request, MPI_STATUSES_IGNORE);

  // Now place the weights back into their original locations
  for ( int i = 0, j = 0; i < num_quads && j < send_ptr[mpi_size]; i++ ){
    if (weights[nweights*i].index == -1){
      memcpy(&weights[nweights*i], &new_weights[nweights*j],
             nweights*sizeof(TMRIndexWeight));
      j++;
    }
  }

  delete [] new_weights;
  delete [] send_request;
  delete [] send_ptr;
  delete [] recv_ptr;
  delete [] dist_weights;

  // The node numbers within the weights are global. Convert them into
  // a local node ordering and create a list of the external node
  // numbers referenced by the weights.

  // Get the node range for the filter design variables
  const int *filter_range;
  filter->getOwnedNodeRange(&filter_range);

  // The number of local nodes
  int num_filter_local = filter_range[mpi_rank+1] - filter_range[mpi_rank];
  // Get the external numbers from the filter itself
  const int *filter_ext;
  int num_filter_ext = filter->getNodeNumbers(&filter_ext);
  // Count up all the external nodes
  int num_ext = 0;
  int max_ext_nodes = nweights*num_quads + num_filter_ext;
  int *ext_nodes = new int[ max_ext_nodes ];

  // Add the external nodes from the filter
  for ( int i = 0; i < num_filter_ext; i++ ){
    int node = filter_ext[i];
    if (node >= 0 &&
        (node < filter_range[mpi_rank] ||
         node >= filter_range[mpi_rank+1])){
      ext_nodes[num_ext] = node;
      num_ext++;
    }
  }

  // Add the external nodes from the element-level connectivity
  for ( int i = 0; i < nweights*num_quads; i++ ){
    int node = weights[i].index;
    if (node < filter_range[mpi_rank] ||
        node >= filter_range[mpi_rank+1]){
      ext_nodes[num_ext] = node;
      num_ext++;
    }
  }

  // Sort the external array of nodes
  qsort(ext_nodes, num_ext, sizeof(int), compare_integers);

  // Remove duplicates from the array
  int len = 0;
  for ( int i = 0; i < num_ext; i++, len++ ){
    while ((i < num_ext-1) && (ext_nodes[i] == ext_nodes[i+1])){
      i++;
    }
    if (i != len){
      ext_nodes[len] = ext_nodes[i];
    }
  }

  // Truncate the array and delete the old array
  int num_ext_nodes = len;
  int *ext_node_nums = new int[ len ];
  memcpy(ext_node_nums, ext_nodes, len*sizeof(int));
  delete [] ext_nodes;

  // Set up the external filter indices for this filter.  The indices
  // objects steals the array for the external nodes.
  filter_indices = new TACSBVecIndices(&ext_node_nums, num_ext_nodes);
  filter_indices->incref();
  filter_indices->setUpInverse();

  // Scan through all of the weights and convert them to the local
  // ordering
  for ( int i = 0; i < nweights*num_quads; i++ ){
    int node = weights[i].index;
    if (node >= filter_range[mpi_rank] && node < filter_range[mpi_rank+1]){
      node = node - filter_range[mpi_rank];
    }
    else {
      node = num_filter_local + filter_indices->findIndex(node);
    }
    weights[i].index = node;
  }

  // Loop over the octants
  quadrants->getArray(&quads, &num_quads);
  for ( int i = 0; i < num_quads; i++ ){
    // Allocate the stiffness object
    elements[i] = createElement(order, &quads[i],
                                nweights, &weights[nweights*i]);
  }
  // Free the weights
  delete [] weights;
}

/*
  Set up a creator class for the given forest problem
*/
TMROctConformTACSTopoCreator::TMROctConformTACSTopoCreator( TMRBoundaryConditions *_bcs,
                                                            int _design_vars_per_node,
                                                            TMROctForest *_forest,
                                                            int order,
                                                            TMRInterpolationType interp_type ){
  // Use the forest as the filter in these cases
  TMROctForest *_filter = NULL;
  if (order < 0 ||
      (order == _forest->getMeshOrder() &&
       interp_type == _forest->getInterpType())){
    _filter = _forest;
  }
  else {
    _filter = _forest->duplicate();
    order = _filter->getMeshOrder()-1;
    if (order < 2){
      order = 2;
    }
    _filter->setMeshOrder(order, interp_type);
  }
  _filter->incref();

  // Create the nodes within the filter
  _filter->createNodes();

  initialize(_bcs, _design_vars_per_node, _filter);
}

/*
  Free the creator object
*/
TMROctConformTACSTopoCreator::~TMROctConformTACSTopoCreator(){}

/*
  Create all of the elements for the topology optimization problem
*/
void TMROctConformTACSTopoCreator::createElements( int order,
                                                   TMROctForest *forest,
                                                   int num_elements,
                                                   TACSElement **elements ){
  // Get the MPI communicator
  int mpi_rank, mpi_size;
  MPI_Comm comm = forest->getMPIComm();
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // Get the quadrants associated with the forest
  TMROctantArray *octants;
  forest->getOctants(&octants);

  // Get the array of quadrants from the forest
  int num_octs;
  TMROctant *octs;
  octants->getArray(&octs, &num_octs);

  // The number of weights/element
  const int filter_order = filter->getMeshOrder();
  const int nweights = filter_order*filter_order*filter_order;

  // Loop over the nodes and convert to the local numbering scheme
  const int *conn;
  filter->getNodeConn(&conn);

  // Loop over the octants
  octants->getArray(&octs, &num_octs);
  for ( int i = 0; i < num_octs; i++ ){
    // Allocate the stiffness object
    elements[i] = createElement(order, &octs[i],
                                nweights, &conn[nweights*i],
                                filter);
  }
}

/*
  Set up a creator class for the given forest problem
*/
TMRQuadConformTACSTopoCreator::TMRQuadConformTACSTopoCreator( TMRBoundaryConditions *_bcs,
                                                              int _design_vars_per_node,
                                                              TMRQuadForest *_forest,
                                                              int order,
                                                              TMRInterpolationType interp_type ){
  // Use the forest as the filter in these cases
  TMRQuadForest *_filter = NULL;
  if (order < 0 ||
      (order == _forest->getMeshOrder() &&
       interp_type == _forest->getInterpType())){
    _filter = _forest;
  }
  else {
    _filter = _forest->duplicate();
    order = _filter->getMeshOrder()-1;
    if (order < 2){
      order = 2;
    }
    _filter->setMeshOrder(order, interp_type);
  }
  _filter->incref();

  // Create the nodes within the filter
  _filter->createNodes();

  initialize(_bcs, _design_vars_per_node, _filter);
}

/*
  Free the creator object
*/
TMRQuadConformTACSTopoCreator::~TMRQuadConformTACSTopoCreator(){}

/*
  Create all of the elements for the topology optimization problem
*/
void TMRQuadConformTACSTopoCreator::createElements( int order,
                                                    TMRQuadForest *forest,
                                                    int num_elements,
                                                    TACSElement **elements ){

  // Get the MPI communicator
  int mpi_rank, mpi_size;
  MPI_Comm comm = forest->getMPIComm();
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // Get the quadrants associated with the forest
  TMRQuadrantArray *quadrants;
  forest->getQuadrants(&quadrants);

  // Get the array of quadrants from the forest
  int num_quads;
  TMRQuadrant *quads;
  quadrants->getArray(&quads, &num_quads);

  // The number of weights/element
  const int filter_order = filter->getMeshOrder();
  const int nweights = filter_order*filter_order;

  // Loop over the nodes and convert to the local numbering scheme
  const int *conn;
  filter->getNodeConn(&conn);

  // Loop over the octants
  quadrants->getArray(&quads, &num_quads);
  for ( int i = 0; i < num_quads; i++ ){
    // Allocate the stiffness object
    elements[i] = createElement(order, &quads[i],
                                nweights, &conn[nweights*i],
                                filter);
  }
}
