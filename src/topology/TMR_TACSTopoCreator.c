#include "TMR_TACSTopoCreator.h"
#include "TMROctStiffness.h"
#include "FElibrary.h"
#include "Solid.h"

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
                                              TMROctForest *_filter,
                                              const char *_shell_attr,
                                              SolidShellWrapper *_shell ):
TMROctTACSCreator(_bcs){
  // Reference the filter
  filter = _filter;
  filter->incref();
  
  int mpi_rank;
  MPI_Comm comm = filter->getMPIComm();
  MPI_Comm_rank(comm, &mpi_rank);

  // Create the nodes within the filter
  filter->createNodes(2);

  // Create the dependent node connectivity
  filter->createDepNodeConn();

  // Get the node range for the filter design variables
  const int *filter_range;
  filter->getOwnedNodeRange(&filter_range);

  // Set up the variable map for the design variable numbers
  int num_filter_local = filter_range[mpi_rank+1] - filter_range[mpi_rank];
  filter_map = new TACSVarMap(comm, num_filter_local);
  filter_map->incref();

  // Set the filter indices to NULL
  filter_indices = NULL;

  // Set shell attributes if they exist 
  shell_attr = NULL;
  shell = NULL;
  if (_shell_attr && _shell){
    shell_attr = new char[ strlen(_shell_attr)+1 ];
    strcpy(shell_attr, _shell_attr);
    shell = _shell;
    shell->incref();
  }
}

/*
  Free the creator object
*/
TMROctTACSTopoCreator::~TMROctTACSTopoCreator(){
  filter->decref();
  filter_map->decref();
  if (filter_indices){ filter_indices->decref(); }
  if (shell_attr){ delete [] shell_attr; }
}

// Get the underlying information about the
void TMROctTACSTopoCreator::getFilter( TMROctForest **_filter ){
  *_filter = filter;
}

void TMROctTACSTopoCreator::getMap( TACSVarMap **_map ){
  *_map = filter_map;
}

void TMROctTACSTopoCreator::getIndices( TACSBVecIndices **_indices ){
  *_indices = filter_indices;
}

/*
  If order != 2, this is not going to work.....
*/
void TMROctTACSTopoCreator::createConnectivity( int order,
                                                TMROctForest *forest,
                                                int **_conn, int **_ptr,
                                                int *_num_elements ){
  // Create the mesh
  int *elem_conn, num_octs;
  forest->createMeshConn(&elem_conn, &num_octs);

  // Set the default attributes
  int num_elements = num_octs;

  // Get the octants from the top and bottom surface
  if (shell){
    TMROctantArray *nodes;
    forest->getNodes(&nodes);

    TMROctantArray *shell_octs = forest->getOctsWithAttribute(shell_attr);
    
    // Get the top and bottom surface octants
    int nshell;
    TMROctant *octs;
    shell_octs->getArray(&octs, &nshell);

    num_elements = num_octs + nshell;
    int *conn = new int[ 8*num_elements ];
    memcpy(conn, elem_conn, 8*num_octs*sizeof(int));
    delete [] elem_conn;

    int index = 8*num_octs;

    // Add the connectivity from the nodes
    for ( int i = 0; i < nshell; i++ ){
      const int32_t h = 1 << (TMR_MAX_LEVEL - octs[i].level);
      for ( int jj = 0; jj < 2; jj++ ){
        for ( int ii = 0; ii < 2; ii++ ){
          int face_index = octs[i].tag;

          TMROctant node;
          node.block = octs[i].block;

          if (face_index < 2){
            node.x = octs[i].x + (face_index % 2)*h;
            node.y = octs[i].y + ii*h;
            node.z = octs[i].z + jj*h;
          }
          else if (face_index < 4){
            node.x = octs[i].x + ii*h;
            node.y = octs[i].y + (face_index % 2)*h;
            node.z = octs[i].z + jj*h;
          }
          else {
            node.x = octs[i].x + ii*h;
            node.y = octs[i].y + jj*h;
            node.z = octs[i].z + (face_index % 2)*h;
          }

          // Transform the node to the global ordering
          forest->transformNode(&node);
        
          const int use_node_search = 1;
          TMROctant *t = nodes->contains(&node, use_node_search);
          conn[index] = t->tag;  index++;
          conn[index] = t->tag+1; index++;
        }
      }
    }

    // Free the bottom octants
    delete shell_octs;
    elem_conn = conn;    
  }

  // Set the element ptr
  int *ptr = new int[ num_elements+1 ];
  for ( int i = 0; i < num_elements+1; i++ ){
    ptr[i] = order*order*order*i;
  }

  *_conn = elem_conn;
  *_ptr = ptr;
  *_num_elements = num_elements;
}

void TMROctTACSTopoCreator::computeWeights( TMROctant *oct,
                                            TMROctant *node,
                                            TMRIndexWeight *welem ){
  // Find the side length of the octant in the filter that contains
  // the element octant
  const int32_t hoct = 1 << (TMR_MAX_LEVEL - oct->level);
    
  // Get the u/v/w values within the filter octant
  double pt[3];
  pt[0] = -1.0 + 2.0*(node->x - oct->x)/hoct;
  pt[1] = -1.0 + 2.0*(node->y - oct->y)/hoct;
  pt[2] = -1.0 + 2.0*(node->z - oct->z)/hoct;
        
  // Get the Lagrange shape functions
  double N[8];
  FElibrary::triLagrangeSF(N, pt, 2);
  
  // Get the dependent node information for this mesh
  const int *dep_ptr, *dep_conn;
  const double *dep_weights;
  filter->getDepNodeConn(&dep_ptr, &dep_conn,
                         &dep_weights);

  // Get the octant array for the nodes
  TMROctantArray *filter_nodes;
  filter->getNodes(&filter_nodes);
  
  // Store the weights for each node
  int nweights = 0;
  TMRIndexWeight weights[32];
    
  // Loop over the adjacent nodes within the filter
  for ( int kk = 0; kk < 2; kk++ ){
    for ( int jj = 0; jj < 2; jj++ ){
      for ( int ii = 0; ii < 2; ii++ ){
        // Set the weights
        double wval = N[ii + 2*jj + 4*kk];
            
        // Compute the location of the node
        TMROctant p;
        p.block = oct->block;
        p.x = oct->x + hoct*ii;
        p.y = oct->y + hoct*jj;
        p.z = oct->z + hoct*kk;

        // Transform the node into the appropriate location
        filter->transformNode(&p);
            
        // Search for the node p within the octant
        const int use_node_search = 1;
        TMROctant *t = filter_nodes->contains(&p, use_node_search);

        // Get the node number
        int node = t->tag;
        if (node >= 0){
          weights[nweights].index = node;
          weights[nweights].weight = wval;
          nweights++;
        }
        else {
          node = -node-1;
          for ( int jp = dep_ptr[node]; jp < dep_ptr[node+1]; jp++ ){
            int dep_node = dep_conn[jp];
            weights[nweights].index = dep_node;
            weights[nweights].weight = wval*dep_weights[jp];
            nweights++;
          }
        }
      }
    }
  }
  
  // Sort and sum the array of weights - there are only 8 nodes
  // per filter point at most
  nweights = TMRIndexWeight::uniqueSort(weights, nweights);
  memcpy(welem, weights, nweights*sizeof(TMRIndexWeight));
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
  const int nweights = 8;

  // Allocate the weights for all of the local elements 
  TMRIndexWeight *weights = new TMRIndexWeight[ nweights*num_octs ];

  for ( int i = 0; i < num_octs; i++ ){
    // Get the original octant from the forest    
    TMROctant node = octs[i];

    // Compute the half the edge length of the octant
    const int32_t h_half = 1 << (TMR_MAX_LEVEL - (octs[i].level + 1));

    // Place the octant node at the element mid-point location
    node.x += h_half;
    node.y += h_half;
    node.z += h_half;
    
    // No need to transform the octant since this will always lie
    // strictly in the interior of the owner octree. Find the
    // enclosing octant within the filter. If it does not exist on
    // this processor, add it to the external queue.
    TMROctant *oct = filter->findEnclosing(&node);

    if (!oct){
      // Push the octant to the external queue. We will handle these
      // cases seperately after a collective communication.
      queue->push(&node);
      weights[nweights*i].index = -1;
    }
    else {
      computeWeights(oct, &node, &weights[nweights*i]);
    }
  }

  // Create a list of octants that are external
  TMROctantArray *nodes = queue->toArray();
  delete queue;
  
  // Distribute the nodes to the processors that own them
  int use_tags = 0;
  int *send_ptr, *recv_ptr;
  TMROctantArray *dist_nodes = filter->distributeOctants(nodes, use_tags,
                                                         &send_ptr, &recv_ptr);
  delete nodes;

  // Get the external nodes that are local to this processor and
  // compute their weights
  int dist_size;
  TMROctant *dist_array;
  dist_nodes->getArray(&dist_array, &dist_size);

  // Create the distributed weights
  TMRIndexWeight *dist_weights = new TMRIndexWeight[ nweights*dist_size ];
  for ( int i = 0; i < dist_size; i++ ){
    TMROctant *oct = filter->findEnclosing(&dist_array[i]);
    if (oct){
      computeWeights(oct, &dist_array[i], &dist_weights[nweights*i]);
    }
  }

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

  // The node numbers within the weights are global. Convert them into
  // a local node ordering and create a list of the external node
  // numbers referenced by the weights.

  // Get the node range for the filter design variables
  const int *filter_range;
  filter->getOwnedNodeRange(&filter_range);

  // The number of local nodes
  int num_filter_local = filter_range[mpi_rank+1] - filter_range[mpi_rank];

  // Get the dependent node information for this mesh
  const int *dep_ptr, *dep_conn;
  const double *dep_weights;
  int num_dep_nodes = filter->getDepNodeConn(&dep_ptr, &dep_conn,
                                             &dep_weights);

  // Get the external numbers from the filter itself
  int *filter_ext;
  int num_filter_ext = filter->getExtNodeNums(&filter_ext);

  // Count up all the external nodes
  int num_ext = 0;
  int max_ext_nodes = nweights*num_octs + 
    dep_ptr[num_dep_nodes] + num_filter_ext;
  int *ext_nodes = new int[ max_ext_nodes ];

  // Add the external nodes from the filter
  for ( int i = 0; i < num_filter_ext; i++ ){
    ext_nodes[num_ext] = filter_ext[i];
    num_ext++;
  }
  delete [] filter_ext;

  // Add the external nodes from the dependent connectivity
  for ( int i = 0; i < dep_ptr[num_dep_nodes]; i++ ){
    int node = dep_conn[i];
    if (node < filter_range[mpi_rank] || 
        node >= filter_range[mpi_rank+1]){
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
                                &weights[nweights*i], nweights);
  }

  for ( int i = num_octs; i < num_elements; i++ ){
    elements[i] = shell;
  }

  delete [] weights;
}
