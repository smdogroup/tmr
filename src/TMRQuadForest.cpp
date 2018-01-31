#include "TMRQuadForest.h"

/*
  Copyright (c) 2016 Graeme Kennedy. All rights reserved. 
*/

/*
  Face to edge node connectivity
*/
const int face_to_edge_nodes[][2] = {{0, 2},
                                     {1, 3},
                                     {0, 1},
                                     {2, 3}};

/*
  Compare integers for sorting
*/
static int compare_integers( const void *a, const void *b ){
  return (*(int*)a - *(int*)b);
}

/*
  Compare tags for sorting
*/
static int compare_quadrant_tags( const void *a, const void *b ){
  const TMRQuadrant *A = static_cast<const TMRQuadrant*>(a);
  const TMRQuadrant *B = static_cast<const TMRQuadrant*>(b);
  return A->tag - B->tag;
}

/*
  Create the TMRQuadForest object
*/
TMRQuadForest::TMRQuadForest( MPI_Comm _comm, int _mesh_order,
                              TMRInterpolationType interp_type ){
  // Initialize the TMR-specific MPI data types
  if (!TMRIsInitialized()){
    TMRInitialize();
  }

  // Set the MPI communicator
  comm = _comm;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // Check that the order falls within allowable bounds
  mesh_order = _mesh_order;
  if (mesh_order < 2){
    mesh_order = 2;
  }

  // Allocate the interpolation knots and set the knot locations
  interp_knots = new double[ mesh_order ];
  if (interp_type == TMR_GAUSS_LOBATTO_POINTS){
    interp_knots[0] = 0.0;
    interp_knots[mesh_order-1] = 1.0;
    for ( int i = 1; i < mesh_order-1; i++ ){
      interp_knots[i] = 0.5*(1.0 - cos(M_PI*i/(mesh_order-1)));
    }
  }
  else {
    // Uniform mesh spacing
    interp_knots[0] = 0.0;
    interp_knots[mesh_order-1] = 1.0;
    for ( int i = 1; i < mesh_order-1; i++ ){
      interp_knots[i] = 1.0*i/(mesh_order-1);
    }
  }

  // Set the topology object to NULL
  topo = NULL;

  // Set the range of nodes
  node_range = NULL;

  // Zero out the nodes/edges/faces and all data
  num_nodes = 0;
  num_edges = 0;
  num_faces = 0;

  // Set all the unallocated pointers to NULL
  face_conn = NULL;
  face_edge_conn = NULL;
  node_face_conn = NULL;
  node_face_ptr = NULL;
  edge_face_conn = NULL;
  edge_face_ptr = NULL;
  edge_face_owners = NULL;
  node_face_owners = NULL;

  // Null the quadrant owners/quadrant list
  owners = NULL;
  quadrants = NULL;
  adjacent = NULL;
  X = NULL;

  // Set data for the number of elements/nodes/dependents
  node_range = NULL;
  num_owned_nodes = 0;
  num_dep_nodes = 0;
  dep_ptr = NULL;
  dep_conn = NULL;
  dep_weights = NULL;
}

/*
  Free the data allocated by the TMRQuadForest object
*/
TMRQuadForest::~TMRQuadForest(){
  // Free the topology object associated with this mesh (if any)
  if (topo){ topo->decref(); }

  freeData();
}

/*
  Free data and prepare for it to be reallocated
*/
void TMRQuadForest::freeData(){
  // Free the connectivity
  if (face_conn){ delete [] face_conn; }
  if (face_edge_conn){ delete [] face_edge_conn; }
  if (node_face_ptr){ delete [] node_face_ptr; }
  if (node_face_conn){ delete [] node_face_conn; }
  if (edge_face_ptr){ delete [] edge_face_ptr; }
  if (edge_face_conn){ delete [] edge_face_conn; }
  
  // Free the ownership data
  if (node_face_owners){ delete [] node_face_owners; }
  if (edge_face_owners){ delete [] edge_face_owners; }

  // Free the quadrants/adjacency
  if (owners){ delete [] owners; }
  if (quadrants){ delete quadrants; }
  if (adjacent){ delete adjacent; }
  if (X){ delete [] X; }

  if (node_range){ delete [] node_range; }
  if (dep_ptr){ delete [] dep_ptr; }
  if (dep_conn){ delete [] dep_conn; }
  if (dep_weights){ delete [] dep_weights; }

  // Set the range of nodes
  node_range = NULL;

  // Zero out the nodes/edges/faces and all data
  num_nodes = 0;
  num_edges = 0;
  num_faces = 0;

  // Set all the unallocated pointers to NULL
  face_conn = NULL;
  face_edge_conn = NULL;
  node_face_ptr = NULL;
  node_face_conn = NULL;
  edge_face_ptr = NULL;
  edge_face_conn = NULL;
  edge_face_owners = NULL;
  node_face_owners = NULL;

  // Null the quadrant owners/quadrant list
  owners = NULL;
  quadrants = NULL;
  adjacent = NULL;
  X = NULL;

  // Set data for the number of elements/nodes/dependents
  node_range = NULL;
  num_owned_nodes = 0;
  num_dep_nodes = 0;
  dep_ptr = NULL;
  dep_conn = NULL;
  dep_weights = NULL;
}

/*
  Copy the connectivity data, but not the quadrants/nodes
*/
void TMRQuadForest::copyData( TMRQuadForest *copy ){
  // Copy over the connectivity data
  copy->num_nodes = num_nodes;
  copy->num_edges = num_edges;
  copy->num_faces = num_faces;

  // Allocate/copy the face connectivities
  copy->face_conn = new int[ 4*num_faces ];
  copy->face_edge_conn = new int[ 4*num_faces ];
  memcpy(copy->face_conn, face_conn, 4*num_faces*sizeof(int));
  memcpy(copy->face_edge_conn, face_edge_conn, 
         4*num_faces*sizeof(int));
    
  // Allocate/copy the inverse relationships
  copy->node_face_ptr = new int[ num_nodes+1 ];
  copy->node_face_conn = new int[ node_face_ptr[num_nodes] ];
  memcpy(copy->node_face_ptr, node_face_ptr, 
         (num_nodes+1)*sizeof(int));
  memcpy(copy->node_face_conn, node_face_conn, 
         node_face_ptr[num_nodes]*sizeof(int));

  copy->edge_face_ptr = new int[ num_edges+1 ];
  copy->edge_face_conn = new int[ edge_face_ptr[num_edges] ];
  memcpy(copy->edge_face_ptr, edge_face_ptr, 
         (num_edges+1)*sizeof(int));
  memcpy(copy->edge_face_conn, edge_face_conn, 
         edge_face_ptr[num_edges]*sizeof(int));

  // Copy the ownership information
  copy->edge_face_owners = new int[ num_edges ];
  copy->node_face_owners = new int[ num_nodes ];
  memcpy(copy->edge_face_owners, edge_face_owners,
         num_edges*sizeof(int));
  memcpy(copy->node_face_owners, node_face_owners,
         num_nodes*sizeof(int));

  // Copy over the topology object
  copy->topo = topo;
  if (copy->topo){
    copy->topo->incref();
  }
}

/*
  Set the mesh topology - this has the effect of resetting the
  data and altering the topology of the mesh.
*/
void TMRQuadForest::setTopology( TMRTopology *_topo ){
  if (_topo){
    // Incref the topology object
    _topo->incref();
    if (topo){ topo->decref(); }
    topo = _topo;

    // Compute the topology and set it internally
    int _num_faces, _num_edges, _num_nodes;
    const int *_face_conn, *_face_edge_conn;
    topo->getConnectivity(&_num_nodes, &_num_edges, &_num_faces,
                          &_face_conn, &_face_edge_conn);

    // Set the full connectivity
    setFullConnectivity(_num_nodes, _num_edges, _num_faces,
                        _face_conn, _face_edge_conn);
  }
}

/*
  Retrieve the topology object from the TMRQuadForest 
  
  Note that this may return NULL if no topology is defined.
*/
TMRTopology *TMRQuadForest::getTopology(){
  return topo;
}

/*
  Set the connectivity of the faces

  This call is collective on all processors. Every processor must make
  a call with the same connectivity information, otherwise the
  inter-quadtree information will be inconsistent.  This code sets the
  face connectivity and generates the following additional data that
  are required:

  1. Face to node connectivity (input)
  2. Node to face connectivity (required for corner balancing)
  3. Unique edge ordering
  4. Face to edge and edge to face connectivity

  This information is required for creating quadtree forests on the
  unstructured mesh.
*/
void TMRQuadForest::setConnectivity( int _num_nodes,
                                     const int *_face_conn,
                                     int _num_faces ){
  // Free any data if it has already been allocated. 
  freeData();

  // Copy over the data locally
  num_nodes = _num_nodes;
  num_edges = 0;
  num_faces = _num_faces;
  
  // Copy over the face connectivity
  face_conn = new int[ 4*num_faces ];
  memcpy(face_conn, _face_conn, 4*num_faces*sizeof(int));

  // Compute the node to face information
  computeNodesToFaces();

  // Compute the edge connectivity from the face data
  computeEdgesFromNodes();
  computeEdgesToFaces();

  // Compute the face owners based on the node, edge and face data
  computeFaceOwners();
}

/*
  Set the full connectivity, specifying the node, edge and face
  numbers independently.
*/
void TMRQuadForest::setFullConnectivity( int _num_nodes, 
                                         int _num_edges,
                                         int _num_faces, 
                                         const int *_face_conn,
                                         const int *_face_edge_conn ){
  // Free any data allocated for other connectivities
  freeData();

  // Copy over the number of geometric entities
  num_nodes = _num_nodes;
  num_edges = _num_edges;
  num_faces = _num_faces;

  // Copy over the face connectivity
  face_conn = new int[ 4*num_faces ];
  memcpy(face_conn, _face_conn, 4*num_faces*sizeof(int));

  // Compute the node to face information
  computeNodesToFaces();

  // Copy over the edge information
  face_edge_conn = new int[ 4*num_faces ];
  memcpy(face_edge_conn, _face_edge_conn, 4*num_faces*sizeof(int));
  
  // Compute the edge to face information
  computeEdgesToFaces();

  // Compute the face owners based on the node, edge and face data
  computeFaceOwners();
}

/*
  Given the face to node connectivity, compute the node to face
  connectivity information
*/
void TMRQuadForest::computeNodesToFaces(){
  // Create the data structure for the node to face connectivity
  node_face_ptr = new int[ num_nodes+1 ];
  memset(node_face_ptr, 0, (num_nodes+1)*sizeof(int));

  // Count the number of times each node is referred to
  for ( int i = 0; i < 4*num_faces; i++ ){
    node_face_ptr[face_conn[i]+1]++;
  }

  // Adjust the counter array so that it points into the full array
  for ( int i = 1; i < num_nodes+1; i++ ){
    node_face_ptr[i] += node_face_ptr[i-1];
  }

  // Allocate the full node to face pointer array
  node_face_conn = new int[ node_face_ptr[num_nodes] ];
  for ( int i = 0; i < num_faces; i++ ){
    for ( int j = 0; j < 4; j++ ){
      int node = face_conn[4*i + j];
      node_face_conn[node_face_ptr[node]] = i;
      node_face_ptr[node]++;
    }
  }

  // Reset the pointer array so that it contains the correct offsets
  for ( int i = num_nodes; i >= 1; i-- ){
    node_face_ptr[i] = node_face_ptr[i-1];
  }
  node_face_ptr[0] = 0;

  // Loop over all the faces and reset node->face connectivity
  // to store both the adjacent block and the corresponding
  // node index into that array
  for ( int node = 0; node < num_nodes; node++ ){
    for ( int ip = node_face_ptr[node];
          ip < node_face_ptr[node+1]; ip++ ){
      int adj = node_face_conn[ip];
      int adj_index = 0;
      for ( ; adj_index < 4; adj_index++ ){
        if (face_conn[4*adj + adj_index] == node){
          break;
        }
      }

      node_face_conn[ip] = 4*adj + adj_index;
    }
  }
}

/*
  Based on the face to node connectivity information alone, compute a
  unique set of edges with associated edge numbers
*/
void TMRQuadForest::computeEdgesToFaces(){
  // Create the data structure for the edge to face connectivity
  edge_face_ptr = new int[ num_edges+1 ];
  memset(edge_face_ptr, 0, (num_edges+1)*sizeof(int));

  // Count the number of times each edge is referred to
  for ( int i = 0; i < 4*num_faces; i++ ){
    edge_face_ptr[face_edge_conn[i]+1]++;
  }

  // Adjust the counter array so that it points into the full array
  for ( int i = 1; i < num_edges+1; i++ ){
    edge_face_ptr[i] += edge_face_ptr[i-1];
  }

  // Allocate the full node to face pointer array
  edge_face_conn = new int[ edge_face_ptr[num_edges] ];
  for ( int face = 0; face < num_faces; face++ ){
    for ( int j = 0; j < 4; j++ ){
      int e = face_edge_conn[4*face + j];
      edge_face_conn[edge_face_ptr[e]] = face;
      edge_face_ptr[e]++;
    }
  }

  // Reset the pointer array so that it contains the correct offsets
  for ( int i = num_edges; i >= 1; i-- ){
    edge_face_ptr[i] = edge_face_ptr[i-1];
  }
  edge_face_ptr[0] = 0;

  // Loop over all edges and determine their relative orientation
  for ( int edge = 0; edge < num_edges; edge++ ){
    int face_owner = num_faces;
    int owner_index = 0;

    // Scan through the faces pointing to this edge to determine
    // the face owner - the face with the lowest index
    for ( int ip = edge_face_ptr[edge]; ip < edge_face_ptr[edge+1]; ip++ ){
      int face = edge_face_conn[ip];
      if (face < face_owner){
        face_owner = face;

        // Find the new owner index
        owner_index = 0;
        for ( int j = 0; j < 4; j++, owner_index++ ){
          if (face_edge_conn[4*face + j] == edge){
            break;
          }
        }
      }
    }

    // Retrieve the first and second node numbers
    int n1 = face_conn[4*face_owner + face_to_edge_nodes[owner_index][0]];
    int n2 = face_conn[4*face_owner + face_to_edge_nodes[owner_index][1]];

    // Now determine the local edge index on each face and adjust
    // connectivity data
    for ( int ip = edge_face_ptr[edge]; ip < edge_face_ptr[edge+1]; ip++ ){
      // Find the face index
      int face = edge_face_conn[ip];
      
      for ( int edge_index = 0; edge_index < 4; edge_index++ ){
        int nn1 = face_conn[4*face + face_to_edge_nodes[edge_index][0]];
        int nn2 = face_conn[4*face + face_to_edge_nodes[edge_index][1]];
        
        // Check if the edges now match up
        if ((n1 == nn1 && n2 == nn2) ||
            (n1 == nn2 && n2 == nn1)){
          edge_face_conn[ip] = 4*face + edge_index;
          break;
        }
      }
    }
  }
}

/*
  Compute the edges from the nodes
*/
void TMRQuadForest::computeEdgesFromNodes(){
  // Now establish a unique ordering of the edges along each face
  face_edge_conn = new int[ 4*num_faces ];
  for ( int i = 0; i < 4*num_faces; i++ ){
    face_edge_conn[i] = -1;
  }

  // Keep track of the edge numbers
  int edge = 0;
  for ( int i = 0; i < num_faces; i++ ){
    // Loop over each edge on this face
    for ( int j = 0; j < 4; j++ ){
      if (face_edge_conn[4*i + j] < 0){
        int n1 = face_conn[4*i + face_to_edge_nodes[j][0]];
        int n2 = face_conn[4*i + face_to_edge_nodes[j][1]];

        // Keep track of the number of edges found
        const int max_nedges = 20;
        int edge_index[max_nedges];
        int nedges = 1;
        edge_index[0] = 4*i + j;

        // Set the edge number - if any is found
        int edge_num = -1;

        // Scan through the faces that share the same
        // node and check if any of the edges are also
        // shared 
        for ( int ip = node_face_ptr[n1];
              ip < node_face_ptr[n1+1]; ip++ ){
          int ii = node_face_conn[ip]/4;
          
          // Loop over each edge in the new face
          for ( int jj = 0; jj < 4; jj++ ){
            int nn1 = face_conn[4*ii + face_to_edge_nodes[jj][0]];
            int nn2 = face_conn[4*ii + face_to_edge_nodes[jj][1]];

            // Check if the face matches
            if ((n1 == nn1 && n2 == nn2) ||
                (n1 == nn2 && n2 == nn1)){
              if (face_edge_conn[4*ii + jj] >= 0){
                // If this edge has been ordered, copy over
                // the edge number
                edge_num = face_edge_conn[4*ii + jj];
              }
              else if (nedges < max_nedges){
                // This edge has not yet been ordered, add it
                // to the unordered list if there is still room
                // if not, we will detect and order it during
                // a future iteration
                edge_index[nedges] = 4*ii + jj;
                nedges++;
              }
            }
          }
        }

        // If this edge does not have an edge number, assign
        // a new one to the list
        if (edge_num < 0){
          edge_num = edge;
          edge++;
        }

        // Set the edge numbers for all the edges that we found
        for ( int ii = 0; ii < nedges; ii++ ){
          face_edge_conn[edge_index[ii]] = edge_num;
        }
      }
    }
  }

  // Set the total number of edges
  num_edges = edge;
}

/*
  Compute the face index that owns the edges/nodes
*/
void TMRQuadForest::computeFaceOwners(){
  // Allocate the edge/node ownership data
  edge_face_owners = new int[ num_edges ];
  node_face_owners = new int[ num_nodes ];

  // The owner is chosen as the connecting face with the lowest face
  // number
  for ( int edge = 0; edge < num_edges; edge++ ){
    edge_face_owners[edge] = num_faces;

    int ipend = edge_face_ptr[edge+1];
    for ( int ip = edge_face_ptr[edge]; ip < ipend; ip++ ){
      int face = edge_face_conn[ip]/4;
      if (face < edge_face_owners[edge]){
        edge_face_owners[edge] = face;
      }
    }
  }

  // Find the node owners
  for ( int node = 0; node < num_nodes; node++ ){
    node_face_owners[node] = num_faces;
    
    int ipend = node_face_ptr[node+1];
    for ( int ip = node_face_ptr[node]; ip < ipend; ip++ ){
      int face = node_face_conn[ip]/4;
      if (face < node_face_owners[node]){
        node_face_owners[node] = face;
      }
    }
  }
}

/*
  Write a representation of the connectivity of the forest out to a
  VTK file.
*/
void TMRQuadForest::writeToVTK( const char *filename ){
  if (mpi_rank == 0 && topo){
    FILE *fp = fopen(filename, "w");
    if (fp){
      fprintf(fp, "# vtk DataFile Version 3.0\n");
      fprintf(fp, "vtk output\nASCII\n");
      fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

      // Write out the points
      fprintf(fp, "POINTS %d float\n", num_nodes);
      for ( int k = 0; k < num_nodes; k++ ){
        // Get the owner of this node
        int face = node_face_owners[k];

        // Check where the node is located
        int corner = 0;
        for ( ; corner < 4; corner++ ){
          if (face_conn[4*face + corner] == k){
            break;
          }
        }

        // Determine the parametric location of the point p
        double u = 0.0, v = 0.0;
        if (corner & 1){ u = 1.0; }
        if (corner & 2){ v = 1.0; }
        
        // Get the surface object and evaluate the point
        TMRFace *surf;
        topo->getFace(face, &surf);

        TMRPoint p;
        if (surf){
          surf->evalPoint(u, v, &p);
        }

        // Write the point to the file
        fprintf(fp, "%e %e %e\n", p.x, p.y, p.z);
      }
        
      // Write out the cells
      fprintf(fp, "\nCELLS %d %d\n", num_faces, 5*num_faces);
      for ( int k = 0; k < num_faces; k++ ){
        fprintf(fp, "4 %d %d %d %d\n", 
                face_conn[4*k], face_conn[4*k+1],
                face_conn[4*k+3], face_conn[4*k+2]);
      }

      // All quadrilaterals
      fprintf(fp, "\nCELL_TYPES %d\n", num_faces);
      for ( int k = 0; k < num_faces; k++ ){
        fprintf(fp, "%d\n", 9);
      }

      // Print out the rest as fields one-by-one
      fprintf(fp, "CELL_DATA %d\n", num_faces);
      fprintf(fp, "SCALARS entity_index float 1\n");
      fprintf(fp, "LOOKUP_TABLE default\n");
      for ( int k = 0; k < num_faces; k++ ){
        fprintf(fp, "%e\n", 1.0*k);
      }
      fclose(fp);
    }
  }
}

/*
  Write a representation of the connectivity of the forest out to a
  VTK file.
*/
void TMRQuadForest::writeToTecplot( const char *filename ){
  if (mpi_rank == 0 && topo){
    FILE *fp = fopen(filename, "w");
    if (fp){
      fprintf(fp, "Variables = X,Y,Z,face\n");
      fprintf(fp, "Zone N = %d E = %d ", num_nodes, num_faces);
      fprintf(fp, "DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL\n");
      fprintf(fp, "VARLOCATION = ([4]=CELLCENTERED)\n");

      // Allocate temporary data for the vertices
      TMRPoint *Xtmp = new TMRPoint[ num_nodes ];

      // Write out the points
      for ( int k = 0; k < num_nodes; k++ ){
        // Get the owner of this node
        int face = node_face_owners[k];

        // Check where the node is located
        int corner = 0;
        for ( ; corner < 4; corner++ ){
          if (face_conn[4*face + corner] == k){
            break;
          }
        }

        // Determine the parametric location of the point p
        double u = 0.0, v = 0.0;
        if (corner & 1){ u = 1.0; }
        if (corner & 2){ v = 1.0; }
        
        // Get the surface object and evaluate the point
        TMRFace *surf;
        topo->getFace(face, &surf);

        TMRPoint p;
        if (surf){
          surf->evalPoint(u, v, &Xtmp[k]);
        }
      }

      // Write out the nodes
      for ( int k = 0; k < num_nodes; k++ ){
        fprintf(fp, "%e\n", Xtmp[k].x);
      }
      for ( int k = 0; k < num_nodes; k++ ){
        fprintf(fp, "%e\n", Xtmp[k].y);
      }
      for ( int k = 0; k < num_nodes; k++ ){
        fprintf(fp, "%e\n", Xtmp[k].z);
      }
      delete [] Xtmp;
        
      // Write out the cell values
      for ( int k = 0; k < num_faces; k++ ){
        fprintf(fp, "%e\n", 1.0*k);
      }

      // Write out the connectivity
      for ( int k = 0; k < num_faces; k++ ){
        fprintf(fp, "%d %d %d %d\n", 
                face_conn[4*k]+1, face_conn[4*k+1]+1,
                face_conn[4*k+3]+1, face_conn[4*k+2]+1);
      }
      fclose(fp);
    }
  }
}

/*
  Write the entire forest to a VTK file
*/
void TMRQuadForest::writeForestToVTK( const char *filename ){
  if (quadrants && topo){
    // Write out the vtk file
    FILE *fp = fopen(filename, "w");
    
    if (fp){
      fprintf(fp, "# vtk DataFile Version 3.0\n");
      fprintf(fp, "vtk output\nASCII\n");
      fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
      
      // Get the quadrants
      int size;
      TMRQuadrant *array;
      quadrants->getArray(&array, &size);
      
      // Write out the points
      fprintf(fp, "POINTS %d float\n", 4*size);

      // Set the edge length
      const int32_t hmax = 1 << TMR_MAX_LEVEL;

      for ( int k = 0; k < size; k++ ){
        const int32_t h = 1 << (TMR_MAX_LEVEL - array[k].level);
                
        // Get the surface object and evaluate the point
        TMRFace *surf;
        topo->getFace(array[k].face, &surf);

        for ( int jj = 0; jj < 2; jj++ ){
          for ( int ii = 0; ii < 2; ii++ ){
            double u = 1.0*(array[k].x + ii*h)/hmax;
            double v = 1.0*(array[k].y + jj*h)/hmax;

            TMRPoint p;
            surf->evalPoint(u, v, &p);
            fprintf(fp, "%e %e %e\n", p.x, p.y, p.z);
          }
        }
      }
  
      // Write out the cell values
      fprintf(fp, "\nCELLS %d %d\n", size, 5*size);
      for ( int k = 0; k < size; k++ ){
        fprintf(fp, "4 %d %d %d %d\n", 4*k, 4*k+1, 4*k+3, 4*k+2);
      }

      // All quadrilaterals
      fprintf(fp, "\nCELL_TYPES %d\n", size);
      for ( int k = 0; k < size; k++ ){
        fprintf(fp, "%d\n", 9);
      }

      // Print out the rest as fields one-by-one
      fprintf(fp, "CELL_DATA %d\n", size);
      fprintf(fp, "SCALARS entity_index float 1\n");
      fprintf(fp, "LOOKUP_TABLE default\n");
      for ( int k = 0; k < size; k++ ){
        fprintf(fp, "%e\n", 1.0*array[k].face);
      }

      fclose(fp);
    }
  }
}

/*
  Write the entire forest to a VTK file
*/
void TMRQuadForest::writeAdjacentToVTK( const char *filename ){
  if (!adjacent){
    computeAdjacentQuadrants();
  }
  if (topo){
    // Write out the vtk file
    FILE *fp = fopen(filename, "w");
    
    if (fp){
      fprintf(fp, "# vtk DataFile Version 3.0\n");
      fprintf(fp, "vtk output\nASCII\n");
      fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
      
      // Get the quadrants
      int size;
      TMRQuadrant *array;
      adjacent->getArray(&array, &size);
      
      // Write out the points
      fprintf(fp, "POINTS %d float\n", 4*size);

      // Set the edge length
      const int32_t hmax = 1 << TMR_MAX_LEVEL;

      for ( int k = 0; k < size; k++ ){
        const int32_t h = 1 << (TMR_MAX_LEVEL - array[k].level);
                
        // Get the surface object and evaluate the point
        TMRFace *surf;
        topo->getFace(array[k].face, &surf);

        for ( int jj = 0; jj < 2; jj++ ){
          for ( int ii = 0; ii < 2; ii++ ){
            double u = 1.0*(array[k].x + ii*h)/hmax;
            double v = 1.0*(array[k].y + jj*h)/hmax;

            TMRPoint p;
            surf->evalPoint(u, v, &p);
            fprintf(fp, "%e %e %e\n", p.x, p.y, p.z);
          }
        }
      }
  
      // Write out the cell values
      fprintf(fp, "\nCELLS %d %d\n", size, 5*size);
      for ( int k = 0; k < size; k++ ){
        fprintf(fp, "4 %d %d %d %d\n", 4*k, 4*k+1, 4*k+3, 4*k+2);
      }

      // All quadrilaterals
      fprintf(fp, "\nCELL_TYPES %d\n", size);
      for ( int k = 0; k < size; k++ ){
        fprintf(fp, "%d\n", 9);
      }

      // Print out the rest as fields one-by-one
      fprintf(fp, "CELL_DATA %d\n", size);
      fprintf(fp, "SCALARS entity_index float 1\n");
      fprintf(fp, "LOOKUP_TABLE default\n");
      for ( int k = 0; k < size; k++ ){
        fprintf(fp, "%e\n", 1.0*array[k].face);
      }

      fclose(fp);
    }
  }
}

/*
  Retrieve the element index
*/
int TMRQuadForest::getElementIndex( TMRQuadrant *element ){
  if (quadrants){
    TMRQuadrant *t = quadrants->contains(element);
    TMRQuadrant *array;
    quadrants->getArray(&array, NULL);
    if (t){
      return t - array;
    }
  }

  return -1;
}

/*
  Retrieve information about the connectivity between faces, edges and
  nodes 
*/
void TMRQuadForest::getConnectivity( int *_nfaces, 
                                     int *_nedges, int *_nnodes, 
                                     const int **_face_conn, 
                                     const int **_face_edge_conn ){
  if (_nfaces){ *_nfaces = num_faces; }
  if (_nedges){ *_nedges = num_edges; }
  if (_nnodes){ *_nnodes = num_nodes; }
  if (_face_conn){ *_face_conn = face_conn; }
  if (_face_edge_conn){ *_face_edge_conn = face_edge_conn; }
}

/*
  Retrieve the inverse of the connectivity
*/
void TMRQuadForest::getInverseConnectivity( const int **_node_face_conn,
                                            const int **_node_face_ptr,
                                            const int **_edge_face_conn,
                                            const int **_edge_face_ptr ){
  if (_node_face_conn){ *_node_face_conn = node_face_conn; }
  if (_node_face_ptr){ *_node_face_ptr = node_face_ptr; }
  if (_edge_face_conn){ *_edge_face_conn = edge_face_conn; }
  if (_edge_face_ptr){ *_edge_face_ptr = edge_face_ptr; }
}

/*
  Create a forest with the specified refinement level
*/
void TMRQuadForest::createTrees( int refine_level ){
  int32_t level = refine_level;
  if (level < 0){
    level = 0;
  }
  else if (level >= TMR_MAX_LEVEL){
    level = TMR_MAX_LEVEL-1;
  }

  // Set who owns what faces
  int nfaces = num_faces/mpi_size;
  int remain = num_faces % mpi_size;
  int start = mpi_rank*nfaces;
  int end = (mpi_rank+1)*nfaces;
  if (mpi_rank < remain){
    nfaces += 1;
    start += mpi_rank;
    end += mpi_rank+1;
  }
  else {
    start += remain;
    end += remain;
  }
  
  // Create an array of the quadrants that will be stored
  int nelems = 1 << level;
  int size = nelems*nelems*nfaces;
  TMRQuadrant *array = new TMRQuadrant[ size ];

  // Generate all of the quadrants on the associated faces
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - level);
  for ( int count = 0, face = start; face < end; face++ ){
    for ( int32_t x = 0; x < hmax; x += h ){
      for ( int32_t y = 0; y < hmax; y += h ){
        array[count].tag = 0;
        array[count].face = face;
        array[count].level = level;
        array[count].x = x;
        array[count].y = y;
        count++;
      }
    }
  }

  // Create the array of quadrants
  quadrants = new TMRQuadrantArray(array, size);
  quadrants->sort();

  // Set the last quadrant
  TMRQuadrant p;
  p.tag = -1;
  p.face = num_faces-1;
  p.x = p.y = hmax;
  if (size > 0){
    p = array[0];
  }

  owners = new TMRQuadrant[ mpi_size ];
  MPI_Allgather(&p, 1, TMRQuadrant_MPI_type, 
                owners, 1, TMRQuadrant_MPI_type, comm);

  // Set the offsets if some of the processors have zero
  // quadrants
  for ( int k = 1; k < mpi_size; k++ ){
    if (owners[k].tag == -1){
      owners[k] = owners[k-1];
    }
  }
}

/*
  Create a forest with the specified refinement level
*/
void TMRQuadForest::createRandomTrees( int nrand,
                                       int min_level, int max_level ){
  // Set who owns what faces
  int nfaces = num_faces/mpi_size;
  int remain = num_faces % mpi_size;
  int start = mpi_rank*nfaces;
  int end = (mpi_rank+1)*nfaces;
  if (mpi_rank < remain){
    nfaces += 1;
    start += mpi_rank;
    end += mpi_rank+1;
  }
  else {
    start += remain;
    end += remain;
  }
  
  // Create an array of the quadrants that will be stored
  int size = nrand*nfaces;
  TMRQuadrant *array = new TMRQuadrant[ size ];

  // Generate a random number of quadrants along random directions
  for ( int count = 0, face = start; face < end; face++ ){
    for ( int i = 0; i < nrand; i++, count++ ){
      int32_t level = min_level + (rand() % (max_level - min_level + 1));

      const int32_t h = 1 << (TMR_MAX_LEVEL - level);
      int32_t x = h*(rand() % (1 << level));
      int32_t y = h*(rand() % (1 << level));
      
      array[count].tag = 0;
      array[count].face = face;
      array[count].level = level;
      array[count].x = x;
      array[count].y = y;
    }
  }

  // Create the array of quadrants
  quadrants = new TMRQuadrantArray(array, size);
  quadrants->sort();

  // Set the last quadrant
  TMRQuadrant p;
  p.tag = -1;
  p.face = num_faces-1;
  p.x = p.y = 1 << TMR_MAX_LEVEL;
  if (size > 0){
    p = array[0];
  }

  owners = new TMRQuadrant[ mpi_size ];
  MPI_Allgather(&p, 1, TMRQuadrant_MPI_type, 
                owners, 1, TMRQuadrant_MPI_type, comm);

  // Set the offsets if some of the processors have zero
  // quadrants
  for ( int k = 1; k < mpi_size; k++ ){
    if (owners[k].tag == -1){
      owners[k] = owners[k-1];
    }
  }
}

/*
  Repartition the quadrants across all processors.

  This does not repartition the nodes. You have to recreate the nodes
  after this call so be careful.
*/
void TMRQuadForest::repartition(){
  // First, this stores the number of elements on quadtrees owned on
  // each processor
  int *ptr = new int[ mpi_size+1 ];
  int size;
  TMRQuadrant *array;
  quadrants->getArray(&array, &size);

  // Gather the sizes from all the arrays
  MPI_Allgather(&size, 1, MPI_INT, &ptr[1], 1, MPI_INT, comm);

  // Set the pointers
  ptr[0] = 0;
  for ( int k = 0; k < mpi_size; k++ ){
    ptr[k+1] += ptr[k];
  }

  // Compute the average size of the new counts
  int average_count = ptr[mpi_size]/mpi_size;
  int remain = ptr[mpi_size] - average_count*mpi_size;
  
  // Figure out what goes where on the new distribution of quadrants
  int *new_ptr = new int[ mpi_size+1 ];
  new_ptr[0] = 0;
  for ( int k = 0; k < mpi_size; k++ ){
    new_ptr[k+1] = new_ptr[k] + average_count;
    if (k < remain){
      new_ptr[k+1] += 1;
    }
  }

  // Allocate the new array of quadrants
  int new_size = new_ptr[mpi_rank+1] - new_ptr[mpi_rank];
  TMRQuadrant *new_array = new TMRQuadrant[ new_size ];

  // Ptr:      |----|---|--------------------|-| 
  // New ptr:  |-------|-------|-------|-------|

  // Count up the number of sends/recvs
  int send_count = 0;
  for ( int i = 0; i < mpi_size; i++ ){
    if (i != mpi_rank && 
        new_ptr[i+1] >= ptr[mpi_rank] &&
        new_ptr[i] < ptr[mpi_rank+1]){
      send_count++;
    }
  }

  // Allocate arrays send/recv arrays
  MPI_Request *send_requests = new MPI_Request[ send_count ];

  // Zero the send/recv counts
  send_count = 0;

  // Match up the intervals for sends
  for ( int i = 0; i < mpi_size; i++ ){
    // If we're on the correct interval, 
    if (new_ptr[i+1] >= ptr[mpi_rank] &&
        new_ptr[i] < ptr[mpi_rank+1]){
      // The start and end locations within the local array
      int start = new_ptr[i] - ptr[mpi_rank];
      if (start < 0){ start = 0; }

      int max_count = ptr[mpi_rank+1] - ptr[mpi_rank];
      int end = new_ptr[i+1] - ptr[mpi_rank];
      if (end > max_count){ end = max_count; }

      // Set the count
      int count = end - start;
      
      if (i == mpi_rank){
        int new_start = ptr[i] - new_ptr[i];
        if (new_start < 0){ new_start = 0; }
        memcpy(&new_array[new_start], &array[start], 
               count*sizeof(TMRQuadrant)); 
      }
      else if (count > 0){
        // Send the element array to the new owner
        MPI_Isend(&array[start], count, TMRQuadrant_MPI_type,
                  i, 0, comm, &send_requests[send_count]);
        send_count++;
      }
    }
  }

  // Match up the intervals for recvs
  for ( int i = 0; i < mpi_size; i++ ){
    // If we're on the correct interval, 
    if (i != mpi_rank && 
        ptr[i+1] >= new_ptr[mpi_rank] && 
        ptr[i] < new_ptr[mpi_rank+1]){
      // The start location within the local array
      int start = ptr[i] - new_ptr[mpi_rank];
      if (start < 0){ start = 0; }

      // Set the end location
      int max_count = new_ptr[mpi_rank+1] - new_ptr[mpi_rank];
      int end = ptr[i+1] - new_ptr[mpi_rank];
      if (end > max_count){ end = max_count; }

      // Set the count
      int count = end - start;

      // Send the element array to the new owner
      if (count > 0){
        MPI_Recv(&new_array[start], count, TMRQuadrant_MPI_type,
                 i, 0, comm, MPI_STATUS_IGNORE);
      }
    }
  }

  // Free the memory
  delete [] ptr;
  delete [] new_ptr;

  // Wait for any remaining sends to complete
  MPI_Waitall(send_count, send_requests, MPI_STATUSES_IGNORE);
  delete [] send_requests;

  // Free the quadrant arrays
  delete quadrants;
  quadrants = new TMRQuadrantArray(new_array, new_size);

  owners = new TMRQuadrant[ mpi_size ];
  MPI_Allgather(&new_array[0], 1, TMRQuadrant_MPI_type, 
                owners, 1, TMRQuadrant_MPI_type, comm);

  // Set the local reordering for the elements
  quadrants->getArray(&new_array, &new_size);
  for ( int i = 0; i < new_size; i++ ){
    new_array[i].tag = i;
  }
}

/*
  Duplicate the forest

  This function creates a duplicate representation of the current
  forest. This function copies the global connectivity of the forest
  and copies each individual tree.
*/
TMRQuadForest *TMRQuadForest::duplicate(){
  TMRQuadForest *dup = new TMRQuadForest(comm);
  if (face_conn){
    copyData(dup);

    // Copy the quadrants
    dup->quadrants = quadrants->duplicate();
    dup->owners = new TMRQuadrant[ mpi_size ];
    memcpy(dup->owners, owners, sizeof(TMRQuadrant)*mpi_size);
  }

  return dup;
}

/*
  Coarsen the entire forest

  This function creates a coarsened representation of the current
  forest. This is done by copying the global connectivity of the
  forest and coarsening each individual tree. Note that the resulting
  forest is not necessarily balanced.  
*/
TMRQuadForest *TMRQuadForest::coarsen(){
  TMRQuadForest *coarse = new TMRQuadForest(comm);
  if (face_conn){
    copyData(coarse);

    // Coarsen the quadrants in the array
    int size;
    TMRQuadrant *array;
    quadrants->getArray(&array, &size);

    // Create a new queue of quadrants
    TMRQuadrantQueue *queue = new TMRQuadrantQueue();
 
    // Scan through the list, if we have offset quadrants which
    // all share the same parent, then coarsen the parent
    for ( int i = 0; i < size; i++ ){
      if (array[i].level > 0){
        if (array[i].childId() == 0){
          TMRQuadrant p;
          array[i].parent(&p);
          queue->push(&p);
        }
      }
      else {
        queue->push(&array[i]);
      }
    }

    // Create the coarse quadrants
    coarse->quadrants = queue->toArray();
    delete queue;

    // Set the owner array
    coarse->quadrants->getArray(&array, &size);
    coarse->owners = new TMRQuadrant[ mpi_size ];
    MPI_Allgather(&array[0], 1, TMRQuadrant_MPI_type, 
                  coarse->owners, 1, TMRQuadrant_MPI_type, comm);
  }

  return coarse;
}

/*
  Refine the quadrant mesh based on the input refinement level
*/
void TMRQuadForest::refine( const int refinement[],
                            int min_level, int max_level ){
  // Adjust the min and max levels to ensure consistency
  if (min_level < 0){ min_level = 0; }
  if (max_level > TMR_MAX_LEVEL){ max_level = TMR_MAX_LEVEL; }

  // This is just a sanity check
  if (min_level > max_level){ min_level = max_level; }

  // Free memory if it has been allocated
  if (adjacent){ delete adjacent; }
  if (X){ delete [] X; }
  adjacent = NULL;
  X = NULL;
  
  // Create a hash table for the refined quadrants and the quadrants
  // that are external (on other processors)
  TMRQuadrantHash *hash = new TMRQuadrantHash();
  TMRQuadrantHash *ext_hash = new TMRQuadrantHash();

  // Get the current array of quadrants
  int size;
  TMRQuadrant *array;
  quadrants->getArray(&array, &size);

  if (refinement){
    for ( int i = 0; i < size; i++ ){
      if (refinement[i] == 0){
        // We know that this quadrant is locally owned
        hash->addQuadrant(&array[i]);
      }
      else if (refinement[i] < 0){
        // Coarsen this quadrant
        if (array[i].level > min_level){
          // Compute the new refinement level
          int new_level = array[i].level + refinement[i];
          if (new_level < min_level){
            new_level = min_level;
          }

          // Copy over the quadrant
          TMRQuadrant q = array[i];
          q.level = new_level;
          
          // Compute the new side-length of the quadrant
          const int32_t h = 1 << (TMR_MAX_LEVEL - q.level);
          q.x = q.x - (q.x % h);
          q.y = q.y - (q.y % h);
          if (mpi_rank == getQuadrantMPIOwner(&q)){
            hash->addQuadrant(&q);
          }
          else {
            ext_hash->addQuadrant(&q);
          }
        }
        else {
          // If it is already at the min level, just add it
          hash->addQuadrant(&array[i]);
        }
      }
      else if (refinement[i] > 0){
        // Refine this quadrant
        if (array[i].level < max_level){
          // Compute the new refinement level
          int new_level = array[i].level + refinement[i];
          if (new_level > max_level){
            new_level = max_level;
          }

          // Compute the relative level of refinement
          int ref = new_level - array[i].level;
          if (ref <= 0){ 
            ref = 1; 
          }
          else {
            ref = 1 << (ref - 1);
          }
          
          // Copy the quadrant and set the new level
          TMRQuadrant q = array[i];
          q.level = new_level;

          // Compute the new side-length of the quadrant
          const int32_t h = 1 << (TMR_MAX_LEVEL - q.level);
          int32_t x = q.x - (q.x % h);
          int32_t y = q.y - (q.y % h);
          for ( int ii = 0; ii < ref; ii++ ){
            for ( int jj = 0; jj < ref; jj++ ){
              q.x = x + 2*ii*h;
              q.y = y + 2*jj*h;
              if (mpi_rank == getQuadrantMPIOwner(&q)){
                hash->addQuadrant(&q);
              }
              else {
                ext_hash->addQuadrant(&q);
              }
            }
          }
        }
        else {
          // If the quadrant is at the max level add it without
          // refinement
          hash->addQuadrant(&array[i]);
        }
      }
    }
  }
  else {
    // No refinement array is provided. Just go ahead and refine
    // everything one level
    for ( int i = 0; i < size; i++ ){
      if (array[i].level < max_level){
        TMRQuadrant q = array[i];
        q.level += 1;
        if (mpi_rank == getQuadrantMPIOwner(&q)){
          hash->addQuadrant(&q);
        }
        else {
          ext_hash->addQuadrant(&q);
        }
      }
      else {
        hash->addQuadrant(&array[i]);
      }
    }
  }
  
  // Free the old quadrants class
  delete quadrants;

  // Sort the list of external quadrants
  TMRQuadrantArray *list = ext_hash->toArray();
  list->sort();
  delete ext_hash;

  // Get the local list quadrants added from other processors
  // and add them to the local hash table
  TMRQuadrantArray *local = distributeQuadrants(list);
  delete list;

  // Get the local quadrants and add them to the hash table
  local->getArray(&array, &size);
  for ( int i = 0; i < size; i++ ){
    hash->addQuadrant(&array[i]);
  }
  delete local;

  // Cover the hash table to a list and uniquely sort it
  quadrants = hash->toArray();
  quadrants->sort();

  delete hash;

  // Get the octants and order their labels
  quadrants->getArray(&array, &size);
  for ( int i = 0; i < size; i++ ){
    array[i].tag = i;
  }
}

/*
  Get the owner of the quadrant
*/
int TMRQuadForest::getQuadrantMPIOwner( TMRQuadrant *quad ){
  int rank = 0;
  for ( ; (rank < mpi_size-1 && 
           owners[rank+1].comparePosition(quad) <= 0); rank++ );

  return rank;
}

/*
  Match the quadrant intervals to the MPI owners. Note that this
  requires that the input array is sorted.
*/
void TMRQuadForest::matchQuadrantIntervals( TMRQuadrant *array,
                                            int size, int *ptr ){
  ptr[0] = 0;

  int index = 0;
  for ( int rank = 0; rank < mpi_size-1; rank++ ){
    while (index < size &&
           owners[rank+1].comparePosition(&array[index]) > 0){
      index++;
    }
    ptr[rank+1] = index;
  }
  ptr[mpi_size] = size;
}

/*
  Match the MPI intervals
*/
void TMRQuadForest::matchTagIntervals( TMRQuadrant *array,
                                       int size, int *ptr ){
  ptr[0] = 0;
  for ( int i = 0, rank = 0; rank < mpi_size; rank++ ){
    while (i < size && array[i].tag <= rank){
      i++;
    }
    ptr[rank+1] = i;
  }
  ptr[mpi_size] = size;
}

/*
  Send a distributed list of quadrants to their owner processors
*/
TMRQuadrantArray *TMRQuadForest::distributeQuadrants( TMRQuadrantArray *list,
                                                      int use_tags,
                                                      int **_quad_ptr, 
                                                      int **_quad_recv_ptr,
                                                      int include_local ){
  // Get the array itself
  int size;
  TMRQuadrant *array;
  list->getArray(&array, &size);

  // The number of quadrants that will be sent from this processor
  // to all other processors in the communicator
  int *quad_ptr = new int[ mpi_size+1 ];
  int *quad_recv_ptr = new int[ mpi_size+1 ];

  // Match the quadrant intervals to determine how mnay quadrants
  // need to be sent to each processor
  if (use_tags){
    matchTagIntervals(array, size, quad_ptr);
  }
  else {
    matchQuadrantIntervals(array, size, quad_ptr);
  }

  // Count up the number of quadrants
  int *quad_counts = new int[ mpi_size ];
  for ( int i = 0; i < mpi_size; i++ ){
    if (!include_local && i == mpi_rank){
      quad_counts[i] = 0;
    }
    else {
      quad_counts[i] = quad_ptr[i+1] - quad_ptr[i];
    }
  }

  // Now distribute the quadrants to their destination quadrees and
  // balance the corresponding quadrees including the new elements.
  int *quad_recv_counts = new int[ mpi_size ];
  MPI_Alltoall(quad_counts, 1, MPI_INT,
               quad_recv_counts, 1, MPI_INT, comm);

  // Now use quad_ptr to point into the recv array
  quad_recv_ptr[0] = 0;
  for ( int i = 0; i < mpi_size; i++ ){
    quad_recv_ptr[i+1] = quad_recv_ptr[i] + quad_recv_counts[i];
  }

  delete [] quad_counts;
  delete [] quad_recv_counts;

  // Create the distributed array
  TMRQuadrantArray *dist = sendQuadrants(list, quad_ptr, quad_recv_ptr);

  // Free other data associated with the parallel communication
  if (_quad_ptr){ 
    *_quad_ptr = quad_ptr;
  }
  else {
    delete [] quad_ptr;
  }
  if (_quad_recv_ptr){
    *_quad_recv_ptr = quad_recv_ptr;
  }
  else {
    delete [] quad_recv_ptr;
  }

  return dist;
}

/*
  Send the quadrants to the processors designated by the pointer arrays
*/
TMRQuadrantArray *TMRQuadForest::sendQuadrants( TMRQuadrantArray *list,
                                                const int *quad_ptr,
                                                const int *quad_recv_ptr ){
  // Get the array itself
  int size;
  TMRQuadrant *array;
  list->getArray(&array, &size);

  // Count up the number of recvs
  int nsends = 0, nrecvs = 0;
  for ( int i = 0; i < mpi_size; i++ ){
    if (i != mpi_rank){
      if (quad_ptr[i+1] - quad_ptr[i] > 0){
        nsends++;
      }
      if (quad_recv_ptr[i+1] - quad_recv_ptr[i] > 0){
        nrecvs++;
      }
    }
  }

  // Allocate the space for the recving array
  int recv_size = quad_recv_ptr[mpi_size];
  TMRQuadrant *recv_array = new TMRQuadrant[ recv_size ];

  // Allocate space for the requests
  MPI_Request *send_request = new MPI_Request[ nsends ];

  // Loop over all the ranks and send 
  for ( int i = 0, j = 0; i < mpi_size; i++ ){
    if (i != mpi_rank && quad_ptr[i+1] - quad_ptr[i] > 0){
      // Post the send to the destination
      int count = quad_ptr[i+1] - quad_ptr[i];
      MPI_Isend(&array[quad_ptr[i]], count, TMRQuadrant_MPI_type, 
                i, 0, comm, &send_request[j]);
      j++;
    }
    else if (i == mpi_rank){
      int count = quad_recv_ptr[i+1] - quad_recv_ptr[i];
      if (count > 0 &&
          (count == quad_ptr[i+1] - quad_ptr[i])){
        memcpy(&recv_array[quad_recv_ptr[i]], 
               &array[quad_ptr[i]], count*sizeof(TMRQuadrant));
      }
    }
  }

  // Loop over the recieve calls
  for ( int i = 0; i < mpi_size; i++ ){
    if (i != mpi_rank && quad_recv_ptr[i+1] > quad_recv_ptr[i]){
      int recv_count = quad_recv_ptr[i+1] - quad_recv_ptr[i];
     
      MPI_Recv(&recv_array[quad_recv_ptr[i]], recv_count, 
               TMRQuadrant_MPI_type,
               i, 0, comm, MPI_STATUS_IGNORE);
    }
  }

  // Wait for any remaining sends to complete
  MPI_Waitall(nsends, send_request, MPI_STATUSES_IGNORE);
  delete [] send_request;

  return new TMRQuadrantArray(recv_array, recv_size);
}

/*
  Add the edge neighbors for adjacent trees

  This function is called to balance the forest across tree edges.
  Given an quadrant p on the specified edge index, this code ensures a
  edge balanced tree, by adding the corresponding edge quadrants to
  all edge-adjacent quadtrees. If the quadrant is added to the hash
  object, it is then appended to the local face queue to ensure that
  it is also balanced locally.

  input:
  edge_index:  index for the given edge
  p:           quadrant to balance
  hash:        the array of hash objects for each processor
  queue:       the array of newly added qudrants for each processor
*/
void TMRQuadForest::addEdgeNeighbors( int edge_index, 
                                      TMRQuadrant p,
                                      TMRQuadrantHash *hash,
                                      TMRQuadrantHash *ext_hash,
                                      TMRQuadrantQueue *queue ){
  // First determine the global edge number
  int face = p.face;
  int edge = face_edge_conn[4*face + edge_index];
  
  // Compute the edge length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - p.level);

  // Store the u coordinate along the edge
  int32_t ucoord = 0;
  if (edge_index < 2){
    ucoord = p.y;
  }
  else {
    ucoord = p.x;
  }

  // Retrieve the first and second node numbers
  int n1 = face_conn[4*face + face_to_edge_nodes[edge_index][0]];
  int n2 = face_conn[4*face + face_to_edge_nodes[edge_index][1]];

  // Now, cycle through all the adjacent faces
  for ( int ip = edge_face_ptr[edge]; ip < edge_face_ptr[edge+1]; ip++ ){
    // Get the face that is adjacent across this edge
    int adj = edge_face_conn[ip]/4;
    if (adj != face){
      // Get the adjacent edge index
      int adj_index = edge_face_conn[ip] % 4;

      // Get the nodes on the adjacent face
      int nn1 = face_conn[4*adj + face_to_edge_nodes[adj_index][0]];
      int nn2 = face_conn[4*adj + face_to_edge_nodes[adj_index][1]];

      // Add the quadrant to the list
      int reverse = (n1 == nn2 && n2 == nn1);
      int32_t u = ucoord;
      if (reverse){
        u = hmax - 2*h - ucoord;
      }
      
      TMRQuadrant neighbor;
      neighbor.face = adj;
      neighbor.level = p.level;
      if (adj_index < 2){
        neighbor.x = (hmax - 2*h)*(adj_index % 2);
        neighbor.y = u;
      }
      else {
        neighbor.x = u;
        neighbor.y = (hmax - 2*h)*(adj_index % 2);
      }

      // Find the quadrant owner and add the quadrant to the hash
      // tables and possibly queue
      int owner = getQuadrantMPIOwner(&neighbor);
      if (owner == mpi_rank){
        if (hash->addQuadrant(&neighbor)){
          queue->push(&neighbor);
        }
      }
      else if (ext_hash && ext_hash->addQuadrant(&neighbor)){
        queue->push(&neighbor);
      }
    }
  }
}

/*
  Add the corner neighbors for a given tree

  This function is called to balance the forest across tree corners.
  Given an quadrant p on the specified corner index, this code ensures a
  corner balanced tree, by adding the corresponding corner quadrants to
  all node-adjacent quadtrees. If the quadrant is added to the hash
  object, it is then appended to the local face queue to ensure that
  it is also balanced locally.

  input:
  corner:  the corner index (p must lie on this corner)
  p:       the local quadrant
  hash:    the array of hash objects for each processor
  queue:   the array of newly added qudrants for each processor
*/
void TMRQuadForest::addCornerNeighbors( int corner,
                                        TMRQuadrant p,
                                        TMRQuadrantHash *hash,
                                        TMRQuadrantHash *ext_hash,
                                        TMRQuadrantQueue *queue ){
  // First determine the global edge number
  int face = p.face;
  int node = face_conn[4*face + corner];
  
  // Compute the edge length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - p.level);

  // Now, cycle through all the adjacent faces
  for ( int ip = node_face_ptr[node];
        ip < node_face_ptr[node+1]; ip++ ){
      
    // Get the faces that are adjacent across this edge
    int adj = node_face_conn[ip]/4;
    if (adj != face){
      int adj_index = node_face_conn[ip] % 4;

      // Compute the quadrant location
      TMRQuadrant neighbor;
      neighbor.face = adj;
      neighbor.level = p.level;
      neighbor.x = (hmax - 2*h)*(adj_index % 2);
      neighbor.y = (hmax - 2*h)*(adj_index/2);
      
      // Find the quadrant owner and add the quadrant to the hash
      // tables and possibly queue
      int owner = getQuadrantMPIOwner(&neighbor);
      if (owner == mpi_rank){
        if (hash->addQuadrant(&neighbor)){
          queue->push(&neighbor);
        }
      }  
      else if (ext_hash && ext_hash->addQuadrant(&neighbor)){
        queue->push(&neighbor);
      }
    }
  }
}

/*
  Balance the quadrant on the entire quadtree

  This code finds the 0-parent of all adjacent quadrants either within
  the current tree or within an adjacent tree and adds those quadrants
  to balance the input quadrant 'quad'.

  input:
  quad:            the quadrant itself
  hash:            the array of hash tables for each face
  ext_hash:        the array of hash tables for each face
  queue:           the quadrant queues for each face
  balance_corner:  balance across corners 
  balance_tree:    balance on the entire tree
*/
void TMRQuadForest::balanceQuadrant( TMRQuadrant *quad,
                                     TMRQuadrantHash *hash,
                                     TMRQuadrantHash *ext_hash,
                                     TMRQuadrantQueue *queue,
                                     const int balance_corner,
                                     const int balance_tree ){
  // Local quadrant data
  TMRQuadrant p, neighbor, q;

  // Get the max level
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  
  // Get the parent of the quadrant, and add the their
  // face-matched quadrants from each face, as long 
  // as they fall within the bounds
  if (quad->level > 1){
    quad->parent(&p);
        
    // Add the edge-adjacent elements
    for ( int edge = 0; edge < 4; edge++ ){
      p.edgeNeighbor(edge, &neighbor);
      neighbor.getSibling(0, &q);
	  
      // If we're in bounds, add the neighbor
      if ((q.x >= 0 && q.x < hmax) &&
          (q.y >= 0 && q.y < hmax)){
        // Get the MPI rank of the quadrant owner
        int owner = getQuadrantMPIOwner(&q);
        if (owner == mpi_rank){
          if (hash->addQuadrant(&q)){
            queue->push(&q);
          }
        }
        else if (ext_hash && ext_hash->addQuadrant(&q)){
          queue->push(&q);
        }
      }
      else if (balance_tree){
        // The node may lie across an edge or face
        int ex = (q.x < 0 || q.x >= hmax);
        int ey = (q.y < 0 || q.y >= hmax);
        
        if (ex || ey){
          // The quadrant lies along a true edge
          addEdgeNeighbors(edge, q, hash, ext_hash, queue);
        }
      }
    }

    // If we're balancing across edges and corners
    if (balance_corner){
      for ( int corner = 0; corner < 4; corner++ ){
        p.cornerNeighbor(corner, &neighbor);
        neighbor.getSibling(0, &q);
        
        if ((q.x >= 0 && q.x < hmax) &&
            (q.y >= 0 && q.y < hmax)){
          // Get the MPI rank of the quadrant owner
          int owner = getQuadrantMPIOwner(&q);
          if (owner == mpi_rank){
            if (hash->addQuadrant(&q)){
              queue->push(&q);
            }
          }
          else if (ext_hash && ext_hash->addQuadrant(&q)){
            queue->push(&q);
          }
        }
        else if (balance_tree){
          // The node may lie across a corner, edge or face
          int ex = (q.x < 0 || q.x >= hmax);
          int ey = (q.y < 0 || q.y >= hmax);

          if (ex && ey){
            // Add the quadrant to the other trees
            addCornerNeighbors(corner, q, hash, ext_hash, queue);
          }
          else {
            // The quadrant lies along a true edge
            int edge = ex*(q.x < 0 ? 0 : 1) + ey*(q.y < 0 ? 2 : 3);
            addEdgeNeighbors(edge, q, hash, ext_hash, queue);
          }
        }
      }
    }
  }
}

/*
  Balance the forest of quadtrees

  This algorithm uses a hash and a queue to balance the forest of
  quadtrees. For each element in the quadtree, we add the neighbors that
  are required to balance to the tree. If the element is not in the
  hash, we add them to a queue, which keeps track of recently added
  elements. After the first pass, the algorithm continues popping
  elements until the all the queues are empty.

  Note that only 0-th siblings are added/popped on the hash/queue.
  Then at the end, all neighboring siblings are added.

  The type of balancing - face/edge balanced or face/edge/corner
  balanced is determined using the balance_corner flag. Face balancing
  is balancing across faces, edge balancing is balancing across edges
  of the elements and corner balances across corners. The code always
  balances faces and edges (so that there is at most one depdent node
  per edge) and balances across corners optionally.
*/
void TMRQuadForest::balance( int balance_corner ){
  // Create a hash table for the balanced tree
  TMRQuadrantHash *hash = new TMRQuadrantHash();
  TMRQuadrantHash *ext_hash = new TMRQuadrantHash();
  TMRQuadrantQueue *queue = new TMRQuadrantQueue();
  
  // Get the array of quadrants
  int quad_size;
  TMRQuadrant *quad_array;
  quadrants->getArray(&quad_array, &quad_size);
    
  // Add all the elements
  for ( int i = 0; i < quad_size; i++ ){
    TMRQuadrant quad;
    quad_array[i].getSibling(0, &quad);

    // Get the quadrant owner
    int owner = getQuadrantMPIOwner(&quad);

    // Add the owner
    if (owner == mpi_rank){
      hash->addQuadrant(&quad);
    }
    else {
      ext_hash->addQuadrant(&quad);
    }
      
    // Balance the quadrants locally
    const int balance_tree = 1;
    balanceQuadrant(&quad, hash, ext_hash, queue, 
                    balance_corner, balance_tree);
  }

  // Free the original quadrant array and set it to NULL
  delete quadrants;

  while (queue->length() > 0){
    // Now continue until the queue of added quadrants is
    // empty. At each iteration, pop an quadrant and add 
    // its neighbours until nothing new is added. This code
    // handles the propagation of quadrants to adjacent quadrants.
    TMRQuadrant quad = queue->pop();
    const int balance_tree = 1;
    balanceQuadrant(&quad, hash, ext_hash, queue, 
                    balance_corner, balance_tree);
  }

  // Now everything is locally balanced - all the elements on the
  // current processor are balanced with all the other elements on the
  // current processor, but nothing is inter-processor balanced yet.
  // Create a sorted list of local the 0-child quadrants. This can be
  // further reduced to limit the amount of memory passed between
  // processors
  TMRQuadrantArray *elems0 = ext_hash->toArray();
  delete ext_hash;
  elems0->sort();

  // Get the array of 0-quadrants 
  int size;
  TMRQuadrant *array;
  elems0->getArray(&array, &size);

  if (size > 0){
    // Get the parent of the quadrant
    TMRQuadrant p, s = array[0];
    s.parent(&p);

    // Loop over all of the local quadrants 
    for ( int i = 0; i < size; i++ ){
      if (!p.contains(&array[i])){
        queue->push(&s);
      }
      // Get the next quadrant and find its parent
      s = array[i];
      s.parent(&p);
    }
    
    // Push the last quadrant onto the queue
    queue->push(&s);
  }
  
  // Free the elements and the hash table
  delete elems0;

  // Send the quadrants to their destination...
  TMRQuadrantArray *list = queue->toArray();
  delete queue;

  // Get the local quadrants from the list
  TMRQuadrantArray *local = distributeQuadrants(list);
  delete list;

  // Allocate a new queue
  queue = new TMRQuadrantQueue();

  // Get the local array of quadrants and add them to the
  // hash table
  local->getArray(&array, &size);
  for ( int i = 0; i < size; i++ ){
    if (hash->addQuadrant(&array[i])){
      queue->push(&array[i]);
    }
  }
  delete local;

  // Now all the received quadrants will balance the tree locally
  // without having to worry about off-processor quadrants.
  while (queue->length() > 0){
    const int balance_tree = 1;
    TMRQuadrant quad = queue->pop();
    balanceQuadrant(&quad, hash, NULL, queue, 
                    balance_corner, balance_tree);
  }

  // Now convert the elements from child-0 elements to
  // elements which cover the full mesh
  TMRQuadrantArray *child0_elems = hash->toArray();
  child0_elems->getArray(&quad_array, &quad_size);

  // Loop over all elements and add their siblings
  for ( int i = 0; i < quad_size; i++ ){
    if (quad_array[i].level > 0){
      for ( int j = 0; j < 4; j++ ){
        TMRQuadrant q;
        quad_array[i].getSibling(j, &q);
        int owner = getQuadrantMPIOwner(&q);
        if (mpi_rank == owner){
          hash->addQuadrant(&q);
        }
        else {
          queue->push(&q);
        }
      }
    }
  }

  // Free the temporary elements
  delete child0_elems;

  // Turn the queue into an array
  list = queue->toArray();
  delete queue;

  // Sort the list before distributing it
  list->sort();

  // Get the local list quadrants added from other processors
  // and add them to the local hash table
  local = distributeQuadrants(list);
  delete list;

  // Get the local quadrants and add them to the hash table
  local->getArray(&array, &size);
  for ( int i = 0; i < size; i++ ){
    hash->addQuadrant(&array[i]);
  }
  delete local;
  
  // Set the elements into the quadtree
  quadrants = hash->toArray();
  quadrants->sort();

  // Get the quadrants and order their labels
  quadrants->getArray(&array, &size);
  for ( int i = 0; i < size; i++ ){
    array[i].tag = i;
  }

  // Free the hash
  delete hash;
}

/*
  Add the quadrant to the processor queue corresponding to the
  non-local faces that touch the given edge
*/
void TMRQuadForest::addAdjacentEdgeToQueue( int edge_index, 
                                            TMRQuadrant p,
                                            TMRQuadrantQueue *queue, 
                                            TMRQuadrant orig ){
  // First determine the global edge number
  int face = p.face;
  int edge = face_edge_conn[4*face + edge_index];
  
  // Compute the edge length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - p.level);

  // Store the u coordinate along the edge
  int32_t ucoord = 0;
  if (edge_index < 2){
    ucoord = p.y;
  }
  else {
    ucoord = p.x;
  }

  // Retrieve the first and second node numbers
  int n1 = face_conn[4*face + face_to_edge_nodes[edge_index][0]];
  int n2 = face_conn[4*face + face_to_edge_nodes[edge_index][1]];

  // Now, cycle through all the adjacent faces
  for ( int ip = edge_face_ptr[edge]; ip < edge_face_ptr[edge+1]; ip++ ){
    // Get the face that is adjacent across this edge
    int adj = edge_face_conn[ip]/4;
    if (adj != face){
      // Get the adjacent edge index
      int adj_index = edge_face_conn[ip] % 4;

      // Get the nodes on the adjacent face
      int nn1 = face_conn[4*adj + face_to_edge_nodes[adj_index][0]];
      int nn2 = face_conn[4*adj + face_to_edge_nodes[adj_index][1]];

      // Add the quadrant to the list
      int reverse = (n1 == nn2 && n2 == nn1);
      int32_t u = ucoord;
      if (reverse){
        u = hmax - h - ucoord;
      }
      
      TMRQuadrant neighbor;
      neighbor.face = adj;
      neighbor.level = p.level;
      if (adj_index < 2){
        neighbor.x = (hmax - h)*(adj_index % 2);
        neighbor.y = u;
      }
      else {
        neighbor.x = u;
        neighbor.y = (hmax - h)*(adj_index % 2);
      }
      
      // Find the quadrant owner
      int owner = getQuadrantMPIOwner(&neighbor);
      if (owner != mpi_rank){
        orig.tag = owner;
        queue->push(&orig);
      }
    }
  }
}

/*
  Add the quadrant to the queue that correspond to the non-local faces
  that touch the corner 
*/
void TMRQuadForest::addAdjacentCornerToQueue( int corner,
                                              TMRQuadrant p,
                                              TMRQuadrantQueue *queue, 
                                              TMRQuadrant orig ){
  // First determine the global edge number
  int face = p.face;
  int node = face_conn[4*face + corner];
  
  // Compute the edge length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - p.level);

  // Now, cycle through all the adjacent faces
  for ( int ip = node_face_ptr[node];
        ip < node_face_ptr[node+1]; ip++ ){
      
    // Get the faces that are adjacent across this edge
    int adj = node_face_conn[ip]/4;
    if (adj != face){
      int adj_index = node_face_conn[ip] % 4;

      // Compute the quadrant location
      TMRQuadrant neighbor;
      neighbor.face = adj;
      neighbor.level = p.level;
      neighbor.x = (hmax - h)*(adj_index % 2);
      neighbor.y = (hmax - h)*(adj_index/2);
      
      // Find the quadrant owner
      int owner = getQuadrantMPIOwner(&neighbor);
      if (owner != mpi_rank){
        orig.tag = owner;
        queue->push(&orig);
      }
    }
  }
}

/*
  The following code exchanges the neighboring quadrants for each
  locally owned quadtree within the forest.

  This code exchanges non-local quadrants across processors so that we
  can locally query quadrants on adjacent quadtrees without having to
  perform parallel communication.

  Note that this code creates partial non-local quadtrees that are
  adjacent to the local quadtrees in the forest. These partial local
  quadtrees should be freed after the nodal ordering has been
  computed. 
*/
void TMRQuadForest::computeAdjacentQuadrants(){
  if (adjacent){ delete adjacent; }

  // Allocate the queue that stores the quadrants destined for each of
  // the processors
  TMRQuadrantQueue *queue = new TMRQuadrantQueue();
  
  // Get the actual quadrant array
  int size;
  TMRQuadrant *array;
  quadrants->getArray(&array, &size);
    
  // Loop over all the elements and check where we need to send
  // the quadrants that are along each edge/face
  for ( int i = 0; i < size; i++ ){
    const int32_t hmax = 1 << TMR_MAX_LEVEL;
    
    // Enumerate the sibling-ids for each edge
    const int edge_ids[][2] =
      {{0, 2}, {1, 3}, {0, 1}, {2, 3}};

    // Add the edge-adjacent quadrant across the boundary
    for ( int edge_index = 0; edge_index < 4; edge_index++ ){
      for ( int k = 0; k < 2; k++ ){
        // Get the quadrant and increase the level
        TMRQuadrant p = array[i];
        p.level += 1;

        // Get the sibling id for each quadrant along the face that
        // we're on right now
        TMRQuadrant q;
        p.getSibling(edge_ids[edge_index][k], &q);

        // Get the edge neighbor
        q.edgeNeighbor(edge_index, &q);
	  
        // If we're in bounds, add the neighbor
        if ((q.x >= 0 && q.x < hmax) &&
            (q.y >= 0 && q.y < hmax)){
          // Get the MPI rank of the quadrant owner
          int owner = getQuadrantMPIOwner(&q);
          if (owner != mpi_rank){
            p = array[i];
            p.tag = owner;
            queue->push(&p);
          }
        }
        else {
          // The node may lie across an edge
          int ex = (q.x < 0 || q.x >= hmax);
          int ey = (q.y < 0 || q.y >= hmax);
        
          if (ex || ey){
            // The quadrant lies along a true edge
            addAdjacentEdgeToQueue(edge_index, q, queue, array[i]);
          }
        }
      }
    }

    // Add corner-adjacent quadrants
    for ( int corner = 0; corner < 4; corner++ ){
      // Get the quadrant and increase the level
      TMRQuadrant p = array[i];
      p.level += 1;

      // Get the sibling
      TMRQuadrant q;
      p.getSibling(corner, &q);
      q.cornerNeighbor(corner, &q);

      if ((q.x >= 0 && q.x < hmax) &&
          (q.y >= 0 && q.y < hmax)){
        // Get the MPI rank of the quadrant owner
        int owner = getQuadrantMPIOwner(&q);
        if (owner != mpi_rank){
          TMRQuadrant p = array[i];
          p.tag = owner;
          queue->push(&p);
        }
      }
      else {
        // The node may lie across a corner, edge or face
        int ex = (q.x < 0 || q.x >= hmax);
        int ey = (q.y < 0 || q.y >= hmax);

        if (ex && ey){
          addAdjacentCornerToQueue(corner, q, queue, array[i]);
        }
        else {
          // Add the quadrant to the other trees
          int edge = ex*(q.x < 0 ? 0 : 1) + ey*(q.y < 0 ? 2 : 3);
          addAdjacentEdgeToQueue(edge, q, queue, array[i]);
        }
      }
    }
  }

  // Convert the local adjacency non-local list of quadrants
  TMRQuadrantArray *list = queue->toArray();
  list->getArray(&array, &size);
  qsort(array, size, sizeof(TMRQuadrant), compare_quadrant_tags);

  // Distribute the quadrants
  int use_tags = 1;
  adjacent = distributeQuadrants(list, use_tags);
  adjacent->sort();

  delete list;
}

/*
  Determine if there is an adjacent quadrant on the connecting edge.

  Return true if an adjacent edge is found across a face-edge and
  false if no quadrant is found.
  
  input:
  edge_index:    the local edge index
  b:             the quadrant
*/
int TMRQuadForest::checkAdjacentDepEdges( int edge_index,
                                          TMRQuadrant *b,
                                          TMRQuadrantArray *adjquads ){
  // Get the side length of the quadrant
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - b->level);

  // Store the u coordinate along the edge
  int32_t ucoord = 0;
  if (edge_index < 2){
    ucoord = b->y;
  }
  else {
    ucoord = b->x;
  }

  // Retrieve the first and second node numbers
  int face_owner = b->face;
  int edge = face_edge_conn[4*face_owner + edge_index];
  int n1 = face_conn[4*face_owner + face_to_edge_nodes[edge_index][0]];
  int n2 = face_conn[4*face_owner + face_to_edge_nodes[edge_index][1]];

  // Now, cycle through all the adjacent edges
  for ( int ip = edge_face_ptr[edge]; ip < edge_face_ptr[edge+1]; ip++ ){
    int face = edge_face_conn[ip]/4;

    if (face_owner != face){
      // Get the adjacent edge index
      int adj_index = edge_face_conn[ip] % 4;

      // Get the nodes on the adjacent face
      int nn1 = face_conn[4*face + face_to_edge_nodes[adj_index][0]];
      int nn2 = face_conn[4*face + face_to_edge_nodes[adj_index][1]];

      // Add the quadrant to the list
      int reverse = (n1 == nn2 && n2 == nn1);
      int32_t u = ucoord;
      if (reverse){
        u = hmax - h - ucoord;
      }

      // Search for the neighboring quadrant
      TMRQuadrant quad;
      quad.face = face;
      quad.level = b->level;
      if (adj_index < 2){
        quad.x = (hmax - h)*(adj_index % 2);
        quad.y = u;
      }
      else {
        quad.x = u;
        quad.y = (hmax - h)*(adj_index % 2);
      }
      
      // If the more-refined element exists then label the
      // corresponding nodes as dependent
      TMRQuadrant *dep = quadrants->contains(&quad);
      if (dep){
        // Set the info flag to the corresponding adjacent index
        dep->info |= 1 << adj_index;
      }
      if (adjquads){
        dep = adjquads->contains(&quad);
        if (dep){
          dep->info |= 1 << adj_index;
        }
      }
    }
  }

  return 0;
}

/*
  Compute the dependent nodes (hanging edge/face nodes) on each face
  and on the interfaces between adjacent faces.

  The hanging face nodes may occur on any face within the face and
  are associated with the 4 parallel hanging edge nodes.  Within this
  code, we only store the corresponding face node and use the
  associated edges. Edge nodes along face interfaces may be hanging
  even if there is no associated hanging face node.

  side effects:
  dep_edges:   a list of dependent edges (aligned with face edges)
*/
void TMRQuadForest::computeDepEdges(){
  if (dep_ptr){ 
    delete [] dep_ptr;  
    delete [] dep_conn;
    delete [] dep_weights; 
    dep_ptr = NULL;  
    dep_conn = NULL; 
    dep_weights = NULL; 
  }

  for ( int iter = 0; iter < 2; iter++ ){
    // Get the elements either in the regular quadrant array or
    // in the adjacent element array
    int size = 0;
    TMRQuadrant *array = NULL;
    TMRQuadrantArray *adjquads = NULL;
    if (iter == 0){
      quadrants->getArray(&array, &size);
      adjquads = adjacent;
    }
    else if (adjacent){
      adjacent->getArray(&array, &size);
      adjquads = NULL;
    }

    for ( int i = 0; i < size; i++ ){
      // Get the side length of the element
      const int32_t hmax = 1 << TMR_MAX_LEVEL;
    
      // Enumerate the sibling-ids for each edge
      const int edge_index_to_children[][2] =
        {{0, 2}, {1, 3}, 
         {0, 1}, {2, 3}};
      const int edge_index_to_adjacent[] = {1, 0, 3, 2};

      // Check whether the next-level refined element exists over an
      // adjacent edge
      for ( int edge_index = 0; edge_index < 4; edge_index++ ){
        for ( int k = 0; k < 2; k++ ){
          // Get the quadrant and increase the level
          TMRQuadrant p = array[i];
          p.level += 1;

          // Get the sibling id for each quadrant along the
          // face that we're on right now
          TMRQuadrant q;
          p.getSibling(edge_index_to_children[edge_index][k], &q);

          // Get the edge neighbor
          q.edgeNeighbor(edge_index, &q);

          // Check if the adjacent quadrant q lies over an quadtree edge
          // or face and if so check for the corresponding quadrant
          int fx0 = (q.x < 0);
          int fy0 = (q.y < 0);
          int fx = (fx0 || q.x >= hmax);
          int fy = (fy0 || q.y >= hmax);
          
          if (fx || fy){
            checkAdjacentDepEdges(edge_index, &q, adjquads);
          }
          else {
            TMRQuadrant *dep = quadrants->contains(&q);
            if (dep){
              // Set the info flag to the corresponding adjacent index
              dep->info |= 1 << edge_index_to_adjacent[edge_index];
            }
            if (adjquads){
              dep = adjquads->contains(&q);
              if (dep){
                dep->info |= 1 << edge_index_to_adjacent[edge_index];
              }
            }
          }
        }
      }
    }
  }
}

/*
  Transform the node from a local coordinate system into the global
  node numbering scheme.

  This transforms the given quadrant to the coordinate system of the
  lowest owner face.
  
  input/output:
  quad:  the quadrant representing a node in the local coordinate system
*/
void TMRQuadForest::transformNode( TMRQuadrant *quad, 
                                   int *edge_reversed ){
  // Get the maximum quadrant length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;

  // Check if this node lies on an quadtree boundary
  int fx0 = (quad->x == 0);
  int fy0 = (quad->y == 0);
  int fx = (fx0 || quad->x == hmax);
  int fy = (fy0 || quad->y == hmax);

  // In most cases, the node coordinates will not 
  // be reversed across an edge
  if (edge_reversed){
    *edge_reversed = 0;
  }

  if (fx || fy){
    // Get the original face index
    int face = quad->face;

    if (fx && fy){
      // This node lies on a corner
      int corner = (fx0 ? 0 : 1) + (fy0 ? 0 : 2);
      int node = face_conn[4*face + corner];

      if (face != node_face_owners[face]){
        // Get the pointer information
        int ptr = node_face_ptr[node];
        int adj = node_face_conn[ptr]/4;
        int adj_index = node_face_conn[ptr] % 4;

        // Copy the coordinates of the adjacent face
        quad->face = adj;
        quad->x = hmax*(adj_index % 2);
        quad->y = hmax*(adj_index / 2);
      }
    }
    else {
      // Which edge index are we dealing with?
      int edge_index = fx*(fx0 ? 0 : 1) + fy*(fy0 ? 2 : 3);
      int edge = face_edge_conn[4*face + edge_index];

      if (face != edge_face_owners[edge]){
        // Get the adjacent edge index on the opposite face
        int ptr = edge_face_ptr[edge]; 
        int adj = edge_face_conn[ptr]/4;
        int adj_index = edge_face_conn[ptr] % 4;

        // Retrieve the first and second node numbers to determine the
        // relative orientation between this edge and each adjacent edge
        int n1 = face_conn[4*face + face_to_edge_nodes[edge_index][0]];
        int n2 = face_conn[4*face + face_to_edge_nodes[edge_index][1]];
          
        // Get the orientation
        int nn1 = face_conn[4*adj + face_to_edge_nodes[adj_index][0]];
        int nn2 = face_conn[4*adj + face_to_edge_nodes[adj_index][1]];
              
        // Determine whether the edges are in the same direction
        // or are reversed
        int reverse = (n1 == nn2 && n2 == nn1);

        if (reverse && edge_reversed){
          *edge_reversed = 1;
        }

        // Get the edge coordinate index
        int32_t u = quad->x;
        if (edge_index < 2){
          u = quad->y;
        }
        
        // Set the u-coordinate along the edge
        int32_t uquad = u;
        if (reverse){
          uquad = hmax - u;
        }
      
        // Transform the quadant to the adjacent coordinate system
        quad->face = adj;
        if (adj_index < 2){
          quad->x = hmax*(adj_index % 2);
          quad->y = uquad;
        }
        else {
          quad->x = uquad;
          quad->y = hmax*(adj_index % 2);
        }
      }
    }
      
    // Truncate the node back into the domain if it is on any of the
    // outer boundaries
    if (quad->x == hmax){ quad->x = hmax-1; }
    if (quad->y == hmax){ quad->y = hmax-1; }
  }
}

/*
  Label the dependent face and edge nodes

  This code is called after all the dependent faces have been
  computed.  Note that this relies on the mesh being edge-balanced
  (which is required).
*/
void TMRQuadForest::labelDependentNodes( int *nodes ){
  int size = 0;
  TMRQuadrant *array = NULL;
  quadrants->getArray(&array, &size);

  // Scan through all the quadrants and label the edges that are
  // dependent
  for ( int i = 0; i < size; i++ ){
    // This element contains dependent node information
    if (array[i].info){
      for ( int edge_index = 0; edge_index < 4; edge_index++ ){
        // Find whether this edge is dependent or not
        if (array[i].info & 1 << edge_index){
          // Get the element number
          int num = array[i].tag;

          // Compute the offset into the array based on the 
          // edge index
          int offset = 0;
          if (edge_index < 2){
            offset = (edge_index % 2)*(mesh_order-1);
          }
          else {
            offset = (edge_index % 2)*(mesh_order-1)*mesh_order;
          }

          // Compute the increment for the local element
          int incr = 1;
          if (edge_index < 2){
            incr = mesh_order;
          }

          // Set the start/end point depending on the childId()
          int id = array[i].childId();
          int start, end;
          if ((edge_index < 2 && id/2 == 0) ||
              (edge_index >= 2 && id % 2 == 0)){
            start = 1;
            end = mesh_order - (mesh_order % 2);
          }
          else {
            start = (mesh_order % 2);
            end = mesh_order-1;
          }

          // Label the dependent node
          if (interp_type == TMR_UNIFORM_POINTS){
            int indx = mesh_order*mesh_order*num + offset;
            for ( int k = start; k < end; k += 2 ){
              nodes[conn[indx + k*incr]] = -1;
            }
          }
          else {
            // Otherwise, all the dependent nodes along the edge,
            // except the first and possibly last are dependent
            int indx = mesh_order*mesh_order*num + offset;
            for ( int k = start; k < end; k++ ){
              nodes[conn[indx + k*incr]] = -1;
            }
          }
        }
      }
    }
  }
}

/*
  Get the corner owner index
*/
int TMRQuadForest::getCornerOwnerIndex( TMRQuadrant *owner, 
                                        TMRQuadrant *node ){
  int32_t h = 1 << (TMR_MAX_LEVEL - owner->level);

  for ( int jj = 0; jj < 2; jj++ ){
    for ( int ii = 0; ii < 2; ii++ ){
      TMRQuadrant n;
      n.face = owner->face;
      n.x = owner->x + ii*h;
      n.y = owner->y + jj*h;
      transformNode(&n);

      // Compare the positions of the transformed node
      if (n.comparePosition(node) == 0){
        return mesh_order*mesh_order*owner->tag + ii + jj*mesh_order;
      }
    }
  }

  return -1;
}

/*
  Get the edge owner index
*/
void TMRQuadForest::getEdgeOwnerIndex( TMRQuadrant *owner, 
                                       TMRQuadrant *edge, 
                                       int *owner_edge_index ){
  int32_t h = 1 << (TMR_MAX_LEVEL - owner->level - 1);

  owner_edge_index[0] = -1;
  owner_edge_index[mesh_order-1] = -1;

  for ( int edge_index = 0; edge_index < 4; edge_index++ ){
    TMRQuadrant e;
    e.face = owner->face;
    if (edge_index < 2){
      e.x = owner->x + 2*h*(edge_index % 2);
      e.y = owner->y + h;
    }
    else {
      e.x = owner->x + h;
      e.y = owner->y + 2*h*(edge_index % 2);
    }
    transformNode(&e);

    // Compare the positions of the transformed node
    if (e.comparePosition(edge) == 0){
      if (edge_index < 2){
        for ( int k = 1; k < mesh_order-1; k++ ){
          owner_edge_index[k] = 
            mesh_order*mesh_order*owner->tag +
            k*mesh_order + (mesh_order-1)*(edge_index % 2);
        }
      }
      else {
        for ( int k = 1; k < mesh_order-1; k++ ){
          owner_edge_index[k] = 
            mesh_order*mesh_order*owner->tag +
            k + mesh_order*(mesh_order-1)*(edge_index % 2);
        }
      }
    }
  }
}

/*
  Convert from the integer coordinate system to a physical coordinate
  with the off-by-one check.
*/
double convertToCoordinate( const int32_t x ){
  static const int32_t hmax = 1 << TMR_MAX_LEVEL;
  if (x == 0){
    return 0.0;
  }
  else if (x == hmax-1){
    return 1.0;
  }
  else {
    return 1.0*x/hmax;
  }
}

/*
  Create the nodes from the element mesh

  Note that the element mesh must be balanced before the nodes can be
  ordered.

  This function first computes the face that owns of each of the
  faces, edges and corners(/nodes) within the super-mesh. Next, the
  the non-local quadrants that border each quadtree are passed back to the
  processors which own the quadtree. This creates a layer of quadrants
  that are temporarily stored in a partial quadtree. The code then
  creates nodes within each quadrant element for all of the quadtrees
  (including the partial quadtrees). Next, the dependent nodes (that are
  hanging on a face or an edge of an element)are labeled according to
  whether they are on an edge or face.
  
  After the dependent nodes are labeled, the nodes are ordered on all
  processors to form a complete global ordering. Next, the nodes are
  communicated locally across quadrant and partial quadrant faces, edges
  and corners. Finally, the new node numbers are returned to the
  processors that border the quadtree owners. And lastly, the non-local
  partial quadtrees are freed.

  input:
  order:   the order of the mesh
*/
void TMRQuadForest::createNodes(){
  // Send/recv the adjacent quadrants
  computeAdjacentQuadrants();

  // Compute the dependent face nodes
  computeDepEdges();
  
  // Allocate the array of elements
  int num_elements;
  TMRQuadrant *quads;
  quadrants->getArray(&quads, &num_elements);

  // Create all the nodes/edges/faces
  int max_nodes = 4*num_elements;
  if (mesh_order >= 3){
    max_nodes = 9*num_elements;
  }
  const int use_node_index = 1;
  TMRQuadrantHash *all_nodes = new TMRQuadrantHash(use_node_index);

  // Set the node, edge and face label
  int node_label = 0, edge_label = 0, face_label = 0;
  if (mesh_order >= 4){
    node_label = 0;
    edge_label = 1;
    face_label = 2;
  }

  // Set the node locations
  int index = 0;
  if (mesh_order == 2){
    for ( int i = 0; i < num_elements; i++ ){
      const int32_t h = 1 << (TMR_MAX_LEVEL - quads[i].level);
      for ( int jj = 0; jj < 2; jj++ ){
        for ( int ii = 0; ii < 2; ii++ ){
          TMRQuadrant node;
          node.face = quads[i].face;
          node.level = 0;
          node.x = quads[i].x + h*ii;
          node.y = quads[i].y + h*jj;
          node.tag = i;
          node.info = node_label;
          transformNode(&node);
          all_nodes->addQuadrant(&node);
        }
      }
    }
  }
  else {
    for ( int i = 0; i < num_elements; i++ ){
      const int32_t h = 1 << (TMR_MAX_LEVEL - quads[i].level - 1);
      for ( int jj = 0; jj < 3; jj++ ){
        for ( int ii = 0; ii < 3; ii++ ){
          TMRQuadrant node;
          node.face = quads[i].face;
          node.level = 0;
          node.x = quads[i].x + h*ii;
          node.y = quads[i].y + h*jj;
          node.tag = i;
          if ((ii == 0 || ii == 2) &&
              (jj == 0 || jj == 2)){
            node.info = node_label;
          }
          else if (ii == 0 || ii == 2 ||
                   jj == 0 || jj == 2){
            node.info = edge_label;
          }
          else {
            node.info = face_label;
          }
          transformNode(&node);
          all_nodes->addQuadrant(&node);
        }
      }
    }
  }

  // Create the array of nodes
  TMRQuadrantArray *nodes = all_nodes->toArray();
  nodes->sort();

  // Allocate the connectivity
  int size = mesh_order*mesh_order*num_elements;
  conn = new int[ size ];
  memset(conn, 0, size*sizeof(int));

  if (mesh_order <= 3){
    TMRQuadrant *array;
    nodes->getArray(&array, &num_local_nodes);

    for ( int i = 0; i < num_elements; i++ ){
      int *c = &conn[mesh_order*mesh_order*i];
      const int32_t h = 1 << (TMR_MAX_LEVEL - quads[i].level - 1);

      // Loop over the element nodes
      for ( int corner_index = 0; corner_index < 4; corner_index++ ){
        TMRQuadrant node;
        node.face = quads[i].face;
        node.x = quads[i].x + 2*h*(corner_index % 2);
        node.y = quads[i].y + 2*h*(corner_index / 2);
        transformNode(&node);
        TMRQuadrant *t = nodes->contains(&node);
        int offset = (mesh_order-1)*(corner_index % 2) +
                     (mesh_order-1)*mesh_order*(corner_index/2);
        c[offset] = t - array;
      }

      if (mesh_order == 3){
        // Loop over the edges and get the owners
        for ( int edge_index = 0; edge_index < 4; edge_index++ ){
          TMRQuadrant node;
          node.face = quads[i].face;
          if (edge_index < 2){
            node.x = quads[i].x + 2*h*(edge_index % 2);
            node.y = quads[i].y + h;
          }
          else {
            node.x = quads[i].x + h;
            node.y = quads[i].y + 2*h*(edge_index % 2);
          }
          transformNode(&node);
          TMRQuadrant *t = nodes->contains(&node);
          if (edge_index < 2){
            int offset = mesh_order + (mesh_order-1)*edge_index;
            c[offset] = t - array;
          }
          else {
            int offset = 1 + (mesh_order-1)*mesh_order*(edge_index % 2);
            c[offset] = t - array;
          }
        }

        TMRQuadrant node;
        node.face = quads[i].face;
        node.x = quads[i].x + h;
        node.y = quads[i].y + h;
        transformNode(&node);
        TMRQuadrant *t = nodes->contains(&node);
        c[4] = t - array;
      }
    }
  }
  else {
    // Get the edge indices
    int *owner_edge_index = new int[ mesh_order ];

    // Loop over all the elements and assign the local index owners
    // for each node
    for ( int i = 0; i < num_elements; i++ ){
      int *c = &conn[mesh_order*mesh_order*i];

      // Compute the half-edge length of the quadrant
      const int32_t h = 1 << (TMR_MAX_LEVEL - quads[i].level - 1);

      // Loop over the nodes and set the
      for ( int corner_index = 0; corner_index < 4; corner_index++ ){
        // Compute the offset to the local node
        int offset = (mesh_order-1)*(corner_index % 2) +
                     (mesh_order-1)*mesh_order*(corner_index/2);

        // Find the node at the corner to determine the owner
        TMRQuadrant node;
        node.face = quads[i].face;
        node.info = node_label;
        node.x = quads[i].x + 2*h*(corner_index % 2);
        node.y = quads[i].y + 2*h*(corner_index / 2);
        transformNode(&node);
        TMRQuadrant *t = nodes->contains(&node);

        // Get the owner of this node and label all the
        // nodes that touch it
        int owner = t->tag;
        int index = 0;
        if (owner == i){
          index = mesh_order*mesh_order*owner + offset;
        }
        else {
          index = getCornerOwnerIndex(&quads[owner], &node);
        }
        c[offset] = -index - 1;
      }

      // Loop over the edges and get the owners
      for ( int edge_index = 0; edge_index < 4; edge_index++ ){
        TMRQuadrant edge;
        edge.face = quads[i].face;
        edge.info = edge_label;
        if (edge_index < 2){
          edge.x = quads[i].x + 2*h*(edge_index % 2);
          edge.y = quads[i].y + h;
        }
        else {
          edge.y = quads[i].y + 2*h*(edge_index % 2);
          edge.x = quads[i].x + h;
        }
        int edge_reversed = 0;
        transformNode(&edge, &edge_reversed);
        TMRQuadrant *t = nodes->contains(&edge);

        // Get the owner of this edge and label all the
        // nodes that touch it
        int owner = t->tag;
        if (owner == i){
          if (edge_index < 2){
            for ( int k = 1; k < mesh_order-1; k++ ){
              int offset = k*mesh_order + (mesh_order-1)*edge_index;
              int index = mesh_order*mesh_order*owner + offset;
              c[offset] = -index - 1;
            }
          }
          else {
            for ( int k = 1; k < mesh_order-1; k++ ){
              int offset = k + (mesh_order-1)*mesh_order*(edge_index % 2);
              int index = mesh_order*mesh_order*owner + offset;
              c[offset] = -index - 1;
            }
          }
        }
        else {
          getEdgeOwnerIndex(&quads[owner], &edge, owner_edge_index);

          // The owner edge is reversed relative to this edge
          if (edge_reversed){
            if (edge_index < 2){
              for ( int k = 1; k < mesh_order-1; k++ ){
                int offset = k*mesh_order + (mesh_order-1)*edge_index;
                c[offset] = -owner_edge_index[mesh_order-1-k] - 1;
              }
            }
            else {
              for ( int k = 1; k < mesh_order-1; k++ ){
                int offset = k + (mesh_order-1)*mesh_order*(edge_index % 2);
                c[offset] = -owner_edge_index[mesh_order-1-k] - 1;
              }
            }
          }
          else {
            if (edge_index < 2){
              for ( int k = 1; k < mesh_order-1; k++ ){
                int offset = k*mesh_order + (mesh_order-1)*edge_index;
                c[offset] = -owner_edge_index[k] - 1;
              }
            }
            else {
              for ( int k = 1; k < mesh_order-1; k++ ){
                int offset = k + (mesh_order-1)*mesh_order*(edge_index % 2);
                c[offset] = -owner_edge_index[k] - 1;
              }
            }
          }
        }
      }

      // Loop over the face owners
      for ( int jj = 1; jj < mesh_order-1; jj++ ){
        for ( int ii = 1; ii < mesh_order-1; ii++ ){
          int offset = ii + jj*mesh_order;
          int index = mesh_order*mesh_order*i + offset;
          c[offset] = -index - 1;
        }
      }
    }

    // Free the mesh indices
    delete [] owner_edge_index;

    // Set the initial number of local nodes/dependent nodes
    num_local_nodes = 0; // Including dependent nodes
    num_owned_nodes = 0;
    num_dep_nodes = 0;

    // Assign the local node numbers
    for ( int i = 0, index = 0; i < num_elements; i++ ){
      int *c = &conn[mesh_order*mesh_order*i];

      // Assign the corner nodes
      for ( int ii = 0; ii < mesh_order; ii++ ){
        for ( int jj = 0; jj < mesh_order; jj++ ){
          // Get the local element node index
          int offset = ii + jj*mesh_order;
          int index = mesh_order*mesh_order*i + offset;

          // If the offset is the owner, then label it right away
          if (c[offset] == -index - 1){
            c[offset] = num_local_nodes;
            num_local_nodes++;
          }
          else if (c[offset] < 0){
            // If the owner has not been ordered, order it now
            int ptr = -c[offset] - 1;
            if (conn[ptr] < 0){
              conn[ptr] = num_local_nodes;
              num_local_nodes++;
            }
            c[offset] = conn[ptr];
          }
        }
      }
    }
  }

  // Set the locally owned numbers
  /*
  // Create an array of all the independent nodes that are
  // owned by other processors and referenced by the
  // elements on this processor.
  TMRQuadrantArray *ext_array = ext_hash->toArray();
  delete ext_hash;
  ext_array->sort();

  // Distribute the non-local nodes back to their owning processors
  // to determine their node numbers
  int use_tags = 0;
  int *send_ptr, *recv_ptr;
  TMRQuadrantArray *dist = distributeQuadrants(ext_array, use_tags,
                                               &send_ptr, &recv_ptr);
  delete ext_array;

  // Loop over the off-processor nodes and search for them in the
  // sorted node list and assign them the correct tag
  int size;
  TMRQuadrant *array;
  dist->getArray(&array, &size);
  for ( int i = 0; i < size; i++ ){
    // Read the data out of the connectivity array
    // const int use_node_search = 1;
    // TMRQuadrant *t = nodes->contains(&array[i], use_node_search);
    // array[i].tag = t->tag;
  }

  // Send the nodes back to the original processors
  TMRQuadrantArray *ext_nodes = sendQuadrants(dist, recv_ptr, send_ptr);
  delete dist;
  delete [] recv_ptr;
  delete [] send_ptr;
  */

  // Retrieve the non-local nodes and set the external node numbers

  // Set the local node numbers

  // Allocate the array of locally owned nodes
  X = new TMRPoint[ num_local_nodes + num_dep_nodes ];
  memset(X, 0, (num_local_nodes + num_dep_nodes)*sizeof(TMRPoint));

  if (topo){
    for ( int i = 0; i < num_elements; i++ ){
      // Get the right surface
      TMRFace *surf;
      topo->getFace(quads[i].face, &surf);

      // Compute the edge length
      const int32_t h = 1 << (TMR_MAX_LEVEL - quads[i].level);

      // Compute the origin of the element in parametric space
      // and the edge length of the element
      double d = convertToCoordinate(h);
      double u = convertToCoordinate(quads[i].x);
      double v = convertToCoordinate(quads[i].y);

      // Look for nodes that are not assigned
      for ( int jj = 0; jj < mesh_order; jj++ ){
        for ( int ii = 0; ii < mesh_order; ii++ ){
          // Compute the mesh index
          int index = conn[mesh_order*mesh_order*i + 
                           ii + jj*mesh_order];
          if (index >= 0 && index < num_local_nodes){
            surf->evalPoint(u + d*interp_knots[ii], 
                            v + d*interp_knots[jj], &X[index]);
          }
          else if (index < 0 && index >= -num_dep_nodes){
            index = num_local_nodes - index-1;
            surf->evalPoint(u + d*interp_knots[ii], 
                            v + d*interp_knots[jj], &X[index]);            
          }
        }
      }
    }
  }

  // Write out the connectivity to a file
  if (mpi_rank == 0 && topo){
    FILE *fp = fopen("output_result.dat", "w");
    if (fp){
      int n = num_local_nodes + num_dep_nodes;
      fprintf(fp, "Variables = X,Y,Z,var\n");
      fprintf(fp, "Zone N = %d E = %d ", n,
              (mesh_order-1)*(mesh_order-1)*num_elements);
      fprintf(fp, "DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL\n");

      // Write out the points
      for ( int k = 0; k < n; k++ ){
        fprintf(fp, "%e\n", X[k].x);
      }
      for ( int k = 0; k < n; k++ ){
        fprintf(fp, "%e\n", X[k].y);
      }
      for ( int k = 0; k < n; k++ ){
        fprintf(fp, "%e\n", X[k].z);
      }
      for ( int k = 0; k < n; k++ ){
        fprintf(fp, "%d\n", k);
      }

      // Write out the connectivity
      for ( int k = 0; k < num_elements; k++ ){
        int *c = &conn[mesh_order*mesh_order*k];
        for ( int jj = 0; jj < mesh_order-1; jj++ ){
          for ( int ii = 0; ii < mesh_order-1; ii++ ){
            int n[4];
            n[0] = c[ii + jj*mesh_order];
            n[1] = c[ii+1 + jj*mesh_order];
            n[2] = c[ii+1 + (jj+1)*mesh_order];
            n[3] = c[ii + (jj+1)*mesh_order];

            for ( int j = 0; j < 4; j++ ){
              if (n[j] < 0 && n[j] >= -num_dep_nodes){
                n[j] = 1; // num_local_nodes - n[j];
              }
              else if (n[j] >= 0){
                n[j] += 1;
              }
              else {
                n[j] = 1;
              }
            }

            fprintf(fp, "%d %d %d %d\n", n[0], n[1], n[2], n[3]);
          }
        }
      }
      fclose(fp);
    }
  }
}

/*
  Get the elements that either lie on a face or curve with a given
  attribute.

  This code loops over all quadrants owned locally on this processor
  and checks if each quadrant lies on a face or boundary. If the face
  attribute matches, the quadrant is added without modification. If
  the quadrant lies on an edge, the quadrant is modified so that the
  tag indicates which edge the quadrant lies on using the regular edge
  ordering scheme.

  input:
  attr:   string attribute associated with the geometric feature

  returns:
  list:   an array of quadrants satisfying the attribute
*/
TMRQuadrantArray* TMRQuadForest::getQuadsWithAttribute( const char *attr ){
  if (!topo){
    return NULL;
  }

  // Create a queue to store the elements that we find
  TMRQuadrantQueue *queue = new TMRQuadrantQueue();

  // Get the quadrants
  int size;
  TMRQuadrant *array;
  quadrants->getArray(&array, &size);

  // Loop over the quadrants and find out whether it touches
  // a face or edge with the prescribed attribute
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  for ( int i = 0; i < size; i++ ){
    const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level);

    // Get the surface quadrant
    TMRFace *surf;
    topo->getFace(array[i].face, &surf);  
    const char *face_attr = surf->getAttribute();
    if ((!attr && !face_attr) ||
        (face_attr && strcmp(face_attr, attr) == 0)){
      queue->push(&array[i]);
    }
    else {
      // If this quadrant was not added from a face
      // attribute, check to see if it should be added
      // as an edge/curve attribute
      TMREdge *edge;
      if (array[i].x == 0){
        int edge_num = face_edge_conn[4*array[i].face];
        topo->getEdge(edge_num,&edge);
        const char *edge_attr = edge->getAttribute();
        if (edge_attr && strcmp(edge_attr, attr) == 0){
          TMRQuadrant p = array[i];
          p.tag = 0;
          queue->push(&p);
        }
      }
      if (array[i].x+h == hmax){
        int edge_num = face_edge_conn[4*array[i].face+1];
        topo->getEdge(edge_num, &edge);
        const char *edge_attr = edge->getAttribute();
        if (edge_attr && strcmp(edge_attr, attr) == 0){
          TMRQuadrant p = array[i];
          p.tag = 1;
          queue->push(&p);
        }
      }
      if (array[i].y == 0){
        int edge_num = face_edge_conn[4*array[i].face+2];
        topo->getEdge(edge_num, &edge);
        const char *edge_attr = edge->getAttribute();
        if (edge_attr && strcmp(edge_attr, attr) == 0){
          TMRQuadrant p = array[i];
          p.tag = 2;
          queue->push(&p);
        }
      }
      if (array[i].y+h == hmax){
        int edge_num = face_edge_conn[4*array[i].face+3];
        topo->getEdge(edge_num, &edge);
        const char *edge_attr = edge->getAttribute();
        if (edge_attr && strcmp(edge_attr, attr) == 0){
          TMRQuadrant p = array[i];
          p.tag = 3;
          queue->push(&p);
        }
      }
    }
  }

  TMRQuadrantArray *list = queue->toArray();
  delete queue;
  return list;
}

/*
  Create an array of the nodes that are lie on a surface, edge or
  corner with a given attribute

  This code loops over all nodes and check whether they lie on a
  geometric entity that has the given attribute. The nodes are not
  unique if they are lie on a shared boundary between processors.

  input:
  attr:   the string of the attribute to search

  returns:
  list:   the nodes matching the specified attribute
*/
int TMRQuadForest::getNodesWithAttribute( const char *attr,
                                          int **_nodes ){
  if (!topo){
    return NULL;
  }

   *_nodes = NULL;
   return 0;
}

/*
  Get the dependent connectivity information (create it if it has not
  been allocated previously).

  output:
  ptr:      pointer for each dependent node number
  conn:     connectivity to each (global) independent node
  weights:  the weight values for each dependent node
*/
int TMRQuadForest::getDepNodeConn( const int **ptr, const int **conn,
                                   const double **weights ){
  if (ptr){ *ptr = dep_ptr; }
  if (conn){ *conn = dep_conn; }
  if (weights){ *weights = dep_weights; }
  return num_dep_nodes;
}

/*
  Create the dependent mesh information for all local dependent
  nodes.

  output:
  ptr:      pointer for each dependent node number
  conn:     connectivity to each (global) independent node
  weights:  the weight values for each dependent node
*/
void TMRQuadForest::createDepNodeConn( int **_ptr, int **_conn,
                                       double **_weights ){}

/*
  Given a node, find the enclosing quadrant

  This code is used to find the quadrant in the quadrant array that
  encloses the given node.
*/
TMRQuadrant* TMRQuadForest::findEnclosing( TMRQuadrant *node ){
  // Retrieve the array of elements
  int size = 0;
  TMRQuadrant *array = NULL;
  quadrants->getArray(&array, &size);

  // Set the lower and upper bounds for the quadrant
  const int32_t face = node->face;
  const int32_t x = node->x;
  const int32_t y = node->y;

  // Set the low and high indices to the first and last
  // element of the element array
  int low = 0;
  int high = size-1;
  int mid = low + (high - low)/2;

  // Maintain values of low/high and mid such that the
  // quadrant is between (elems[low], elems[high]).
  // Note that if high-low=1, then mid = high
  while (high != mid){
    // Check if array[mid] contains the provided quadrant
    const int32_t h = 1 << (TMR_MAX_LEVEL - array[mid].level);
    if ((array[mid].face == face) &&
        (array[mid].x <= x && x <= array[mid].x+h) &&
        (array[mid].y <= y && y <= array[mid].y+h)){
      return &array[mid];
    }
    
    // Compare the ordering of the two quadrants - if the
    // quadrant is less than the other, then adjust the mid point 
    if (node->comparePosition(&array[mid]) < 0){
      high = mid-1;
    } 
    else {
      low = mid+1;
    }
    
    // Re compute the mid-point and repeat
    mid = high - (int)((high - low)/2);
  }

  // Check if array[mid] contains the provided quadrant
  const int32_t h1 = 1 << (TMR_MAX_LEVEL - array[mid].level);
  if ((array[mid].face == face) &&
      (array[mid].x <= x && x <= array[mid].x+h1) &&
      (array[mid].y <= y && y <= array[mid].y+h1)){
    return &array[mid];
  }

  // Check if elems[low] contains the provided quadrant
  const int32_t h2 = 1 << (TMR_MAX_LEVEL - array[low].level);
  if ((array[low].face == face) &&
      (array[low].x <= x && x <= array[low].x+h2) &&
      (array[low].y <= y && y <= array[low].y+h2)){
   return &array[low];
 }

 // No quadrant was found, return NULL
 return NULL;
}

/*
  Create the interpolation operator from the coarse to the fine mesh.

  Each processor builds the interpolation for the locally owned nodes
  in the mesh. This interpolation will refer to other nodes, but the
  output is local.
  
  input:
  coarse:   the coarse quadtree forest that has the same layout as this
  
  output:
  ptr:      the pointer into the local rows
  conn:     the connectivity using global numbers
  weights:  the interpolation weights for each point
*/
void TMRQuadForest::createInterpolation( TMRQuadForest *coarse,
                                         TACSBVecInterp *interp ){
  /*
  // Create the dependent node connectivity on the coarse mesh
  coarse->createDepNodeConn();

  // Get the dependent node information
  const int *cdep_ptr, *cdep_conn;
  const double *cdep_weights;
  coarse->getDepNodeConn(&cdep_ptr, &cdep_conn, &cdep_weights);

  // Create the array of local nodes on the fine mesh. The original
  // node array contains duplicates along processor boundaries.
  int node_size;
  TMRQuadrant *node_array;
  nodes->getArray(&node_array, &node_size);
  TMRQuadrant *local_array = new TMRQuadrant[ node_size ];

  // Copy only the nodes that are locally owned
  int count = 0;
  for ( int i = 0; i < node_size; i++ ){
    if (node_array[i].tag >= node_range[mpi_rank] &&
        node_array[i].tag < node_range[mpi_rank+1]){
      local_array[count] = node_array[i];
      count++;
    }
  }

  // Allocate the quadrant array and distribute it across all 
  // processors on the coarse mesh
  TMRQuadrantArray *local = new TMRQuadrantArray(local_array, count);

  // Distribute the quadrants to the owners - include the local
  // quadrants in the new array since everything has to be
  // interpolated
  int use_tags = 0; // Use the quadrant ownership to distribute (not the tags)
  int include_local = 1; // Include the locally owned quadrants
  TMRQuadrantArray *fine_nodes = 
    coarse->distributeQuadrants(local, use_tags, NULL, NULL, include_local);
  delete local;

  // Get the number of locally owned nodes on this processor
  int fine_size;
  TMRQuadrant *fine;
  fine_nodes->getArray(&fine, &fine_size);

  // Loop over the nodes in the fine mesh that are owned by quadrants
  // on the coarse mesh stored on this processor
  for ( int i = 0; i < fine_size; i++ ){
    // Find the quadrant that encloses the node - this is not unique,
    // but does produce a unique interpolation (since
    // edges/face/corners will be treated the same when adjacent
    // quadrants that both share a common node location touch)
    TMRQuadrant *quad = coarse->findEnclosing(&fine[i]);
    
    if (quad){
      // The maximum possible size of the array of weights. Note that
      // this is found if every node is a dependent node (which is
      // impossible) which points to a dependent face node (also
      // impossible). It is an upper bound.
      const int max_size = (4*4*4)*(4*4);
      TMRIndexWeight weights[max_size];

      // Get the element size for coarse element
      const int32_t h = 1 << (TMR_MAX_LEVEL - quad->level);

      // Compute the parametric location within the element. Note that
      // this modifies the end-point locations so that they are flush
      // with the upper limit
      int32_t u = (fine[i].x == hmax-1 ? hmax : fine[i].x) - quad->x;
      int32_t v = (fine[i].y == hmax-1 ? hmax : fine[i].y) - quad->y;

      // Set the base node location
      int32_t x = quad->x + (u == h ? h : 0);
      int32_t y = quad->y + (v == h ? h : 0);


      // Compute the interpolation weights
      double Nu[4], Nv[4];
      int nu = computeInterpWeights(coarse->mesh_order, u, h, Nu);
      int nv = computeInterpWeights(coarse->mesh_order, v, h, Nv);
    
      // Loop over the nodes that are within this quadrant
      int nweights = 0;
      for ( int jj = 0; jj < nv; jj++ ){
        for ( int ii = 0; ii < nu; ii++ ){
          // Compute the interpolation weight
          double weight = Nu[ii]*Nv[jj];
          
          // Set the node locations
          TMRQuadrant node;
          node.face = quad->face;
          node.x = x + hc*ii;
          node.y = y + hc*jj;
            
          // Transform the node using the coarse transform
          coarse->transformNode(&node);
          
          // Find the coarse mesh
          int use_node_search = 1;
          TMRQuadrant *t = coarse->nodes->contains(&node, use_node_search);
          
          if (t->tag >= 0){
            weights[nweights].index = t->tag;
            weights[nweights].weight = weight;
            nweights++;
          }
          else {
            // Unravel the dependent node connectivity
            int node = -t->tag-1;
            for ( int jp = cdep_ptr[node]; jp < cdep_ptr[node+1]; jp++ ){
              weights[nweights].index = cdep_conn[jp];
              weights[nweights].weight = weight*cdep_weights[jp];
              nweights++;
            }
          }
        }
      }

      // Sort the dependent weight values
      nweights = TMRIndexWeight::uniqueSort(weights, nweights);

      // The interpolation variables/weights on the coarse mesh
      int vars[49];
      double wvals[49];
      for ( int k = 0; k < nweights; k++ ){
        vars[k] = weights[k].index;
        wvals[k] = weights[k].weight;
      }

      // Add the weights/indices to the interpolation object
      interp->addInterp(fine[i].tag, wvals, vars, nweights);
    }
  }

  delete fine_nodes;
  */
}

/*
  Create a sorted, unique array of the external node numbers that are
  referenced on this processor, but are not local.

  This function can be used to determine a local order for the nodes
  on this processor.

  in/out:
  ext_nodes:   the external nodes
  
  returns:     the number of external nodes
*/
int TMRQuadForest::getExtNodeNums( int **_ext_nodes ){
  /*
  // Determine the number of fine nodes
  int node_size;
  TMRQuadrant *node_array;
  nodes->getArray(&node_array, &node_size);

  // Number of local nodes
  int local_size = node_range[mpi_rank+1] - node_range[mpi_rank];
    
  // The maximum number of external nodes
  int max_ext_nodes = node_size - local_size;
  int num_ext = 0;
  int *ext_nodes = new int[ max_ext_nodes ];

  // Scan through and add any external, independent nodes
  for ( int i = 0; i < node_size; i++ ){
    if (node_array[i].tag >= 0 &&
        (node_array[i].tag < node_range[mpi_rank] ||
         node_array[i].tag >= node_range[mpi_rank+1])){
      ext_nodes[num_ext] = node_array[i].tag;
      num_ext++;
    }
  }

  // Sort the array using quicksort
  qsort(ext_nodes, num_ext, sizeof(int), compare_integers);

  // Create a smaller array to store the result
  *_ext_nodes = new int[ num_ext ];
  memcpy(*_ext_nodes, ext_nodes, num_ext*sizeof(int));

  // Free the larger array
  delete [] ext_nodes;

  return num_ext;
  */
}

