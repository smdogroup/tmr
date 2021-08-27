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

#include "TMRQuadForest.h"
#include "TMRInterpolation.h"
#include "tmrlapack.h"
#include <stdlib.h>

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
  Convert from the integer coordinate system to a physical coordinate
  with the off-by-one check.
*/
static double convert_to_coordinate( const int32_t x ){
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
  Create the TMRQuadForest object
*/
TMRQuadForest::TMRQuadForest( MPI_Comm _comm, int _mesh_order,
                              TMRInterpolationType _interp_type ){
  // Initialize the TMR-specific MPI data types
  if (!TMRIsInitialized()){
    TMRInitialize();
  }

  // Set the MPI communicator
  comm = _comm;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // Set default mesh data
  mesh_order = 2;
  interp_knots = NULL;

  // Set the topology object to NULL
  topo = NULL;

  // Null out the face data
  fdata = NULL;

  // Null the quadrant owners/quadrant list
  owners = NULL;
  quadrants = NULL;
  adjacent = NULL;
  X = NULL;

  // Set data for the number of elements/nodes/dependents
  conn = NULL;
  node_numbers = NULL;
  node_range = NULL;
  num_local_nodes = 0;
  num_owned_nodes = 0;
  num_dep_nodes = 0;
  ext_pre_offset = 0;

  dep_ptr = NULL;
  dep_conn = NULL;
  dep_weights = NULL;

  // Set the mesh order
  setMeshOrder(_mesh_order, _interp_type);
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
  if (fdata){
    fdata->decref();
  }
  fdata = NULL;

  // Free the quadrants/adjacency
  if (owners){ delete [] owners; }
  if (quadrants){ delete quadrants; }
  if (adjacent){ delete adjacent; }
  if (X){ delete [] X; }

  if (conn){ delete [] conn; }
  if (node_numbers){ delete [] node_numbers; }
  if (node_range){ delete [] node_range; }
  if (dep_ptr){ delete [] dep_ptr; }
  if (dep_conn){ delete [] dep_conn; }
  if (dep_weights){ delete [] dep_weights; }

  // Null the quadrant owners/quadrant list
  owners = NULL;
  quadrants = NULL;
  adjacent = NULL;
  X = NULL;

  // Set data for the number of elements/nodes/dependents
  conn = NULL;
  node_numbers = NULL;
  node_range = NULL;
  num_local_nodes = 0;
  num_owned_nodes = 0;
  num_dep_nodes = 0;
  ext_pre_offset = 0;

  dep_ptr = NULL;
  dep_conn = NULL;
  dep_weights = NULL;
}

/*
  Free the data associated with the mesh
*/
void TMRQuadForest::freeMeshData( int free_quads,
                                  int free_owners ){
  if (free_quads){
    if (quadrants){ delete quadrants; }
    quadrants = NULL;
  }
  if (free_owners){
    if (owners){ delete [] owners; }
    owners = NULL;
  }

  // Free any data associated with the mesh
  if (adjacent){ delete adjacent; }
  if (X){ delete [] X; }

  if (conn){ delete [] conn; }
  if (node_numbers){ delete [] node_numbers; }
  if (node_range){ delete [] node_range; }
  if (dep_ptr){ delete [] dep_ptr; }
  if (dep_conn){ delete [] dep_conn; }
  if (dep_weights){ delete [] dep_weights; }

  // Reset the data
  adjacent = NULL;
  X = NULL;

  // Set data for the number of elements/nodes/dependents
  conn = NULL;
  node_numbers = NULL;
  node_range = NULL;
  dep_ptr = NULL;
  dep_conn = NULL;
  dep_weights = NULL;

  // Set the data to NULL
  num_local_nodes = 0;
  num_owned_nodes = 0;
  num_dep_nodes = 0;
  ext_pre_offset = 0;
}

/*
  Copy the connectivity data, but not the quadrants/nodes
*/
void TMRQuadForest::copyData( TMRQuadForest *copy ){
  // Copy over the connectivity data
  fdata->incref();
  if (copy->fdata){
    copy->fdata->decref();
  }
  copy->fdata = fdata;

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
  freeData();

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

  // Create the new face conn data
  fdata = new TMRFaceConn();
  fdata->incref();

  // Copy over the data locally
  fdata->num_nodes = _num_nodes;
  fdata->num_edges = 0;
  fdata->num_faces = _num_faces;

  // Copy over the face connectivity
  fdata->face_conn = new int[ 4*_num_faces ];
  memcpy(fdata->face_conn, _face_conn, 4*_num_faces*sizeof(int));

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

  // Create the new face conn data
  fdata = new TMRFaceConn();
  fdata->incref();

  // Copy over the number of geometric entities
  fdata->num_nodes = _num_nodes;
  fdata->num_edges = _num_edges;
  fdata->num_faces = _num_faces;

  // Copy over the face connectivity
  fdata->face_conn = new int[ 4*_num_faces ];
  memcpy(fdata->face_conn, _face_conn, 4*_num_faces*sizeof(int));

  // Compute the node to face information
  computeNodesToFaces();

  // Copy over the edge information
  fdata->face_edge_conn = new int[ 4*_num_faces ];
  memcpy(fdata->face_edge_conn, _face_edge_conn, 4*_num_faces*sizeof(int));

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
  const int num_faces = fdata->num_faces;
  const int num_nodes = fdata->num_nodes;

  // Create the data structure for the node to face connectivity
  fdata->node_face_ptr = new int[ num_nodes+1 ];
  memset(fdata->node_face_ptr, 0, (num_nodes+1)*sizeof(int));

  // Count the number of times each node is referred to
  for ( int i = 0; i < 4*num_faces; i++ ){
    fdata->node_face_ptr[fdata->face_conn[i]+1]++;
  }

  // Adjust the counter array so that it points into the full array
  for ( int i = 1; i < num_nodes+1; i++ ){
    fdata->node_face_ptr[i] += fdata->node_face_ptr[i-1];
  }

  // Allocate the full node to face pointer array
  fdata->node_face_conn = new int[ fdata->node_face_ptr[num_nodes] ];
  for ( int i = 0; i < num_faces; i++ ){
    for ( int j = 0; j < 4; j++ ){
      int node = fdata->face_conn[4*i + j];
      fdata->node_face_conn[fdata->node_face_ptr[node]] = i;
      fdata->node_face_ptr[node]++;
    }
  }

  // Reset the pointer array so that it contains the correct offsets
  for ( int i = num_nodes; i >= 1; i-- ){
    fdata->node_face_ptr[i] = fdata->node_face_ptr[i-1];
  }
  fdata->node_face_ptr[0] = 0;

  // Loop over all the faces and reset node->face connectivity
  // to store both the adjacent block and the corresponding
  // node index into that array
  for ( int node = 0; node < num_nodes; node++ ){
    for ( int ip = fdata->node_face_ptr[node];
          ip < fdata->node_face_ptr[node+1]; ip++ ){
      int adj = fdata->node_face_conn[ip];
      int adj_index = 0;
      for ( ; adj_index < 4; adj_index++ ){
        if (fdata->face_conn[4*adj + adj_index] == node){
          break;
        }
      }
      fdata->node_face_conn[ip] = 4*adj + adj_index;
    }
  }
}

/*
  Based on the face to node connectivity information alone, compute a
  unique set of edges with associated edge numbers
*/
void TMRQuadForest::computeEdgesToFaces(){
  const int num_faces = fdata->num_faces;
  const int num_edges = fdata->num_edges;

  // Create the data structure for the edge to face connectivity
  fdata->edge_face_ptr = new int[ num_edges+1 ];
  memset(fdata->edge_face_ptr, 0, (num_edges+1)*sizeof(int));

  // Count the number of times each edge is referred to
  for ( int i = 0; i < 4*num_faces; i++ ){
    fdata->edge_face_ptr[fdata->face_edge_conn[i]+1]++;
  }

  // Adjust the counter array so that it points into the full array
  for ( int i = 1; i < num_edges+1; i++ ){
    fdata->edge_face_ptr[i] += fdata->edge_face_ptr[i-1];
  }

  // Allocate the full node to face pointer array
  fdata->edge_face_conn = new int[ fdata->edge_face_ptr[num_edges] ];
  for ( int face = 0; face < num_faces; face++ ){
    for ( int j = 0; j < 4; j++ ){
      int e = fdata->face_edge_conn[4*face + j];
      fdata->edge_face_conn[fdata->edge_face_ptr[e]] = face;
      fdata->edge_face_ptr[e]++;
    }
  }

  // Reset the pointer array so that it contains the correct offsets
  for ( int i = num_edges; i >= 1; i-- ){
    fdata->edge_face_ptr[i] = fdata->edge_face_ptr[i-1];
  }
  fdata->edge_face_ptr[0] = 0;

  // Loop over all edges and determine their relative orientation
  for ( int edge = 0; edge < num_edges; edge++ ){
    int face_owner = num_faces;
    int owner_index = 0;

    // Scan through the faces pointing to this edge to determine
    // the face owner - the face with the lowest index
    for ( int ip = fdata->edge_face_ptr[edge];
          ip < fdata->edge_face_ptr[edge+1]; ip++ ){
      int face = fdata->edge_face_conn[ip];
      if (face < face_owner){
        face_owner = face;

        // Find the new owner index
        owner_index = 0;
        for ( int j = 0; j < 4; j++, owner_index++ ){
          if (fdata->face_edge_conn[4*face + j] == edge){
            break;
          }
        }
      }
    }

    // Retrieve the first and second node numbers
    int n1 = fdata->face_conn[4*face_owner +
                              face_to_edge_nodes[owner_index][0]];
    int n2 = fdata->face_conn[4*face_owner +
                              face_to_edge_nodes[owner_index][1]];

    // Now determine the local edge index on each face and adjust
    // connectivity data
    for ( int ip = fdata->edge_face_ptr[edge];
          ip < fdata->edge_face_ptr[edge+1]; ip++ ){
      // Find the face index
      int face = fdata->edge_face_conn[ip];

      for ( int edge_index = 0; edge_index < 4; edge_index++ ){
        int nn1 = fdata->face_conn[4*face +
                                   face_to_edge_nodes[edge_index][0]];
        int nn2 = fdata->face_conn[4*face +
                                   face_to_edge_nodes[edge_index][1]];

        // Check if the edges now match up
        if ((n1 == nn1 && n2 == nn2) ||
            (n1 == nn2 && n2 == nn1)){
          fdata->edge_face_conn[ip] = 4*face + edge_index;
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
  const int num_faces = fdata->num_faces;

  // Now establish a unique ordering of the edges along each face
  fdata->face_edge_conn = new int[ 4*num_faces ];
  for ( int i = 0; i < 4*num_faces; i++ ){
    fdata->face_edge_conn[i] = -1;
  }

  // Keep track of the edge numbers
  int edge = 0;
  for ( int i = 0; i < num_faces; i++ ){
    // Loop over each edge on this face
    for ( int j = 0; j < 4; j++ ){
      if (fdata->face_edge_conn[4*i + j] < 0){
        int n1 = fdata->face_conn[4*i + face_to_edge_nodes[j][0]];
        int n2 = fdata->face_conn[4*i + face_to_edge_nodes[j][1]];

        // Keep track of the number of edges found
        const int max_nedges = 20;
        int edge_index[max_nedges];
        int nedges = 1;
        edge_index[0] = 4*i + j;

        // Set the edge number - if any is found
        int edge_num = -1;

        // Scan through the faces that share the same node and check
        // if any of the edges are also shared
        for ( int ip = fdata->node_face_ptr[n1];
              ip < fdata->node_face_ptr[n1+1]; ip++ ){
          int ii = fdata->node_face_conn[ip]/4;

          // Loop over each edge in the new face
          for ( int jj = 0; jj < 4; jj++ ){
            int nn1 = fdata->face_conn[4*ii + face_to_edge_nodes[jj][0]];
            int nn2 = fdata->face_conn[4*ii + face_to_edge_nodes[jj][1]];

            // Check if the face matches
            if ((n1 == nn1 && n2 == nn2) ||
                (n1 == nn2 && n2 == nn1)){
              if (fdata->face_edge_conn[4*ii + jj] >= 0){
                // If this edge has been ordered, copy over the edge
                // number
                edge_num = fdata->face_edge_conn[4*ii + jj];
              }
              else if (nedges < max_nedges){
                // This edge has not yet been ordered, add it to the
                // unordered list if there is still room if not, we
                // will detect and order it during a future iteration
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
          fdata->face_edge_conn[edge_index[ii]] = edge_num;
        }
      }
    }
  }

  // Set the total number of edges
  fdata->num_edges = edge;
}

/*
  Compute the face index that owns the edges/nodes
*/
void TMRQuadForest::computeFaceOwners(){
  const int num_faces = fdata->num_faces;
  const int num_edges = fdata->num_edges;
  const int num_nodes = fdata->num_nodes;

  // Allocate the edge/node ownership data
  fdata->edge_face_owners = new int[ num_edges ];
  fdata->node_face_owners = new int[ num_nodes ];

  // The owner is chosen as the connecting face with the lowest face
  // number
  for ( int edge = 0; edge < num_edges; edge++ ){
    fdata->edge_face_owners[edge] = num_faces;

    int ipend = fdata->edge_face_ptr[edge+1];
    for ( int ip = fdata->edge_face_ptr[edge]; ip < ipend; ip++ ){
      int face = fdata->edge_face_conn[ip]/4;
      if (face < fdata->edge_face_owners[edge]){
        fdata->edge_face_owners[edge] = face;
      }
    }
  }

  // Find the node owners
  for ( int node = 0; node < num_nodes; node++ ){
    fdata->node_face_owners[node] = num_faces;

    int ipend = fdata->node_face_ptr[node+1];
    for ( int ip = fdata->node_face_ptr[node]; ip < ipend; ip++ ){
      int face = fdata->node_face_conn[ip]/4;
      if (face < fdata->node_face_owners[node]){
        fdata->node_face_owners[node] = face;
      }
    }
  }
}

/*
  Write a representation of the connectivity of the forest out to a
  VTK file.
*/
void TMRQuadForest::writeToVTK( const char *filename ){
  const int num_faces = fdata->num_faces;
  const int num_nodes = fdata->num_nodes;
  const int *face_conn = fdata->face_conn;
  const int *node_face_owners = fdata->node_face_owners;

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
  const int num_faces = fdata->num_faces;
  const int num_nodes = fdata->num_nodes;
  const int *face_conn = fdata->face_conn;
  const int *node_face_owners = fdata->node_face_owners;

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
  Set the mesh order
*/
void TMRQuadForest::setMeshOrder( int _mesh_order,
                                  TMRInterpolationType _interp_type ){

  // Don't free the octants/owner information if it exists,
  // but free the connectivity and node data
  freeMeshData(0, 0);

  // Free the interpolation knots
  if (interp_knots){
    delete [] interp_knots;
  }

  // Check that the order falls within allowable bounds
  mesh_order = _mesh_order;
  if (mesh_order < 2){
    mesh_order = 2;
  }
  else if (mesh_order > MAX_ORDER){
    mesh_order = MAX_ORDER;
  }

  // Allocate the interpolation knots and set the knot locations
  interp_type = _interp_type;
  interp_knots = new double[ 2*(mesh_order) ];
  memset(interp_knots, 0.0, 2*(mesh_order)*sizeof(double));
  if (interp_type == TMR_GAUSS_LOBATTO_POINTS ||
      interp_type == TMR_BERNSTEIN_POINTS){
    interp_knots[0] = -1.0;
    interp_knots[mesh_order-1] = 1.0;
    for ( int i = 1; i < mesh_order-1; i++ ){
      interp_knots[i] = -cos(M_PI*i/(mesh_order-1));
    }
  }
  else {
    // Uniform mesh spacing
    interp_knots[0] = -1.0;
    interp_knots[mesh_order-1] = 1.0;
    for ( int i = 1; i < mesh_order-1; i++ ){
      interp_knots[i] = -1.0 + 2.0*i/(mesh_order-1);
    }
  }
}

/*
  Retrieve the mesh order
*/
int TMRQuadForest::getMeshOrder(){
  return mesh_order;
}

/*
  Retrieve the interpolation type for this mesh
*/
TMRInterpolationType TMRQuadForest::getInterpType(){
  return interp_type;
}

/*
  Get the node-processor ownership range
*/
int TMRQuadForest::getOwnedNodeRange( const int **_node_range ){
  if (_node_range){
    *_node_range = node_range;
  }
  return mpi_size;
}

/*
  Get the quadrants and the nodes
*/
void TMRQuadForest::getQuadrants( TMRQuadrantArray **_quadrants ){
  if (_quadrants){
    *_quadrants = quadrants;
  }
}

/*
  Get the node numbers (note that this may be NULL)
*/
int TMRQuadForest::getNodeNumbers( const int **_node_numbers ){
  if (_node_numbers){
    *_node_numbers = node_numbers;
  }
  return num_local_nodes;
}

/*
  Get the node locations
*/
int TMRQuadForest::getPoints( TMRPoint **_X ){
  if (_X){
    *_X = X;
  }
  return num_local_nodes;
}

/*
  Retrieve the local node number
*/
int TMRQuadForest::getLocalNodeNumber( int node ){
  if (node_numbers){
    int *item = (int*)bsearch(&node, node_numbers, num_local_nodes,
                              sizeof(int), compare_integers);
    if (item){
      return item - node_numbers;
    }
  }
  return -1;
}

/*
  Get the knot points for the interpolation
*/
int TMRQuadForest::getInterpKnots( const double **_knots ){
  if (_knots){
    *_knots = interp_knots;
  }
  return mesh_order;
}

/*
  Evaluate the interpolant at the given parametric point
*/
void TMRQuadForest::evalInterp( const double pt[], double N[] ){
  double Nu[MAX_ORDER], Nv[MAX_ORDER];
  if (interp_type == TMR_BERNSTEIN_POINTS){
    // Evaluate the bernstein shape functions
    bernstein_shape_functions(mesh_order, pt[0], Nu);
    bernstein_shape_functions(mesh_order, pt[1], Nv);
  }
  else{
    // Evaluate the lagrange shape functions
    lagrange_shape_functions(mesh_order, pt[0], interp_knots, Nu);
    lagrange_shape_functions(mesh_order, pt[1], interp_knots, Nv);
  }

  for ( int j = 0; j < mesh_order; j++ ){
    for ( int i = 0; i < mesh_order; i++ ){
      N[0] = Nu[i]*Nv[j];
      N++;
    }
  }
}

/*
  Evaluate the interpolant at the given parametric point
*/
void TMRQuadForest::evalInterp( const double pt[], double N[],
                                double N1[], double N2[] ){
  double Nu[MAX_ORDER], Nv[MAX_ORDER];
  double Nud[MAX_ORDER], Nvd[MAX_ORDER];

  if (interp_type == TMR_BERNSTEIN_POINTS){
    // Evaluate the bernstein shape functions
    bernstein_shape_func_derivative(mesh_order, pt[0], Nu, Nud);
    bernstein_shape_func_derivative(mesh_order, pt[1], Nv, Nvd);
  }
  else {
    // Evaluate the shape functions
    lagrange_shape_func_derivative(mesh_order, pt[0], interp_knots, Nu, Nud);
    lagrange_shape_func_derivative(mesh_order, pt[1], interp_knots, Nv, Nvd);
  }

  for ( int j = 0; j < mesh_order; j++ ){
    for ( int i = 0; i < mesh_order; i++ ){
      N[0] = Nu[i]*Nv[j];
      N1[0] = Nud[i]*Nv[j];
      N2[0] = Nu[i]*Nvd[j];
      N++;
      N1++;
      N2++;
    }
  }
}

/*
  Evaluate the interpolant at the given parametric point
*/
void TMRQuadForest::evalInterp( const double pt[], double N[],
                                double N1[], double N2[],
                                double N11[], double N22[], double N12[] ){
  double Nu[MAX_ORDER], Nud[MAX_ORDER], Nudd[MAX_ORDER];
  double Nv[MAX_ORDER], Nvd[MAX_ORDER], Nvdd[MAX_ORDER];

  if (interp_type == TMR_BERNSTEIN_POINTS){
    // Evaluate the shape functions
    lagrange_shape_func_second_derivative(mesh_order, pt[0], interp_knots,
                                          Nu, Nud, Nudd);
    lagrange_shape_func_second_derivative(mesh_order, pt[1], interp_knots,
                                          Nv, Nvd, Nvdd);
  }
  else {
    // Evaluate the bernstein shape functions
    bernstein_shape_func_second_derivative(mesh_order, pt[0], Nu, Nud, Nudd);
    bernstein_shape_func_second_derivative(mesh_order, pt[0], Nv, Nvd, Nvdd);
  }

  for ( int j = 0; j < mesh_order; j++ ){
    for ( int i = 0; i < mesh_order; i++ ){
      N[0] = Nu[i]*Nv[j];
      N1[0] = Nud[i]*Nv[j];
      N2[0] = Nu[i]*Nvd[j];
      N11[0] = Nudd[i]*Nv[j];
      N22[0] = Nu[i]*Nvdd[j];
      N12[0] = Nud[i]*Nvd[j];
      N++;
      N1++;
      N2++;
      N11++;
      N22++;
      N12++;
    }
  }
}

/*
  Retrieve information about the connectivity between faces, edges and
  nodes
*/
void TMRQuadForest::getConnectivity( int *_nfaces,
                                     int *_nedges, int *_nnodes,
                                     const int **_face_conn,
                                     const int **_face_edge_conn ){
  if (fdata){
    if (_nfaces){ *_nfaces = fdata->num_faces; }
    if (_nedges){ *_nedges = fdata->num_edges; }
    if (_nnodes){ *_nnodes = fdata->num_nodes; }
    if (_face_conn){ *_face_conn = fdata->face_conn; }
    if (_face_edge_conn){ *_face_edge_conn = fdata->face_edge_conn; }
  }
  else {
    if (_nfaces){ *_nfaces = 0; }
    if (_nedges){ *_nedges = 0; }
    if (_nnodes){ *_nnodes = 0; }
    if (_face_conn){ *_face_conn = NULL; }
    if (_face_edge_conn){ *_face_edge_conn = NULL; }
  }
}

/*
  Retrieve the inverse of the connectivity
*/
void TMRQuadForest::getInverseConnectivity( const int **_node_face_conn,
                                            const int **_node_face_ptr,
                                            const int **_edge_face_conn,
                                            const int **_edge_face_ptr ){
  if (fdata){
    if (_node_face_conn){ *_node_face_conn = fdata->node_face_conn; }
    if (_node_face_ptr){ *_node_face_ptr = fdata->node_face_ptr; }
    if (_edge_face_conn){ *_edge_face_conn = fdata->edge_face_conn; }
    if (_edge_face_ptr){ *_edge_face_ptr = fdata->edge_face_ptr; }
  }
  else {
    if (_node_face_conn){ *_node_face_conn = NULL; }
    if (_node_face_ptr){ *_node_face_ptr = NULL; }
    if (_edge_face_conn){ *_edge_face_conn = NULL; }
    if (_edge_face_ptr){ *_edge_face_ptr = NULL; }
  }
}

/*
  Create a forest with the specified refinement level
*/
void TMRQuadForest::createTrees( int refine_level ){
  const int num_faces = fdata->num_faces;

  // Free all the mesh-specific data
  freeMeshData();

  // Set the max/min refinement level
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
        array[count].info = 0;
        array[count].x = x;
        array[count].y = y;
        count++;
      }
    }
  }

  // Create the array of quadrants
  quadrants = new TMRQuadrantArray(array, size);
  quadrants->sort();

  // Set the local ordering for the elements
  int quad_size;
  TMRQuadrant *quads;
  quadrants->getArray(&quads, &quad_size);
  for ( int i = 0; i < quad_size; i++ ){
    quads[i].tag = i;
  }

  // Set the last quadrant
  TMRQuadrant p;
  p.tag = -1;
  p.face = num_faces-1;
  p.x = p.y = hmax;
  p.info = 0;
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
  const int num_faces = fdata->num_faces;

  // Free all the mesh-specific data
  freeMeshData();

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
      array[count].info = 0;
      array[count].x = x;
      array[count].y = y;
    }
  }

  // Create the array of quadrants
  quadrants = new TMRQuadrantArray(array, size);
  quadrants->sort();

  // Set the local ordering for the elements
  int quad_size;
  TMRQuadrant *quads;
  quadrants->getArray(&quads, &quad_size);
  for ( int i = 0; i < quad_size; i++ ){
    quads[i].tag = i;
  }

  // Set the last quadrant
  TMRQuadrant p;
  p.tag = -1;
  p.face = num_faces-1;
  p.x = p.y = 1 << TMR_MAX_LEVEL;
  p.info = 0;
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
  const int num_faces = fdata->num_faces;

  // Free everything but the quadrants
  freeMeshData(0);

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

  // Set the last quadrant
  TMRQuadrant q;
  q.face = num_faces-1;
  q.tag = -1;
  q.level = 0;
  q.info = 0;
  q.x = q.y = 1 << TMR_MAX_LEVEL;
  if (new_size > 0){
    q = new_array[0];
  }

  // Free the quadrant arrays
  delete quadrants;
  quadrants = new TMRQuadrantArray(new_array, new_size);

  owners = new TMRQuadrant[ mpi_size ];
  MPI_Allgather(&q, 1, TMRQuadrant_MPI_type,
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
  TMRQuadForest *dup = new TMRQuadForest(comm, mesh_order, interp_type);
  if (fdata){
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
  TMRQuadForest *coarse = new TMRQuadForest(comm, mesh_order, interp_type);
  if (fdata){
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

    // Order the labels of the coarse octants
    coarse->quadrants->getArray(&array, &size);
    for ( int i = 0; i < size; i++ ){
      array[i].tag = i;
    }

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
  // Free the data associated with the mesh but not the quadrants/owners
  freeMeshData(0, 0);

  // Adjust the min and max levels to ensure consistency
  if (min_level < 0){ min_level = 0; }
  if (max_level > TMR_MAX_LEVEL){ max_level = TMR_MAX_LEVEL; }

  // This is just a sanity check
  if (min_level > max_level){ min_level = max_level; }

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
          q.info = 0;

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
          q.info = 0;

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
        q.info = 0;
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

  // Get the quadrants and order their labels
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
TMRQuadrantArray
  *TMRQuadForest::distributeQuadrants( TMRQuadrantArray *list,
                                       int use_tags,
                                       int **_quad_ptr,
                                       int **_quad_recv_ptr,
                                       int include_local,
                                       int use_node_index ){
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
  TMRQuadrantArray *dist = sendQuadrants(list, quad_ptr,
                                         quad_recv_ptr, use_node_index);

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
                                                const int *quad_recv_ptr,
                                                int use_node_index  ){
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

  return new TMRQuadrantArray(recv_array, recv_size, use_node_index);
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
  int edge = fdata->face_edge_conn[4*face + edge_index];

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
  int n1 = fdata->face_conn[4*face + face_to_edge_nodes[edge_index][0]];
  int n2 = fdata->face_conn[4*face + face_to_edge_nodes[edge_index][1]];

  // Now, cycle through all the adjacent faces
  for ( int ip = fdata->edge_face_ptr[edge];
        ip < fdata->edge_face_ptr[edge+1]; ip++ ){
    // Get the face that is adjacent across this edge
    int adj = fdata->edge_face_conn[ip]/4;
    if (adj != face){
      // Get the adjacent edge index
      int adj_index = fdata->edge_face_conn[ip] % 4;

      // Get the nodes on the adjacent face
      int nn1 = fdata->face_conn[4*adj + face_to_edge_nodes[adj_index][0]];
      int nn2 = fdata->face_conn[4*adj + face_to_edge_nodes[adj_index][1]];

      // Add the quadrant to the list
      int reverse = (n1 == nn2 && n2 == nn1);
      int32_t u = ucoord;
      if (reverse){
        u = hmax - 2*h - ucoord;
      }

      TMRQuadrant neighbor;
      neighbor.face = adj;
      neighbor.level = p.level;
      neighbor.info = 0;
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
  int node = fdata->face_conn[4*face + corner];

  // Compute the edge length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - p.level);

  // Now, cycle through all the adjacent faces
  for ( int ip = fdata->node_face_ptr[node];
        ip < fdata->node_face_ptr[node+1]; ip++ ){

    // Get the faces that are adjacent across this edge
    int adj = fdata->node_face_conn[ip]/4;
    if (adj != face){
      int adj_index = fdata->node_face_conn[ip] % 4;

      // Compute the quadrant location
      TMRQuadrant neighbor;
      neighbor.face = adj;
      neighbor.level = p.level;
      neighbor.info = 0;
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
  if (!quadrants){
    fprintf(stderr, "TMRQuadForest Error: Cannot call balance(), "
            "no quadrants have been created\n");
    return;
  }

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
  int edge = fdata->face_edge_conn[4*face + edge_index];

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
  int n1 = fdata->face_conn[4*face + face_to_edge_nodes[edge_index][0]];
  int n2 = fdata->face_conn[4*face + face_to_edge_nodes[edge_index][1]];

  // Now, cycle through all the adjacent faces
  for ( int ip = fdata->edge_face_ptr[edge];
        ip < fdata->edge_face_ptr[edge+1]; ip++ ){
    // Get the face that is adjacent across this edge
    int adj = fdata->edge_face_conn[ip]/4;
    if (adj != face){
      // Get the adjacent edge index
      int adj_index = fdata->edge_face_conn[ip] % 4;

      // Get the nodes on the adjacent face
      int nn1 = fdata->face_conn[4*adj + face_to_edge_nodes[adj_index][0]];
      int nn2 = fdata->face_conn[4*adj + face_to_edge_nodes[adj_index][1]];

      // Add the quadrant to the list
      int reverse = (n1 == nn2 && n2 == nn1);
      int32_t u = ucoord;
      if (reverse){
        u = hmax - h - ucoord;
      }

      TMRQuadrant neighbor;
      neighbor.face = adj;
      neighbor.level = p.level;
      neighbor.info = 0;
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
  int node = fdata->face_conn[4*face + corner];

  // Compute the edge length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - p.level);

  // Now, cycle through all the adjacent faces
  for ( int ip = fdata->node_face_ptr[node];
        ip < fdata->node_face_ptr[node+1]; ip++ ){

    // Get the faces that are adjacent across this edge
    int adj = fdata->node_face_conn[ip]/4;
    if (adj != face){
      int adj_index = fdata->node_face_conn[ip] % 4;

      // Compute the quadrant location
      TMRQuadrant neighbor;
      neighbor.face = adj;
      neighbor.level = p.level;
      neighbor.info = 0;
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
  if (adjacent){
    delete adjacent;
  }

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
  delete list;
  adjacent->sort();

  // Set the local quadrant tags to be their local index
  adjacent->getArray(&array, &size);
  for ( int i = 0; i < size; i++ ){
    array[i].tag = i;
  }
}

/*
  Compute the dependent nodes (hanging edge) on each face and on
  the interfaces between adjacent quadtrees.

  The dependent edge nodes may occur on any face within the quadtree.
  Each quadrant stores a flag in its info member that contains the
  bits indicating which edge is a dependent edge. Note that the
  mid-edge node may also be dependent, depending on the order of the mesh.

  side effects:
  info flags set in all quads in the quadrants list and adjacent list
*/
void TMRQuadForest::computeDepEdges(){
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
            computeAdjacentDepEdges(edge_index, &q, adjquads);
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
  Determine if there is an adjacent quadrant on the connecting edge.
*/
void TMRQuadForest::computeAdjacentDepEdges( int edge_index,
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
  int edge = fdata->face_edge_conn[4*face_owner + edge_index];
  int n1 = fdata->face_conn[4*face_owner +
                            face_to_edge_nodes[edge_index][0]];
  int n2 = fdata->face_conn[4*face_owner +
                            face_to_edge_nodes[edge_index][1]];

  // Now, cycle through all the adjacent edges
  for ( int ip = fdata->edge_face_ptr[edge];
        ip < fdata->edge_face_ptr[edge+1]; ip++ ){
    int face = fdata->edge_face_conn[ip]/4;

    if (face_owner != face){
      // Get the adjacent edge index
      int adj_index = fdata->edge_face_conn[ip] % 4;

      // Get the nodes on the adjacent face
      int nn1 = fdata->face_conn[4*face + face_to_edge_nodes[adj_index][0]];
      int nn2 = fdata->face_conn[4*face + face_to_edge_nodes[adj_index][1]];

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
      quad.info = 0;
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
}

/*
  Label the dependent face and edge nodes

  This code is called after all the dependent faces have been
  computed.  Note that this relies on the mesh being edge-balanced
  (which is required). It also relies on the connectivity still
  being in a local state (such that conn[] refers to the local nodes)
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

          // Compute the offset into the array based on the edge index
          int offset = (edge_index % 2)*(mesh_order-1)*mesh_order;
          if (edge_index < 2){
            offset = (edge_index % 2)*(mesh_order-1);
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
            end = mesh_order;
            if (mesh_order == 3 && interp_type != TMR_BERNSTEIN_POINTS){
              end = mesh_order-1;
            }
          }
          else {
            start = 0;
            if (mesh_order == 3 && interp_type != TMR_BERNSTEIN_POINTS){
              start = 1;
            }
            end = mesh_order-1;
          }

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
      int node = fdata->face_conn[4*face + corner];

      if (face != fdata->node_face_owners[node]){
        // Get the pointer information
        int ptr = fdata->node_face_ptr[node];
        int adj = fdata->node_face_conn[ptr]/4;
        int adj_index = fdata->node_face_conn[ptr] % 4;

        // Copy the coordinates of the adjacent face
        quad->face = adj;
        quad->x = hmax*(adj_index % 2);
        quad->y = hmax*(adj_index / 2);
      }
    }
    else {
      // Which edge index are we dealing with?
      int edge_index = fx*(fx0 ? 0 : 1) + fy*(fy0 ? 2 : 3);
      int edge = fdata->face_edge_conn[4*face + edge_index];

      if (face != fdata->edge_face_owners[edge]){
        // Get the adjacent edge index on the opposite face
        int ptr = fdata->edge_face_ptr[edge];
        int adj = fdata->edge_face_conn[ptr]/4;
        int adj_index = fdata->edge_face_conn[ptr] % 4;

        // Retrieve the first and second node numbers to determine the
        // relative orientation between this edge and each adjacent edge
        int n1 = fdata->face_conn[4*face + face_to_edge_nodes[edge_index][0]];
        int n2 = fdata->face_conn[4*face + face_to_edge_nodes[edge_index][1]];

        // Get the orientation
        int nn1 = fdata->face_conn[4*adj + face_to_edge_nodes[adj_index][0]];
        int nn2 = fdata->face_conn[4*adj + face_to_edge_nodes[adj_index][1]];

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
*/
void TMRQuadForest::createNodes(){
  if (!quadrants){
    fprintf(stderr, "TMRQuadForest Error: Cannot call createNodes(), "
            "no quadrants have been created\n");
    return;
  }
  if (conn){
    // The connectivity has already been created and not deleted so
    // there is no need to create it a second time.
    return;
  }

  // Send/recv the adjacent quadrants
  computeAdjacentQuadrants();

  // Compute the dependent face nodes
  computeDepEdges();

  // Create and assign the ownership for the local node numbers
  TMRQuadrantArray *nodes = createLocalNodes();

  // Retrieve the size of the node array and count up the offsets for
  // each node. When mesh_order <= 3, the offset array will be equal
  // to the index since each quadrant in the node array will represent
  // only one node. When mesh_order >= 4 the quadrants may represent
  // more than one node.
  int node_size;
  TMRQuadrant *node_array;
  nodes->getArray(&node_array, &node_size);

  int *node_offset = new int[ node_size ];
  num_local_nodes = 0;
  for ( int i = 0; i < node_size; i++ ){
    node_offset[i] = num_local_nodes;
    num_local_nodes += node_array[i].level;
  }

  // Create the connectivity based on the node array
  createLocalConn(nodes, node_offset);

  // Allocate an array that will store the new node numbers
  node_numbers = new int[ num_local_nodes ];
  memset(node_numbers, 0, num_local_nodes*sizeof(int));

  // Label any node that is dependent as a negative value
  labelDependentNodes(node_numbers);

  // Count up and set the dependent node numbers
  num_dep_nodes = 0;
  for ( int i = 0; i < num_local_nodes; i++ ){
    if (node_numbers[i] == -1){
      num_dep_nodes++;
      node_numbers[i] = -num_dep_nodes;
    }
  }

  // Create the local connectivyt based on the node array
  createDependentConn(node_numbers, nodes, node_offset);

  // Loop over all the nodes, check whether they are local (all
  // dependent nodes are dependent)
  const int use_node_index = 1;
  TMRQuadrantHash *ext_nodes = new TMRQuadrantHash(use_node_index);

  // Add nodes that are externally owned
  nodes->getArray(&node_array, &node_size);

  // Get the MPI owner for each node
  int index = 0;
  for ( int i = 0; i < node_size; i++ ){
    // Check if this is a dependent node or not...
    if (node_numbers[index] >= 0){
      // Send it to the owner processor
      int mpi_owner = node_array[i].tag;
      if (mpi_owner != mpi_rank){
        TMRQuadrant node = node_array[i];
        node.tag = mpi_owner;
        ext_nodes->addQuadrant(&node);

        // Label the nodes here as owned by another processor
        for ( int k = 0; k < node_array[i].level; k++ ){
          node_numbers[index + k] = -num_dep_nodes-1;
        }
      }
    }
    index += node_array[i].level;
  }

  // Now all the external and dependent nodes are labeled, any remaining
  // nodes that have a non-negative value are independent and must
  // be ordered. These are the locally owned nodes.
  num_owned_nodes = 0;
  for ( int i = 0; i < num_local_nodes; i++ ){
    if (node_numbers[i] >= 0){
      num_owned_nodes++;
    }
  }

  // Gather the owned node counts from each processor
  node_range = new int[ mpi_size+1 ];
  memset(node_range, 0, (mpi_size+1)*sizeof(int));
  MPI_Allgather(&num_owned_nodes, 1, MPI_INT,
                &node_range[1], 1, MPI_INT, comm);

  // Set the offsets to each node
  for ( int i = 0; i < mpi_size; i++ ){
    node_range[i+1] += node_range[i];
  }

  // Set the global node numbers for the owned nodes
  num_owned_nodes = 0;
  for ( int i = 0; i < num_local_nodes; i++ ){
    if (node_numbers[i] >= 0){
      node_numbers[i] = node_range[mpi_rank] + num_owned_nodes;
      num_owned_nodes++;
    }
  }

  // Create an array of all the independent nodes that are
  // owned by other processors and referenced by the
  // elements on this processor.
  TMRQuadrantArray *ext_array = ext_nodes->toArray();
  delete ext_nodes;

    // Sort based on the tags
  int ext_size;
  TMRQuadrant *ext_quads;
  ext_array->getArray(&ext_quads, &ext_size);
  qsort(ext_quads, ext_size, sizeof(TMRQuadrant), compare_quadrant_tags);

  // Distribute the non-local nodes back to their owning processors
  // to determine their node numbers
  int use_tags = 1;
  int *send_ptr, *recv_ptr;
  TMRQuadrantArray *dist_nodes = distributeQuadrants(ext_array, use_tags,
                                                     &send_ptr, &recv_ptr);
  delete ext_array;

  // Loop over the off-processor nodes and search for them in the
  // sorted node list and assign them the correct tag
  int dist_size;
  TMRQuadrant *dist_quads;
  dist_nodes->getArray(&dist_quads, &dist_size);
  for ( int i = 0; i < dist_size; i++ ){
    TMRQuadrant *t = nodes->contains(&dist_quads[i]);
    if (t){
      // Compute the node number
      int index = t - node_array;
      dist_quads[i].tag = node_numbers[node_offset[index]];
    }
  }

  // Send the nodes back to the original processors
  TMRQuadrantArray *return_nodes = sendQuadrants(dist_nodes,
                                                 recv_ptr, send_ptr);
  delete dist_nodes;
  delete [] recv_ptr;
  delete [] send_ptr;

  // Now go back through and set the external node numbers
  int return_size;
  TMRQuadrant *return_quads;
  return_nodes->getArray(&return_quads, &return_size);
  for ( int i = 0; i < return_size; i++ ){
    TMRQuadrant *t = nodes->contains(&return_quads[i]);
    for ( int k = 0; k < t->level; k++ ){
      int index = t - node_array;
      node_numbers[node_offset[index] + k] = return_quads[i].tag + k;
    }
  }
  delete return_nodes;

  // Free the local node array
  delete nodes;
  delete [] node_offset;

  // Apply the node numbering scheme to the local connectivity to
  // give us a global numbering scheme
  int num_elements;
  quadrants->getArray(NULL, &num_elements);
  for ( int i = 0; i < mesh_order*mesh_order*num_elements; i++ ){
    conn[i] = node_numbers[conn[i]];
  }

  // Apply the node numbering scheme to the local dependent node
  // connectivity
  for ( int i = 0; i < dep_ptr[num_dep_nodes]; i++ ){
    dep_conn[i] = node_numbers[dep_conn[i]];
  }

  // Now, sort the global numbers
  qsort(node_numbers, num_local_nodes, sizeof(int), compare_integers);

  // Compute num_ext_pre_nodes -- the number of external pre nodes
  int *item = (int*)bsearch(&node_range[mpi_rank], node_numbers,
                            num_local_nodes, sizeof(int), compare_integers);
  ext_pre_offset = item - node_numbers;

  // Evaluate the node locations
  evaluateNodeLocations();
}

/*
  Create the local nodes and assign their owners

  The ownership rules are as follows:

  1. Dependent nodes are locally owned (this ownership property is
  assigned in a second step since it can be decided locally)

  2. Independent nodes are owned on processors where the node is
  created from an element (not just because it is required
  for a dependent node)

  3. If multiple elements create the same node, its owner is the proc
  with the lower processor rank

  The returns an array of quadrants with the following information:
  1. tag represents the MPI owner
  2. info represents the label (node/edge/face)
  3. level represents the number of nodes represented by the quad
*/
TMRQuadrantArray *TMRQuadForest::createLocalNodes(){
  // Allocate the array of elements
  int num_elements;
  TMRQuadrant *quads;
  quadrants->getArray(&quads, &num_elements);

  // Create all the nodes/edges/faces
  const int use_node_index = 1;
  TMRQuadrantHash *local_nodes = new TMRQuadrantHash(use_node_index);

  // Set the node, edge and face label
  int node_label, edge_label, face_label;
  node_label = edge_label = face_label = TMR_QUAD_NODE_LABEL;

  int label_type[3];
  // If the mesh order is high enough, we will have multiple nodes
  // per edge/face
  initLabel(mesh_order, interp_type, label_type);

  node_label = label_type[0];
  edge_label = label_type[1];
  face_label = label_type[2];

  // Set the node locations
  if (mesh_order == 2){
    // First of all, add all the nodes from the local elements
    // on this processor
    for ( int i = 0; i < num_elements; i++ ){
      const int32_t h = 1 << (TMR_MAX_LEVEL - quads[i].level);
      for ( int jj = 0; jj < 2; jj++ ){
        for ( int ii = 0; ii < 2; ii++ ){
          TMRQuadrant node;
          node.face = quads[i].face;
          node.level = 1;
          node.x = quads[i].x + h*ii;
          node.y = quads[i].y + h*jj;
          node.tag = mpi_rank;
          node.info = node_label;
          transformNode(&node);
          local_nodes->addQuadrant(&node);
        }
      }
    }

    // Add the nodes that the dependent nodes depend on
    for ( int i = 0; i < num_elements; i++ ){
      // Add the external nodes from dependent edges
      if (quads[i].info){
        for ( int edge_index = 0; edge_index < 4; edge_index++ ){
          if (quads[i].info & 1 << edge_index){
            TMRQuadrant parent;
            quads[i].parent(&parent);

            const int32_t hp = 1 << (TMR_MAX_LEVEL - parent.level);
            for ( int ii = 0; ii < 2; ii++ ){
              TMRQuadrant node;
              node.face = parent.face;
              node.level = 1;
              if (edge_index < 2){
                node.x = parent.x + hp*(edge_index % 2);
                node.y = parent.y + hp*ii;
              }
              else {
                node.x = parent.x + hp*ii;
                node.y = parent.y + hp*(edge_index % 2);
              }
              // Assign a negative rank index for now...
              node.tag = -1;
              node.info = node_label;
              transformNode(&node);
              local_nodes->addQuadrant(&node);
            }
          }
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
          node.x = quads[i].x + h*ii;
          node.y = quads[i].y + h*jj;
          if ((ii == 0 || ii == 2) &&
              (jj == 0 || jj == 2)){
            node.level = 1;
            node.info = node_label;
          }
          else if (ii == 0 || ii == 2 ||
                   jj == 0 || jj == 2){
            node.level = mesh_order-2;
            node.info = edge_label;
          }
          else {
            node.level = (mesh_order-2)*(mesh_order-2);
            node.info = face_label;
          }
          node.tag = mpi_rank;
          transformNode(&node);
          local_nodes->addQuadrant(&node);
        }
      }
    }

    // Add the nodes from the dependent edges (these will only be
    // overwritten if required)
    for ( int i = 0; i < num_elements; i++ ){
      if (quads[i].info){
        const int32_t h = 1 << (TMR_MAX_LEVEL - quads[i].level);

        for ( int edge_index = 0; edge_index < 4; edge_index++ ){
          if (quads[i].info & 1 << edge_index){
            TMRQuadrant parent;
            quads[i].parent(&parent);

            for ( int ii = 0; ii < 3; ii++ ){
              TMRQuadrant node;
              node.face = parent.face;
              if (edge_index < 2){
                node.x = parent.x + 2*h*(edge_index % 2);
                node.y = parent.y + h*ii;
              }
              else {
                node.x = parent.x + h*ii;
                node.y = parent.y + 2*h*(edge_index % 2);
              }
              if (ii == 0 || ii == 2){
                node.level = 1;
                node.info = node_label;
              }
              else {
                node.level = mesh_order-2;
                node.info = edge_label;
              }
              // Assign the negative rank to this processor
              node.tag = -1;
              transformNode(&node);
              local_nodes->addQuadrant(&node);
            }
          }
        }
      }
    }
  }

  // Now the local_nodes hash table contains all of the nodes
  // (dependent, indepdnent and non-local) that are referenced by
  // this processor
  TMRQuadrantArray *nodes = local_nodes->toArray();
  delete local_nodes;
  nodes->sort();

  // Now, determine the node ownership - if nodes that are not
  // dependent on this processor
  int use_tags = 0, include_local = 0;
  int *send_ptr, *recv_ptr;
  TMRQuadrantArray *recv_nodes =
    distributeQuadrants(nodes, use_tags, &send_ptr, &recv_ptr,
                        include_local, use_node_index);

  // Create a unique list of the nodes sent to this processor
  TMRQuadrantArray *recv_sorted = recv_nodes->duplicate();
  recv_sorted->sort();

  // Now loop over nodes sent from other processors and decide
  // which processor owns the node.
  int recv_size;
  TMRQuadrant *recv_array;
  recv_nodes->getArray(&recv_array, &recv_size);

  // Loop over all the nodes and see if they have a donor
  // element from another processor that is not from a dependent
  // node relationship
  for ( int i = 0; i < recv_size; i++ ){
    // This is the processor that donates
    if (recv_array[i].tag >= 0){
      TMRQuadrant *t = recv_sorted->contains(&recv_array[i]);
      if (t->tag < 0){
        // t is not the owner, it is defined from a dependent edge
        t->tag = recv_array[i].tag;
      }
      else { // t->tag >= 0
        // *t is not the owner since it has a higher rank than the
        // other, equivalent node -- re-assign the node number
        if (recv_array[i].tag < t->tag){
          t->tag = recv_array[i].tag;
        }
      }
    }
  }

  // Now search and make consistent the owners and the internal
  // nodes on this processor
  int sorted_size;
  TMRQuadrant *sorted_array;
  recv_sorted->getArray(&sorted_array, &sorted_size);
  for ( int i = 0; i < sorted_size; i++ ){
    TMRQuadrant *t = nodes->contains(&sorted_array[i]);
    // Note that even though these nodes are mapped to this processor,
    // they may not be defined on it for some corner cases...
    if (t){
      if (t->tag < 0){
        t->tag = sorted_array[i].tag;
      }
      else if (sorted_array[i].tag < 0){
        sorted_array[i].tag = t->tag;
      }
      else if (t->tag < sorted_array[i].tag){
        sorted_array[i].tag = t->tag;
      }
      else {
        t->tag = sorted_array[i].tag;
      }
    }
  }

  // Make the return nodes consistent with the sorted list that
  // is unique
  for ( int i = 0; i < recv_size; i++ ){
    // This is the processor that donates from an owner
    TMRQuadrant *t = recv_sorted->contains(&recv_array[i]);
    recv_array[i].tag = t->tag;
  }

  delete recv_sorted;

  // Adjust the send_ptr array since there will be a gap
  // for the quadrants that are processor-local
  int offset = send_ptr[mpi_rank+1] - send_ptr[mpi_rank];
  for ( int i = mpi_rank+1; i <= mpi_size; i++ ){
    send_ptr[i] -= offset;
  }

  // Return the nodes back to the senders with the new
  // owner information attached
  TMRQuadrantArray *owner_nodes =
    sendQuadrants(recv_nodes, recv_ptr, send_ptr, use_node_index);
  delete recv_nodes;
  delete [] recv_ptr;
  delete [] send_ptr;

  // Go trhough the owner nodes and assign the MPI owner
  int owner_size;
  TMRQuadrant *owner_array;
  owner_nodes->getArray(&owner_array, &owner_size);

  for ( int i = 0; i < owner_size; i++ ){
    // Get the owner of the node on this processor
    TMRQuadrant *t = nodes->contains(&owner_array[i]);

    // Assign the MPI owner rank
    t->tag = owner_array[i].tag;
  }
  delete owner_nodes;

  // Return the owners for each node
  return nodes;
}

/*
  Create the local connectivity based on the ordering in the node
  array and the offset arrays

  This local connectivity is based on the local ordering of the nodes
  on this processor. This local ordering is overwritten in a second
  step once the global order of the nodes is finalized.

  The nodes are represented by quadrants with the additional info:
  1. The tag member is the MPI owner of the node
  2. The info member is the node/edge/face label
  3. The level member contains the number of nodes per node object
  which depends on the order of the mesh.

  input:
  nodes:        the array of quadrants that represent nodes
  node_offset:  the array of offsets for each node
*/
void TMRQuadForest::createLocalConn( TMRQuadrantArray *nodes,
                                     const int *node_offset ){
  // Retrieve the quadrants on this processor
  int num_elements;
  TMRQuadrant *quads;
  quadrants->getArray(&quads, &num_elements);

  // Retrieve the nodes
  int node_size;
  TMRQuadrant *node_array;
  nodes->getArray(&node_array, &node_size);

  // Set the node, edge and face label
  int node_label, edge_label, face_label;
  node_label = edge_label = face_label = TMR_QUAD_NODE_LABEL;

  int label_type[3];
  // If the mesh order is high enough, we will have multiple nodes
  // per edge/face
  initLabel(mesh_order, interp_type, label_type);

  node_label = label_type[0];
  edge_label = label_type[1];
  face_label = label_type[2];

  // Allocate the connectivity
  int size = mesh_order*mesh_order*num_elements;
  conn = new int[ size ];
  memset(conn, 0, size*sizeof(int));

  if (mesh_order <= 3){
    for ( int i = 0; i < num_elements; i++ ){
      int *c = &conn[mesh_order*mesh_order*i];
      const int32_t h = 1 << (TMR_MAX_LEVEL - quads[i].level - 1);

      // Loop over the element nodes
      for ( int corner_index = 0; corner_index < 4; corner_index++ ){
        TMRQuadrant node;
        node.face = quads[i].face;
        node.x = quads[i].x + 2*h*(corner_index % 2);
        node.y = quads[i].y + 2*h*(corner_index / 2);
        node.info = node_label;
        transformNode(&node);
        TMRQuadrant *t = nodes->contains(&node);
        int index = t - node_array;
        int offset = (mesh_order-1)*(corner_index % 2) +
                     (mesh_order-1)*mesh_order*(corner_index/2);
        c[offset] = node_offset[index];
      }

      if (mesh_order == 3){
        // Loop over the edges and get the owners
        for ( int edge_index = 0; edge_index < 4; edge_index++ ){
          TMRQuadrant node;
          node.face = quads[i].face;
          node.info = edge_label;
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
          int index = t - node_array;
          if (edge_index < 2){
            int offset = mesh_order + (mesh_order-1)*edge_index;
            c[offset] = node_offset[index];
          }
          else {
            int offset = 1 + (mesh_order-1)*mesh_order*(edge_index % 2);
            c[offset] = node_offset[index];
          }
        }

        TMRQuadrant node;
        node.face = quads[i].face;
        node.x = quads[i].x + h;
        node.y = quads[i].y + h;
        node.info = face_label;
        transformNode(&node);
        TMRQuadrant *t = nodes->contains(&node);
        int index = t - node_array;
        c[4] = node_offset[index];
      }
    }
  }
  else {
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
        int index = t - node_array;
        c[offset] = node_offset[index];
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
        int index = t - node_array;

        // The owner edge is reversed relative to this edge
        if (edge_reversed){
          if (edge_index < 2){
            for ( int k = 1; k < mesh_order-1; k++ ){
              int offset = k*mesh_order + (mesh_order-1)*edge_index;
              c[offset] = node_offset[index] + mesh_order-2-k;
            }
          }
          else {
            for ( int k = 1; k < mesh_order-1; k++ ){
              int offset = k + (mesh_order-1)*mesh_order*(edge_index % 2);
              c[offset] = node_offset[index] + mesh_order-2-k;
            }
          }
        }
        else {
          if (edge_index < 2){
            for ( int k = 1; k < mesh_order-1; k++ ){
              int offset = k*mesh_order + (mesh_order-1)*edge_index;
              c[offset] = node_offset[index] + k-1;
            }
          }
          else {
            for ( int k = 1; k < mesh_order-1; k++ ){
              int offset = k + (mesh_order-1)*mesh_order*(edge_index % 2);
              c[offset] = node_offset[index] + k-1;
            }
          }
        }
      }

      // Loop over the face owners
      TMRQuadrant face;
      face.face = quads[i].face;
      face.info = face_label;
      face.x = quads[i].x + h;
      face.y = quads[i].y + h;
      transformNode(&face);
      TMRQuadrant *t = nodes->contains(&face);
      int index = t - node_array;

      for ( int jj = 1; jj < mesh_order-1; jj++ ){
        for ( int ii = 1; ii < mesh_order-1; ii++ ){
          int offset = ii + jj*mesh_order;
          c[offset] = node_offset[index] + (ii-1) + (jj-1)*(mesh_order-2);
        }
      }
    }
  }
}

/*
  Create the dependent mesh information for all local dependent
  nodes.

  output:
  ptr:      pointer for each dependent node number
  conn:     connectivity to each (global) independent node
  weights:  the weight values for each dependent node
*/
void TMRQuadForest::createDependentConn( const int *node_nums,
                                         TMRQuadrantArray *nodes,
                                         const int *node_offset ){
  // Allocate space for the connectivity
  dep_ptr = new int[ num_dep_nodes+1 ];
  for ( int k = 0; k < num_dep_nodes+1; k++ ){
    dep_ptr[k] = k*mesh_order;
  }

  // Allocate the space for the node numbers
  dep_conn = new int[ mesh_order*num_dep_nodes ];
  dep_weights = new double[ mesh_order*num_dep_nodes ];

  // Get the quadrants
  int num_elements;
  TMRQuadrant *quads;
  quadrants->getArray(&quads, &num_elements);

  // Get the quadrants
  int node_size;
  TMRQuadrant *node_array;
  nodes->getArray(&node_array, &node_size);

  // Set the node, edge and face label
  int node_label, edge_label, face_label;
  node_label = edge_label = face_label = TMR_QUAD_NODE_LABEL;

  int label_type[3];
  // If the mesh order is high enough, we will have multiple nodes
  // per edge/face
  initLabel(mesh_order, interp_type, label_type);

  node_label = label_type[0];
  edge_label = label_type[1];
  face_label = label_type[2];

  // Allocate space to store the free node variables
  int *edge_nodes = new int[ mesh_order ];

  for ( int i = 0; i < num_elements; i++ ){
    if (quads[i].info){
      // Set the edge length based on the quadrant
      const int32_t h = 1 << (TMR_MAX_LEVEL - quads[i].level);

      for ( int edge_index = 0; edge_index < 4; edge_index++ ){
        if (quads[i].info & 1 << edge_index){
          TMRQuadrant parent;
          quads[i].parent(&parent);

          // Get the local edge nodes
          if (mesh_order == 2){
            // Transform from the global ordering to the local ordering
            for ( int ii = 0; ii < 2; ii++ ){
              TMRQuadrant node;
              node.face = parent.face;
              if (edge_index < 2){
                node.x = parent.x + 2*h*(edge_index % 2);
                node.y = parent.y + 2*h*ii;
              }
              else {
                node.x = parent.x + 2*h*ii;
                node.y = parent.y + 2*h*(edge_index % 2);
              }
              node.info = node_label;
              transformNode(&node);
              TMRQuadrant *t = nodes->contains(&node);
              int index = t - node_array;
              edge_nodes[ii] = node_offset[index];
            }
          }
          else {
            // Transform from the global ordering to the local ordering
            for ( int ii = 0; ii < 3; ii++ ){
              TMRQuadrant node;
              node.face = parent.face;
              if (edge_index < 2){
                node.x = parent.x + 2*h*(edge_index % 2);
                node.y = parent.y + h*ii;
              }
              else {
                node.x = parent.x + h*ii;
                node.y = parent.y + 2*h*(edge_index % 2);
              }
              if (ii == 0 || ii == 2){
                node.info = node_label;
                transformNode(&node);
                TMRQuadrant *t = nodes->contains(&node);
                int index = t - node_array;
                edge_nodes[(ii/2)*(mesh_order-1)] = node_offset[index];
              }
              else {
                int reversed = 0;
                node.info = edge_label;
                transformNode(&node, &reversed);
                TMRQuadrant *t = nodes->contains(&node);
                int index = t - node_array;
                if (reversed){
                  for ( int k = 1; k < mesh_order-1; k++ ){
                    edge_nodes[k] = node_offset[index] + mesh_order-2-k;
                  }
                }
                else {
                  for ( int k = 1; k < mesh_order-1; k++ ){
                    edge_nodes[k] = node_offset[index] + k-1;
                  }
                }
              }
            }
          }

          // Set the offset into the local connectivity array
          const int *c = &conn[mesh_order*mesh_order*i];
          if (interp_type == TMR_BERNSTEIN_POINTS){
            for ( int k = 0; k < mesh_order; k++ ){
              // Compute the offset to the local edge
              int offset = 0;
              if (edge_index < 2){
                offset = k*mesh_order + (mesh_order-1)*edge_index;
              }
              else {
                offset = k + (mesh_order-1)*mesh_order*(edge_index % 2);
              }

              // If it's a negative number, it's a dependent node
              // whose interpolation must be set
              int index = node_nums[c[offset]];
              if (index < 0){
                index = -index-1;

                // Compute the shape functions
                int ptr = dep_ptr[index];
                for ( int j = 0; j < mesh_order; j++ ){
                  dep_conn[ptr + j] = edge_nodes[j];
                }

                // Compute parametric location along the edge
                int u = 0.0;
                if (edge_index < 2){
                  u = (mesh_order-1)*(-1 + (quads[i].childId()/2)) + k;
                }
                else {
                  u = (mesh_order-1)*(-1 + (quads[i].childId() % 2)) + k;
                }

                // Evaluate dependent weights
                eval_bernstein_weights(mesh_order, u, &dep_weights[ptr]);
              }
            }
          }
          else {
            for ( int k = 0; k < mesh_order; k++ ){
              // Compute the offset to the local edge
              int offset = 0;
              if (edge_index < 2){
                offset = k*mesh_order + (mesh_order-1)*edge_index;
              }
              else {
                offset = k + (mesh_order-1)*mesh_order*(edge_index % 2);
              }

              // If it's a negative number, it's a dependent node
              // whose interpolation must be set
              int index = node_nums[c[offset]];
              if (index < 0){
                index = -index-1;

                // Compute the shape functions
                int ptr = dep_ptr[index];
                for ( int j = 0; j < mesh_order; j++ ){
                  dep_conn[ptr + j] = edge_nodes[j];
                }

                // Compute parametric location along the edge
                double u = 0.0;
                if (edge_index < 2){
                  u = 1.0*(-1 + (quads[i].childId()/2)) +
                    0.5*(1.0 + interp_knots[k]);
                }
                else {
                  u = 1.0*(-1 + (quads[i].childId() % 2)) +
                    0.5*(1.0 + interp_knots[k]);
                }

                // Evaluate the shape functions
                lagrange_shape_functions(mesh_order, u, interp_knots,
                                         &dep_weights[ptr]);
              }
            }
          }
        }
      }
    }
  }

  delete [] edge_nodes;
}

/*
  Evaluate the node locations based on the parametric locations
*/
void TMRQuadForest::evaluateNodeLocations(){
  // Allocate the array of locally owned nodes
  X = new TMRPoint[ num_local_nodes ];
  memset(X, 0, num_local_nodes*sizeof(TMRPoint));

  int *flags = new int[ num_local_nodes ];
  memset(flags, 0, num_local_nodes*sizeof(int));

  int num_elements;
  TMRQuadrant *quads;
  quadrants->getArray(&quads, &num_elements);

  // Set the knots to use in the interpolation
  const double *knots = interp_knots;

  if (topo && (interp_type == TMR_BERNSTEIN_POINTS && mesh_order > 2)){
    // Form the interpolating matrix
    int size = mesh_order*mesh_order;
    double *interp = new double[ size*size ];

    // For each point within the mesh, evaluate the shape functions
    for ( int i = 0; i < size; i++ ){
      double pt[2];
      pt[0] = knots[i % mesh_order];
      pt[1] = knots[i / mesh_order];

      // Evaluate the interpolant
      evalInterp(pt, &interp[size*i]);
    }

    // Compute the inverse for interpolation
    int info = 0;
    int *ipiv = new int[ size ];
    TmrLAPACKdgetrf(&size, &size, interp, &size, ipiv, &info);

    // Apply the factoriziation to the inverse
    double *inverse = new double[ size*size ];
    memset(inverse, 0, size*size*sizeof(double));
    for ( int i = 0; i < size; i++ ){
      inverse[(size+1)*i] = 1.0;
    }

    // Compute the inverse -- note that the transpose is used here
    // since the interpolation matrix is stored in a row-major order
    // not column-major order.
    TmrLAPACKdgetrs("N", &size, &size, interp, &size, ipiv,
                    inverse, &size, &info);

    // Free the rest of the data
    delete [] ipiv;
    delete [] interp;

    // Allocate storage for the element node locations
    TMRPoint *Xtmp = new TMRPoint[ size ];

    for ( int i = 0; i < num_elements; i++ ){
      // Get the right surface
      TMRFace *surf;
      topo->getFace(quads[i].face, &surf);

      // Compute the edge length
      const int32_t h = 1 << (TMR_MAX_LEVEL - quads[i].level);

      // Compute the origin of the element in parametric space
      // and the edge length of the element
      double d = convert_to_coordinate(h);
      double u = convert_to_coordinate(quads[i].x);
      double v = convert_to_coordinate(quads[i].y);

      for ( int jj = 0; jj < mesh_order; jj++ ){
        for ( int ii = 0; ii < mesh_order; ii++ ){
          int local_index = ii + jj*mesh_order;
          surf->evalPoint(u + 0.5*d*(1.0 + knots[ii]),
                          v + 0.5*d*(1.0 + knots[jj]),
                          &Xtmp[local_index]);
        }
      }

      for ( int jj = 0; jj < mesh_order; jj++ ){
        for ( int ii = 0; ii < mesh_order; ii++ ){
          // Compute the mesh index
          int local_index = ii + jj*mesh_order;
          int node = conn[mesh_order*mesh_order*i + local_index];
          int index = getLocalNodeNumber(node);
          if (!flags[index]){
            flags[index] = 1;

            X[index].x = X[index].y = X[index].z = 0.0;
            for ( int j = 0; j < size; j++ ){
              X[index].x += inverse[local_index*size + j]*Xtmp[j].x;
              X[index].y += inverse[local_index*size + j]*Xtmp[j].y;
              X[index].z += inverse[local_index*size + j]*Xtmp[j].z;
            }
          }
        }
      }
    }

    delete [] inverse;
    delete [] Xtmp;
  }
  else if (topo){
    for ( int i = 0; i < num_elements; i++ ){
      // Get the right surface
      TMRFace *surf;
      topo->getFace(quads[i].face, &surf);

      // Compute the edge length
      const int32_t h = 1 << (TMR_MAX_LEVEL - quads[i].level);

      // Compute the origin of the element in parametric space
      // and the edge length of the element
      double d = convert_to_coordinate(h);
      double u = convert_to_coordinate(quads[i].x);
      double v = convert_to_coordinate(quads[i].y);

      // Look for nodes that are not assigned
      for ( int jj = 0; jj < mesh_order; jj++ ){
        for ( int ii = 0; ii < mesh_order; ii++ ){
          // Compute the mesh index
          int node = conn[mesh_order*mesh_order*i +
                          ii + jj*mesh_order];
          int index = getLocalNodeNumber(node);
          if (!flags[index]){
            flags[index] = 1;
            surf->evalPoint(u + 0.5*d*(1.0 + knots[ii]),
                            v + 0.5*d*(1.0 + knots[jj]),
                            &X[index]);
          }
        }
      }
    }
  }

  delete [] flags;
}

/*
  Get the nodal connectivity. This can only be called after the nodes
  have been created.

  output:
  conn:             the connectivity
  num_elements:     the number of elements
  num_owned_nodes:  the number of owned nodes on this proc
*/
void TMRQuadForest::getNodeConn( const int **_conn,
                                 int *_num_elements,
                                 int *_num_owned_nodes,
                                 int *_num_local_nodes ){
  int num_elements = 0;
  if (quadrants){
    quadrants->getArray(NULL, &num_elements);
  }
  if (_conn){ *_conn = conn; }
  if (_num_elements){ *_num_elements = num_elements; }
  if (_num_owned_nodes){ *_num_owned_nodes = num_owned_nodes; }
  if (_num_owned_nodes){ *_num_owned_nodes = num_owned_nodes; }
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
  Get the elements that either lie on a face or curve with a given
  name.

  This code loops over all quadrants owned locally on this processor
  and checks if each quadrant lies on a face or boundary. If the face
  name matches, the quadrant is added without modification. If
  the quadrant lies on an edge, the quadrant is modified so that the
  info indicates which edge the quadrant lies on using the regular edge
  ordering scheme.

  input:
  name:   string name associated with the geometric feature

  returns:
  list:   an array of quadrants satisfying the name
*/
TMRQuadrantArray* TMRQuadForest::getQuadsWithName( const char *name ){
  if (!topo){
    fprintf(stderr, "TMRQuadForest Error: Must define topology to use "
            "getQuadsWithName()\n");
    return NULL;
  }
  if (!quadrants){
    fprintf(stderr, "TMRQuadForest Error: Must create quadrants to use "
            "getQuadsWithName()\n");
    return NULL;
  }

  // Create a queue to store the elements that we find
  TMRQuadrantQueue *queue = new TMRQuadrantQueue();

  // Get the quadrants
  int size;
  TMRQuadrant *array;
  quadrants->getArray(&array, &size);

  // Loop over the quadrants and find out whether it touches
  // a face or edge with the prescribed name
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  for ( int i = 0; i < size; i++ ){
    const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level);

    // Get the surface quadrant
    TMRFace *surf;
    topo->getFace(array[i].face, &surf);
    const char *face_name = surf->getName();
    if ((!name && !face_name) ||
        (face_name && strcmp(face_name, name) == 0)){
      queue->push(&array[i]);
    }
    else {
      // If this quadrant was not added from a face
      // name, check to see if it should be added
      // as an edge/curve name
      TMREdge *edge;
      if (array[i].x == 0){
        int edge_num = fdata->face_edge_conn[4*array[i].face];
        topo->getEdge(edge_num,&edge);
        const char *edge_name = edge->getName();
        if (edge_name && strcmp(edge_name, name) == 0){
          TMRQuadrant p = array[i];
          p.info = 0;
          queue->push(&p);
        }
      }
      if (array[i].x+h == hmax){
        int edge_num = fdata->face_edge_conn[4*array[i].face+1];
        topo->getEdge(edge_num, &edge);
        const char *edge_name = edge->getName();
        if (edge_name && strcmp(edge_name, name) == 0){
          TMRQuadrant p = array[i];
          p.info = 1;
          queue->push(&p);
        }
      }
      if (array[i].y == 0){
        int edge_num = fdata->face_edge_conn[4*array[i].face+2];
        topo->getEdge(edge_num, &edge);
        const char *edge_name = edge->getName();
        if (edge_name && strcmp(edge_name, name) == 0){
          TMRQuadrant p = array[i];
          p.info = 2;
          queue->push(&p);
        }
      }
      if (array[i].y+h == hmax){
        int edge_num = fdata->face_edge_conn[4*array[i].face+3];
        topo->getEdge(edge_num, &edge);
        const char *edge_name = edge->getName();
        if (edge_name && strcmp(edge_name, name) == 0){
          TMRQuadrant p = array[i];
          p.info = 3;
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
  corner with a given name

  This code loops over all nodes and check whether they lie on a
  geometric entity that has the given name. The nodes are not
  unique if they are lie on a shared boundary between processors.

  input:
  name:   the string of the name to search

  returns:
  list:   the nodes matching the specified name
*/
int TMRQuadForest::getNodesWithName( const char *name,
                                     int **_nodes ){
  if (!topo){
    fprintf(stderr, "TMRQuadForest Error: Must define topology to use "
            "getNodesWithName()\n");
    *_nodes = NULL;
    return 0;
  }
  if (!conn){
    fprintf(stderr, "TMRQuadForest Error: Nodes must be created before "
            "calling getNodesWithName()\n");
    *_nodes = NULL;
    return 0;
  }

  // The maximum quadrant edge length
  const int32_t hmax = 1 << TMR_MAX_LEVEL;

  // Get the local nodal quadrants
  int size;
  TMRQuadrant *quads;
  quadrants->getArray(&quads, &size);

  int count = 0; // Current node count
  int max_len = 1024; // max length of the node list
  int *node_list = new int[ max_len ]; // Nodes touching this name

  // Max number of nodes added by one quadrant
  const int max_node_incr = 4 + 4*mesh_order + mesh_order*mesh_order;

  // Loop over the quadrants and find out whether it touches a face or
  // edge with the prescribed name
  for ( int i = 0; i < size; i++ ){
    if (count + max_node_incr > max_len){
      // Extend the length of the array
      max_len = 2*max_len + max_node_incr;
      int *tmp = new int[ max_len ];
      memcpy(tmp, node_list, count*sizeof(int));
      delete [] node_list;
      node_list = tmp;
    }

    // Compute the quadrant edge length
    const int32_t h = 1 << (TMR_MAX_LEVEL - quads[i].level);

    // Check if this node is on a corner, edge or face, and whether it
    // shares the appropriate name
    int fx0 = (quads[i].x == 0);
    int fy0 = (quads[i].y == 0);
    int fx = (fx0 || quads[i].x + h == hmax);
    int fy = (fy0 || quads[i].y + h == hmax);

    if (fx && fy){
      // Keep track of which corners this element touches
      int ncorners = 0;
      int corner_index[4];
      if (quads[i].x == 0 && quads[i].y == 0){
        corner_index[ncorners] = 0; ncorners++;
      }
      if (quads[i].x + h == hmax && quads[i].y == 0){
        corner_index[ncorners] = 1; ncorners++;
      }
      if (quads[i].x == 0 && quads[i].y + h == hmax){
        corner_index[ncorners] = 2; ncorners++;
      }
      if (quads[i].x + h == hmax && quads[i].y + h == hmax){
        corner_index[ncorners] = 3; ncorners++;
      }

      for ( int ii = 0; ii < ncorners; ii++ ){
        TMRVertex *vert;
        int vert_num = fdata->face_conn[4*quads[i].face + corner_index[ii]];
        topo->getVertex(vert_num, &vert);
        const char *vert_name = vert->getName();
        if (vert_name && strcmp(vert_name, name) == 0){
          int offset = ((mesh_order-1)*(corner_index[ii] % 2) +
                        (mesh_order-1)*mesh_order*(corner_index[ii]/2));
          node_list[count] =
            conn[mesh_order*mesh_order*quads[i].tag + offset];
          count++;
        }
      }
    }
    if (fx || fy){
      // Keep track of which edges this element touches
      int nedges = 0;
      int edge_index[4];
      if (quads[i].x == 0){
        edge_index[nedges] = 0; nedges++;
      }
      if (quads[i].x + h == hmax){
        edge_index[nedges] = 1; nedges++;
      }
      if (quads[i].y == 0){
        edge_index[nedges] = 2; nedges++;
      }
      if (quads[i].y + h == hmax){
        edge_index[nedges] = 3; nedges++;
      }

      for ( int ii = 0; ii < nedges; ii++ ){
        TMREdge *edge;
        int edge_num = fdata->face_edge_conn[4*quads[i].face + edge_index[ii]];
        topo->getEdge(edge_num, &edge);
        const char *edge_name = edge->getName();
        if (edge_name && strcmp(edge_name, name) == 0){
          for ( int k = 0; k < mesh_order; k++ ){
            int offset = 0;
            if (edge_index[ii] < 2){
              offset = k*mesh_order + (mesh_order-1)*edge_index[ii];
            }
            else {
              offset = k + (mesh_order-1)*mesh_order*(edge_index[ii] % 2);
            }
            node_list[count] =
              conn[mesh_order*mesh_order*quads[i].tag + offset];
            count++;
          }
        }
      }
    }

    // This node lies on the face
    TMRFace *face;
    topo->getFace(quads[i].face, &face);
    const char *face_name = face->getName();
    if (face_name && strcmp(face_name, name) == 0){
      for ( int jj = 0; jj < mesh_order; jj++ ){
        for ( int ii = 0; ii < mesh_order; ii++ ){
          int offset = ii + jj*mesh_order;
          node_list[count] =
            conn[mesh_order*mesh_order*quads[i].tag + offset];
          count++;
        }
      }
    }
  }

  // Now, sort the node numbers and remove duplicates
  qsort(node_list, count, sizeof(int), compare_integers);

  // Remove duplicates from the array
  int len = 0;
  for ( int ptr = 0; ptr < count; ptr++, len++ ){
    while ((ptr < count-1) && (node_list[ptr] == node_list[ptr+1])){
      ptr++;
    }

    if (ptr != len){
      node_list[len] = node_list[ptr];
    }
  }

  *_nodes = node_list;
  return len;
}

/*
  Given a node, find the enclosing quadrant

  This code is used to find the quadrant in the quadrant array that
  encloses the given node.
*/
TMRQuadrant* TMRQuadForest::findEnclosing( const int order,
                                           const double *knots,
                                           TMRQuadrant *node,
                                           int *mpi_owner ){
  // Assume that we'll find octant on this processor for now..
  if (mpi_owner){
    *mpi_owner = mpi_rank;
  }

  // Retrieve the array of elements
  int size = 0;
  TMRQuadrant *array = NULL;
  quadrants->getArray(&array, &size);

  // Set the lower and upper bounds for the quadrant
  const int32_t face = node->face;
  const int32_t h = 1 << (TMR_MAX_LEVEL - node->level);

  // Compute the ii/jj locations
  const int ii = node->info % order;
  const int jj = node->info/order;

  // Compute the integer locations for the x/y nodes if they lie
  // exactly along a coordinate line. These will take precidence over
  // the real value parametric locations since comparisons will be
  // exact.
  int32_t xi = -1, yi = -1;
  if (ii == 0 || ii == order-1){
    xi = node->x + (ii/(order-1))*h;
  }
  else if (order % 2 == 1 && ii == order/2){
    xi = node->x + h/2;
  }
  if (jj == 0 || jj == order-1){
    yi = node->y + (jj/(order-1))*h;
  }
  else if (order % 2 == 1 && jj == order/2){
    yi = node->y + h/2;
  }

  // Compute the parametric node location on this block
  const double xd = node->x + 0.5*h*(1.0 + knots[ii]);
  const double yd = node->y + 0.5*h*(1.0 + knots[jj]);

  // Set the low and high indices to the first and last
  // element of the element array
  int low = 0;
  int high = size-1;
  int mid = low + (high - low)/2;

  // Maintain values of low/high and mid such that the octant is
  // between (elems[low], elems[high]).  Note that if high-low=1, then
  // mid = low
  while (mid != low){
    // Check if the node is contained by the mid octant
    if (array[mid].contains(node)){
      break;
    }

    // Compare the ordering of the two octants - if the octant is less
    // than the other, then adjust the mid point
    int stat = array[mid].comparePosition(node);

    // array[mid] ? node
    if (stat == 0){
      break;
    }
    else if (stat < 0){
      low = mid+1;
    }
    else {
      high = mid-1;
    }

    // Re compute the mid-point and repeat
    mid = low + (int)((high - low)/2);
  }

  // Compute the bounding quadrant. Quadrants greater than this quad
  // cannot own the node so a further search is futile.
  TMRQuadrant quad;
  quad.face = face;
  quad.x = node->x + h;
  quad.y = node->y + h;

  while (mid < size && array[mid].comparePosition(&quad) <= 0){
    // Check if array[mid] contains the provided octant
    const int32_t hm = 1 << (TMR_MAX_LEVEL - array[mid].level);

    // First, make sure that we're on the right block
    if (array[mid].face == face){
      // Check the intervals. If the integers are non-negative, use
      // the integer comparison, otherwise use the double values.
      int xinterval = 0, yinterval = 0;
      if (xi >= 0){
        xinterval = (array[mid].x <= xi && xi <= array[mid].x+hm);
      }
      else {
        xinterval = (array[mid].x <= xd && xd <= array[mid].x+hm);
      }
      if (yi >= 0){
        yinterval = (array[mid].y <= yi && yi <= array[mid].y+hm);
      }
      else {
        yinterval = (array[mid].y <= yd && yd <= array[mid].y+hm);
      }

      // If all the intervals are satisfied, return the array
      if (xinterval && yinterval){
        return &array[mid];
      }
    }
    mid++;
  }

  if (mpi_owner){
    const int32_t hmax = 1 << TMR_MAX_LEVEL;
    TMRQuadrant n;
    n.face = face;
    n.x = (xi < 0 ? (int)xd : xi);
    n.y = (yi < 0 ? (int)yd : yi);

    if (n.x == 0){ n.x += 1; }
    else if (n.x == hmax){ n.x -= 1; }
    if (n.y == 0){ n.y += 1; }
    else if (n.y == hmax){ n.y -= 1; }
    *mpi_owner = getQuadrantMPIOwner(&n);
  }

  // No quadrant was found, return NULL
  return NULL;
}

/*
  Compute the interpolant from the given fine node to the
  coarse mesh quadrants.

  Note that the node is defined as a TMRQuadrant class with the
  quadrant information for the node where the info member is
  the local index of the node on the element.

  input:
  node:     the node (element quad with info = element node index)
  coarse:   the coarse TMRQuadForest object
  quad:     an enclosing quadrant on the coarse mesh
  tmp:      temporary array (must be of size 3*coarse->mesh_order)

  output:
  weights:  the index/weight pairs for the mesh
*/
int TMRQuadForest::computeElemInterp( TMRQuadrant *node,
                                      TMRQuadForest *coarse,
                                      TMRQuadrant *quad,
                                      TMRIndexWeight *weights,
                                      double *tmp ){
  // Loop over the array of nodes
  const int coarse_nodes_per_element =
    coarse->mesh_order*coarse->mesh_order;

  // Compute the i, j location of the fine mesh node on the element
  const int i = node->info % mesh_order;
  const int j = node->info/mesh_order;

  // Get the element size for coarse element
  const int32_t h = 1 << (TMR_MAX_LEVEL - node->level);
  const int32_t hc = 1 << (TMR_MAX_LEVEL - quad->level);

  // Set pointers to create the interpolation
  int istart = 0, iend = coarse->mesh_order;
  int jstart = 0, jend = coarse->mesh_order;
  double *Nu = &tmp[0];
  double *Nv = &tmp[coarse->mesh_order];

  // Check that the interpolation type between the meshes are identical
  if (interp_type != coarse->interp_type){
    fprintf(stderr, "TMRQuadForest Error: Interpolation types between "
            "meshes are not identical\n");
  }

  if (interp_type == TMR_BERNSTEIN_POINTS &&
      mesh_order - coarse->mesh_order > 1){
    fprintf(stderr, "TMRQuadForest Error: mesh order difference across "
            "grids should be 1\n");
  }

  if (interp_type == TMR_BERNSTEIN_POINTS){
    // Create the evenly spaced "Bern knots"
    const double *bern_knots = interp_knots;

    // Check whether the node is on a coarse mesh surface in either
    // the x,y,z directions
    if ((i == 0 && quad->x == node->x) ||
        (i == mesh_order-1 && quad->x == node->x + h)){
      istart = 0;
      iend = 1;
      Nu[istart] = 1.0;
    }
    else if ((i == 0 && quad->x + hc == node->x) ||
             (i == mesh_order-1 && quad->x + hc == node->x + h)){
      istart = coarse->mesh_order-1;
      iend = coarse->mesh_order;
      Nu[istart] = 1.0;
    }
    else {
      if (mesh_order == coarse->mesh_order){
        double u = -1.0 + 2.0*(node->x + 0.5*h*(1.0 + bern_knots[i]) -
                               quad->x)/hc;
        bernstein_shape_functions(mesh_order, u, Nu);
      }
      else {
        eval_bernstein_interp_weights(mesh_order, coarse->mesh_order, i, Nu);
      }
    }
    if ((j == 0 && quad->y == node->y) ||
        (j == mesh_order-1 && quad->y == node->y + h)){
      jstart = 0;
      jend = 1;
      Nv[jstart] = 1.0;
    }
    else if ((j == 0 && quad->y + hc == node->y) ||
             (j == mesh_order-1 && quad->y + hc == node->y + h)){
      jstart = coarse->mesh_order-1;
      jend = coarse->mesh_order;
      Nv[jstart] = 1.0;
    }
    else {
      if (mesh_order == coarse->mesh_order){
        double v = -1.0 + 2.0*(node->y + 0.5*h*(1.0 + bern_knots[j]) -
                               quad->y)/hc;
        bernstein_shape_functions(mesh_order, v, Nv);
      }
      else{
        eval_bernstein_interp_weights(mesh_order, coarse->mesh_order, j, Nv);
      }
    }
  }
  else {
    // Check whether the node is on a coarse mesh surface in either
    // the x,y,z directions
    if ((i == 0 && quad->x == node->x) ||
        (i == mesh_order-1 && quad->x == node->x + h)){
      istart = 0;
      iend = 1;
      Nu[istart] = 1.0;
    }
    else if ((i == 0 && quad->x + hc == node->x) ||
             (i == mesh_order-1 && quad->x + hc == node->x + h)){
      istart = coarse->mesh_order-1;
      iend = coarse->mesh_order;
      Nu[istart] = 1.0;
    }
    else {
      double u = -1.0 + 2.0*(node->x + 0.5*h*(1.0 + interp_knots[i]) -
                             quad->x)/hc;
      lagrange_shape_functions(coarse->mesh_order, u,
                               coarse->interp_knots, Nu);
    }
    if ((j == 0 && quad->y == node->y) ||
        (j == mesh_order-1 && quad->y == node->y + h)){
      jstart = 0;
      jend = 1;
      Nv[jstart] = 1.0;
    }
    else if ((j == 0 && quad->y + hc == node->y) ||
             (j == mesh_order-1 && quad->y + hc == node->y + h)){
      jstart = coarse->mesh_order-1;
      jend = coarse->mesh_order;
      Nv[jstart] = 1.0;
    }
    else {
      double v = -1.0 + 2.0*(node->y + 0.5*h*(1.0 + interp_knots[j]) -
                             quad->y)/hc;
      lagrange_shape_functions(coarse->mesh_order, v,
                               coarse->interp_knots, Nv);
    }
  }

  // Get the coarse grid information
  const int *cdep_ptr;
  const int *cdep_conn;
  const double *cdep_weights;
  coarse->getDepNodeConn(&cdep_ptr, &cdep_conn, &cdep_weights);

  // Get the coarse connectivity array
  const int num = quad->tag;
  const int *c = &(coarse->conn[coarse_nodes_per_element*num]);

  // Loop over the nodes that are within this octant
  int nweights = 0;
  for ( int jj = jstart; jj < jend; jj++ ){
    for ( int ii = istart; ii < iend; ii++ ){
      // Compute the offset into the coarse mesh
      int offset = ii + jj*coarse->mesh_order;

      // Compute the interpolation weight
      double weight = Nu[ii]*Nv[jj];

      // Get the tag number
      if (c[offset] >= 0){
        weights[nweights].index = c[offset];
        weights[nweights].weight = weight;
        nweights++;
      }
      else {
        int node = -c[offset]-1;
        for ( int jp = cdep_ptr[node]; jp < cdep_ptr[node+1]; jp++ ){
          weights[nweights].index = cdep_conn[jp];
          weights[nweights].weight = weight*cdep_weights[jp];
          nweights++;
        }
      }
    }
  }

  // Sort and add up the weight values
  nweights = TMRIndexWeight::uniqueSort(weights, nweights);

  return nweights;
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
  // Ensure that the nodes are allocated on both octree forests
  createNodes();
  coarse->createNodes();

  // Get the dependent node information
  const int *cdep_ptr, *cdep_conn;
  const double *cdep_weights;
  coarse->getDepNodeConn(&cdep_ptr, &cdep_conn, &cdep_weights);

  // First, loop over the local list
  int local_size = node_range[mpi_rank+1] - node_range[mpi_rank];
  int *flags = new int[ local_size ];
  memset(flags, 0, local_size*sizeof(int));

  // Allocate additional space for the interpolation
  double *tmp = new double[ 2*coarse->mesh_order ];

  // The interpolation variables/weights on the coarse mesh
  const int order = coarse->mesh_order;
  int max_nodes = order*order;
  int *vars = new int[ max_nodes ];
  double *wvals = new double[ max_nodes ];

  // Maximum number of weights
  int max_weights = order*order*order*order;
  TMRIndexWeight *weights = new TMRIndexWeight[ max_weights ];

  // Loop over the array of nodes
  const int nodes_per_element = mesh_order*mesh_order;

  // Get the quadrants
  int num_elements;
  TMRQuadrant *quads;
  quadrants->getArray(&quads, &num_elements);

  // Allocate a queue to store the nodes that are on other procs
  TMRQuadrantQueue *ext_queue = new TMRQuadrantQueue();

  // Set the knots to use in the interpolation
  const double *knots = interp_knots;

  for ( int i = 0; i < num_elements; i++ ){
    const int *c = &conn[nodes_per_element*i];
    for ( int j = 0; j < nodes_per_element; j++ ){
      // Check if the fine node is owned by this processor
      if (c[j] >= node_range[mpi_rank] &&
          c[j] < node_range[mpi_rank+1]){
        int index = c[j] - node_range[mpi_rank];
        if (!flags[index]){
          // We're going to handle this node now, mark it as done
          flags[index] = 1;

          // Find the enclosing coarse quad on this
          // processor if it exits
          TMRQuadrant node = quads[i];
          node.info = j;

          // Find the MPI owner or the
          int mpi_owner = mpi_rank;
          TMRQuadrant *t = coarse->findEnclosing(mesh_order, knots,
                                                 &node, &mpi_owner);

          // The node is owned a coarse element on this processor
          if (t){
            // Compute the element interpolation
            int nweights = computeElemInterp(&node, coarse, t, weights, tmp);

            for ( int k = 0; k < nweights; k++ ){
              vars[k] = weights[k].index;
              wvals[k] = weights[k].weight;
            }
            interp->addInterp(c[j], wvals, vars, nweights);
          }
          else {
            // We've got to transfer the node to the processor that
            // owns an enclosing element. To do that, add the quad to
            // the list of externals and store its mpi owner
            node.tag = mpi_owner;
            ext_queue->push(&node);
          }
        }
      }
    }
  }

  // Free the data
  delete [] flags;

  // Sort the sending quadrants by MPI rank
  TMRQuadrantArray *ext_array = ext_queue->toArray();
  delete ext_queue;

  // Sort the node
  int size;
  TMRQuadrant *array;
  ext_array->getArray(&array, &size);
  qsort(array, size, sizeof(TMRQuadrant), compare_quadrant_tags);

  // The number of quadrants that will be sent from this processor
  // to all other processors in the communicator
  int *quad_ptr = new int[ mpi_size+1 ];
  int *quad_recv_ptr = new int[ mpi_size+1 ];

  // Match the quad intervals to determine how mnay quad
  // need to be sent to each processor
  matchTagIntervals(array, size, quad_ptr);

  // Now convert the node tags nodes back to node numbers
  // from the connectivity
  for ( int i = 0; i < size; i++ ){
    // Search for the quad in the quadrants array
    TMRQuadrant *t = quadrants->contains(&array[i]);

    // Set the tag value as the global node number
    array[i].tag = conn[nodes_per_element*t->tag + array[i].info];
  }

  // Count up the number of quad destined for other procs
  int *quad_counts = new int[ mpi_size ];
  for ( int i = 0; i < mpi_size; i++ ){
    if (i == mpi_rank){
      quad_counts[i] = 0;
    }
    else {
      quad_counts[i] = quad_ptr[i+1] - quad_ptr[i];
    }
  }

  // Now distribute the quad to their destination processors
  int *quad_recv_counts = new int[ mpi_size ];
  MPI_Alltoall(quad_counts, 1, MPI_INT,
               quad_recv_counts, 1, MPI_INT, comm);

  // Now use oct_ptr to point into the recv array
  quad_recv_ptr[0] = 0;
  for ( int i = 0; i < mpi_size; i++ ){
    quad_recv_ptr[i+1] = quad_recv_ptr[i] + quad_recv_counts[i];
  }

  delete [] quad_counts;
  delete [] quad_recv_counts;

  // Distribute the quad based on the oct_ptr/oct_recv_ptr arrays
  TMRQuadrantArray *recv_array =
    sendQuadrants(ext_array, quad_ptr, quad_recv_ptr);
  delete [] quad_ptr;
  delete [] quad_recv_ptr;
  delete ext_array;

  // Get the nodes recv'd from other processors
  int recv_size;
  TMRQuadrant *recv_nodes;
  recv_array->getArray(&recv_nodes, &recv_size);

  // Recv the nodes and loop over the connectivity
  for ( int i = 0; i < recv_size; i++ ){
    TMRQuadrant *t = coarse->findEnclosing(mesh_order, knots,
                                           &recv_nodes[i]);
    if (t){
      // Compute the element interpolation
      int nweights = computeElemInterp(&recv_nodes[i], coarse, t,
                                       weights, tmp);

      for ( int k = 0; k < nweights; k++ ){
        vars[k] = weights[k].index;
        wvals[k] = weights[k].weight;
      }
      interp->addInterp(recv_nodes[i].tag, wvals, vars, nweights);
    }
    else {
      // This should not happen. Print out an error message here.
      fprintf(stderr, "[%d] TMRQuadForest Error: Destination processor does "
              "not own node\n", mpi_rank);
    }
  }

  // Free the recv array
  delete recv_array;

  // Free the temporary arrays
  delete [] tmp;
  delete [] vars;
  delete [] wvals;
  delete [] weights;
}

/*
  Initialize the node label
*/
void TMRQuadForest::initLabel( int mesh_order, TMRInterpolationType interp_type,
                               int label_type[] ){
  if (mesh_order > 3 ||
      (mesh_order >= 3 && interp_type == TMR_BERNSTEIN_POINTS)){
    label_type[0] = TMR_QUAD_NODE_LABEL;
    label_type[1] = TMR_QUAD_EDGE_LABEL;
    label_type[2] = TMR_QUAD_FACE_LABEL;
  }
  else {
    for (int i = 0; i < 3; i++){
      label_type[i] = TMR_QUAD_NODE_LABEL;
    }
  }
}
