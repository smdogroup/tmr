#include "TMROctForest.h"
#include "TACSMeshLoader.h"

/*
  Bottom surface      Top surface
  12-------- 14       13 ------- 15
  | \      / |        | \      / |
  |  2 -- 3  |        |  6 -- 7  |
  |  |    |  |        |  |    |  |
  |  0 -- 1  |        |  4 -- 5  |
  | /      \ |        | /      \ |
  8 -------- 10       9 -------- 11
*/

const double box_xpts[] = 
  {-.5, -.5, -.5,
   .5, -.5, -.5,
   -.5, .5, -.5,
   .5, .5, -.5,
   -.5, -.5, .5,
   .5, -.5, .5,
   -.5, .5, .5,
   .5, .5, .5,
   -1, -1, -1,
   -1, -1, 1,
   1, -1, -1,
   1, -1, 1,
   -1, 1, -1,
   -1, 1, 1,
   1, 1, -1,
   1, 1, 1};

const int box_conn[] =
  {0, 1, 2, 3, 4, 5, 6, 7,
   8, 10, 0, 1, 9, 11, 4, 5,
   5, 11, 1, 10, 7, 15, 3, 14,
   7, 15, 3, 14, 6, 13, 2, 12,
   9, 13, 4, 6, 8, 12, 0, 2, 
   10, 14, 8, 12, 1, 3, 0, 2,
   4, 5, 6, 7, 9, 11, 13, 15};

/*
  Interpoalte from the connectivity/node locations
*/
void getLocation( int i, const int *elem_node_conn, 
                  const double *Xpts,
                  const TMROctant *oct, TMRPoint *pt ){
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  double u = (1.0*oct->x)/hmax;
  double v = (1.0*oct->y)/hmax;
  double w = (1.0*oct->z)/hmax;

  double N[8];
  N[0] = (1.0 - u)*(1.0 - v)*(1.0 - w);
  N[1] = u*(1.0 - v)*(1.0 - w);
  N[2] = (1.0 - u)*v*(1.0 - w);
  N[3] = u*v*(1.0 - w);
  N[4] = (1.0 - u)*(1.0 - v)*w;
  N[5] = u*(1.0 - v)*w;
  N[6] = (1.0 - u)*v*w;
  N[7] = u*v*w;

  pt->x = pt->y = pt->z = 0.0;
  for ( int k = 0; k < 8; k++ ){
    int node = elem_node_conn[8*i + k];
    
    pt->x += Xpts[3*node]*N[k];
    pt->y += Xpts[3*node+1]*N[k];
    pt->z += Xpts[3*node+2]*N[k];
  }
}

void computeShapeDeriv( double u, double v, double w,
                        double Na[], double Nb[], double Nc[] ){
  Na[0] = -(1.0 - v)*(1.0 - w);
  Na[1] = (1.0 - v)*(1.0 - w);
  Na[2] = -v*(1.0 - w);
  Na[3] = v*(1.0 - w);
  Na[4] = -(1.0 - v)*w;
  Na[5] = (1.0 - v)*w;
  Na[6] = -v*w;
  Na[7] = v*w;

  Nb[0] = -(1.0 - u)*(1.0 - w);
  Nb[1] = -u*(1.0 - w);
  Nb[2] = (1.0 - u)*(1.0 - w);
  Nb[3] = u*(1.0 - w);
  Nb[4] = -(1.0 - u)*w;
  Nb[5] = -u*w;
  Nb[6] = (1.0 - u)*w;
  Nb[7] = u*w;

  Nc[0] = -(1.0 - u)*(1.0 - v);
  Nc[1] = -u*(1.0 - v);
  Nc[2] = -(1.0 - u)*v;
  Nc[3] = -u*v;
  Nc[4] = (1.0 - u)*(1.0 - v);
  Nc[5] = u*(1.0 - v);
  Nc[6] = (1.0 - u)*v;
  Nc[7] = u*v;
}

/*
  Check volume of the mesh to see if it is valid - is there a better
  way to do this?
*/
double computeVolume( int i, const int *elem_node_conn, 
                      const double *Xpts ){
  const double pt = 1.0/sqrt(3.0);
  double V = 0.0;

  for ( int kk = 0; kk < 2; kk++ ){
    for ( int jj = 0; jj < 2; jj++ ){
      for ( int ii = 0; ii < 2; ii++ ){
        double u = 0.5 + (ii-0.5)*pt;
        double v = 0.5 + (jj-0.5)*pt;
        double w = 0.5 + (kk-0.5)*pt;

        double Na[8], Nb[8], Nc[8];
        computeShapeDeriv(u, v, w, Na, Nb, Nc);

        double Xd[9];
        memset(Xd, 0, 9*sizeof(double));
        for ( int k = 0; k < 8; k++ ){
          int node = elem_node_conn[8*i + k];
          Xd[0] += Xpts[3*node]*Na[k];
          Xd[3] += Xpts[3*node+1]*Na[k];
          Xd[6] += Xpts[3*node+2]*Na[k];
          
          Xd[1] += Xpts[3*node]*Nb[k];
          Xd[4] += Xpts[3*node+1]*Nb[k];
          Xd[7] += Xpts[3*node+2]*Nb[k];
          
          Xd[2] += Xpts[3*node]*Nc[k];
          Xd[5] += Xpts[3*node+1]*Nc[k];
          Xd[8] += Xpts[3*node+2]*Nc[k];
        }

        V += 0.125*(Xd[8]*(Xd[0]*Xd[4] - Xd[3]*Xd[1]) - 
                    Xd[7]*(Xd[0]*Xd[5] - Xd[3]*Xd[2]) +
                    Xd[6]*(Xd[1]*Xd[5] - Xd[2]*Xd[4]));
      }
    }
  }

  return V;
}

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  int partition = 0;
  int box_problem = 0;
  for ( int k = 0; k < argc; k++ ){
    if (strcmp(argv[k], "partition") == 0){
      partition = 1;
    }
    if (strcmp(argv[k], "box") == 0){
      box_problem = 1;
    }
  }

  // Define the different forest levels
  MPI_Comm comm = MPI_COMM_WORLD;
  const int MAX_NUM_MESH = 5;
  TMROctForest *forest[MAX_NUM_MESH];

  // The dependent node information and the interpolation
  int *cdep_ptr[MAX_NUM_MESH], *cdep_conn[MAX_NUM_MESH];
  double *cdep_weights[MAX_NUM_MESH];

  // Create the forests
  forest[0] = new TMROctForest(comm);

  // Get the MPI rank
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  // The "super-node" locations
  const double *Xpts = NULL;
  const int *elem_node_conn = NULL;
  if (box_problem){
    int nnodes = 16;
    int nelems = 7;
    Xpts = box_xpts;
    elem_node_conn = box_conn;

    forest[0]->setConnectivity(nnodes, box_conn,
                               nelems, partition);
    forest[0]->createRandomTrees(50, 0, 10);
  }
  else {
    // Create the TACSMeshLoader class
    TACSMeshLoader *mesh = new TACSMeshLoader(MPI_COMM_SELF);
    mesh->incref();
    mesh->scanBDFFile("uCRM_3D_box_mesh.bdf");
    
    // Extract the connectivity
    int nnodes, nelems;
    mesh->getConnectivity(&nnodes, &nelems, NULL, 
                          &elem_node_conn, &Xpts);
    forest[0]->setConnectivity(nnodes, elem_node_conn, 
                               nelems, partition);
    
    // Set the refinement increasing out the wing
    int *refine = new int[ nelems ];
    memset(refine, 0, nelems*sizeof(int));
    
    int max_refine = 5;
    int min_refine = 2;
    double y_max = 30.0;
    for ( int k = 0; k < nelems; k++ ){
      double y_ref = Xpts[3*elem_node_conn[8*k]+1];
      refine[k] = min_refine + 
        (max_refine - min_refine)*(1.0 - (y_ref/y_max));
    }
  
    forest[0]->createTrees(refine);
    delete [] refine;
  }

  if (mpi_rank == 0){
    // Get the face ids and count up how many of each we have
    int nblocks, nfaces, nedges, nnodes;
    const int *face_ids;
    const int *block_faces;
    forest[0]->getConnectivity(&nblocks, &nfaces, &nedges, &nnodes,
                               NULL, &block_faces, NULL, &face_ids);

    // Check if any of the blocks have element
    for ( int i = 0; i < nblocks; i++ ){
      double V = computeVolume(i, elem_node_conn, Xpts);
      if (V < 0.0){
        printf("Negative volume in element %d\n", i);
      }
    }

    // Count up the face ids
    int face_id_count[8] = {0, 0, 0, 0, 
                            0, 0, 0, 0};
    for ( int k = 0; k < 6*nblocks; k++ ){
      if (face_ids[k] >= 0){
        face_id_count[face_ids[k]]++;
      }
    }

    printf("nblocks = %d\nnfaces = %d\nnedges = %d\nnnodes = %d\n", 
           nblocks, nfaces, nedges, nnodes);
    // Print out the face connectivity
    for ( int k = 0; k < nblocks; k++ ){
      printf("block[%d] = %d %d  %d %d  %d %d\n",
             k, block_faces[6*k], block_faces[6*k+1], 
             block_faces[6*k+2], block_faces[6*k+3], 
             block_faces[6*k+4], block_faces[6*k+5]);
    }
    // Print out the face id counts
    for ( int k = 0; k < 8; k++ ){
      printf("face_id_count[%d] = %d\n", k, face_id_count[k]);
    }
  }

  // Repartition the octrees
  printf("[%d] Repartition\n", mpi_rank);
  forest[0]->repartition();

  for ( int level = 0; level < MAX_NUM_MESH; level++ ){  
    printf("[%d] Balance\n", mpi_rank);
    double tbal = MPI_Wtime();
    forest[level]->balance((level == 0));
    tbal = MPI_Wtime() - tbal;

    printf("[%d] Create nodes\n", mpi_rank);
    double tnodes = MPI_Wtime();
    forest[level]->createNodes(2);
    tnodes = MPI_Wtime() - tnodes;

    // Create the mesh
    double tmesh = MPI_Wtime();
    int *conn, nfe = 0;
    forest[level]->createMeshConn(&conn, &nfe);
    tmesh = MPI_Wtime() - tmesh;

    // Create the local dependent node information
    forest[level]->createDepNodeConn(&cdep_ptr[level], &cdep_conn[level],
                                     &cdep_weights[level]);


    if (level > 0){
      int *ptr, *conn;
      double *weights;
      forest[level-1]->createInterpolation(forest[level],
                                           cdep_ptr[level], 
                                           cdep_conn[level],
                                           cdep_weights[level],
                                           &ptr, &conn, &weights);
      delete [] ptr;
      delete [] conn;
      delete [] weights;
    }


    int ntotal = 0;
    MPI_Allreduce(&nfe, &ntotal, 1, MPI_INT, MPI_SUM, comm);

    // Get the rank
    if (mpi_rank == 0){
      printf("balance:  %15.5f s\n", tbal);
      printf("nodes:    %15.5f s\n", tnodes);
      printf("mesh:     %15.5f s\n", tmesh);
      printf("nelems:   %15d\n", ntotal);
    }
  
    // Get the octrees within the forest
    TMROctree **octrees;
    int ntrees = forest[level]->getOctrees(&octrees);

    const int *owned;
    int nowned = forest[level]->getOwnedOctrees(&owned);
    for ( int k = 0; k < nowned; k++ ){
      int block = owned[k];

      // The list of points
      TMRPoint *X;
      octrees[block]->getPoints(&X);
    
      // Get the octant nodes
      TMROctantArray *nodes;
      octrees[block]->getNodes(&nodes);
      
      // Get the array
      int size;
      TMROctant *array;
      nodes->getArray(&array, &size);

      // Loop over all the nodes
      for ( int i = 0; i < size; i++ ){
        getLocation(block, elem_node_conn, Xpts,
                    &array[i], &X[i]);
      }
    }
    if (level+1 < MAX_NUM_MESH){
      forest[level+1] = forest[level]->coarsen();
    }
  }

  // Get the octrees within the forest
  int print_level = MAX_NUM_MESH-1;
  TMROctree **octrees;
  int ntrees = forest[print_level]->getOctrees(&octrees);
  
  // Get the owned trees
  const int *owned;
  int nowned = forest[print_level]->getOwnedOctrees(&owned);

  // Write out a file for each processor - bad practice!
  char filename[128];
  sprintf(filename, "parallel%d.dat", mpi_rank);
  FILE *fp = fopen(filename, "w");

  // Write the tecplot header
  fprintf(fp, "Variables = X, Y, Z, dv\n");

  for ( int k = 0; k < nowned; k++ ){
    int block = owned[k];

    // Retrieve the node/element arrays
    TMROctantArray *elements, *nodes;
    octrees[block]->getNodes(&nodes);
    octrees[block]->getElements(&elements);

    // Get the points
    TMRPoint *X;
    int nnodes = octrees[block]->getPoints(&X);
    int nelems = octrees[block]->getNumElements();

    // Get the node array
    TMROctant *array;
    nodes->getArray(&array, NULL);

    fprintf(fp, "ZONE T=TMR%d N=%d E=%d ", block, nnodes, nelems);
    fprintf(fp, "DATAPACKING=POINT ZONETYPE=FEBRICK\n");

    // Write out this portion of the forrest
    for ( int i = 0; i < nnodes; i++ ){
      fprintf(fp, "%e %e %e %d\n", X[i].x, X[i].y, X[i].z, 
              array[i].tag);
    }

    // Get the elements
    elements->getArray(&array, &nelems);

    TMROctant *node_array;
    nodes->getArray(&node_array, NULL);

    for ( int i = 0; i < nelems; i++ ){
      int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level);
      int index[8];

      for ( int kk = 0; kk < 2; kk++ ){
        for ( int jj = 0; jj < 2; jj++ ){
          for ( int ii = 0; ii < 2; ii++ ){
            TMROctant oct;
            oct.x = array[i].x + ii*h;
            oct.y = array[i].y + jj*h;
            oct.z = array[i].z + kk*h;
            
            const int use_nodes = 1;
            TMROctant *t = nodes->contains(&oct, use_nodes);
            index[ii + 2*jj + 4*kk] = (t - node_array) + 1;
          }
        }
      }

      fprintf(fp, "%d %d %d %d %d %d %d %d\n",
              index[0], index[1], index[3], index[2],
              index[4], index[5], index[7], index[6]);
    }
  }
  
  fclose(fp);
  for ( int level = 0; level < MAX_NUM_MESH; level++ ){
    delete forest[level];
  }

  TMRFinalize();
  MPI_Finalize();
  return (0);
}
