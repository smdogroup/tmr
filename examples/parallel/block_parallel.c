#include "TMROctForest.h"
#include "TACSMeshLoader.h"
#include "TACSAssembler.h"
#include "TACSMg.h"
#include "BVecInterp.h"
#include "Solid.h"

/*
  The box problem

  Bottom surface      Top surface
  12-------- 14       13 ------- 15
  | \      / |        | \      / |
  |  2 -- 3  |        |  6 -- 7  |
  |  |    |  |        |  |    |  |
  |  0 -- 1  |        |  4 -- 5  |
  | /      \ |        | /      \ |
  8 -------- 10       9 -------- 11
*/
const int box_npts = 16;
const int box_nelems = 7;

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
  Definitions for the connector problem
*/
const int connector_npts = 52;
const int connector_nelems = 15;

const double connector_xpts[] = 
  {-0.375, -0.375, -0.125,
   0.375, -0.375, -0.125,
   -0.125, -0.125, -0.125,
   0.125, -0.125, -0.125,
   -0.125, 0.125, -0.125,
   0.125, 0.125, -0.125,
   -0.075, 0.25, -0.125,
   0.075, 0.25, -0.125,
   -0.375, 0.375, -0.125,
   0.375, 0.375, -0.125,
   -0.25, 0.475, -0.125,
   0.25, 0.475, -0.125,
   -0.25, 1.475, -0.125,
   0.25, 1.475, -0.125,
   -0.45, 1.675, -0.125,
   0.45, 1.675, -0.125,
   -0.3125, 1.875, -0.125,
   0.3125, 1.875, -0.125,
   -0.175, 1.825, -0.125,
   0.175, 1.825, -0.125,
   -0.45, 2.425, -0.125,
   0.45, 2.425, -0.125,
   -0.3125, 2.425, -0.125,
   0.3125, 2.425, -0.125,
   -0.175, 2.425, -0.125,
   0.175, 2.425, -0.125,
   -0.375, -0.375, 0.125,
   0.375, -0.375, 0.125,
   -0.125, -0.125, 0.125,
   0.125, -0.125, 0.125,
   -0.125, 0.125, 0.125,
   0.125, 0.125, 0.125,
   -0.075, 0.25, 0.125,
   0.075, 0.25, 0.125,
   -0.375, 0.375, 0.125,
   0.375, 0.375, 0.125,
   -0.25, 0.475, 0.125,
   0.25, 0.475, 0.125,
   -0.25, 1.475, 0.125,
   0.25, 1.475, 0.125,
   -0.45, 1.675, 0.125,
   0.45, 1.675, 0.125,
   -0.3125, 1.875, 0.125,
   0.3125, 1.875, 0.125,
   -0.175, 1.825, 0.125,
   0.175, 1.825, 0.125,
   -0.45, 2.425, 0.125,
   0.45, 2.425, 0.125,
   -0.3125, 2.425, 0.125,
   0.3125, 2.425, 0.125,
   -0.175, 2.425, 0.125,
   0.175, 2.425, 0.125};

const int connector_conn[] = 
  {0, 1, 2, 3, 26, 27, 28, 29,
   0, 2, 8, 4, 26, 28, 34, 30,
   3, 1, 5, 9, 29, 27, 31, 35,
   4, 5, 6, 7, 30, 31, 32, 33,
   6, 7, 10, 11, 32, 33, 36, 37,
   8, 4, 10, 6, 34, 30, 36, 32,
   7, 5, 11, 9, 33, 31, 37, 35,
   10, 11, 12, 13, 36, 37, 38, 39,
   12, 13, 18, 19, 38, 39, 44, 45,
   14, 12, 16, 18, 40, 38, 42, 44,
   13, 15, 19, 17, 39, 41, 45, 43,
   14, 16, 20, 22, 40, 42, 46, 48,
   16, 18, 22, 24, 42, 44, 48, 50,
   19, 17, 25, 23, 45, 43, 51, 49,
   17, 15, 23, 21, 43, 41, 49, 47};

/*
  Interpoalte from the connectivity/node locations
*/
void getLocation( int i, const int *elem_node_conn, 
                  const double *Xpts,
                  const TMROctant *oct, TacsScalar x[] ){
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

  x[0] = x[1] = x[2] = 0.0;
  for ( int k = 0; k < 8; k++ ){
    int node = elem_node_conn[8*i + k];
    
    x[0] += Xpts[3*node]*N[k];
    x[1] += Xpts[3*node+1]*N[k];
    x[2] += Xpts[3*node+2]*N[k];
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
  
  // The "super-node" locations
  double omega = 1.0;
  int order = 2;
  int npts = 0;
  int nelems = 0;
  const double *Xpts = NULL;
  const int *elem_node_conn = NULL;

  for ( int k = 0; k < argc; k++ ){
    if (strcmp(argv[k], "partition") == 0){
      partition = 1;
    }
    if (strcmp(argv[k], "box") == 0){
      npts = box_npts;
      nelems = box_nelems;
      Xpts = box_xpts;
      elem_node_conn = box_conn;
    }
    if (strcmp(argv[k], "connector") == 0){
      npts = connector_npts;
      nelems = connector_nelems;
      Xpts = connector_xpts;
      elem_node_conn = connector_conn;
    }
    if (sscanf(argv[k], "order=%d", &order) == 1){
      if (order < 2){ order = 2; }
      if (order > 3){ order = 3; }
    }
    if (sscanf(argv[k], "omega=%lf", &omega) == 1){
      if (omega < 0.05){ omega = 0.05; }
    }
  }

  // Define the different forest levels
  MPI_Comm comm = MPI_COMM_WORLD;
  const int MAX_NUM_MESH = 4;
  TMROctForest *forest[MAX_NUM_MESH];
  TACSAssembler *tacs[MAX_NUM_MESH];
  TACSBVecInterp *interp[MAX_NUM_MESH-1];

  // Create the forests
  forest[0] = new TMROctForest(comm);

  // Get the MPI rank
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  if (mpi_rank == 0){
    printf("order = %d\n", order);
    printf("omega = %g\n", omega);
  }

  if (Xpts && elem_node_conn){
    forest[0]->setConnectivity(npts, elem_node_conn,
                               nelems, partition);
    // forest[0]->createRandomTrees(100, 0, 10);
    forest[0]->createTrees(5);
  }
  else {
    // Create the TACSMeshLoader class
    TACSMeshLoader *mesh = new TACSMeshLoader(MPI_COMM_SELF);
    mesh->incref();
    mesh->scanBDFFile("uCRM_3D_box_mesh.bdf");
    
    // Extract the connectivity
    mesh->getConnectivity(&npts, &nelems, NULL, 
                          &elem_node_conn, &Xpts);
    forest[0]->setConnectivity(npts, elem_node_conn, 
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

    // Print out the number of blocks, faces, edges and nodes
    printf("nblocks = %d\nnfaces = %d\nnedges = %d\nnnodes = %d\n", 
           nblocks, nfaces, nedges, nnodes);

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
    forest[level]->createNodes(order);
    tnodes = MPI_Wtime() - tnodes;

    // Get the octrees within the forest
    TMROctree **octrees;
    int ntrees = forest[level]->getOctrees(&octrees);

    // Create the mesh
    double tmesh = MPI_Wtime();
    int *conn, num_elements = 0;
    forest[level]->createMeshConn(&conn, &num_elements);
    tmesh = MPI_Wtime() - tmesh;

    // Find the number of nodes for this processor
    const int *range;
    forest[level]->getOwnedNodeRange(&range);
    int num_nodes = range[mpi_rank+1] - range[mpi_rank];

    // Get the dependent node information
    const int *dep_ptr, *dep_conn;
    const double *dep_weights;
    int num_dep_nodes = 
      forest[level]->getDepNodeConn(&dep_ptr, &dep_conn,
                                    &dep_weights);

    // Create the associated TACSAssembler object
    int vars_per_node = 3;
    tacs[level] =
      new TACSAssembler(comm, vars_per_node,
                        num_nodes, num_elements,
                        num_dep_nodes);
    tacs[level]->incref();

    // Set the element ptr
    int *ptr = new int[ order*order*order*num_elements ];
    for ( int i = 0; i < num_elements+1; i++ ){
      ptr[i] = order*order*order*i;
    }
    
    // Set the element connectivity into TACSAssembler
    tacs[level]->setElementConnectivity(conn, ptr);
    delete [] conn;
    delete [] ptr;
    
    // Set the dependent node information
    tacs[level]->setDependentNodes(dep_ptr, dep_conn,
                                   dep_weights);

    // Add nodes associated with the boundary conditions
    if (octrees[0]){
      // Get the octant nodes
      TMROctantArray *nodes;
      octrees[0]->getNodes(&nodes);
      
      // Get the array
      int size;
      TMROctant *array;
      nodes->getArray(&array, &size);

      for ( int i = 0; i < size; i++ ){
        if (array[i].x == 0 && 
            (array[i].tag >= range[mpi_rank] && 
             array[i].tag < range[mpi_rank+1])){
          int node = array[i].tag;
          tacs[level]->addBCs(1, &node);
        }
      }
    }

    TacsScalar rho = 2550.0, E = 70e9, nu = 0.3;
    SolidStiffness *stiff = new SolidStiffness(rho, E, nu);
    TACSElement *solid = NULL;
    if (order == 2){
      solid = new Solid<2>(stiff);
    }
    else if (order == 3){
      solid = new Solid<3>(stiff);
    }

    // Create the elements
    TACSElement **elements = new TACSElement*[ num_elements ];
    for ( int i = 0; i < num_elements; i++ ){
      elements[i] = solid;
    }
    tacs[level]->setElements(elements);
    delete [] elements;

    // Initialize
    tacs[level]->initialize();

    if (level > 0){
      // Create the interpolation on the TMR side
      int *ptr, *conn;
      double *weights;
      forest[level-1]->createInterpolation(forest[level],
                                           &ptr, &conn, &weights);

      // Create the interpolation object
      interp[level-1] = new TACSBVecInterp(tacs[level]->getVarMap(),
                                           tacs[level-1]->getVarMap(),
                                           tacs[level]->getVarsPerNode());
      interp[level-1]->incref();

      // Get the range of nodes
      const int *node_range;
      forest[level-1]->getOwnedNodeRange(&node_range);

      // Add all the values in the interpolation
      for ( int node = node_range[mpi_rank], k = 0;
            node < node_range[mpi_rank+1]; node++, k++ ){
        int len = ptr[k+1] - ptr[k];
        interp[level-1]->addInterp(node, &weights[ptr[k]], 
                                   &conn[ptr[k]], len);
      }

      interp[level-1]->initialize();

      delete [] ptr;
      delete [] conn;
      delete [] weights;
    }

    // Get the rank
    if (mpi_rank == 0){
      printf("balance:  %15.5f s\n", tbal);
      printf("nodes:    %15.5f s\n", tnodes);
      printf("mesh:     %15.5f s\n", tmesh);
    }
  
    // Create the nodal vector
    TACSBVec *Xvec = tacs[level]->createNodeVec();
    Xvec->incref();

    TacsScalar *X = NULL;
    Xvec->getArray(&X);

    // Scan through the nodes
    const int *owned;
    int nowned = forest[level]->getOwnedOctrees(&owned);
    for ( int k = 0; k < nowned; k++ ){
      int block = owned[k];
    
      // Get the octant nodes
      TMROctantArray *nodes;
      octrees[block]->getNodes(&nodes);
      
      // Get the array
      int size;
      TMROctant *array;
      nodes->getArray(&array, &size);

      // Loop over all the nodes
      for ( int i = 0; i < size; i++ ){
        if (array[i].tag >= 0 && 
            (array[i].tag >= range[mpi_rank] && 
             array[i].tag < range[mpi_rank+1])){
          int index = array[i].tag - range[mpi_rank];
          getLocation(block, elem_node_conn, Xpts,
                      &array[i], &X[3*index]);
        }
      }
    }

    // Set the nodal vector
    tacs[level]->setNodes(Xvec);
    Xvec->decref();

    if (level+1 < MAX_NUM_MESH){
      if (order == 3){
        // Duplicate the forest for a lower-order mesh
        forest[level+1] = forest[level]->duplicate();
        order = 2;
      }
      else {
        forest[level+1] = forest[level]->coarsen();
      }
    }
  }

  // Create the multigrid object
  int sor_iters = 1;
  int sor_symm = 0;
  TACSMg *mg = new TACSMg(comm, MAX_NUM_MESH, 
                          omega, sor_iters, sor_symm);
  mg->incref();

  for ( int level = 0; level < MAX_NUM_MESH; level++ ){
    if (level < MAX_NUM_MESH-1){
      mg->setLevel(level, tacs[level], interp[level], 2);
    }
    else {
      mg->setLevel(level, tacs[level], NULL);
    }
  }

  // Create the residual and solution vectors on the finest TACS mesh
  TACSBVec *res = tacs[0]->createVec();  res->incref();
  TACSBVec *ans = tacs[0]->createVec();  ans->incref();

  // Allocate the GMRES solution method
  int gmres_iters = 100;
  int nrestart = 1;
  int is_flexible = 0;
  GMRES *gmres = new GMRES(mg->getMat(0), mg,
                           gmres_iters, nrestart, is_flexible);
  gmres->incref();

  // Set a monitor to check on solution progress
  int freq = 1;
  gmres->setMonitor(new KSMPrintStdout("GMRES", mpi_rank, freq));

  // Assemble the Jacobian matrix for each level
  mg->assembleJacobian(1.0, 0.0, 0.0, res);
  tacs[0]->getNodes(res);
  res->scale(1e3);
  res->applyBCs();

  // "Factor" the preconditioner
  mg->factor();

  // Compute the solution using GMRES
  gmres->solve(res, ans);
  
  // Set the variables into TACS
  ans->scale(-1.0);

  // Set the variables on all levels
  mg->setVariables(ans);
  
  for ( int level = 0; level < MAX_NUM_MESH; level++ ){
    // Output for visualization
    unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                               TACSElement::OUTPUT_DISPLACEMENTS);
    TACSToFH5 *f5 = new TACSToFH5(tacs[level], SOLID, write_flag);
    f5->incref();
    char filename[128];
    sprintf(filename, "parallel_output%d.f5", level);
    f5->writeToFile(filename);
    f5->decref();
  }

  mg->decref();

  for ( int level = 0; level < MAX_NUM_MESH; level++ ){
    delete forest[level];
    tacs[level]->decref();
  }

  for ( int level = 0; level < MAX_NUM_MESH-1; level++ ){
    interp[level]->decref();
  }

  TMRFinalize();
  MPI_Finalize();
  return (0);
}
