#include "TMROctForest.h"
#include "TMROctStiffness.h"
#include "TACSMeshLoader.h"
#include "TACSAssembler.h"
#include "TACSMg.h"
#include "BVecInterp.h"
#include "Solid.h"
#include "TMR_STLTools.h"
#include "TMRTopoProblem.h"
#include "ParOpt.h"

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

/*
  Compute the shape functions and their derivatives
*/
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
  
  // Set the communicator
  MPI_Comm comm = MPI_COMM_WORLD;

  // Get the MPI rank
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // Repartition the mesh immediately
  int partition = 0;

  // The material properties
  TacsScalar rho = 2550.0, E = 70e9, nu = 0.3;

  // The "super-node" locations
  double omega = 1.0;
  int order = 2;
  int npts = 0;
  int nelems = 0;
  const double *Xpts = NULL;
  const int *elem_node_conn = NULL;

  // Set up the model dimensions
  int nx = 80, ny = 8, nz = 8;
  nelems = nx*ny*nz;
  npts = (nx+1)*(ny+1)*(nz+1);
  TacsScalar *X = new TacsScalar[ 3*npts ];
  int *conn = new int[ 8*nelems ];
  
  // The length of the model
  double xlen = (1.0*nx)/ny;
  double ylen = 1.0;
  double zlen = 1.0;

  // Set the arrays
  Xpts = X;
  elem_node_conn = conn;
  
  for ( int iz = 0; iz < nz+1; iz++ ){
    for ( int iy = 0; iy < ny+1; iy++ ){
      for ( int ix = 0; ix < nx+1; ix++ ){
        X[0] = xlen*ix/nx;  X++;
        X[0] = ylen*iy/ny;  X++;
        X[0] = zlen*iz/nz;  X++;
      }
    }
  }
  
  // Loop over all the elements in the mesh
  for ( int iz = 0; iz < nz; iz++ ){
    for ( int iy = 0; iy < ny; iy++ ){
      for ( int ix = 0; ix < nx; ix++ ){
        // Set the element-level connectivity
        for ( int iiz = 0; iiz < 2; iiz++ ){
          for ( int iiy = 0; iiy < 2; iiy++ ){
            for ( int iix = 0; iix < 2; iix++ ){
              conn[0] = 
                ((ix+iix) + (iy+iiy)*(nx+1) + (iz+iiz)*(nx+1)*(ny+1));
              conn++;
            }
          }
        }
      }
    }
  }

  // Set the initial prefix for the output files
  char prefix[256];
  sprintf(prefix, "results/");

  // Take command-line parameters as input
  double mass_fraction = 0.05;
  int scale_objective = 0;
  int use_inverse_vars = 0;
  int use_linear_method = 0;
  int max_num_bfgs = 20;
  int tree_depth = 3;
  for ( int k = 0; k < argc; k++ ){
    if (strcmp(argv[k], "partition") == 0){
      partition = 1;
    }
    if (sscanf(argv[k], "tree_depth=%d", &tree_depth) == 1){
      if (tree_depth < 0){ tree_depth = 1; }
    }
    if (sscanf(argv[k], "order=%d", &order) == 1){
      if (order < 2){ order = 2; }
      if (order > 3){ order = 3; }
    }
    if (sscanf(argv[k], "omega=%lf", &omega) == 1){
      if (omega < 0.05){ omega = 0.05; }
    }
    if (sscanf(argv[k], "prefix=%s", prefix) == 1){
      if (mpi_rank == 0){
        printf("Using prefix = %s\n", prefix);
      }
    }
    if (sscanf(argv[k], "max_num_bfgs=%d", &max_num_bfgs) == 1){
      if (max_num_bfgs < 0){
	max_num_bfgs = 1;
      }
    }
    if (strcmp(argv[k], "use_inverse_vars") == 0){
      use_inverse_vars = 1;
    }
    if (strcmp(argv[k], "use_linear_method") == 0){
      use_linear_method = 1;
    }
    if (strcmp(argv[k], "scale_objective") == 0){
      scale_objective = 1;
    }
    if (sscanf(argv[k], "mass_fraction=%lf", &mass_fraction) == 1){
      if (mass_fraction < 0.05){
	mass_fraction = 0.05;
      }
      else if (mass_fraction > 0.5){
	mass_fraction = 0.5;
      }
    }
  }

  // Define the different forest levels
  const int MAX_NUM_MESH = 4;
  TMROctForest *forest[MAX_NUM_MESH];
  TMROctForest *filter[MAX_NUM_MESH];

  // The TACSAssembler models at each level within the mesh
  TACSAssembler *tacs[MAX_NUM_MESH];

  // The interpolation/restriction operator between solution levels
  TACSBVecInterp *interp[MAX_NUM_MESH-1];

  // Set up the variable map for the design variable numbers
  TACSVarMap *filter_maps[MAX_NUM_MESH];
  TACSBVecIndices *filter_indices[MAX_NUM_MESH];

  // Create the forests
  forest[0] = new TMROctForest(comm);

  if (mpi_rank == 0){
    printf("order = %d\n", order);
    printf("omega = %g\n", omega);
    printf("prefix = %s\n", prefix);
    printf("mass_fraction = %f\n", mass_fraction);
    printf("use_inverse_vars = %d\n", use_inverse_vars);
    printf("use_linear_method = %d\n", use_linear_method);
    printf("scale_objective = %d\n", scale_objective);
    printf("tree_depth = %d\n", tree_depth);
    printf("max_num_bfgs = %d\n", max_num_bfgs);
  }

  if (Xpts && elem_node_conn){
    forest[0]->setConnectivity(npts, elem_node_conn,
                               nelems, partition);
    if (order == 3){
      forest[0]->createTrees(tree_depth);
    }
    else {
      forest[0]->createTrees(tree_depth);
    }
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

    // Print out the number of blocks, faces, edges and nodes
    printf("nblocks = %d\nnfaces = %d\nnedges = %d\nnnodes = %d\n", 
           nblocks, nfaces, nedges, nnodes);
  }

  // Repartition the octrees
  printf("[%d] Repartition\n", mpi_rank);
  forest[0]->repartition();

  for ( int level = 0; level < MAX_NUM_MESH; level++ ){
    // Balance the
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

    // Set the node locations
    const int *owned;
    int nowned = forest[level]->getOwnedOctrees(&owned);
    for ( int k = 0; k < nowned; k++ ){
      int block = owned[k];
    
      // Get the octant nodes
      TMROctantArray *nodes;
      octrees[block]->getNodes(&nodes);
      
      // Get the nodal array
      int size;
      TMROctant *array;
      nodes->getArray(&array, &size);

      // Get the nodes
      TMRPoint *X;
      octrees[block]->getPoints(&X);

      // Loop over all the nodes
      for ( int i = 0; i < size; i++ ){
        TacsScalar x[3];
        getLocation(block, elem_node_conn, Xpts, &array[i], x);
        X[i].x = x[0];
        X[i].y = x[1];
        X[i].z = x[2];
      }
    }

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

    for ( int iz = 0; iz < nz; iz++ ){
      for ( int iy = 0; iy < ny; iy++ ){
        int block = nx*iy + nx*ny*iz;

        // Add nodes associated with the boundary conditions
        if (octrees[block]){
          // Get the octant nodes
          TMROctantArray *nodes;
          octrees[block]->getNodes(&nodes);
      
          // Get the array
          int size;
          TMROctant *array;
          nodes->getArray(&array, &size);
          
          // Set the boundary conditions
          int nbc = 0;
          for ( int i = 0; i < size; i++ ){
            if (array[i].x == 0 && 
                (array[i].tag >= range[mpi_rank] && 
                 array[i].tag < range[mpi_rank+1])){
              nbc++;
              int node = array[i].tag;
              tacs[level]->addBCs(1, &node);
            }
          }
        }
      }
    }

    // Allocate the element array
    TACSElement **elems = new TACSElement*[ num_elements ];

    // Set up the filter for the local ordering
    if (level == 0){
      filter[level] = forest[level]->coarsen();
    }
    else {
      filter[level] = filter[level-1]->coarsen();
    }
    filter[level]->balance();
    filter[level]->createNodes(2);
    filter[level]->createMeshConn(NULL, NULL);

    // Get the node range for the filter
    const int *filter_range;
    filter[level]->getOwnedNodeRange(&filter_range);

    // Get a sorted list of the external node numbers
    int *ext_node_nums;
    int num_ext_nodes = filter[level]->getExtNodeNums(&ext_node_nums);

    // Set up the variable map for the design variable numbers
    int num_filter_local = filter_range[mpi_rank+1] - filter_range[mpi_rank];
    filter_maps[level] = new TACSVarMap(comm, num_filter_local);

    // Set the external filter indices
    filter_indices[level] = new TACSBVecIndices(&ext_node_nums, num_ext_nodes);
    ext_node_nums = NULL;
    filter_indices[level]->setUpInverse();

    // Get the filter octrees
    TMROctree **filter_octrees;
    filter[level]->getOctrees(&filter_octrees);
    filter[level]->getDepNodeConn(&dep_ptr, &dep_conn,
                                  &dep_weights);

    // Set the node locations in the filter
    const int *filter_owned;
    int filter_nowned = filter[level]->getOwnedOctrees(&filter_owned);
    for ( int k = 0; k < filter_nowned; k++ ){
      int block = filter_owned[k];
    
      // Get the octant nodes
      TMROctantArray *nodes;
      filter_octrees[block]->getNodes(&nodes);
      
      // Get the nodal array
      int size;
      TMROctant *array;
      nodes->getArray(&array, &size);

      // Get the nodes
      TMRPoint *X;
      filter_octrees[block]->getPoints(&X);

      // Loop over all the nodes
      for ( int i = 0; i < size; i++ ){
        TacsScalar x[3];
        getLocation(block, elem_node_conn, Xpts, &array[i], x);
        X[i].x = x[0];
        X[i].y = x[1];
        X[i].z = x[2];
      }
    }

    // Set the material properties to use
    double qval = 5.0;

    // Loop over all of the elements within the array
    for ( int k = 0, num = 0; k < nowned; k++ ){
      int block = owned[k];
    
      // Get the filter nodes
      TMROctantArray *filter_nodes;
      filter_octrees[block]->getNodes(&filter_nodes);

      // Get the elements from the mesh
      TMROctantArray *elements;
      octrees[block]->getElements(&elements);
      
      // Get the nodal array
      int size;
      TMROctant *array;
      elements->getArray(&array, &size);

      // Get the nodes
      TMRPoint *X;
      octrees[block]->getPoints(&X);

      // Loop over all the elements
      for ( int i = 0; i < size; i++, num++ ){
        // Get the enclosing element within the filter
        TMROctant *oct = filter_octrees[block]->findEnclosing(&array[i]);

        // Get the side-length of the container
        const int32_t h = 1 << (TMR_MAX_LEVEL - array[i].level);
        const int32_t hoct = 1 << (TMR_MAX_LEVEL - oct->level);
    
        // Get the u/v/w values within the filter octant
        double pt[3];
        pt[0] = -1.0 + 2.0*(array[i].x + 0.5*h - oct->x)/hoct;
        pt[1] = -1.0 + 2.0*(array[i].y + 0.5*h - oct->y)/hoct;
        pt[2] = -1.0 + 2.0*(array[i].z + 0.5*h - oct->z)/hoct;
        
        // Get the Lagrange shape functions
        double N[8];
        FElibrary::triLagrangeSF(N, pt, 2);
  
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
              p.x = oct->x + hoct*ii;
              p.y = oct->y + hoct*jj;
              p.z = oct->z + hoct*kk;

              const int use_node_search = 1;
              TMROctant *t = filter_nodes->contains(&p, use_node_search);

              int node = t->tag;
              if (node >= 0){
                if (node >= filter_range[mpi_rank] && 
                    node < filter_range[mpi_rank+1]){
                  node = node - filter_range[mpi_rank];
                }
                else {
                  node = num_filter_local + 
                    filter_indices[level]->findIndex(node);
                }
                weights[nweights].index = node;
                weights[nweights].weight = wval;
                nweights++;
              }
              else {
                node = -node-1;
                for ( int jp = dep_ptr[node]; jp < dep_ptr[node+1]; jp++ ){
                  int dep_node = dep_conn[jp];
                  if (dep_node >= filter_range[mpi_rank] && 
                      dep_node < filter_range[mpi_rank+1]){
                    dep_node = dep_node - filter_range[mpi_rank];
                  }
                  else { 
                    dep_node = num_filter_local + 
                      filter_indices[level]->findIndex(dep_node);
                  }
                  weights[nweights].index = dep_node;
                  weights[nweights].weight = wval*dep_weights[jp];
                  nweights++;
                }
              }
            }
          }
        }

        // Sort and sum the array of weights
        nweights = TMRIndexWeight::uniqueSort(weights, nweights);
        
        // Allocate the stiffness object
        /*
          SolidStiffness *stiff = new TMROctStiffness(weights, nweights,
          rho, E, nu, qval);
        */
        SolidStiffness *stiff = 
          new TMRLinearOctStiffness(weights, nweights, rho, E, nu, qval);

        TACSElement *solid = NULL;
        if (order == 2){
          solid = new Solid<2>(stiff);
        }
        else if (order == 3){
          solid = new Solid<3>(stiff);
        }

        // Set the element
        elems[num] = solid;
      }      
    }

    // Set the element array
    tacs[level]->setElements(elems);
    delete [] elems;

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

    TacsScalar *Xpts = NULL;
    Xvec->getArray(&Xpts);

    // Scan through the nodes
    for ( int k = 0; k < nowned; k++ ){
      int block = owned[k];
    
      // Get the octant nodes
      TMROctantArray *nodes;
      octrees[block]->getNodes(&nodes);
      
      // Get the array
      int size;
      TMROctant *array;
      nodes->getArray(&array, &size);

      // Get the nodes
      TMRPoint *X;
      octrees[block]->getPoints(&X);

      // Loop over all the nodes
      for ( int i = 0; i < size; i++ ){
        if (array[i].tag >= 0 && 
            (array[i].tag >= range[mpi_rank] && 
             array[i].tag < range[mpi_rank+1])){
          int index = array[i].tag - range[mpi_rank];
          Xpts[3*index] = X[i].x;
          Xpts[3*index+1] = X[i].y;
          Xpts[3*index+2] = X[i].z;
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
  int sor_iters = 2;
  int sor_symm = 1;
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

  // Create a force vector
  TACSBVec *force = tacs[0]->createVec();

  // Get the octrees within the forest
  TMROctree **octrees;
  int ntrees = forest[0]->getOctrees(&octrees);
  int block = nx-1;

  // Find the far node associated with the block
  if (octrees[block]){
    // Get the octant nodes
    TMROctantArray *nodes;
    octrees[block]->getNodes(&nodes);
    
    const int *range;
    forest[0]->getOwnedNodeRange(&range);
    
    // Get the array
    int size;
    TMROctant *array;
    nodes->getArray(&array, &size);
          
    // Set the boundary conditions
    int nbc = 0;
    const int hmax = 1 << TMR_MAX_LEVEL;
    for ( int i = 0; i < size; i++ ){
      if (array[i].x == hmax &&
          array[i].y == 0 &&
          array[i].z == 0){
        int index = 3*(array[i].tag - range[mpi_rank]);

        TacsScalar *f;
        force->getArray(&f);
        f[index+1] = 1e4;
        f[index+2] = 1e4;
      }
    }
  }

  // Set the log/output file
  char outfile[256];
  sprintf(outfile, "%s//paropt_output.out", prefix);

  // Set the target mass
  double target_mass = rho*mass_fraction*xlen*ylen*zlen;

  // Create the ParOpt problem class
  TMRTopoProblem *prob = 
    new TMRTopoProblem(MAX_NUM_MESH, tacs, force, filter,
                       filter_maps, filter_indices, mg,
                       target_mass, prefix);
  
  // Set the problem to use reciprocal variables
  if (use_inverse_vars){
    prob->setUseReciprocalVariables();
  }

  if (scale_objective){
    // Allocate space for the design variables
    ParOptVec *x = prob->createDesignVec();
    ParOptVec *g = prob->createDesignVec();
    ParOptVec *A = prob->createDesignVec();
    
    // Evaluate the constrain and the constraint gradient
    ParOptScalar fobj, con;
    prob->getVarsAndBounds(x, NULL, NULL);
    prob->evalObjCon(x, &fobj, &con);
    prob->evalObjConGradient(x, g, &A);  
    
    // Evaluate the average gradient component and reset the scaling
    ParOptScalar fobj_scale = prob->getObjectiveScaling();
    
    // Scale the objective by the l1 norm of the gradient
    fobj_scale *= 1.0/g->maxabs();
    
    // Set the objective scaling
    prob->setObjectiveScaling(fobj_scale);

    // Free the data that is not needed
    delete x;
    delete g;
    delete A;
  }

  // Create the topology optimization object
  ParOpt *opt = new ParOpt(prob, max_num_bfgs);

  // Use a sequential linear method
  if (use_linear_method){
    opt->setSequentialLinearMethod(1);
  }

  // Set the optimization parameters
  opt->setMaxMajorIterations(2000);
  opt->setOutputFrequency(1);
  opt->setAbsOptimalityTol(1e-5);
  opt->setOutputFile(outfile);

  // Set the problem up to use the Hessian-vector products
  opt->setUseLineSearch(0);
  opt->setUseHvecProduct(1);
  opt->setGMRESSubspaceSize(50);
  opt->setNKSwitchTolerance(1.0);
  opt->setEisenstatWalkerParameters(0.5, 0.0);
  opt->setGMRESTolerances(1.0, 1e-30);

  // Set the Hessian reset frequency
  opt->setBFGSUpdateType(LBFGS::DAMPED_UPDATE);

  // Set the barrier parameter information
  opt->setBarrierPower(1.5);
  opt->setBarrierFraction(0.25);

  // Check the gradients
  opt->checkGradients(1e-6);

  // Set the history/restart file
  char restartfile[256];
  sprintf(restartfile, "%s//paropt_restart.bin", prefix);
  opt->optimize(restartfile);

  delete opt;
  delete prob;

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
