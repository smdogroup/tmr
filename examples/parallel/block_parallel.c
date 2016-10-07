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

  // The material properties
  TacsScalar rho = 2550.0, E = 70e9, nu = 0.3;

  // The "super-node" locations
  int order = 2;

  // The BDF file
  const char *bdf_file = "uCRM_3D_box_mesh.bdf";

  // Set the initial prefix for the output files
  char prefix[256];
  sprintf(prefix, "results/");

  // Take command-line parameters as input
  double omega = 1.0;
  double mass_fraction = 0.05;
  int scale_objective = 0;
  int use_inverse_vars = 0;
  int use_linear_method = 0;
  int max_num_bfgs = 20;
  int tree_depth = 3;
  
  // Set the mesh
  int use_mesh = 0;

  // Read out the command line arguments
  for ( int k = 0; k < argc; k++ ){
    if (strcmp(argv[k], "use_mesh") == 0){
      use_mesh = 1;
    }
    if (sscanf(argv[k], "tree_depth=%d", &tree_depth) == 1){
      if (tree_depth < 0){ tree_depth = 0; }
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

  // Define the mesh information
  int npts = 0;
  int nelems = 0;
  int num_bc_blocks = 0;
  int *bc_blocks = NULL;
  int *refine = NULL;
  double *Xpts = NULL;
  int *elem_node_conn = NULL;
  double *nodal_forces = 0;

  if (use_mesh){
    // Create the TACSMeshLoader class
    TACSMeshLoader *mesh = new TACSMeshLoader(MPI_COMM_SELF);
    mesh->incref();
    mesh->scanBDFFile(bdf_file);
    
    // Extract the connectivity
    const int *conn;
    const double *X;
    mesh->getConnectivity(&npts, &nelems, NULL, 
                          &conn, &X);

    // Copy the element data
    elem_node_conn = new int[ 8*nelems ];
    Xpts = new double[ 3*npts ];
    memcpy(elem_node_conn, conn, 8*nelems*sizeof(int));
    memcpy(Xpts, X, 3*npts*sizeof(double));    

    // Set the forces at each node
    nodal_forces = new double[ 3*npts ];
    memset(nodal_forces, 0, 3*npts*sizeof(double));
    for ( int i = 0; i < npts; i++ ){
      nodal_forces[3*i+2] = 1e3;
    }

    // Set the refinement increasing out the wing
    refine = new int[ nelems ];
    memset(refine, 0, nelems*sizeof(int));
    
    int max_refine = 2 + tree_depth;
    int min_refine = tree_depth;
    double y_max = 30.0;
    for ( int k = 0; k < nelems; k++ ){
      double y_ref = Xpts[3*elem_node_conn[8*k]+1];
      refine[k] = min_refine + 
        (max_refine - min_refine)*(1.0 - (y_ref/y_max));
    }

    // Deallocate the mesh
    mesh->decref();
  }
  else {
    // Allocate the data that is required to define the problem
    int nx = 80, ny = 8, nz = 8;
    nelems = nx*ny*nz;
    npts = (nx+1)*(ny+1)*(nz+1);

    refine = new int[ nelems ];
    elem_node_conn = new int[ 8*nelems ];
    Xpts = new double[ 3*npts ];
    nodal_forces = new double[ 3*npts ];
    
    // The length of the model
    double xlen = (1.0*nx)/ny;
    double ylen = 1.0;
    double zlen = 1.0;

    // Set the boundary conditions
    num_bc_blocks = ny*nz;
    bc_blocks = new int[ num_bc_blocks ];
    for ( int count = 0, iz = 0; iz < nz; iz++ ){
      for ( int iy = 0; iy < ny; iy++, count++ ){
        int block = nx*iy + nx*ny*iz;
        bc_blocks[count] = 6*block;
      }
    }
    
    // Set the octree refinement
    for ( int i = 0; i < nelems; i++ ){
      refine[i] = tree_depth;
    }
    
    // Set the node locations
    TacsScalar *X = Xpts;
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
    int *conn = elem_node_conn;
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

    // Set the nodal forces
    memset(nodal_forces, 0, 3*npts*sizeof(double));
    nodal_forces[3*(nx+1)+1] = 1e4;
    nodal_forces[3*(nx+1)+2] = 1e4;
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
  forest[0]->setConnectivity(npts, elem_node_conn, nelems);
  forest[0]->createTrees(tree_depth);

  // Compute the volume of the super mesh
  double volume = 0.0;

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
      volume += V;
      if (V < 0.0){
        printf("Negative volume in element %d\n", i);
      }
    }
  }

  // Broadcast the volume
  MPI_Bcast(&volume, 1, MPI_DOUBLE, 0, comm);

  // Repartition the octrees
  if (mpi_rank == 0){ printf("[%d] Repartition\n", mpi_rank); }
  forest[0]->repartition();

  for ( int level = 0; level < MAX_NUM_MESH; level++ ){
    // Balance the
    if (mpi_rank == 0){ printf("[%d] Balance\n", mpi_rank); }
    double tbal = MPI_Wtime();
    forest[level]->balance((level == 0));
    tbal = MPI_Wtime() - tbal;

    if (mpi_rank == 0){ printf("[%d] Create nodes\n", mpi_rank); }
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

    for ( int k = 0; k < num_bc_blocks; k++ ){
      // Get the boundary condition information
      int block = bc_blocks[k]/6;
      int face = bc_blocks[k] % 6;

      if (octrees[block]){
        // Get the octant nodes
        TMROctantArray *nodes;
        octrees[block]->getNodes(&nodes);
      
        // Get the array
        int size;
        TMROctant *array;
        nodes->getArray(&array, &size);
        
        // Set the boundary conditions
        for ( int i = 0; i < size; i++ ){
          // Label the boundary conditions on the face
          int node = -1;
          const int32_t hmax = 1 << TMR_MAX_LEVEL;
          if (face < 2 && array[i].x == (face % 2)*hmax){
            node = array[i].tag;            
          }
          else if (face < 4 && array[i].y == (face % 2)*hmax){
            node = array[i].tag;
          }
          else if (array[i].z == (face % 2)*hmax){
            node = array[i].tag;
          }

          // Add the boundary condition if it is in range
          if (array[i].tag >= range[mpi_rank] && 
              array[i].tag < range[mpi_rank+1]){
            tacs[level]->addBCs(1, &node);
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
        SolidStiffness *stiff = 
          new TMRLinearOctStiffness(weights, nweights, mass_fraction,
				    rho, E, nu, qval);

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
  int sor_iters = 1;
  int sor_symm = 1;
  TACSMg *mg = new TACSMg(comm, MAX_NUM_MESH, 
                          omega, sor_iters, sor_symm);
  mg->incref();

  for ( int level = 0; level < MAX_NUM_MESH; level++ ){
    if (level < MAX_NUM_MESH-1){
      mg->setLevel(level, tacs[level], interp[level], 1);
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

  // Get the connectivity and the inverse connectivity to enable easy
  // node->block look up
  const int *block_conn;
  const int *node_block_conn, *node_block_ptr;
  forest[0]->getConnectivity(NULL, NULL, NULL, NULL,
                             &block_conn, NULL, NULL, NULL);
  forest[0]->getInverseConnectivity(&node_block_conn, &node_block_ptr,
                                    NULL, NULL, NULL, NULL);

  // Get the inverse connectivity for each node
  for ( int node = 0; node < npts; node++ ){
    int ip = node_block_ptr[node];
    int owner = node_block_conn[ip];
    ip++;

    // Find the block owner
    for ( ; ip < node_block_ptr[node+1]; ip++ ){
      int block = node_block_conn[ip];
      if (block < owner){ owner = block; }
    }

    if (octrees[owner]){
      // Find the corner that matches
      int corner = 0;
      for ( ; corner < 8; corner++ ){
        if (block_conn[8*owner + corner] == node){
          break;
        }
      }

      // Set the node location within the octant
      const int32_t hmax = 1 << TMR_MAX_LEVEL;
      TMROctant n;
      n.x = hmax*(corner % 2);
      n.y = hmax*((corner % 4)/2);
      n.z = hmax*(corner/4);

      // Get the octant nodes
      TMROctantArray *nodes;
      octrees[owner]->getNodes(&nodes);
      
      // Use the node search
      const int use_node_search = 1;
      TMROctant *t = nodes->contains(&n, use_node_search);

      // Compute the index within the local array
      const int *range;
      forest[0]->getOwnedNodeRange(&range);    
      int index = 3*(t->tag - range[mpi_rank]);

      // Get the local force array and apply the nodal forces
      TacsScalar *f;
      force->getArray(&f);
      f[index] = nodal_forces[3*node];
      f[index+1] = nodal_forces[3*node+1];
      f[index+2] = nodal_forces[3*node+2];
    }
  }
  
  // Set the target mass
  double target_mass = rho*mass_fraction*volume;

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
    ParOptScalar mass_scale = prob->getMassScaling();

    // Scale the objective by the l1 norm of the gradient
    fobj_scale *= 1.0/g->maxabs();
    mass_scale *= 1.0/A->maxabs();
    
    // Set the objective scaling
    prob->setObjectiveScaling(fobj_scale);
    prob->setMassScaling(mass_scale);

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

  // Set the problem up to use the Hessian-vector products
  opt->setUseLineSearch(1);
  
  // Set the Hessian reset frequency
  opt->setBFGSUpdateType(LBFGS::DAMPED_UPDATE);
  
  // Set the barrier parameter information
  opt->setBarrierPower(1.0);
  opt->setBarrierFraction(0.25);

  // The penality parameter
  double q = 0.0;
  double tol = 1e-3;

  int max_iters = 100;
  for ( int k = 0; k < max_iters; k++ ){
    // Print out the results at the current iterate
    if (mpi_rank == 0){ 
      printf("[%d] New iteration %d with q = %g\n", mpi_rank, k, q); 
    }

    // Set the barrier parameter if this is the second or greater
    // time throught
    opt->setAbsOptimalityTol(tol);    
    if (k > 0){
      opt->setInitBarrierParameter(tol);
    }

    // Scale the tolerance
    tol *= 0.25;
    if (tol < 1e-6){
      tol = 1e-6;
    }

    // Set the log/output file
    char outfile[256];
    sprintf(outfile, "%s//paropt_output%d.out", prefix, k);
    opt->setOutputFile(outfile);
  
    // Set the history/restart file
    char restartfile[256];
    sprintf(restartfile, "%s//paropt_restart%d.bin", prefix, k);
    opt->optimize(restartfile);

    // Set the new value of the penalization using an incremental approach.
    // Note that this will destroy the value of thelinearization
    q = 1.0 + 1.0*(k/2);
    if (q > 5.0){
      q = 5.0;
    }

    // Get the optimized design point and set the new linearized point
    ParOptVec *x;
    opt->getOptimizedPoint(&x, NULL, NULL, NULL, NULL);
    prob->setLinearization(q, x);

    // Do not re-initialize the starting point: Use the design variable
    // and dual variable estimates.
    opt->setInitStartingPoint(0);

    // Create and write out an fh5 file
    unsigned int write_flag = (TACSElement::OUTPUT_NODES |
			       TACSElement::OUTPUT_DISPLACEMENTS |
			       TACSElement::OUTPUT_EXTRAS);
    TACSToFH5 *f5 = new TACSToFH5(tacs[0], SOLID, write_flag);
    f5->incref();
    
    // Write out the solution
    sprintf(outfile, "%s//tacs_output%d.f5", prefix, k);
    f5->writeToFile(outfile);
    f5->decref();
  }

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
