#include "TMRQuadForest.h"
#include "TACSAssembler.h"
#include "isoFSDTStiffness.h"
#include "PlaneStressQuad.h"
#include "TACSMg.h"
#include "TACSToFH5.h"

/*
  Create the following L-bracket connectivity

  0----1----2
  |    |    |
  |    |    |
  3----4----5----9----10
  |    |    |    |     |
  |    |    |    |     |
  6----7----8----11---12
*/

const int nfaces = 6;
const int npts = 13;
const int conn[] = 
  {3, 4, 0, 1,
   4, 5, 1, 2,
   6, 7, 3, 4,
   7, 8, 4, 5,
   8, 11, 5, 9,
   9, 11, 10, 12};

const double Xpts[] =
  {0.0, 2.0, 0.0,
   1.0, 2.0, 0.0,
   2.0, 2.0, 0.0,
   0.0, 1.0, 0.0,
   1.0, 1.0, 0.0,
   2.0, 1.0, 0.0,
   0.0, 0.0, 0.0,
   1.0, 0.0, 0.0,
   2.0, 0.0, 0.0,
   3.0, 1.0, 0.0,
   4.0, 1.0, 0.0,
   3.0, 0.0, 0.0,
   4.0, 0.0, 0.0};

/*
  Interpoalte from the connectivity/node locations
*/
void getLocation( const int *cn,
                  const double *x,
                  const TMRQuadrant *quad, TacsScalar X[] ){
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  double u = 0.0, v = 0.0;
  if (quad->x == hmax-1){ u = 1.0; }
  else { u = 1.0*quad->x/hmax; }
  if (quad->y == hmax-1){ v = 1.0; }
  else { v = 1.0*quad->y/hmax; }

  double N[4];
  N[0] = (1.0 - u)*(1.0 - v);
  N[1] = u*(1.0 - v);
  N[2] = (1.0 - u)*v;
  N[3] = u*v;

  X[0] = X[1] = X[2] = 0.0;
  for ( int k = 0; k < 4; k++ ){
    int node = cn[4*quad->face + k];
    X[0] += x[3*node]*N[k];
    X[1] += x[3*node+1]*N[k];
    X[2] += x[3*node+2]*N[k];
  }
}

void writeSerialMesh( TMRQuadForest *forest ){
  // Print out the full mesh
  TMRQuadrantArray *nodes;
  forest->getNodes(&nodes);

  // Get the quadrants associated with the nodes
  int size;
  TMRQuadrant *array;
  nodes->getArray(&array, &size);
  
  // Set the node locations
  double *X = new double[ 3*size ];
  memset(X, 0, 3*size*sizeof(double));
  for ( int i = 0; i < size; i++ ){
    int node = array[i].tag;
    if (node >= 0){
      getLocation(conn, Xpts, &array[i], &X[3*node]);
    }
    else {
      getLocation(conn, Xpts, &array[i], &X[3*(size+node)]);
    }
  }

  // Get the mesh connectivity
  int *elem_conn, num_elements = 0;
  forest->createMeshConn(&elem_conn, &num_elements);

  // Write out the vtk file
  FILE *fp = fopen("mesh.vtk", "w");
  fprintf(fp, "# vtk DataFile Version 3.0\n");
  fprintf(fp, "vtk output\nASCII\n");
  fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
    
  // Write out the points
  fprintf(fp, "POINTS %d float\n", size);
  for ( int k = 0; k < size; k++ ){
    fprintf(fp, "%e %e %e\n", X[3*k], X[3*k+1], X[3*k+2]);
  }
  delete [] X;
  
  // Write out the cell values
  fprintf(fp, "\nCELLS %d %d\n", 
          num_elements, 5*num_elements);
  for ( int k = 0; k < num_elements; k++ ){
    fprintf(fp, "4 ");
    const int order[] = {0, 1, 3, 2};
    for ( int j = 0; j < 4; j++ ){
      int node = elem_conn[4*k+order[j]];
      if (node < 0){
        node = size+node;
      }
      fprintf(fp, "%d ", node);
    }
    fprintf(fp, "\n");
  }
  delete [] elem_conn;

  // All quadrilaterals
  fprintf(fp, "\nCELL_TYPES %d\n", num_elements);
  for ( int k = 0; k < num_elements; k++ ){
    fprintf(fp, "%d\n", 9);
  }

  // Print out the rest as fields one-by-one
  fprintf(fp, "POINT_DATA %d\n", size);

  fprintf(fp, "SCALARS var_nums float 1\n");
  fprintf(fp, "LOOKUP_TABLE default\n");
    
  for ( int k = 0; k < size; k++ ){
    fprintf(fp, "%d\n", k);
  }
  fclose(fp);
}

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  // Define the different forest levels
  MPI_Comm comm = MPI_COMM_WORLD;
  
  // Get the MPI rank
  int mpi_size, mpi_rank;
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);

  // Create the forests
  TMRQuadForest *forest[3];
  forest[0] = new TMRQuadForest(comm);
  
  forest[0]->setConnectivity(npts, conn, nfaces);
  forest[0]->createRandomTrees(250, 5, 15);
  // forest[0]->createTrees(6);
  forest[0]->repartition();

  for ( int level = 0; level < 3; level++ ){
    double tbal = MPI_Wtime();
    forest[level]->balance(0);
    tbal = MPI_Wtime() - tbal;
    printf("[%d] Balance: %f\n", mpi_rank, tbal);
    forest[level]->repartition();

    // Create the nodes
    double tnodes = MPI_Wtime();
    forest[level]->createNodes(2);
    tnodes = MPI_Wtime() - tnodes;
    printf("[%d] Nodes: %f\n", mpi_rank, tnodes);
  
    // Create the coarse mesh
    if (level < 2){
      forest[level+1] = forest[level]->coarsen();
    }
  }

  if (mpi_size == 1){
    writeSerialMesh(forest[0]);
  }

  // Allocate the stiffness object
  TacsScalar rho = 2570.0, E = 70e9, nu = 0.3;
  PlaneStressStiffness *stiff = 
    new PlaneStressStiffness(rho, E, nu);

  // Allocate the solid element class
  TACSElement *elem = new PlaneStressQuad<2>(stiff, LINEAR, mpi_rank);
  
  // Create the TACSAssembler objects
  TACSAssembler *tacs[3];

  for ( int level = 0; level < 3; level++ ){
    // Find the number of nodes for this processor
    const int *range;
    forest[level]->getOwnedNodeRange(&range);
    int num_nodes = range[mpi_rank+1] - range[mpi_rank];

    // Create the mesh
    double tmesh = MPI_Wtime();
    int *elem_conn, num_elements = 0;
    forest[level]->createMeshConn(&elem_conn, &num_elements);
    tmesh = MPI_Wtime() - tmesh;
    printf("[%d] Mesh: %f\n", mpi_rank, tmesh);

    // Get the dependent node information
    const int *dep_ptr, *dep_conn;
    const double *dep_weights;

    // Create/retrieve the dependent node information
    double tdep = MPI_Wtime();
    int num_dep_nodes = 
      forest[level]->getDepNodeConn(&dep_ptr, &dep_conn,
                             &dep_weights);
    tdep = MPI_Wtime() - tdep;
    printf("[%d] Dependent nodes: %f\n", mpi_rank, tdep);

    // Create the associated TACSAssembler object
    int vars_per_node = 2;
    tacs[level] = new TACSAssembler(comm, vars_per_node,
                                    num_nodes, num_elements,
                                    num_dep_nodes);
    tacs[level]->incref();

    // Set the element ptr
    int *ptr = new int[ num_elements+1 ];
    for ( int i = 0; i < num_elements+1; i++ ){
      ptr[i] = 4*i;
    }
    
    // Set the element connectivity into TACSAssembler
    tacs[level]->setElementConnectivity(elem_conn, ptr);
    delete [] elem_conn;
    delete [] ptr;
    
    // Set the dependent node information
    tacs[level]->setDependentNodes(dep_ptr, dep_conn, dep_weights);

    // Set the elements
    TACSElement **elems = new TACSElement*[ num_elements ];
    for ( int k = 0; k < num_elements; k++ ){
      elems[k] = elem;
    }
    
    // Set the element array
    tacs[level]->setElements(elems);
    delete [] elems;

    // Get the nodes
    TMRQuadrantArray *nodes;
    forest[level]->getNodes(&nodes);

    // Get the quadrants associated with the nodes
    int size;
    TMRQuadrant *array;
    nodes->getArray(&array, &size);

    // Loop over all the nodes
    for ( int i = 0; i < size; i++ ){
      // Evaluate the point
      TacsScalar Xpoint[3];
      getLocation(conn, Xpts, &array[i], Xpoint);

      if (Xpoint[0] < 1e-3){
        int node = array[i].tag;
        if (node >= 0){
          tacs[level]->addBCs(1, &node);
        }
      }
    }

    // Initialize the TACSAssembler object
    tacs[level]->initialize();

    // Create the node vector
    TacsScalar *Xn;
    TACSBVec *X = tacs[level]->createNodeVec();
    X->getArray(&Xn);

    // Get the points
    TMRPoint *Xp;
    forest[level]->getPoints(&Xp);

    // Loop over all the nodes
    for ( int i = 0; i < size; i++ ){
      // Evaluate the point
      TacsScalar Xpoint[3];
      getLocation(conn, Xpts, &array[i], Xpoint);

      // Set the point
      Xp[i].x = Xpoint[0];
      Xp[i].y = Xpoint[1];
      Xp[i].z = Xpoint[2];

      if (array[i].tag >= range[mpi_rank] &&
          array[i].tag < range[mpi_rank+1]){
        int loc = array[i].tag - range[mpi_rank];
        Xn[3*loc] = Xpoint[0];
        Xn[3*loc+1] = Xpoint[1];
        Xn[3*loc+2] = Xpoint[2];
      }
    }
    
    // Set the node locations into TACSAssembler
    tacs[level]->setNodes(X);
  }

  // Create the interpolation
  TACSBVecInterp *interp[2];

  for ( int level = 0; level < 2; level++ ){
    // Create the interpolation object
    interp[level] = new TACSBVecInterp(tacs[level+1]->getVarMap(),
                                       tacs[level]->getVarMap(),
                                       tacs[level]->getVarsPerNode());
    interp[level]->incref();

    // Set the interpolation
    forest[level]->createInterpolation(forest[level+1], interp[level]);
    
    // Initialize the interpolation
    interp[level]->initialize();
  }

  // Create the multigrid object
  double omega = 1.0;
  int mg_sor_iters = 1;
  int mg_sor_symm = 1;
  int mg_iters_per_level = 1;
  TACSMg *mg = new TACSMg(comm, 3, omega, mg_sor_iters, mg_sor_symm);
  mg->incref();

  for ( int level = 0; level < 3; level++ ){
    if (level < 2){
      mg->setLevel(level, tacs[level], interp[level], mg_iters_per_level);
    }
    else {
      mg->setLevel(level, tacs[level], NULL);
    }
  }

  // Assemble the matrix
  mg->assembleJacobian(1.0, 0.0, 0.0, NULL);
  mg->factor();

  // Create a force vector
  TACSBVec *force = tacs[0]->createVec();
  force->incref();
  force->set(1.0);
  tacs[0]->applyBCs(force);

  TACSBVec *ans = tacs[0]->createVec();
  ans->incref();
  
  // Set up the solver
  int gmres_iters = 100; 
  int nrestart = 2;
  int is_flexible = 1;
  GMRES *gmres = new GMRES(mg->getMat(0), mg, 
                           gmres_iters, nrestart, is_flexible);
  gmres->incref();
  gmres->setMonitor(new KSMPrintStdout("GMRES", mpi_rank, 10));
  gmres->setTolerances(1e-10, 1e-30);

  gmres->solve(force, ans);
  ans->scale(-1.0);
  tacs[0]->setVariables(ans);

  // Create and write out an fh5 file
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS |
                             TACSElement::OUTPUT_STRESSES |
                             TACSElement::OUTPUT_EXTRAS);
  TACSToFH5 *f5 = new TACSToFH5(tacs[0], PLANE_STRESS, write_flag);
  f5->incref();
    
  // Write out the solution
  f5->writeToFile("output.f5");
  f5->decref();

  // Create the level
  for ( int level = 0; level < 3; level++ ){
    tacs[level]->decref();
    delete forest[level];
  }

  TMRFinalize();
  MPI_Finalize();
  return (0);
}
