#include "TMROctForest.h"
#include "TACSMeshLoader.h"
#include "TACSAssembler.h"
#include "Solid.h"
#include "TACSMg.h"
#include "TMR_STLTools.h"

/*
  Bottom
  3---4---5
  |   |   |
  0---1---2

  Top
  9--10--11
  |   |   |
  6---7---8
*/
const int rectangle_npts = 12;
const int rectangle_nelems = 2;

const double rectangle_xpts[] =
  {0.0, 0.0, 0.0,
   1.0, 0.0, 0.0,
   2.0, 0.0, 0.0,
   0.0, 1.0, 0.0,
   1.0, 1.0, 0.0,
   2.0, 1.0, 0.0,
   0.0, 0.0, 1.0,
   1.0, 0.0, 1.0,
   2.0, 0.0, 1.0,
   0.0, 1.0, 1.0,
   1.0, 1.0, 1.0,
   2.0, 1.0, 1.0};

const int rectangle_conn[] =
  {0, 1, 3, 4, 6, 7, 9, 10,
   8, 11, 2, 5, 7, 10, 1, 4};

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
void getLocation( const int *elem_node_conn, 
                  const double *Xpts,
                  const TMROctant *oct,
                  const int order, int index,
                  const double knots[],
                  TacsScalar X[] ){
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - oct->level);
  int ii = index % order;
  int jj = (index % (order*order))/order;
  int kk = index/(order*order);
  double u = (oct->x + h*knots[ii])/hmax;
  double v = (oct->y + h*knots[jj])/hmax;
  double w = (oct->z + h*knots[kk])/hmax;

  double N[8];
  N[0] = (1.0 - u)*(1.0 - v)*(1.0 - w);
  N[1] = u*(1.0 - v)*(1.0 - w);
  N[2] = (1.0 - u)*v*(1.0 - w);
  N[3] = u*v*(1.0 - w);
  N[4] = (1.0 - u)*(1.0 - v)*w;
  N[5] = u*(1.0 - v)*w;
  N[6] = (1.0 - u)*v*w;
  N[7] = u*v*w;

  X[0] = X[1] = X[2] = 0.0;
  for ( int k = 0; k < 8; k++ ){
    int node = elem_node_conn[8*oct->block + k];
    X[0] += Xpts[3*node]*N[k];
    X[1] += Xpts[3*node+1]*N[k];
    X[2] += Xpts[3*node+2]*N[k];
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
  Check volume of the mesh to see if it is valid
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

  const int NUM_LEVELS = 4;

  // Define the different forest levels
  MPI_Comm comm = MPI_COMM_WORLD;
  
  // Get the MPI rank
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  // The "super-node" locations
  int npts = box_npts;
  int nelems = box_nelems;
  const double *Xpts = box_xpts;
  const int *conn = box_conn;

  for ( int k = 0; k < argc; k++ ){
    if (strcmp(argv[k], "connector") == 0){
      npts = connector_npts;
      nelems = connector_nelems;
      Xpts = connector_xpts;
      conn = connector_conn;
    }
    else if (strcmp(argv[k], "rectangle") == 0){
      npts = rectangle_npts;
      nelems = rectangle_nelems;
      Xpts = rectangle_xpts;
      conn = rectangle_conn;
    }
  }

  double knots4[4] = {0.0, 0.25, 0.75, 1.0};
  double knots3[3] = {0.0, 0.5, 1.0};
  double knots2[2] = {0.0, 1.0};
  double *knots[NUM_LEVELS];
  knots[0] = knots4;
  knots[1] = knots3;
  knots[2] = knots2;
  knots[3] = knots2;
  
  // Create the forests
  int order = 4;
  TMROctForest *forest[NUM_LEVELS];
  forest[0] = new TMROctForest(comm, order, TMR_GAUSS_LOBATTO_POINTS);
  forest[0]->incref();
  forest[0]->setConnectivity(npts, conn, nelems);
  forest[0]->createRandomTrees(15, 0, 10);  
  forest[0]->repartition();

  for ( int level = 0; level < NUM_LEVELS; level++ ){
    double tbal = MPI_Wtime();
    forest[level]->balance(0);
    tbal = MPI_Wtime() - tbal;
    printf("[%d] Balance: %f\n", mpi_rank, tbal);
    forest[level]->repartition();

    // Create the nodes
    double tnodes = MPI_Wtime();
    forest[level]->createNodes();
    tnodes = MPI_Wtime() - tnodes;
    printf("[%d] Nodes: %f\n", mpi_rank, tnodes);
  
    // Create the coarse mesh
    if (level < NUM_LEVELS-1){
      if (order > 2){
        forest[level+1] = forest[level]->duplicate();
        order = order-1;
        forest[level+1]->setMeshOrder(order, TMR_GAUSS_LOBATTO_POINTS);
      }
      else {
        forest[level+1] = forest[level]->coarsen();
      }
      forest[level+1]->incref();
    }
  }

  // Allocate the stiffness object
  TacsScalar rho = 2570.0, E = 70e9, nu = 0.3, ys = 1.0;
  SolidStiffness *stiff = new SolidStiffness(rho, E, nu, ys);

  // Allocate the solid element class
  TACSElement *solid[NUM_LEVELS];
  solid[0] = new Solid<4>(stiff, LINEAR, mpi_rank);
  solid[1] = new Solid<3>(stiff, LINEAR, mpi_rank);
  solid[2] = new Solid<2>(stiff, LINEAR, mpi_rank);
  solid[3] = new Solid<2>(stiff, LINEAR, mpi_rank);
  
  // Create the TACSAssembler objects
  TACSAssembler *tacs[NUM_LEVELS];

  for ( int level = 0; level < NUM_LEVELS; level++ ){
    // Find the number of nodes for this processor
    const int *range;
    forest[level]->getOwnedNodeRange(&range);
    int num_nodes = range[mpi_rank+1] - range[mpi_rank];

    // Create the mesh
    const int *elem_conn;
    int num_elements = 0;
    forest[level]->getNodeConn(&elem_conn, &num_elements);

    // Get the dependent node information
    const int *dep_ptr, *dep_conn;
    const double *dep_weights;
    int num_dep_nodes = 
      forest[level]->getDepNodeConn(&dep_ptr, &dep_conn,
                                    &dep_weights);

    // Create the associated TACSAssembler object
    int vars_per_node = 3;
    tacs[level] = new TACSAssembler(comm, vars_per_node,
                                    num_nodes, num_elements,
                                    num_dep_nodes);
    tacs[level]->incref();

    // Set the element ptr
    int order = forest[level]->getMeshOrder();
    int *ptr = new int[ order*order*order*num_elements ];
    for ( int i = 0; i < num_elements+1; i++ ){
      ptr[i] = order*order*order*i;
    }
    
    // Set the element connectivity into TACSAssembler
    tacs[level]->setElementConnectivity(elem_conn, ptr);
    delete [] ptr;
    
    // Set the dependent node information
    tacs[level]->setDependentNodes(dep_ptr, dep_conn, dep_weights);

    // Set the elements
    TACSElement **elems = new TACSElement*[ num_elements ];
    for ( int k = 0; k < num_elements; k++ ){
      elems[k] = solid[level];
      elems[k]->incref();
    }
    
    // Set the element array
    tacs[level]->setElements(elems);
    delete [] elems;

    // Initialize the TACSAssembler object
    tacs[level]->initialize();

    // Get the octant locations
    TMROctantArray *octants;
    forest[level]->getOctants(&octants);
    
    // Get the octants associated with the nodes
    int oct_size;
    TMROctant *octs;
    octants->getArray(&octs, &oct_size);

    // Create the node vector
    TacsScalar *Xn;
    TACSBVec *X = tacs[level]->createNodeVec();
    X->getArray(&Xn);

    // Get the points
    TMRPoint *Xp;
    forest[level]->getPoints(&Xp);

    // Loop over all the nodes
    for ( int i = 0; i < oct_size; i++ ){
      const int *c = &elem_conn[order*order*order*i];
      for ( int j = 0; j < order*order*order; j++ ){
        if (c[j] >= range[mpi_rank] &&
            c[j] < range[mpi_rank+1]){
          int index = c[j] - range[mpi_rank];

          // Evaluate the point
          getLocation(conn, Xpts, &octs[i], order, j, 
                      knots[level], &Xn[3*index]);
        }
      }
    }
    
    // Set the node locations into TACSAssembler
    tacs[level]->setNodes(X);
  }

  // Create the interpolation
  TACSBVecInterp *interp[NUM_LEVELS-1];

  for ( int level = 0; level < NUM_LEVELS-1; level++ ){
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

  // Create a vector on the finest level
  TACSBVec *x[NUM_LEVELS];
  x[0] = tacs[0]->createVec();
  x[0]->incref();
  tacs[0]->getNodes(x[0]);
  for ( int level = 0; level < NUM_LEVELS-1; level++ ){
    x[level+1] = tacs[level+1]->createVec();
    x[level+1]->incref();
    interp[level]->multWeightTranspose(x[level], x[level+1]);
  }

  // Create and write out an fh5 file
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS |
                             TACSElement::OUTPUT_EXTRAS);

  for ( int level = 0; level < NUM_LEVELS; level++ ){
    tacs[level]->setVariables(x[level]);
    TACSToFH5 *f5 = new TACSToFH5(tacs[level], TACS_SOLID, write_flag);
    f5->incref();
    
    // Write out the solution
    char filename[128];
    sprintf(filename, "output%d.f5", level);
    f5->writeToFile(filename);
    f5->decref();
  }

  // Create the level
  for ( int level = 0; level < NUM_LEVELS; level++ ){
    x[level]->decref();
    tacs[level]->decref();
    forest[level]->decref();
  }

  TMRFinalize();
  MPI_Finalize();
  return (0);
}
