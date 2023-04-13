#include "PlaneStressQuad.h"
#include "TACSAssembler.h"
#include "TACSMg.h"
#include "TACSToFH5.h"
#include "TMRQuadForest.h"
#include "isoFSDTStiffness.h"

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
const int conn[] = {3, 4, 0, 1, 4, 5,  1, 2, 6, 7,  3,  4,
                    7, 8, 4, 5, 8, 11, 5, 9, 9, 11, 10, 12};

const double Xpts[] = {0.0, 2.0, 0.0, 1.0, 2.0, 0.0, 2.0, 2.0, 0.0, 0.0,
                       1.0, 0.0, 1.0, 1.0, 0.0, 2.0, 1.0, 0.0, 0.0, 0.0,
                       0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.0, 1.0, 0.0,
                       4.0, 1.0, 0.0, 3.0, 0.0, 0.0, 4.0, 0.0, 0.0};

/*
  Interpoalte from the connectivity/node locations
*/
void getLocation(const int *elem_node_conn, const double *Xpts,
                 const TMRQuadrant *quad, const int order, int index,
                 const double knots[], TacsScalar X[]) {
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  const int32_t h = 1 << (TMR_MAX_LEVEL - quad->level);
  int ii = index % order;
  int jj = index / order;
  double u = (quad->x + h * knots[ii]) / hmax;
  double v = (quad->y + h * knots[jj]) / hmax;

  double N[4];
  N[0] = (1.0 - u) * (1.0 - v);
  N[1] = u * (1.0 - v);
  N[2] = (1.0 - u) * v;
  N[3] = u * v;

  X[0] = X[1] = X[2] = 0.0;
  for (int k = 0; k < 4; k++) {
    int node = elem_node_conn[4 * quad->face + k];
    X[0] += Xpts[3 * node] * N[k];
    X[1] += Xpts[3 * node + 1] * N[k];
    X[2] += Xpts[3 * node + 2] * N[k];
  }
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  TMRInitialize();

  // Set the number of levels to use
  const int NUM_LEVELS = 4;

  // Define the different forest levels
  MPI_Comm comm = MPI_COMM_WORLD;

  // Get the MPI rank
  int mpi_size, mpi_rank;
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);

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
  TMRQuadForest *forest[NUM_LEVELS];
  forest[0] = new TMRQuadForest(comm, order, TMR_GAUSS_LOBATTO_POINTS);
  forest[0]->incref();
  forest[0]->setConnectivity(npts, conn, nfaces);
  forest[0]->createRandomTrees(15, 0, 10);
  forest[0]->repartition();

  for (int level = 0; level < NUM_LEVELS; level++) {
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
    if (level < NUM_LEVELS - 1) {
      if (order > 2) {
        forest[level + 1] = forest[level]->duplicate();
        order = order - 1;
        forest[level + 1]->setMeshOrder(order, TMR_GAUSS_LOBATTO_POINTS);
      } else {
        forest[level + 1] = forest[level]->coarsen();
      }
      forest[level + 1]->incref();
    }
  }

  // Allocate the stiffness object
  TacsScalar rho = 2570.0, E = 70e9, nu = 0.3;
  PlaneStressStiffness *stiff = new PlaneStressStiffness(rho, E, nu);

  // Allocate the solid element class
  TACSElement *elem[NUM_LEVELS];
  elem[0] = new PlaneStressQuad<4>(stiff, LINEAR, mpi_rank);
  elem[1] = new PlaneStressQuad<3>(stiff, LINEAR, mpi_rank);
  elem[2] = new PlaneStressQuad<2>(stiff, LINEAR, mpi_rank);
  elem[3] = new PlaneStressQuad<2>(stiff, LINEAR, mpi_rank);

  // Create the TACSAssembler objects
  TACSAssembler *tacs[NUM_LEVELS];

  for (int level = 0; level < NUM_LEVELS; level++) {
    // Find the number of nodes for this processor
    const int *range;
    forest[level]->getOwnedNodeRange(&range);
    int num_nodes = range[mpi_rank + 1] - range[mpi_rank];

    // Create the mesh
    const int *elem_conn;
    int num_elements = 0;
    forest[level]->getNodeConn(&elem_conn, &num_elements);

    // Get the dependent node information
    const int *dep_ptr, *dep_conn;
    const double *dep_weights;
    int num_dep_nodes =
        forest[level]->getDepNodeConn(&dep_ptr, &dep_conn, &dep_weights);

    // Create the associated TACSAssembler object
    int vars_per_node = 2;
    tacs[level] = new TACSAssembler(comm, vars_per_node, num_nodes,
                                    num_elements, num_dep_nodes);
    tacs[level]->incref();

    // Set the element ptr
    int order = forest[level]->getMeshOrder();
    int *ptr = new int[order * order * num_elements];
    for (int i = 0; i < num_elements + 1; i++) {
      ptr[i] = order * order * i;
    }

    // Set the element connectivity into TACSAssembler
    tacs[level]->setElementConnectivity(elem_conn, ptr);
    delete[] ptr;

    // Set the dependent node information
    tacs[level]->setDependentNodes(dep_ptr, dep_conn, dep_weights);

    // Set the elements
    TACSElement **elems = new TACSElement *[num_elements];
    for (int k = 0; k < num_elements; k++) {
      elems[k] = elem[level];
      elems[k]->incref();
    }

    // Set the element array
    tacs[level]->setElements(elems);
    delete[] elems;

    // Initialize the TACSAssembler object
    tacs[level]->initialize();

    // Get the quad locations
    TMRQuadrantArray *quadrants;
    forest[level]->getQuadrants(&quadrants);

    // Get the quadrants associated with the nodes
    int quad_size;
    TMRQuadrant *quads;
    quadrants->getArray(&quads, &quad_size);

    // Create the node vector
    TacsScalar *Xn;
    TACSBVec *X = tacs[level]->createNodeVec();
    X->getArray(&Xn);

    // Get the points
    TMRPoint *Xp;
    forest[level]->getPoints(&Xp);

    // Loop over all the nodes
    for (int i = 0; i < quad_size; i++) {
      const int *c = &elem_conn[order * order * i];
      for (int j = 0; j < order * order; j++) {
        if (c[j] >= range[mpi_rank] && c[j] < range[mpi_rank + 1]) {
          int index = c[j] - range[mpi_rank];

          // Evaluate the point
          getLocation(conn, Xpts, &quads[i], order, j, knots[level],
                      &Xn[3 * index]);
        }
      }
    }

    // Set the node locations into TACSAssembler
    tacs[level]->setNodes(X);
  }

  // Create the interpolation
  TACSBVecInterp *interp[NUM_LEVELS - 1];

  for (int level = 0; level < NUM_LEVELS - 1; level++) {
    // Create the interpolation object
    interp[level] = new TACSBVecInterp(tacs[level + 1]->getVarMap(),
                                       tacs[level]->getVarMap(),
                                       tacs[level]->getVarsPerNode());
    interp[level]->incref();

    // Set the interpolation
    forest[level]->createInterpolation(forest[level + 1], interp[level]);

    // Initialize the interpolation
    interp[level]->initialize();
  }

  // Create a vector on the finest level
  TACSBVec *x[NUM_LEVELS];
  x[0] = tacs[0]->createVec();
  x[0]->incref();
  x[0]->setRand(-1.0, 1.0);
  for (int level = 0; level < NUM_LEVELS - 1; level++) {
    x[level + 1] = tacs[level + 1]->createVec();
    x[level + 1]->incref();
    interp[level]->multWeightTranspose(x[level], x[level + 1]);
  }

  // Create and write out an fh5 file
  unsigned int write_flag =
      (TACSElement::OUTPUT_NODES | TACSElement::OUTPUT_DISPLACEMENTS |
       TACSElement::OUTPUT_EXTRAS);

  for (int level = 0; level < NUM_LEVELS; level++) {
    tacs[level]->setVariables(x[level]);
    TACSToFH5 *f5 = new TACSToFH5(tacs[level], TACS_PLANE_STRESS, write_flag);
    f5->incref();

    // Write out the solution
    char filename[128];
    sprintf(filename, "output%d.f5", level);
    f5->writeToFile(filename);
    f5->decref();
  }

  // Create the level
  for (int level = 0; level < 3; level++) {
    x[level]->decref();
    tacs[level]->decref();
    delete forest[level];
  }

  TMRFinalize();
  MPI_Finalize();
  return (0);
}
