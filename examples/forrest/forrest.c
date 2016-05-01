#include "TMRForrest.h"
#include "isoFSDTStiffness.h"
#include "MITCShell.h"
#include "TACSCreator.h"
#include "TACSAssembler.h"
#include "PlaneStressQuad.h"
#include "TACSToFH5.h"

// Create the nodal locations for the mesh
const double test_Xpts[] = {0.0, 0.0, 0.0,
                            1.0, 0.0, 0.0,
                            0.3, 0.7, 0.0, 
                            0.8, 0.2, 0.0, 
                            0.25, 0.21, 0.0,
                            0.75, 0.6, 0.0,
                            0.0, 1.0, 0.0,
                            1.0, 1.0, 0.0};

const int test_conn[] = {0, 1, 4, 3,
                         2, 4, 5, 3, 
                         6, 0, 2, 4, 
                         2, 5, 6, 7,
                         3, 1, 5, 7};

/*
  Retrieve the x/y location on the face based on the u/v coordinates
*/
void get_location( int face, double u, double v, 
                   double *X ){
  double N[4];
  N[0] = (1.0 - u)*(1.0 - v);
  N[1] = u*(1.0 - v);
  N[2] = (1.0 - u)*v;
  N[3] = u*v;

  X[0] = (N[0]*test_Xpts[3*test_conn[4*face]] +
          N[1]*test_Xpts[3*test_conn[4*face+1]] +
          N[2]*test_Xpts[3*test_conn[4*face+2]] +
          N[3]*test_Xpts[3*test_conn[4*face+3]]);
  X[1] = (N[0]*test_Xpts[3*test_conn[4*face]+1] +
          N[1]*test_Xpts[3*test_conn[4*face+1]+1] +
          N[2]*test_Xpts[3*test_conn[4*face+2]+1] +
          N[3]*test_Xpts[3*test_conn[4*face+3]+1]);
  X[2] = (N[0]*test_Xpts[3*test_conn[4*face]+2] +
          N[1]*test_Xpts[3*test_conn[4*face+1]+2] +
          N[2]*test_Xpts[3*test_conn[4*face+2]+2] +
          N[3]*test_Xpts[3*test_conn[4*face+3]+2]);
}

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);

  srand(time(NULL));

  const int ORDER = 3;

  // Set the communicator
  MPI_Comm comm = MPI_COMM_WORLD;

  // Get the MPI communicator rank
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Create the TACSCreator object
  int vars_per_node = 6;
  TACSCreator *creator = new TACSCreator(comm, vars_per_node);
  creator->incref();

  if (rank == 0){
    TMRQuadForrest *forrest = new TMRQuadForrest(MPI_COMM_SELF);
    
    // Set the connectivity
    int num_nodes = 8;
    int num_faces = 5;
    forrest->setConnectivity(num_nodes, 
                             test_conn, num_faces);
    
    // Allocate the trees (random trees for now)
    int nrand = 25;
    int min_refine = 0;
    int max_refine = 5;
    forrest->createRandomTrees(nrand, min_refine, max_refine);
    
    // Balance the forrest so that we can use it!
    forrest->balance(1);
    
    // Create the nodes
    forrest->createNodes(ORDER);
    
    // Extract the mesh
    int nnodes, nelems;
    int *elem_ptr, *elem_conn;
    forrest->getMesh(&nnodes, &nelems,
                     &elem_ptr, &elem_conn);
    
    // Set the element id numbers
    int *elem_id_nums = new int[ nelems ];
    memset(elem_id_nums, 0, nelems*sizeof(int));
    
    // Create the Xpts array
    double *Xpts = new double[ 3*nnodes ];
    
    // Fill in the positions
    TMRQuadtree **quad;
    forrest->getQuadtrees(&quad);
    
    int num_bcs = 0;
    int max_bcs = 10000;
    int *bc_nodes = new int[ max_bcs ];
    
    int elem = 0;
    for ( int face = 0; face < num_faces; face++ ){
      // Retrieve the node quadrants
      TMRQuadrantArray *nodes;
      quad[face]->getNodes(&nodes);
      
      // Get the array of node quadrants from this face
      int size;
      TMRQuadrant *array;
      nodes->getArray(&array, &size);
      
      // Iterate through and evaluate the x/y/z locations
      const double dh = 1.0/(1 << TMR_MAX_LEVEL);
      for ( int i = 0; i < size; i++ ){
        int node = array[i].tag;
        if (node >= 0){
          double u = dh*array[i].x;
          double v = dh*array[i].y;         
          get_location(face, u, v, &Xpts[3*node]);
          
          // Set the boundary conditions based on the spatial location
          if (((Xpts[3*node] < 1e-6 || Xpts[3*node] > 0.999999) ||
               (Xpts[3*node+1] < 1e-6 || Xpts[3*node+1] > 0.999999))
              && num_bcs < max_bcs){
            bc_nodes[num_bcs] = node;
            num_bcs++;
          }
        }
      }
    }
    
    // Set the connectivity
    creator->setGlobalConnectivity(nnodes, nelems,
                                   elem_ptr, elem_conn,
                                   elem_id_nums);
    
    // Set the boundary conditions
    creator->setBoundaryConditions(num_bcs, bc_nodes, NULL, NULL);
    
    // Set the nodal locations
    creator->setNodes(Xpts);
    
    // Set the dependent nodes
    int num_dep_nodes;
    const int *dep_ptr;
    const int *dep_conn;
    const double *dep_weights;
    forrest->getDependentNodes(&num_dep_nodes, &dep_ptr,
                               &dep_conn, &dep_weights);
    creator->setDependentNodes(num_dep_nodes, dep_ptr,
                               dep_conn, dep_weights);
    
    // Partition the mesh
    creator->partitionMesh();

    delete forrest;
  }

  // Create the element
  // PlaneStressStiffness *stiff = new PlaneStressStiffness(1.0, 70e3, 0.3);
  // TACSElement *element = new PlaneStressQuad<ORDER>(stiff);
  isoFSDTStiffness *stiff = new isoFSDTStiffness(1.0, 70e9, 0.3,
                                                 5.0/6.0, 350e6, 0.01);
  TACSElement *element = new MITCShell<ORDER>(stiff);

  TestElement *test = new TestElement(element);
  test->setPrintLevel(2);
  test->testResidual();
  test->testJacobian();


  creator->setElements(&element, 1);
  
  // Create the TACSAssembler object
  TACSAssembler *tacs = creator->createTACS();
  tacs->incref();

  // Create the preconditioner
  BVec *res = tacs->createVec();
  BVec *ans = tacs->createVec();
  BVec *tmp = tacs->createVec();
  FEMat *mat = tacs->createFEMat();

  // Increment the reference count to the matrix/vectors
  res->incref();
  ans->incref();
  mat->incref();

  // Allocate the factorization
  int lev = 4500;
  double fill = 10.0;
  int reorder_schur = 1;
  PcScMat *pc = new PcScMat(mat, lev, fill, reorder_schur);
  pc->incref();

  // Assemble and factor the stiffness/Jacobian matrix
  double alpha = 1.0, beta = 0.0, gamma = 0.0;
  tacs->assembleJacobian(res, mat, alpha, beta, gamma);
  mat->applyBCs();
  pc->factor();

  TacsScalar *array;
  int size = res->getArray(&array);
  for ( int k = 2; k < size; k += 6){
    array[k] = 1.0;
  }
  res->applyBCs();
  pc->applyFactor(res, ans);
  tacs->setVariables(ans);

  mat->mult(ans, tmp);
  tmp->axpy(-1.0, res);

  // Print the solution norm
  TacsScalar norm = tmp->norm()/res->norm();
  if (rank == 0){
    printf("Solution residual norm: %15.5e\n", norm);
  }

  // Create an TACSToFH5 object for writing output to files
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS |
                             TACSElement::OUTPUT_STRAINS |
                             TACSElement::OUTPUT_STRESSES |
                             TACSElement::OUTPUT_EXTRAS);
  // TACSToFH5 * f5 = new TACSToFH5(tacs, PLANE_STRESS, write_flag);
  TACSToFH5 * f5 = new TACSToFH5(tacs, SHELL, write_flag);
  f5->incref();
  f5->writeToFile("forrest.f5");

  // Free everything
  f5->decref();

  tacs->decref();
  creator->decref();

  MPI_Finalize();
  return (0);
}
