#include "TMROctree.h"
#include "TACSCreator.h"
#include "TACSAssembler.h"
#include "Solid.h"
#include "TACSToFH5.h"

int main( int argc, char * argv[] ){
  MPI_Init(&argc, &argv);

  // Allocate the TACS creator
  int vars_per_node = 3;
  TACSCreator *creator = new TACSCreator(MPI_COMM_WORLD, vars_per_node);
  creator->incref();

  // Create the octree
  TMROctree *tree = new TMROctree(500, 3, 6);
  tree->balance(7);

  // Create the mesh
  int order = 2;
  int num_nodes, num_dep_nodes, num_elems;
  int *elem_ptr, *elem_conn;
  int *dep_ptr, *dep_conn;
  double *dep_weights;
  
  // Create the mesh
  tree->createMesh(order, &num_nodes, 
		   &num_dep_nodes, &num_elems,
		   &elem_ptr, &elem_conn, 
		   &dep_ptr, &dep_conn, &dep_weights);

  // Allocate the node locations
  double *Xpts = new double[ 3*num_nodes ];
  
  // Get the array of octant nodes
  int size;
  TMROctant *array;
  TMROctantArray *nodes;
  tree->getNodes(&nodes);
  nodes->getArray(&array, &size);

  // Set the length of the cube
  double dh = 1.0/(1 << TMR_MAX_LEVEL);

  int num_bcs = 0;
  int *bc_nodes = new int[ num_nodes ];
  for ( int i = 0; i < size; i++ ){
    // Get the node location
    int node = array[i].tag;
    if (node > 0){
      Xpts[3*node] = dh*array[i].x;
      Xpts[3*node+1] = dh*array[i].y;
      Xpts[3*node+2] = dh*array[i].z;
    
      // Add the boundary condition
      if (array[i].z == 0){
	bc_nodes[num_bcs] = node;
	num_bcs++;
      }
    }
  }

  // Set all the element ids to 0
  int *elem_id_nums = new int[ num_elems ];
  memset(elem_id_nums, 0, num_elems*sizeof(int));

  // Create the Solid element in TACS
  SolidStiffness *stiff = new SolidStiffness(1.0, 70e3, 0.3);
  Solid<2> *elem = new Solid<2>(stiff);
  elem->incref();

  // Set the connectivity
  creator->setGlobalConnectivity(num_nodes, num_elems,
				 elem_ptr, elem_conn,
				 elem_id_nums);
  
  // Set the boundary conditions
  creator->setBoundaryConditions(num_bcs, bc_nodes, NULL, NULL);
  
  // Set the nodal locations
  creator->setNodes(Xpts);

  // Set the dependent nodes
  creator->setDependentNodes(num_dep_nodes, dep_ptr,
			     dep_conn, dep_weights);

  TACSElement *elements = elem;

  // This call must occur on all processor
  creator->setElements(&elements, 1);

  elem->decref();
  
  // Create the TACSAssembler object
  TACSAssembler *tacs = creator->createTACS();
  tacs->incref();

  // Create the preconditioner
  BVec *res = tacs->createVec();
  BVec *ans = tacs->createVec();
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

  res->set(1.0);
  res->applyBCs();
  pc->applyFactor(res, ans);
  tacs->setVariables(ans);

  // Create an TACSToFH5 object for writing output to files
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS |
                             TACSElement::OUTPUT_STRAINS |
                             TACSElement::OUTPUT_STRESSES |
                             TACSElement::OUTPUT_EXTRAS);
  TACSToFH5 * f5 = new TACSToFH5(tacs, SOLID, write_flag);
  f5->incref();
  f5->writeToFile("octree.f5");

  // Free everything
  f5->decref();


  /*
  TMROctree *next = NULL;
  for ( int k = 0; k < 4; k++ ){
    TMROctree *tree = NULL;
    if (k == 0){
      tree = new TMROctree(500, 3, 4);
      double t1 = MPI_Wtime();
      tree->balance(7);
      t1 = MPI_Wtime() - t1;

      if (next){ delete next; }
      next = tree;

      printf("Balance time: %15.8f\n", t1);
    }
    else {
      double t1 = MPI_Wtime();
      tree = next->coarsen();
      double t2 = MPI_Wtime();
      tree->balance(1);
      double t3 = MPI_Wtime();      
      next = tree;

      printf("Coarsen time: %15.8f\n", t2-t1);
      printf("Balance time: %15.8f\n", t3-t2);
    }

    // Get the element connectivity
    int order = 2;
    int num_nodes, num_dep_nodes, num_elems;
    int *elem_ptr, *elem_conn;
    int *dep_ptr, *dep_conn;
    double *dep_weights;
    
    // Create the mesh
    double t1 = MPI_Wtime();
    tree->createMesh(order, &num_nodes, 
                     &num_dep_nodes, &num_elems,
                     &elem_ptr, &elem_conn, 
		     &dep_ptr, &dep_conn, &dep_weights);
    t1 = MPI_Wtime() - t1;
    printf("Nodes time:   %15.8f\n", t1);

    int size;
    TMROctant *array;
    TMROctantArray *nodes;
    tree->getNodes(&nodes);
    nodes->getArray(&array, &size);

    char filename[128];
    sprintf(filename, "octree%d.dat", k);
    FILE * fp = fopen(filename, "w");
    if (fp){    
      fprintf(fp, "Variables = X, Y, Z\n");
      fprintf(fp, "ZONE T=TMR N=%d E=%d ", 
              num_nodes + num_dep_nodes, num_elems);
      fprintf(fp, "DATAPACKING=POINT ZONETYPE=FEBRICK\n");

      // Set the length of the cube
      double dh = 1.0/(1 << TMR_MAX_LEVEL);
      
      for ( int i = 0; i < size; i++ ){
        if (array[i].tag >= 0){
          fprintf(fp, "%e %e %e\n", 
                  dh*array[i].x, dh*array[i].y, dh*array[i].z);
        }
      }
      for ( int i = 0; i < size; i++ ){
        if (array[i].tag < 0){
          fprintf(fp, "%e %e %e\n", 
                  dh*array[i].x, dh*array[i].y, dh*array[i].z);
        }
      }
      
      const int node_to_tec[8] = 
        {0, 1, 3, 2, 4, 5, 7, 6};
      for ( int i = 0; i < num_elems; i++ ){
        for ( int jp = 0; jp < 8; jp++ ){
          int node = elem_conn[elem_ptr[i] + node_to_tec[jp]];
          if (node < 0){ node = num_nodes - node - 1; }
          fprintf(fp, "%d ", node+1);
        }
        fprintf(fp, "\n");
      }
      fclose(fp);
    }

    delete [] elem_ptr;
    delete [] elem_conn;
  }

  if (next){ delete next; }
  */

  MPI_Finalize();
  return (0);
}
