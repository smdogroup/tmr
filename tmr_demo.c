#include "TMROctree.h"
#include "TACSCreator.h"
#include "TACSAssembler.h"
#include "Solid.h"
#include "BVecInterp.h"
#include "TACSMg.h"
#include "TACSToFH5.h"

/*
  Set up the TACSCreator object and set the internal data required to
  create the TACSAssembler object
*/
void setUpTACSCreator( int order, TACSCreator *creator, 
                       TMROctree *tree ){  
  // Create the mesh
  double t1 = MPI_Wtime();
  tree->createMesh(order);
  t1 = MPI_Wtime() - t1;
  printf("Time spent creating the mesh %f [s]\n", t1);

  // Get the mesh connectivity information
  int num_nodes, num_dep_nodes, num_elems;
  const int *elem_ptr, *elem_conn;
  const int *dep_ptr, *dep_conn;
  const double *dep_weights;
  tree->getMesh(&num_nodes, &num_elems, 
                &elem_ptr, &elem_conn);
  tree->getDependentMesh(&num_dep_nodes,
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
    if (node >= 0){
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
}

/*
  Run the main function
*/
int main( int argc, char * argv[] ){
  MPI_Init(&argc, &argv);

  // Get the rank 
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Create the Solid element in TACS
  SolidStiffness *stiff = new SolidStiffness(1.0, 70e3, 0.3);
  Solid<2> *elem = new Solid<2>(stiff);
  elem->incref();

  // The order of the mesh
  int order = 2;

  // The number of variables per node
  int vars_per_node = 3;

  // Set the levels
  const int NUM_LEVELS = 4;

  TMROctree *root_tree = NULL;
  if (rank == 0){
    // Create a uniformly refined octree
    root_tree = new TMROctree(3);
  }

  for ( int iter = 0; iter < 4; iter++ ){
    // Keep pointers to the octree, creator and assembler objects
    TMROctree *tree = root_tree;
    TACSCreator *creator[NUM_LEVELS];
    TACSAssembler *tacs[NUM_LEVELS];
    BVecInterp *restrct[NUM_LEVELS-1], *interp[NUM_LEVELS-1];

    for ( int k = 0; k < NUM_LEVELS-1; k++ ){
      // Create the finest octant tree
      if (k == 0){
	if (rank == 0){
	  tree->balance(1);
	}

	// Allocate the TACS creator
	creator[k] = new TACSCreator(MPI_COMM_WORLD, vars_per_node);
	creator[k]->incref();
	
	// Set up the TACSCreator object
	if (rank == 0){
	  setUpTACSCreator(order, creator[k], tree);
	}

	// This call must occur on all processor
	TACSElement *elements = elem;
	creator[k]->setElements(&elements, 1);
  
	// Create the TACSAssembler object
	tacs[k] = creator[k]->createTACS(TACSAssembler::NATURAL_ORDER);
	tacs[k]->incref();
      }
      
      // Create the coarser octree
      TMROctree *coarse = NULL;
      if (rank == 0){
	coarse = tree->coarsen();
	coarse->balance(0);
      }
      
      // Allocate the TACS creator
      creator[k+1] = new TACSCreator(MPI_COMM_WORLD, vars_per_node);
      creator[k+1]->incref();
      
      // Set up the TACSCreator object
      if (rank == 0){
	setUpTACSCreator(order, creator[k+1], coarse);
      }
      
      // This call must occur on all processor
      TACSElement *elements = elem;
      creator[k+1]->setElements(&elements, 1);
      
      // Create the TACSAssembler object
      tacs[k+1] = creator[k+1]->createTACS(TACSAssembler::NATURAL_ORDER);
      tacs[k+1]->incref();
      
      // Create the interpolation/restriction objects
      restrct[k] = new BVecInterp(tacs[k], tacs[k+1]);
      restrct[k]->incref();      
      interp[k] = new BVecInterp(tacs[k+1], tacs[k]);
      interp[k]->incref();

      if (rank == 0){      
	// Get the node numbers for the coarse and fine operators
	const int *nodes, *coarse_nodes;
	int num_nodes = creator[k]->getNodeNums(&nodes);
	int num_coarse_nodes = creator[k+1]->getNodeNums(&coarse_nodes);
	
	// Create the restriction operator
	int *restrct_ptr, *restrct_conn;
	double *restrct_weights;
	tree->createRestriction(coarse, &restrct_ptr,
				&restrct_conn, &restrct_weights);
	
	for ( int i = 0; i < num_coarse_nodes; i++ ){
	  // Convert the node numbers
	  for ( int jp = restrct_ptr[i]; jp < restrct_ptr[i+1]; jp++ ){
	    restrct_conn[jp] = nodes[restrct_conn[jp]]; 
	  }
	  
	  // Set the values
	  int len = restrct_ptr[i+1] - restrct_ptr[i];
	  restrct[k]->addInterp(coarse_nodes[i],
				&restrct_weights[restrct_ptr[i]],
				&restrct_conn[restrct_ptr[i]], len);
	}
	
	// Free the data
	delete [] restrct_ptr;
	delete [] restrct_conn;
	delete [] restrct_weights;
	
	// Create the interpolation
	int *interp_ptr, *interp_conn;
	double *interp_weights;
	tree->createInterpolation(coarse, &interp_ptr, 
				  &interp_conn, &interp_weights);
	
	for ( int i = 0; i < num_nodes; i++ ){
	  // Convert the node numbers
	  for ( int jp = interp_ptr[i]; jp < interp_ptr[i+1]; jp++ ){
	    interp_conn[jp] = coarse_nodes[interp_conn[jp]]; 
	  }
	  
	  // Set the values
	  int len = interp_ptr[i+1] - interp_ptr[i];
	  interp[k]->addInterp(nodes[i],
			       &interp_weights[interp_ptr[i]],
			       &interp_conn[interp_ptr[i]], len);
	}
	
	// Free the data
	delete [] interp_ptr;
	delete [] interp_conn;
	delete [] interp_weights;
      }
      
      // Finalize the values
      restrct[k]->finalize();
      interp[k]->finalize();
      
      // Free the tree
      if (tree != root_tree){
	delete tree;
      }
      
      // Assign the coarser mesh for the next level
      tree = coarse;
    }
    
    // Set up the multigrid method
    double sor_omega = 1.2;
    int sor_iters = 2;
    int sor_symm = 1;
    TACSMg *mg = new TACSMg(MPI_COMM_WORLD, NUM_LEVELS, 
			    sor_omega, sor_iters, sor_symm);
    mg->incref();
    
    for ( int k = 0; k < NUM_LEVELS; k++ ){
      if (k < NUM_LEVELS-1){
	mg->setLevel(k, tacs[k], restrct[k], interp[k]);
      }
      else {
	mg->setLevel(k, tacs[k], NULL, NULL);
      }
    }
    
    // Assemble the matrix
    mg->assembleMatType();
    mg->factor();
    
    // Create the preconditioner
    BVec *res = tacs[0]->createVec();
    BVec *ans = tacs[0]->createVec();
    res->incref();
    ans->incref();

    int freq = 1;
    mg->setMonitor(new KSMPrintStdout("MG", rank, freq));
    
    res->set(1.0);
    res->applyBCs();
    
    // Now, set up the solver
    int gmres_iters = 100; 
    int nrestart = 2;
    int is_flexible = 1;
    GMRES * ksm = new GMRES(mg->getMat(0), mg, 
			    gmres_iters, nrestart, is_flexible);
    ksm->setTolerances(1e-8, 1e-30);
    ksm->incref();
    
    ksm->setMonitor(new KSMPrintStdout("MG", rank, freq));
    ksm->solve(res, ans);
    
    tacs[0]->setVariables(ans);
    
    // Create an TACSToFH5 object for writing output to files
    unsigned int write_flag = (TACSElement::OUTPUT_NODES |
			       TACSElement::OUTPUT_DISPLACEMENTS |
			       TACSElement::OUTPUT_STRAINS |
			       TACSElement::OUTPUT_STRESSES |
			       TACSElement::OUTPUT_EXTRAS);
    TACSToFH5 * f5 = new TACSToFH5(tacs[0], SOLID, write_flag);
    f5->incref();
    
    char filename[128];
    sprintf(filename, "octree_iter%d.f5", iter);
    f5->writeToFile(filename);

    // Based on the values of the local strain energy, perform
    // a local refinement
    int nelems = tacs[0]->getNumElements();
    TacsScalar *SE = new TacsScalar[ nelems ];
    TacsScalar SE_total = 0.0;    

    for ( int i = 0; i < nelems; i++ ){
      // Get the node locations and variables
      TacsScalar Xpts[3*8], vars[3*8], dvars[3*8], ddvars[3*8];
      TACSElement *elem = tacs[0]->getElement(i, Xpts, vars, dvars, ddvars);

      // Evaluate the element residual
      TacsScalar res[3*8];
      elem->getResidual(res, Xpts, vars, dvars, ddvars);

      // Take the inner product to get the strain energy
      SE[i] = 0.0;
      for ( int j = 0; j < 24; j++ ){
	SE[i] += res[j]*vars[j];
      }
      SE_total += SE[i];
    }

    SE_total = 0.5*SE_total/nelems;

    // Set the refinement
    int *refinement = new int[ nelems ];

    for ( int i = 0; i < nelems; i++ ){
      refinement[i] = -1;
      if (SE[i] > SE_total){
	refinement[i] = 1;
      }
    }

    // Refine the root tree
    root_tree->refine(refinement);

    // Free everything 
    mg->decref();
    f5->decref();
    res->decref();
    ans->decref();
    ksm->decref();

    // Decrease the reference count to things
    for ( int k = 0; k < NUM_LEVELS; k++ ){
      tacs[k]->decref();
      creator[k]->decref();
    }

    // Deallocate the interpolation and restriction operators
    for ( int k = 0; k < NUM_LEVELS-1; k++ ){
      interp[k]->decref();
      restrct[k]->decref();
    }
  }

  elem->decref();
  if (root_tree){ delete root_tree; } 

  MPI_Finalize();
  return (0);
}
