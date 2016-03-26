#include "TMROctree.h"
#include "TACSCreator.h"
#include "TACSAssembler.h"
#include "Solid.h"
#include "BVecInterp.h"
#include "TACSMg.h"
#include "TACSToFH5.h"
#include "OctStiffness.h"
#include "TMRTopoOpt.h"

/*
  TODOs:

  Filter:
  1) Create nested filter objects for each grid level 
  2) Create the interpolation/restrictions between these levels

  Design:
  1) Test the design derivatives with the adjoint equations

  Refinement: 
  1) Choose a rational refinement criteria which depends on some
  combination of the design and the solution.

  - If compliance is the objective, refinement with a fixed design
  should produce a design that is more compliant. This must be true
  from finte-element theory. However, if the objective is to minimize
  the compliance, and we want a more accurate result, some
  adjoint-based mesh refinement criteria could be used.

  Long term:
  1) Implement parallel design variable operations. Need to think
  about how these will be handled in TACS.
*/

/*
  Create the elements that are required based on the global filter
*/
static TMROctree *global_filter;
static TMROctant *global_octant_list;

static TACSElement* create_element( int local_num, int id ){
  TacsScalar density = 1.0;
  TacsScalar E = 70e3;
  TacsScalar nu = 0.3;
  TacsScalar q = 5.0;  
  SolidStiffness *stiff = new OctStiffness(global_filter, 
                                           &global_octant_list[local_num],
                                           density, E, nu, q);
  Solid<2> *elem = new Solid<2>(stiff);
  return elem;
}

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

  // Set all the element ids
  int *elem_id_nums = new int[ num_elems ];
  for ( int i = 0; i < num_elems; i++ ){
    elem_id_nums[i] = i;
  }

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

  // Partition the mesh
  creator->partitionMesh();
}

/*
  Run the main function
*/
int main( int argc, char * argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  MPI_Comm comm = MPI_COMM_WORLD;

  // Get the rank 
  int mpi_size, mpi_rank;
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);

  // The order of the mesh
  int order = 2;

  // The number of variables per node
  int vars_per_node = 3;

  // Set the levels
  const int NUM_LEVELS = 3;

  // Set the minimum and maximum refinement levels within the
  // octree meshes
  const int min_refine = 4;
  const int max_refine = 10;

  // Create the filter octrees
  TMROctree *filters[NUM_LEVELS];

  // Create the finest octree for the design variables
  filters[0] = new TMROctree(min_refine);
  filters[0]->createMesh(order);

  // Coarsen the octree filter
  for ( int k = 1; k < NUM_LEVELS; k++ ){
    filters[k] = filters[k-1]->coarsen();
    filters[k]->balance(0);
    filters[k]->createMesh(order);
  }

  // Set up the octree for the analysis
  TMROctree *root_tree = NULL;
  if (mpi_rank == 0){
    // Create a uniformly refined octree
    root_tree = new TMROctree(min_refine);
  }

  for ( int iter = 0; iter < 3; iter++ ){
    // Keep pointers to the octree, creator and assembler objects
    TMROctree *tree = root_tree;
    TACSCreator *creator[NUM_LEVELS];
    TACSAssembler *tacs[NUM_LEVELS];
    BVecInterp *restrct[NUM_LEVELS-1], *interp[NUM_LEVELS-1];

    for ( int k = 0; k < NUM_LEVELS; k++ ){
      TMROctree *fine = NULL;
      // Create the finest octant tree
      if (mpi_rank == 0){
        if (k == 0){
          tree->balance(1);
        }
        else {
          fine = tree;
          tree = fine->coarsen();
          tree->balance(0);
        }
      }

      // Allocate the TACS creator
      creator[k] = new TACSCreator(comm, vars_per_node);
      creator[k]->incref();
	
      // Set up the TACSCreator object
      if (mpi_rank == 0){
        setUpTACSCreator(order, creator[k], tree);
      
        // Scatter the number of locally owned elements
        int num_owned = 0;
        int *owned_elements;
        creator[k]->getNumOwnedElements(&owned_elements);
        MPI_Scatter(owned_elements, 1, MPI_INT,
                    &num_owned, 1, MPI_INT, 0, comm);

        // Compute the maximum of locally owned
        int max_owned = 0;
        for ( int i = 0; i < mpi_size; i++ ){
          if (owned_elements[i] > max_owned){
            max_owned = owned_elements[i];
          }
        }
        
        // Allocate the local partition
        TMROctant *local = new TMROctant[ max_owned ];

        // Get the local element partition
        const int *partition = NULL;
        creator[k]->getElementPartition(&partition);;
        
        TMROctantArray *elems;
        tree->getElements(&elems);
        TMROctant *elements;
        elems->getArray(&elements, NULL);

        // Scan over each of the non-local processors
        for ( int i = mpi_size-1; i >= 0; i-- ){
          int j = 0, index = 0;
          while (index < owned_elements[i]){
            if (partition[j] == i){
              local[index] = elements[j];
              index++;
            }
            j++;
          }

          // Send the octants for the given partition to the
          // external processors
          if (i != 0){
            int tag = 0;
            MPI_Send(local, owned_elements[i], TMROctant_MPI_type,
                     i, tag, comm);
          }
        }
        
        // Set the elements and create TACS
        global_filter = filters[k];
        global_octant_list = local;
        creator[k]->setElementCreator(create_element);

	// Create the TACSAssembler object
	tacs[k] = creator[k]->createTACS();
	tacs[k]->incref();

        // Null the arrays so they are not use accidentally
        global_filter = NULL;
        global_octant_list = NULL;

        // Free the local octant array
        delete [] local;
      }
      else {
        int num_owned = 0;
        MPI_Scatter(NULL, 1, MPI_INT,
                    &num_owned, 1, MPI_INT, 0, comm);

        TMROctant *local = new TMROctant[ num_owned ];
        int tag = 0;
        MPI_Recv(local, num_owned, TMROctant_MPI_type, 0, tag, comm, 
                 MPI_STATUS_IGNORE);
        
        // Set the elements and the callback object
        global_filter = filters[k];
        global_octant_list = local;
        creator[k]->setElementCreator(create_element);

	// Create the TACSAssembler object
	tacs[k] = creator[k]->createTACS();
	tacs[k]->incref();

        // Null the arrays so they are not use accidentally
        global_filter = NULL;
        global_octant_list = NULL;

        // Free the local octant array
        delete [] local;
      }

      if (k > 0){
        // Create the interpolation/restriction objects
        restrct[k-1] = new BVecInterp(tacs[k-1], tacs[k]);
        restrct[k-1]->incref();      
        interp[k-1] = new BVecInterp(tacs[k], tacs[k-1]);
        interp[k-1]->incref();

        if (mpi_rank == 0){      
          // Get the node numbers for the coarse and fine operators
          const int *nodes, *coarse_nodes;
          int num_nodes = creator[k-1]->getNodeNums(&nodes);
          int num_coarse_nodes = creator[k]->getNodeNums(&coarse_nodes);
          
          // Create the restriction operator
          int *restrct_ptr, *restrct_conn;
          double *restrct_weights;
          fine->createRestriction(tree, &restrct_ptr,
                                  &restrct_conn, &restrct_weights);
	
          for ( int i = 0; i < num_coarse_nodes; i++ ){
            // Convert the node numbers
            for ( int jp = restrct_ptr[i]; jp < restrct_ptr[i+1]; jp++ ){
              restrct_conn[jp] = nodes[restrct_conn[jp]]; 
            }
            
            // Set the values
            int len = restrct_ptr[i+1] - restrct_ptr[i];
            restrct[k-1]->addInterp(coarse_nodes[i],
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
          fine->createInterpolation(tree, &interp_ptr, 
                                    &interp_conn, &interp_weights);
          
          for ( int i = 0; i < num_nodes; i++ ){
            // Convert the node numbers
            for ( int jp = interp_ptr[i]; jp < interp_ptr[i+1]; jp++ ){
              interp_conn[jp] = coarse_nodes[interp_conn[jp]]; 
            }
            
            // Set the values
            int len = interp_ptr[i+1] - interp_ptr[i];
            interp[k-1]->addInterp(nodes[i],
                                   &interp_weights[interp_ptr[i]],
                                   &interp_conn[interp_ptr[i]], len);
          }
          
          // Free the data
          delete [] interp_ptr;
          delete [] interp_conn;
          delete [] interp_weights;
        }
        
        // Finalize the values
        restrct[k-1]->finalize();
        interp[k-1]->finalize();
      }

      // Free the tree
      if (fine != root_tree){
	delete fine;
      }
    }
    
    // Set up the multigrid method
    double sor_omega = 1.2;
    int sor_iters = 1;
    int sor_symm = 1;
    TACSMg *mg = new TACSMg(comm, NUM_LEVELS, 
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

    // Set the KSM monitor
    int freq = 10;
    mg->setMonitor(new KSMPrintStdout("MG", mpi_rank, freq));

    // Create the force vector
    BVec *force = tacs[0]->createVec();

    // Get the ownership range
    const int *owners;
    force->getOwnership(NULL, NULL, &owners);

    // Determine the force node
    int force_node = 0;
    if (mpi_rank == 0){
      // Set the point in the octant where we want to
      // set the force vector
      TMROctant p;
      p.x = 1 << TMR_MAX_LEVEL;
      p.y = 1 << TMR_MAX_LEVEL;
      p.z = 1 << TMR_MAX_LEVEL;
      
      // Find the vector in the tree
      const int use_nodes = 1;
      TMROctantArray *nodes;
      root_tree->getNodes(&nodes);
      TMROctant *t = nodes->contains(&p, use_nodes);

      // Set the force node
      force_node = t->tag;
    }

    // Broadcast the node
    MPI_Bcast(&force_node, 1, MPI_INT, 0, comm);

    // Find the ownership range
    if (force_node >= owners[mpi_rank] &&
        force_node < owners[mpi_rank+1]){
      TacsScalar *force_vals;
      force->getArray(&force_vals);
      int loc = force_node - owners[mpi_rank];
      force_vals[3*loc+1] = 1e3;
      force_vals[3*loc+2] = 1e3;
    }

    // Create the topology optimization problem
    TMRTopoOpt *problem = new TMRTopoOpt(NUM_LEVELS, tacs, filters,
                                         mg, force);

    // Allocate the ParOpt object
    int max_num_bfgs = 40;
    ParOpt *opt = new ParOpt(problem, max_num_bfgs);
    opt->checkGradients(1e-6);
    opt->setMaxMajorIterations(20);

    opt->optimize();

    delete opt;
    delete problem;

    // Perform a local refinement of the nodes based on the strain energy
    // within each element
    int nelems = tacs[0]->getNumElements();
    TacsScalar *SE = new TacsScalar[ nelems ];
    TacsScalar SE_proc_sum = 0.0;

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
      SE_proc_sum += SE[i];
    }

    // Count up the total strain energy 
    TacsScalar SE_sum = 0.0;
    MPI_Reduce(&SE_proc_sum, &SE_sum, 1, TACS_MPI_TYPE, MPI_SUM, 0, comm);

    // Count up the total number of elements
    int nelems_total = 0;
    MPI_Reduce(&nelems, &nelems_total, 1, MPI_INT, MPI_SUM, 0, comm);

    // Gather the number of elements from each processor
    int *elems_per_proc = NULL;
    int *elem_proc_offset = NULL;
    TacsScalar *all_SE = NULL;

    if (mpi_rank == 0){
      elems_per_proc = new int[ mpi_size ];
      elem_proc_offset = new int[ mpi_size ];
      all_SE = new TacsScalar[ nelems_total ];
    }

    // Gather the number of elements from each processor to the root
    // processor
    MPI_Gather(&nelems, 1, MPI_INT, elems_per_proc, 1, MPI_INT, 0, comm);

    if (mpi_rank == 0){
      elem_proc_offset[0] = 0;
      for ( int i = 0; i < mpi_size-1; i++ ){
	elem_proc_offset[i+1] = elem_proc_offset[i] + elems_per_proc[i];
      }
    }
    
    // Gather the element strain energy values to the root processor
    MPI_Gatherv(SE, nelems, TACS_MPI_TYPE, 
	       all_SE, elems_per_proc, elem_proc_offset, TACS_MPI_TYPE, 
                0, comm);

    // Free the element strain energy values for this processor
    delete [] SE;

    if (mpi_rank == 0){
      // Set the cut-off strain energy value
      TacsScalar SE_cut_off = 0.5*SE_sum/nelems_total;

      // Set the refinement
      int *refinement = new int[ nelems_total ];
      memset(refinement, 0, nelems_total*sizeof(int));

      // Get the element partition from the TACSCreator object
      const int *partition;
      creator[0]->getElementPartition(&partition);

      for ( int i = 0; i < nelems_total; i++ ){
	// Get the partition that this element is on
	int proc = partition[i];

	// Retrieve the offset from the 
	int index = elem_proc_offset[proc];
	elem_proc_offset[proc]++;

	if (all_SE[index] > SE_cut_off){
	  refinement[i] = 1;
	}
	else {
	  refinement[i] = -1;
	}
      }

      // Free the strain energy information and offset values
      delete [] elems_per_proc;
      delete [] elem_proc_offset;
      delete [] all_SE;

      // Refine the root tree
      root_tree->refine(refinement, min_refine, max_refine);
      
      delete [] refinement;
    }

    // Free everything 
    mg->decref();
  
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

  if (root_tree){ delete root_tree; } 

  TMRFinalize();
  MPI_Finalize();
  return (0);
}
