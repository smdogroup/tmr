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

void strain_energy_refinement( TACSAssembler *tacs, 
                               TACSCreator *creator,
                               TMROctree *root_tree,
                               int min_refine, int max_refine ){
  // Get the communicator
  MPI_Comm comm = tacs->getMPIComm();
  int mpi_size, mpi_rank;
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);

  // Perform a local refinement of the nodes based on the strain energy
  // within each element
  int nelems = tacs->getNumElements();
  TacsScalar *SE = new TacsScalar[ nelems ];
  TacsScalar SE_proc_sum = 0.0;
  
  for ( int i = 0; i < nelems; i++ ){
    // Get the node locations and variables
    TacsScalar Xpts[3*8], vars[3*8], dvars[3*8], ddvars[3*8];
    TACSElement *elem = tacs->getElement(i, Xpts, vars, dvars, ddvars);
    
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
    creator->getElementPartition(&partition);
    
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
}

/*
  The following code sets up an adaptive refinement technique for
  topology optimization problems on octree meshes.

  The code automatically refines the mesh around the location where
  the structure is present. The underlying filter is fixed (much like
  a level-set type method). Nested meshes are formed by repeated
  coarsening of the mesh.

  Local length-scale control could be implemented by refining the
  filter mesh in areas of interest - for instance around corners that
  might have more-detailed features. The code uses a compliance
  objective with a mass constraint.

  Useage:

  ./tmr_demo [min_refine=%d] [target_mass=%f] [problem_type]

  where problem_type is one of: block, beam, bracket

  which produce different problem domains.
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
  const int NUM_LEVELS = 4;

  // Set the minimum and maximum refinement levels within the
  // octree meshes
  int min_refine = 4;
  int max_refine = 8;

  // Create the filter octrees
  TMROctree *filters[NUM_LEVELS];

  // The ParOpt problem instance
  ParOpt *opt = NULL;

  // Set up the domain - if any
  int ndomain = 0;
  TMROctant *domain = NULL;

  // Set the default target mass
  double target_mass = 0.15;

  // Set the problem instance
  int problem_index = 0;
  const char *problem_type[] = 
    {"block", "beam", "bracket"};

  // Set up the domain for the beam
  TMROctant beam[4];
  for ( int i = 0; i < 4; i++ ){
    beam[i].level = 2;
    int32_t hbeam = 1 << (TMR_MAX_LEVEL - beam[i].level);
    beam[i].x = 0;
    beam[i].y = 0;
    beam[i].z = i*hbeam;
  }

  // Set up the domain for the bracket
  TMROctant bracket[4];
  const int bracket_corner[][3] = 
    {{0, 0, 0}, {0, 0, 1}, {0, 1, 1}, {1, 0, 1}};

  for ( int i = 0; i < 4; i++ ){
    bracket[i].level = 1;
    int32_t hbracket = 1 << (TMR_MAX_LEVEL - bracket[0].level);
    bracket[i].x = bracket_corner[i][0]*hbracket;
    bracket[i].y = bracket_corner[i][1]*hbracket;
    bracket[i].z = bracket_corner[i][2]*hbracket;
  }

  // Scan the input arguments to find the minimum refinement
  for ( int k = 0; k < argc; k++ ){
    int tmp = 0;
    if (sscanf(argv[k], "min_refine=%d", &tmp) == 1){
      // Let's not get crazy setting a really high refinement
      if (tmp <= 6){ 
        max_refine += (tmp - min_refine);
        min_refine = tmp;
      }
    }
    double tmass;
    if (sscanf(argv[k], "target_mass=%lf", &tmass) == 1){
      if (tmass > 0.05 && tmass < 1.0){
        target_mass = tmass;
      }
    }
  }

  // Set the domain for the
  for ( int k = 0; k < argc; k++ ){
    if (strcmp("block", argv[k]) == 0){ break; }
    if (strcmp("beam", argv[k]) == 0){
      // Set the problem index
      problem_index = 1;

      // Set the domain
      ndomain = 4;
      domain = beam;

      // Set the new target mass
      target_mass *= (0.25*0.25);

      // Change the minimum refinement level to match the
      // new domain
      min_refine += 2;
      max_refine += 2;
      break;
    }
    if (strcmp("bracket", argv[k]) == 0){
      // Set the problem index
      problem_index = 2;

      // Set the domain
      ndomain = 4;
      domain = bracket;

      // Set the new target mass
      target_mass *= 0.25;

      // Change the minimum refinement level to match the
      // new domain
      min_refine += 1;
      max_refine += 1;
      break;      
    }
  }

  // Write out the problem properties
  if (mpi_rank == 0){
    printf("Problem type:         %s\n", problem_type[problem_index]);
    printf("Min refinement:       %d\n", min_refine);
    printf("Target mass fraction: %f\n", target_mass);
  }

  // Create the finest octree for the design variables
  filters[0] = new TMROctree(min_refine, domain, ndomain);
  filters[0]->createMesh(order);

  // Coarsen the octree filter
  for ( int k = 1; k < NUM_LEVELS; k++ ){
    filters[k] = filters[k-1]->coarsen();
    filters[k]->balance(0);
    filters[k]->createMesh(order);
  }

  // Get the number of design variables from the filter class
  int num_design_vars = 0;
  filters[0]->getMesh(&num_design_vars, NULL, NULL, NULL);
  
  // Keep track of the design variables between iterations
  TacsScalar *x_design = new TacsScalar[ num_design_vars ];
  for ( int i = 0; i < num_design_vars; i++ ){
    x_design[i] = 0.95;
  }

  // Set up the octree for the analysis
  TMROctree *root_tree = NULL;
  if (mpi_rank == 0){
    // Create a uniformly refined octree
    root_tree = new TMROctree(min_refine, domain, ndomain);
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
    
    // Create the force vector
    BVec *force = tacs[0]->createVec();

    // Get the ownership range
    const int *owners;
    force->getOwnership(NULL, NULL, &owners);

    // Determine the nodes at which to apply forces
    int force_nodes[2] = {0, 0};
    if (mpi_rank == 0){     
      // Get the nodes from the octree
      const int use_nodes = 1;
      TMROctantArray *nodes;
      root_tree->getNodes(&nodes);

      // Get the node corresponding to the numbering
      const int *node_nums;
      creator[0]->getNodeNums(&node_nums);

      // Set the level depending on the problem type
      int level = 0;
      if (problem_index == 1){ // Beam-type problem
        level = 2;
      }

      // Set the point in the octant where we want to
      // set the force vector
      TMROctant p, *t;
      p.x = 0;
      p.y = 1 << (TMR_MAX_LEVEL - level);
      p.z = 1 << TMR_MAX_LEVEL;
      t = nodes->contains(&p, use_nodes);
      force_nodes[0] = node_nums[t->tag];

      p.x = 1 << (TMR_MAX_LEVEL - level);
      p.y = 0;
      p.z = 1 << TMR_MAX_LEVEL;
      t = nodes->contains(&p, use_nodes);
      force_nodes[1] = node_nums[t->tag];
    }

    // Broadcast the node
    MPI_Bcast(force_nodes, 2, MPI_INT, 0, comm);

    // Find the ownership range
    for ( int k = 0; k < 2; k++ ){
      if (force_nodes[k] >= owners[mpi_rank] &&
          force_nodes[k] < owners[mpi_rank+1]){
        TacsScalar *force_vals;
        force->getArray(&force_vals);
        int loc = force_nodes[k] - owners[mpi_rank];
        if (k == 0){
          force_vals[3*loc] = -10.0;
        }
        else if (k == 1){
          force_vals[3*loc+1] = 10.0;
        }
      }
    }

    // Decrease the reference count to the creator objects - they are
    // not required anymore
    for ( int k = 0; k < NUM_LEVELS; k++ ){
      creator[k]->decref();
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

    // Create an TACSToFH5 object for writing output to files
    unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                               TACSElement::OUTPUT_DISPLACEMENTS |
                               TACSElement::OUTPUT_EXTRAS);
    TACSToFH5 *f5 = new TACSToFH5(tacs[0], SOLID, write_flag);

    // Create the topology optimization problem
    char prefix[128];
    sprintf(prefix, "results/%s_refine%d_", 
            problem_type[problem_index], iter);
    TMRTopoOpt *problem = new TMRTopoOpt(NUM_LEVELS, tacs, filters,
                                         mg, force, target_mass, f5, prefix);

    // Set the initial design variable values
    problem->setInitDesignVars(x_design);

    // Allocate the ParOpt object
    if (!opt){
      int max_num_bfgs = 10;
      opt = new ParOpt(problem, max_num_bfgs);
      opt->setMaxMajorIterations(401);
    }
    else {
      opt->resetProblemInstance(problem);
      opt->setInitStartingPoint(0);
      opt->setMaxMajorIterations(201);
    }
    opt->checkGradients(1e-6);
    opt->setAbsOptimalityTol(1e-5);
    opt->optimize();

    // Get the final optimized design variable values
    ParOptVec *x_final;
    opt->getOptimizedPoint(&x_final, NULL, NULL, NULL, NULL);
    problem->setGlobalComponents(x_final, x_design);

    // Free the optimization problem - now you can't touch ParOpt again
    // until we reset the problem type at the next iteration
    delete problem;

    // Free all the analysis objects
    mg->decref();
  
    // Decrease the reference count to things
    for ( int k = 0; k < NUM_LEVELS; k++ ){
      tacs[k]->decref();
    }

    // Deallocate the interpolation and restriction operators
    for ( int k = 0; k < NUM_LEVELS-1; k++ ){
      interp[k]->decref();
      restrct[k]->decref();
    }

    // Determine the refinement based on the design variables alone.
    // Refine any element when it is contained within a region that
    // has one design variable exceeding a cut-off value
    if (mpi_rank == 0){
      double x_cutoff = 0.75;

      // Scan over the filter elements and label those that need to be
      // refiend
      int nelem_filter;
      const int *filter_ptr, *filter_conn;
      filters[0]->getMesh(NULL, &nelem_filter, &filter_ptr, &filter_conn);

      int *filter_refine = new int[ nelem_filter ];
      memset(filter_refine, 0, nelem_filter*sizeof(int));

      for ( int i = 0; i < nelem_filter; i++ ){
        for ( int jp = filter_ptr[i]; jp < filter_ptr[i+1]; jp++ ){
          int node = filter_conn[jp];
          
          // Label the element for refinement if any one of its nodes
          // is greater than the refinement cutoff
          if (node >= 0 && x_design[node] > x_cutoff){
            filter_refine[i] = 1;
          }
        }
      }

      // Get the array of filter octants
      TMROctantArray *elements;
      filters[0]->getElements(&elements);

      // Now get the array of filter octants from the element list
      TMROctant *filter_elems;
      elements->getArray(&filter_elems, NULL);

      // Create the refinement array for the analysis mesh
      int nelems = 0;
      root_tree->getMesh(NULL, &nelems, NULL, NULL);
      int *refine = new int[ nelems ];
      memset(refine, 0, nelems*sizeof(int));

      for ( int i = 0; i < nelem_filter; i++ ){
        // If the filter element has been labeled for refinement,
        // refine the enclosing range of elements within the mesh.
        if (filter_refine[i] > 0){
          int low, high;
          root_tree->findEnclosingRange(&filter_elems[i], &low, &high);
          for ( int k = low; k < high; k++ ){
            refine[k] = 1;
          }
        }
        else if (filter_refine[i] < 0){
          int low, high;
          root_tree->findEnclosingRange(&filter_elems[i], &low, &high);
          for ( int k = low; k < high; k++ ){
            refine[k] = -1;
          }
        }
      }

      // Refine the root tree
      root_tree->refine(refine, min_refine, max_refine);
    
      delete [] refine;
      delete [] filter_refine;
    }
  }

  delete opt;
  if (root_tree){ delete root_tree; } 

  TMRFinalize();
  MPI_Finalize();
  return (0);
}
