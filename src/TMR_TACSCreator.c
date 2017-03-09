#include "TMR_TACSCreator.h"

/*
  Initialize the TMRQuadTACSCreator data in the abstract base class
*/
TMRQuadTACSCreator::TMRQuadTACSCreator(){
  bc_root = NULL;
  bc_current = NULL;
}

/*
  Free the data - including the boundary condition information
*/
TMRQuadTACSCreator::~TMRQuadTACSCreator(){
  // Free the boundary condition information
  while (bc_root){
    bc_current = bc_root->next;
    delete bc_root;
    bc_root = bc_current;
  }
}

/*
  Constructor for the BCNode sub-class
*/
TMRQuadTACSCreator::BCNode::BCNode( const char *_attr,
                                    int _num_bcs, 
                                    const int *_bc_nums ){
  next = NULL;
  attr = new char[ strlen(_attr)+1 ];
  strcpy(attr, _attr);
  num_bcs = _num_bcs;
  bc_nums = new int[ num_bcs ];
  memcpy(bc_nums, _bc_nums, num_bcs*sizeof(int));
}

/*
  Destructor for the BCNode sub-class
*/
TMRQuadTACSCreator::BCNode::~BCNode(){
  delete [] attr;
  delete [] bc_nums;
}

/*
  Add the boundary conditions that will be associated with the
  specified attribute to the boundary condition linked list.
*/
void TMRQuadTACSCreator::addBoundaryCondition( const char *attr, 
                                               int num_bcs, 
                                               const int bc_nums[] ){
  BCNode *node = new BCNode(attr, num_bcs, bc_nums);
  if (!bc_root){
    bc_root = node;
    bc_current = node;
  }
  else {
    bc_current->next = node;
    bc_current = bc_current->next;
  }
}

/*
  Create the TACSAssembler object
*/
TACSAssembler* TMRQuadTACSCreator::createTACS( int order, 
                                               TMRQuadForest *forest ){
  // Get the communicator and the rank
  MPI_Comm comm = forest->getMPIComm();
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);
    
  // Ensure that the forest is balanced, and create the nodes for
  // the underlying finite-element mesh
  forest->balance();
  forest->createNodes(order);

  // Find the number of nodes for this processor
  const int *range;
  forest->getOwnedNodeRange(&range);
  int num_nodes = range[mpi_rank+1] - range[mpi_rank];

  // Create the mesh
  int *elem_conn, num_elements = 0;
  forest->createMeshConn(&elem_conn, &num_elements);

  // Create/retrieve the dependent node information
  const int *dep_ptr, *dep_conn;
  const double *dep_weights;
  int num_dep_nodes = forest->getDepNodeConn(&dep_ptr, &dep_conn,
                                             &dep_weights);

  // Get the quadrant array associated with this element
  TMRQuadrantArray *quadrants;
  forest->getQuadrants(&quadrants);
  int size;
  TMRQuadrant *array;
  quadrants->getArray(&array, &size);
  
  // Create the first element - and read out the number of
  // variables-per-node
  TACSElement *first = createElement(order, forest, array[0]);
  int vars_per_node = first->numDisplacements();

  // Create the associated TACSAssembler object
  TACSAssembler *tacs = 
    new TACSAssembler(forest->getMPIComm(), vars_per_node,
                      num_nodes, num_elements, num_dep_nodes);

  // Set the element ptr
  int *ptr = new int[ num_elements+1 ];
  for ( int i = 0; i < num_elements+1; i++ ){
    ptr[i] = order*order*i;
  }
    
  // Set the element connectivity into TACSAssembler
  tacs->setElementConnectivity(elem_conn, ptr);
  delete [] elem_conn;
  delete [] ptr;
    
  // Set the dependent node information
  tacs->setDependentNodes(dep_ptr, dep_conn, dep_weights);

    
  // Set the elements
  TACSAuxElements *aux = new TACSAuxElements(num_elements);
  TACSElement **elements = new TACSElement*[ num_elements ];
  for ( int k = 0; k < num_elements; k++ ){
    if (k == 0){
      elements[0] = first;
    }
    else {
      elements[k] = createElement(order, forest, array[k]);
    }

    // Create the auxiliary elements - if any
    TACSElement *elem = createAuxElement(order, forest, array[k]);
    if (elem){
      aux->addElement(k, elem);
    }
  }
    
  // Set the element array
  tacs->setElements(elements);
  delete [] elements;

  // Specify the boundary conditions
  setBoundaryConditions(forest, tacs);
    
  // Initialize the TACSAssembler object
  tacs->initialize();

  // Set the auxiliary elements
  tacs->setAuxElements(aux);

  // Set the node locations
  setNodeLocations(forest, tacs);

  return tacs;
}

/*
  Set the boundary condtiions into the TACSAssembler object

  This must be called before the TACSAssembler object is initialized
*/
void TMRQuadTACSCreator::setBoundaryConditions( TMRQuadForest *forest,
                                                TACSAssembler *tacs ){
  // Get the communicator and the rank
  MPI_Comm comm = forest->getMPIComm();
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);
    
  // Find the number of nodes for this processor
  const int *range;
  forest->getOwnedNodeRange(&range);

  BCNode *bc = bc_root;   
  while (bc){
    // Retrieve the nodes associated with the specified attribute
    TMRQuadrantArray *nodes = forest->getNodesWithAttribute(bc->attr);
    int size;
    TMRQuadrant *array;
    nodes->getArray(&array, &size);

    // Allocate the array of the node numbers associated with the BC
    int num = 0;
    int *vars = new int[ size ];
    for ( int i = 0; i < size; i++ ){
      if (array[i].tag >= range[mpi_rank] &&
          array[i].tag < range[mpi_rank+1]){
        vars[num] = array[i].tag;
        num++;
      }
    }

    // Add the boundary conditions to TACSAssembler
    tacs->addBCs(num, vars, bc->num_bcs, bc->bc_nums);
    delete [] vars;

    // Go to the next boundary condition
    bc = bc->next;
  }
}

/*
  Set the node locations from the TMRQuadForest into the TACSAssembler
  object
*/
void TMRQuadTACSCreator::setNodeLocations( TMRQuadForest *forest, 
                                           TACSAssembler *tacs ){
  // Get the communicator and the rank
  MPI_Comm comm = forest->getMPIComm();
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);
    
  // Find the number of nodes for this processor
  const int *range;
  forest->getOwnedNodeRange(&range);

  // Get the nodes
  TMRQuadrantArray *nodes;
  forest->getNodes(&nodes);

  // Get the quadrants associated with the nodes
  int size;
  TMRQuadrant *array;
  nodes->getArray(&array, &size);

  // Get the points
  TMRPoint *Xp;
  forest->getPoints(&Xp);

  TACSBVec *X = tacs->createNodeVec();
  X->incref();

  // Get the node array from the TACSBVec object
  TacsScalar *Xn;
  X->getArray(&Xn);

  // Loop over all the nodes
  for ( int i = 0; i < size; i++ ){
    if (array[i].tag >= range[mpi_rank] &&
        array[i].tag < range[mpi_rank+1]){
      int loc = array[i].tag - range[mpi_rank];
      Xn[3*loc] = Xp[i].x;
      Xn[3*loc+1] = Xp[i].y;
      Xn[3*loc+2] = Xp[i].z;
    }
  }

  tacs->setNodes(X);
  X->decref();
}

