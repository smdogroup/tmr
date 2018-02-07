#include "TMR_TACSCreator.h"

/*
  Allocate a new boundary condition set
*/
TMRBoundaryConditions::TMRBoundaryConditions(){
  num_bcs = 0;
  bc_root = NULL;
  bc_current = NULL;
}

/*
  Free the data - including the boundary condition information
*/
TMRBoundaryConditions::~TMRBoundaryConditions(){
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
TMRBoundaryConditions::BCNode::BCNode( const char *_attr,
                                       int _num_bcs, 
                                       const int *_bc_nums,
                                       const TacsScalar *_bc_vals ){
  next = NULL;
  attr = new char[ strlen(_attr)+1 ];
  strcpy(attr, _attr);
  num_bcs = _num_bcs;

  if (num_bcs > 0){
    bc_nums = new int[ num_bcs ];
    bc_vals = new TacsScalar[ num_bcs ];
    memcpy(bc_nums, _bc_nums, num_bcs*sizeof(int));
    if (_bc_vals){
      memcpy(bc_vals, _bc_vals,num_bcs*sizeof(TacsScalar));
    }
    else {
      memset(bc_vals, 0.0, num_bcs*sizeof(TacsScalar));
    }
  }
  else {
    num_bcs = -1;
    bc_nums = NULL;
    bc_vals = NULL;  
  }
}

/*
  Destructor for the BCNode sub-class
*/
TMRBoundaryConditions::BCNode::~BCNode(){
  delete [] attr;
  if (bc_nums){ delete [] bc_nums; }
  if (bc_vals){ delete [] bc_vals; }
}

/*
  Add the boundary conditions that will be associated with the
  specified attribute to the boundary condition linked list.
*/
void TMRBoundaryConditions::addBoundaryCondition( const char *attribute, 
                                                  int num_bc_nums,
                                                  const int bc_nums[],
                                                  const TacsScalar *bc_vals ){
  num_bcs++;
  BCNode *node = new BCNode(attribute, num_bc_nums, bc_nums, bc_vals);
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
  Get the number of boundary conditions
*/
int TMRBoundaryConditions::getNumBoundaryConditions(){
  return num_bcs;
}

/*
  Retrieve the boundary conditions
*/
void TMRBoundaryConditions::getBoundaryCondition( int bc, const char **_attr, 
                                                  int *_num_bcs,
                                                  const int **_bc_nums, 
                                                  const TacsScalar **_bc_vals ){
  *_attr = NULL;
  *_num_bcs = -1;
  *_bc_nums = NULL;
  *_bc_vals = NULL;
  
  // Increment until we've found the boundary condition or the end of
  // the linked list
  int count = 0;
  BCNode *node = bc_root;
  while (node && count < bc){
    node = node->next;
    count++;
  }

  // If the node exists, write out the node
  if (node){
    *_attr = node->attr;
    *_num_bcs = node->num_bcs;
    *_bc_nums = node->bc_nums;
    *_bc_vals = node->bc_vals;
  }
}

/*
  Initialize the TMRQuadTACSCreator data in the abstract base class
*/
TMRQuadTACSCreator::TMRQuadTACSCreator( TMRBoundaryConditions *_bcs ){
  bcs = _bcs;
  if (bcs){ bcs->incref(); }
}

/*
  Free the data - including the boundary condition information
*/
TMRQuadTACSCreator::~TMRQuadTACSCreator(){
  if (bcs){ bcs->decref(); }
}

/*
  Create the TACSAssembler object
*/
TACSAssembler* 
  TMRQuadTACSCreator::createTACS( TMRQuadForest *forest,
                                  TACSAssembler::OrderingType ordering,
                                  TacsScalar _scale ){
  // Get the communicator and the rank
  MPI_Comm comm = forest->getMPIComm();
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);
    
  // Create the nodes for the underlying finite-element mesh if they
  // don't yet exist
  forest->createNodes();

  // Get the element order
  int order = forest->getMeshOrder();

  // Get the local part of the connectivity
  const int *conn;
  int num_elements = 0, num_owned_nodes = 0;
  forest->getNodeConn(&conn, &num_elements, &num_owned_nodes);

  // Allocate the pointer array
  int *ptr = new int[ num_elements+1 ];
  for ( int i = 0; i <= num_elements; i++ ){
    ptr[i] = order*order*i;
  }

  // Create/retrieve the dependent node information
  const int *dep_ptr, *dep_conn;
  const double *dep_weights;
  int num_dep_nodes = forest->getDepNodeConn(&dep_ptr, &dep_conn,
                                             &dep_weights);

  // Create the elements using the virtual call
  TACSElement **elements = new TACSElement*[ num_elements ];
  createElements(order, forest, num_elements, elements);

  // Create the first element - and read out the number of
  // variables-per-node
  int vars_per_node = elements[0]->numDisplacements();

  // Create the associated TACSAssembler object
  TACSAssembler *tacs = 
    new TACSAssembler(forest->getMPIComm(), vars_per_node,
                      num_owned_nodes, num_elements, num_dep_nodes);

  // Set the element connectivity into TACSAssembler
  tacs->setElementConnectivity(conn, ptr);
  delete [] ptr;
    
  // Set the dependent node information
  tacs->setDependentNodes(dep_ptr, dep_conn, dep_weights);

  // Set the element array
  tacs->setElements(elements);
  delete [] elements;

  // Specify the boundary conditions
  setBoundaryConditions(forest, tacs);

  // Reordering everything - if needed
  if (ordering != TACSAssembler::NATURAL_ORDER){
    tacs->computeReordering(ordering, TACSAssembler::GAUSS_SEIDEL);
  }

  // Initialize the TACSAssembler object
  tacs->initialize();

  // Create the auxiliary elements
  TACSAuxElements *aux = createAuxElements(order, forest);
  if (aux){
    tacs->setAuxElements(aux);
  }

  // Set the node locations
  setNodeLocations(forest, tacs,_scale);

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

  if (bcs){
    for ( int k = 0; k < bcs->getNumBoundaryConditions(); k++ ){
      // Retrieve the boundary condition
      const char *attribute;
      int num_bcs;
      const int *bc_nums;
      const TacsScalar *bc_vals;
      bcs->getBoundaryCondition(k, &attribute, &num_bcs, 
                                &bc_nums, &bc_vals);

      if (attribute){
        // Retrieve the nodes associated with the specified attribute
        int *nodes;
        int num_nodes = forest->getNodesWithAttribute(attribute, &nodes);

        // Add the boundary conditions to TACSAssembler
        tacs->addBCs(num_nodes, nodes, num_bcs, bc_nums, bc_vals);

        delete [] nodes;
      }
    }
  }
}

/*
  Set the node locations from the TMRQuadForest into the TACSAssembler
  object
*/
void TMRQuadTACSCreator::setNodeLocations( TMRQuadForest *forest, 
                                           TACSAssembler *tacs,
                                           TacsScalar _scale ){
  TacsScalar scale = _scale;

  // Get the communicator and the rank
  int mpi_rank;
  MPI_Comm comm = forest->getMPIComm();
  MPI_Comm_rank(comm, &mpi_rank);
    
  // Find the number of nodes for this processor
  const int *range;
  forest->getOwnedNodeRange(&range);

  // Get the node numbers
  const int *nodes;
  int num_local_nodes = forest->getNodeNumbers(&nodes);

  // Get the points
  TMRPoint *Xp;
  forest->getPoints(&Xp);

  TACSBVec *X = tacs->createNodeVec();
  X->incref();

  // Get the node array from the TACSBVec object
  TacsScalar *Xn;
  X->getArray(&Xn);

  // Loop over all the nodes
  for ( int i = 0; i < num_local_nodes; i++ ){
    if (nodes[i] >= range[mpi_rank] && 
        nodes[i] < range[mpi_rank+1]){
      int loc = nodes[i] - range[mpi_rank];
      Xn[3*loc] = scale*Xp[i].x;
      Xn[3*loc+1] = scale*Xp[i].y;
      Xn[3*loc+2] = scale*Xp[i].z;
    }
  }
  
  // Reorder the vector if needed
  tacs->reorderVec(X);
  tacs->setNodes(X);
  X->decref();
}

/*
  Initialize the TMROctTACSCreator data in the abstract base class
*/
TMROctTACSCreator::TMROctTACSCreator( TMRBoundaryConditions *_bcs ){
  bcs = _bcs;
  if (bcs){ bcs->incref(); }
}

/*
  Free the data - including the boundary condition information
*/
TMROctTACSCreator::~TMROctTACSCreator(){
  if (bcs){ bcs->decref(); }
}

/*
  Create the TACSAssembler object
*/
TACSAssembler* 
  TMROctTACSCreator::createTACS( TMROctForest *forest,
                                 TACSAssembler::OrderingType ordering,
                                 TacsScalar _scale ){
  // Get the communicator and the rank
  MPI_Comm comm = forest->getMPIComm();
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);
    
  // Create the nodes for the underlying finite-element mesh if they
  // don't yet exist
  forest->createNodes();

  // Get the element order
  int order = forest->getMeshOrder();

  // Get the local part of the connectivity
  const int *conn;
  int num_elements = 0, num_owned_nodes = 0;
  forest->getNodeConn(&conn, &num_elements, &num_owned_nodes);

  // Allocate the pointer array
  int *ptr = new int[ num_elements+1 ];
  for ( int i = 0; i <= num_elements; i++ ){
    ptr[i] = order*order*order*i;
  }

  // Create/retrieve the dependent node information
  const int *dep_ptr, *dep_conn;
  const double *dep_weights;
  int num_dep_nodes = forest->getDepNodeConn(&dep_ptr, &dep_conn,
                                             &dep_weights);

  // Create the elements using the virtual call
  TACSElement **elements = new TACSElement*[ num_elements ];
  createElements(order, forest, num_elements, elements);

  // Create the first element - and read out the number of
  // variables-per-node
  int vars_per_node = elements[0]->numDisplacements();

  // Create the associated TACSAssembler object
  TACSAssembler *tacs = 
    new TACSAssembler(forest->getMPIComm(), vars_per_node,
                      num_owned_nodes, num_elements, num_dep_nodes);

  // Set the element connectivity into TACSAssembler
  tacs->setElementConnectivity(conn, ptr);
  delete [] ptr;
    
  // Set the dependent node information
  tacs->setDependentNodes(dep_ptr, dep_conn, dep_weights);

  // Set the element array
  tacs->setElements(elements);
  delete [] elements;

  // Specify the boundary conditions
  setBoundaryConditions(forest, tacs);

  // Reordering everything - if needed
  if (ordering != TACSAssembler::NATURAL_ORDER){
    tacs->computeReordering(ordering, TACSAssembler::GAUSS_SEIDEL);
  }

  // Initialize the TACSAssembler object
  tacs->initialize();

  // Create the auxiliary elements
  TACSAuxElements *aux = createAuxElements(order, forest);
  if (aux){
    tacs->setAuxElements(aux);
  }

  // Set the node locations
  setNodeLocations(forest, tacs,_scale);

  return tacs;
}

/*
  Set the boundary condtiions into the TACSAssembler object

  This must be called before the TACSAssembler object is initialized
*/
void TMROctTACSCreator::setBoundaryConditions( TMROctForest *forest,
                                               TACSAssembler *tacs ){
  // Get the communicator and the rank
  MPI_Comm comm = forest->getMPIComm();
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);
    
  // Find the number of nodes for this processor
  const int *range;
  forest->getOwnedNodeRange(&range);

  if (bcs){
    for ( int k = 0; k < bcs->getNumBoundaryConditions(); k++ ){
      // Retrieve the boundary condition
      const char *attribute;
      int num_bcs;
      const int *bc_nums;
      const TacsScalar *bc_vals;
      bcs->getBoundaryCondition(k, &attribute, &num_bcs, 
                                &bc_nums, &bc_vals);

      if (attribute){
        // Retrieve the nodes associated with the specified attribute
        int *nodes;
        int num_nodes = forest->getNodesWithAttribute(attribute, &nodes);

        // Add the boundary conditions to TACSAssembler
        tacs->addBCs(num_nodes, nodes, num_bcs, bc_nums, bc_vals);

        delete [] nodes;
      }
    }
  }
}

/*
  Set the node locations from the TMROctForest into the TACSAssembler
  object
*/
void TMROctTACSCreator::setNodeLocations( TMROctForest *forest, 
                                          TACSAssembler *tacs,
                                          TacsScalar _scale ){
  TacsScalar scale = _scale;

  // Get the communicator and the rank
  int mpi_rank;
  MPI_Comm comm = forest->getMPIComm();
  MPI_Comm_rank(comm, &mpi_rank);
    
  // Find the number of nodes for this processor
  const int *range;
  forest->getOwnedNodeRange(&range);

  // Get the node numbers
  const int *nodes;
  int num_local_nodes = forest->getNodeNumbers(&nodes);

  // Get the points
  TMRPoint *Xp;
  forest->getPoints(&Xp);

  TACSBVec *X = tacs->createNodeVec();
  X->incref();

  // Get the node array from the TACSBVec object
  TacsScalar *Xn;
  X->getArray(&Xn);

  // Loop over all the nodes
  for ( int i = 0; i < num_local_nodes; i++ ){
    if (nodes[i] >= range[mpi_rank] && 
        nodes[i] < range[mpi_rank+1]){
      int loc = nodes[i] - range[mpi_rank];
      Xn[3*loc] = scale*Xp[i].x;
      Xn[3*loc+1] = scale*Xp[i].y;
      Xn[3*loc+2] = scale*Xp[i].z;
    }
  }
  
  // Reorder the vector if needed
  tacs->reorderVec(X);
  tacs->setNodes(X);
  X->decref();
}

