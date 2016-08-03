/*
  The following code is used to generate geometry for a problem
*/





int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();


  





  MPI_Comm comm = MPI_COMM_WORLD;
  TMROctForest *forest = new TMROctForest(comm);

  // Create the TACSMeshLoader class
  TACSMeshLoader *mesh = new TACSMeshLoader(MPI_COMM_SELF);
  mesh->incref();
  mesh->scanBDFFile("uCRM_3D_box_mesh.bdf");

  // Extract the connectivity
  int nnodes, nelems;
  const int *elem_node_conn;
  const double *Xpts;
  mesh->getConnectivity(&nnodes, &nelems, NULL, 
                        &elem_node_conn, &Xpts);
  forest->setConnectivity(nnodes, elem_node_conn, nelems);


  // Create the random trees
  forest->createRandomTrees(25, 0, 15);

  double tbal = MPI_Wtime();
  forest->balance(1);
  tbal = MPI_Wtime() - tbal;

  double tnodes = MPI_Wtime();
  forest->createNodes(2);
  tnodes = MPI_Wtime() - tnodes;

  // Get the octrees within the forest
  TMROctree **trees;
  int ntrees = forest->getOctrees(&trees);

  





  TMRFinalize();
  MPI_Finalize();
  return (0);
}

