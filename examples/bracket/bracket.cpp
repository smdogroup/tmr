#ifdef TMR_HAS_OPENCASCADE

#include "TMROpenCascade.h"
#include "TMRMesh.h"
#include "TACSMeshLoader.h"
#include "Solid.h"
#include "TMROctForest.h"
#include "TMR_TACSTopoCreator.h"
#include "TMR_STLTools.h"

class CreateMe : public TMROctTACSTopoCreator {
 public:
  CreateMe( TMRBoundaryConditions *_bcs, TMROctForest *_filter,
            TMRStiffnessProperties *_props ):
    TMROctTACSTopoCreator(_bcs, _filter){
    props = _props;
  }

  TACSElement *createElement( int order, TMROctant *oct,
                              TMRIndexWeight *weights, int nweights ){
    double q = 5.0;
    TMROctStiffness *stiff = new TMROctStiffness(weights, nweights,
                                                 props, q);
    if (order == 2){
      return new Solid<2>(stiff);
    }
    else if (order == 3){
      return new Solid<3>(stiff);
    }
    else if (order == 4){
      return new Solid<4>(stiff);
    }
    else if (order == 5){
      return new Solid<5>(stiff);
    }    
    return NULL;
  }

 private:
  TMRStiffnessProperties *props;
};

/*
  Test the STL output generator using the bracket example
*/
void test_stl_output( const char *filename, TMROctForest *forest ){
  // Create the forest and balance it
  TMROctForest *filter = forest->coarsen();
  filter->incref();
  filter->balance();
  filter->repartition();

  // Create an empty set of boundary conditions and default material
  // properties
  TMRBoundaryConditions *bcs = new TMRBoundaryConditions();

  // Create the stiffness properties object
  TacsScalar rho = 2700.0;
  TacsScalar E = 70e9;
  TacsScalar nu = 0.3;
  TMRStiffnessProperties *props =
    new TMRStiffnessProperties(1, &rho, &E, &nu);
  
  // Allocate a creator object
  CreateMe *creator = new CreateMe(bcs, filter, props);
  creator->incref();
  
  // Create TACS
  TACSAssembler *tacs = creator->createTACS(forest);
  tacs->incref();

  // Get the underlying variable objects
  TACSVarMap *var_map;
  TACSBVecIndices *indices;
  creator->getMap(&var_map);
  creator->getIndices(&indices);
  TACSBVecDistribute *dist = new TACSBVecDistribute(var_map, indices);
  
  // Create the design vector
  TACSBVec *vars = new TACSBVec(var_map, 1, dist);
  vars->incref();

  // Get the range of the nodes
  const int *range;
  filter->getOwnedNodeRange(&range);

  // Get the communicator rank
  int mpi_rank;
  MPI_Comm_rank(forest->getMPIComm(), &mpi_rank);
  
  // Get the filter points
  TMRPoint *X;
  filter->getPoints(&X);

  // Get the connectivity
  int num_elements;
  const int *conn;
  filter->getNodeConn(&conn, &num_elements);

  // Get the mesh order
  const int order = filter->getMeshOrder(); 
  const int size = order*order*order*num_elements;

  // Get the elements of the array
  TacsScalar *x;
  vars->getArray(&x);
  
  for ( int i = 0; i < size; i++ ){
    if (conn[i] >= range[mpi_rank] &&
        conn[i] < range[mpi_rank+1]){
      int var = conn[i] - range[mpi_rank];
      int index = filter->getLocalNodeNumber(conn[i]);
      double xval = X[index].x - 30.0;
      double yval = X[index].y - 40.0;
      double zval = X[index].z;
      double d = sqrt(xval*xval + yval*yval + 4.0*zval*zval);
      x[var] = cos(0.5*d);
    }
  }

  double cutoff = 0.5;
  int var_offset = 0;
  TMR_GenerateBinFile(filename, filter, vars, var_offset, cutoff);

  // Free the objects
  vars->decref();
  filter->decref();
  creator->decref();
  tacs->decref();
}

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  // Don't write anything to a file, unless a flag is 
  // set on the command line
  int write_faces_to_vtk = 0;
  int test_bdf_file = 1;
  int test_stl_file = 1;
  for ( int k = 0; k < argc; k++ ){
    if (strcmp(argv[k], "--write_faces") == 0){
      write_faces_to_vtk = 1;
    }
    if (strcmp(argv[k], "--test_bdf") == 0){
      test_bdf_file = 1;
    }
    if (strcmp(argv[k], "--test_stl") == 0){
      test_stl_file = 1;
    }
  }

  // Get the communicator
  MPI_Comm comm = MPI_COMM_WORLD;

  // This is all tied to this STEP file
  const char *filename = "bracket_solid.stp";
  double htarget = 10.0;

  // Load in the geometry file
  TMRModel *geo = TMR_LoadModelFromSTEPFile(filename);
  if (geo){
    geo->incref();

    // Get the volume
    TMRVolume **volume;
    geo->getVolumes(NULL, &volume);

    // Get the faces from the volume
    int num_faces;
    TMRFace **faces;
    const int *dir;
    volume[0]->getFaces(&num_faces, &faces, &dir);

    // Write the surfaces files out, if needed
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0 && write_faces_to_vtk){
      for ( int i = 0; i < num_faces; i++ ){
        char output[128];
        sprintf(output, "faces%d.vtk", i);
        faces[i]->writeToVTK(output);
      }
    }

    // Set the upper/lower face numbers. These are based on the
    // ordering in the STEP file. This will hopefully be preserved
    // independent of how the STEP file is loaded.
    int lower_face_num = 1;
    int upper_face_num = 4;
    TMRFace *lower = faces[lower_face_num];
    TMRFace *upper = faces[upper_face_num];
    upper->setSource(volume[0], lower);

    TMRMesh *mesh = new TMRMesh(comm, geo);
    mesh->incref();

    // Mesh the geometry
    TMRMeshOptions options;
    options.num_smoothing_steps = 5;
    mesh->mesh(options, htarget);

    TMRModel *model = mesh->createModelFromMesh();
    model->incref();

    // Create the topology object from the geo-mesh
    TMRTopology *topo = new TMRTopology(comm, model);
    topo->incref();
    
    // Set the topology
    TMROctForest *forest = new TMROctForest(comm);
    forest->incref();

    // Create the random trees
    forest->setTopology(topo);
    forest->createRandomTrees(10, 0, 4);
    forest->balance();
    forest->createNodes();

    // Test the output file
    if (test_stl_file){
      test_stl_output("level_set_test.bstl", forest);

      char fname[128];
      sprintf(fname, "full_forest%d.vtk", rank);
      forest->writeForestToVTK(fname);
    }

    // Write the volume mesh
    if (test_bdf_file){
      mesh->writeToBDF("volume-mesh.bdf", TMRMesh::TMR_HEX);
    }

    topo->decref();
    forest->decref();
    model->decref();
    mesh->decref();
    geo->decref();
  }

  if (test_bdf_file){
    TACSMeshLoader *loader = new TACSMeshLoader(comm);
    loader->incref();

    loader->scanBDFFile("volume-mesh.bdf");

    // Create the solid stiffness object
    SolidStiffness *stiff = new SolidStiffness(1.0, 1.0, 0.3);
    Solid<2> *elem = new Solid<2>(stiff);
    loader->setElement(0, elem);
    
    // Create the TACSAssembler object
    int vars_per_node = 3;
    TACSAssembler *tacs = loader->createTACS(vars_per_node);
    tacs->incref();

    // Create the f5 visualization object
    TACSToFH5 *f5 = loader->createTACSToFH5(tacs, TACS_SOLID,
                                            TACSElement::OUTPUT_NODES);
    f5->incref();
    f5->writeToFile("volume-mesh.f5");
  }

  TMRFinalize();
  MPI_Finalize();

  return 0;
}

#endif // TMR_HAS_OPENCASCADE
