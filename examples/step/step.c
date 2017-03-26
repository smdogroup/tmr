#ifdef TMR_HAS_OPENCASCADE

#include "TMROpenCascade.h"
#include "TMRMesh.h"

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  // Get the communicator
  MPI_Comm comm = MPI_COMM_WORLD;

  char filename[256];
  sprintf(filename, "misc1.step");

  double htarget = 4.0;
  for ( int k = 0; k < argc; k++ ){
    if (sscanf(argv[k], "h=%lf", &htarget) == 1){
      if (htarget < 0.1){ htarget = 0.1; }
    }
    if (sscanf(argv[k], "file=%s", filename) == 1){
      printf("file=%s\n", filename);
    }
  }

  // Load in the geometry file
  TMRModel *geo = TMR_LoadModelFromSTEPFile(filename);
  if (geo){
    geo->incref();

    // Get the vertices
    int num_verts;
    TMRVertex **verts;
    geo->getVertices(&num_verts, &verts);

    // Get the edges
    int num_edges;
    TMREdge **edges;
    geo->getEdges(&num_edges, &edges);

    // Separate and plot the different surfaces
    int num_faces;
    TMRFace **faces;
    geo->getFaces(&num_faces, &faces);

    TMRModel *model = new TMRModel(num_verts, verts,
                                   num_edges, edges, 
                                   num_faces, faces);

    // Allocate the new mesh
    TMRMesh *mesh = new TMRMesh(MPI_COMM_WORLD, model);
    mesh->incref();

    // Adjust the quality factor
    TMRMeshOptions options;
    options.frontal_quality_factor = 1.5;
    options.num_smoothing_steps = 20;

    // Mesh the object of interest
    mesh->mesh(options, htarget);
    mesh->writeToVTK("surface-mesh.vtk");
    mesh->writeToBDF("surface-mesh.bdf");

    // Decref the objects
    mesh->decref();
    geo->decref();
  }
  
  TMRFinalize();
  MPI_Finalize();

  return 0;
}

#endif // TMR_HAS_OPENCASCADE
