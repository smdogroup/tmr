#ifdef TMR_HAS_OPENCASCADE

#include "TMROpenCascade.h"
#include "TMRMesh.h"

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  char filename[256];
  sprintf(filename, "misc1.step");

  double htarget = 4.0;
  for ( int k = 0; k < argc; k++ ){
    if (sscanf(argv[k], "h=%lf", &htarget) == 1){
      if (htarget < 0.1){ htarget = 0.1; }
      if (htarget > 10.0){ htarget = 10.0; }
    }
    if (sscanf(argv[k], "file=%s", filename) == 1){
      printf("file=%s\n", filename);
    }
  }

  // Load in the geometry file
  TMRModel *geo = TMR_LoadModelFromSTEPFile(filename);
  if (geo){
    geo->incref();

    // Allocate the new mesh
    TMRMesh *mesh = new TMRMesh(MPI_COMM_WORLD, geo);
    mesh->incref();

    // Adjust the quality factor
    TMRMeshOptions options;
    options.frontal_quality_factor = 1.25;

    // Mesh the object of interest
    mesh->mesh(options, htarget);
    mesh->writeToVTK("surface-mesh.vtk");

    mesh->decref();
    geo->decref();
  }
  
  TMRFinalize();
  MPI_Finalize();

  return 0;
}

#endif // TMR_HAS_OPENCASCADE
