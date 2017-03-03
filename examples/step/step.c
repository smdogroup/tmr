#ifdef TMR_HAS_OPENCASCADE

#include "TMROpenCascade.h"
#include "TMRMesh.h"

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  const char *filename = "misc1.step";

  // Load in the geometry file
  TMRModel *geo = TMR_LoadModelFromSTEPFile(filename);
  if (geo){
    geo->incref();

    // Get the faces that have been created - if any
    // and write them all to different VTK files
    int num_faces;
    TMRFace **faces;
    geo->getFaces(&num_faces, &faces);

    // Allocate the new mesh
    double htarget = 2.0;

    TMRMesh *mesh = new TMRMesh(geo);
    mesh->incref();
    mesh->mesh(htarget);
    
    for ( int k = 0; k < num_faces; k++ ){
      TMRFaceMesh *fmesh;
      faces[k]->getMesh(&fmesh);

      char fname[128];
      sprintf(fname, "face_mesh%02d.vtk", k);
      fmesh->writeToVTK(fname);
    }

    geo->decref();
  }
  
  TMRFinalize();
  MPI_Finalize();

  return 0;
}

#endif // TMR_HAS_OPENCASCADE
