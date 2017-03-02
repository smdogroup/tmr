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
    printf("Successful load\n");

    // Get the faces that have been created - if any
    // and write them all to different VTK files
    int num_faces;
    TMRFace **faces;
    geo->getFaces(&num_faces, &faces);

    // Allocate the new mesh
    double htarget = 5.0;

    TMRMesh *mesh = new TMRMesh(geo);
    mesh->incref();
    mesh->mesh(htarget);
    
    TMRFaceMesh *fmesh;
    faces[0]->getMesh(&fmesh);
    fmesh->writeToVTK("face_mesh.vtk");

    /*
    TMRFaceMesh *fmesh = new TMRFaceMesh(faces[0]);
    fmesh->incref();

    // Mesh the geometry with the target spacing
    fmesh->mesh(htarget);
    fmesh->writeToVTK("face_mesh.vtk");

    fmesh->decref();
    */

    geo->decref();
  }
  
  TMRFinalize();
  MPI_Finalize();

  return 0;
}

#endif // TMR_HAS_OPENCASCADE
