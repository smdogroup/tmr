#ifdef TMR_HAS_OPENCASCADE

#include "TMROpenCascade.h"
#include "TMRMesh.h"

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  // Don't write anything to a file, unless a flag is 
  // set on the command line
  int write_faces_to_vtk = 0;
  for ( int k = 0; k < argc; k++ ){
    if (strcmp(argv[k], "--write_faces") == 0){
      write_faces_to_vtk = 1;
    }
  }

  // Get the communicator
  MPI_Comm comm = MPI_COMM_WORLD;

  // This is all tied to this STEP file
  const char *filename = "bracket_solid.stp";
  double htarget = 4.0;

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
    upper->setMaster(lower);

    TMRMesh *mesh = new TMRMesh(comm, geo);
    mesh->incref();

    TMRMeshOptions options;
    mesh->mesh(options, htarget);

    // Write the volume mesh
    mesh->writeToVTK("volume-mesh.vtk");
    mesh->writeToBDF("volume-mesh.bdf");

    // TMRVolumeMesh* volmesh;
    // volume[0]->getMesh(&volmesh);
    // volmesh->writeToVTK("volume-mesh.vtk");

    mesh->decref();
    geo->decref();
  }
  
  TMRFinalize();
  MPI_Finalize();

  return 0;
}

#endif // TMR_HAS_OPENCASCADE
