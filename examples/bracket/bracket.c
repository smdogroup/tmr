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
  const char *filename = "bracket_shell.stp";
  double htarget = 4.0;

  // Load in the geometry file
  TMRModel *geo = TMR_LoadModelFromSTEPFile(filename);
  if (geo){
    geo->incref();

    // Now, separate and plot the different surfaces
    int num_faces;
    TMRFace **faces;
    geo->getFaces(&num_faces, &faces);

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
    int lower_face_num = 4;
    int upper_face_num = 1;
    TMRFace *lower = faces[lower_face_num];
    TMRFace *upper = faces[upper_face_num];

    int topo_equivalent = 1;

    // Check that the topology of the two surfaces are equivalent
    if (lower[i]->getNumEdgeLoops() !=
        upper[i]->getNumEdgeLoops()){
      topo_equivalent = 0;
    }


    // Mesh the lower surface
    




    geo->decref();
  }
  
  TMRFinalize();
  MPI_Finalize();

  return 0;
}

#endif // TMR_HAS_OPENCASCADE
