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

    // Get the vertices
    int num_verts;
    TMRVertex **verts;
    geo->getVertices(&num_verts, &verts);

    // Get the edges
    int num_edges;
    TMREdge **edges;
    geo->getEdges(&num_edges, &edges);

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
    upper->setMaster(lower);

    // Create the volumes
    int num_vols = 1;
    TMRVolume* vol = new TMRVolume(num_faces, faces);

    TMRModel *model = new TMRModel(num_verts, verts,
                                   num_edges, edges, 
                                   num_faces, faces, 
                                   num_vols, &vol);
    model->incref();

    TMRMesh *mesh = new TMRMesh(comm, model);
    mesh->incref();

    TMRMeshOptions options;
    mesh->mesh(options, htarget);

    // Write the volume mesh
    mesh->writeToVTK("volume-mesh.vtk");
    mesh->writeToBDF("volume-mesh.bdf");

    mesh->decref();
    model->decref();
    geo->decref();
  }
  
  TMRFinalize();
  MPI_Finalize();

  return 0;
}

#endif // TMR_HAS_OPENCASCADE
