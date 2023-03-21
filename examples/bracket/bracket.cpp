#ifdef TMR_HAS_OPENCASCADE

#include "TACSMeshLoader.h"
#include "TMRMesh.h"
#include "TMROctForest.h"
#include "TMROpenCascade.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  TMRInitialize();

  // Don't write anything to a file, unless a flag is
  // set on the command line
  int write_faces_to_vtk = 0;
  int test_surf_mesh = 0;
  int order = 2;
  for (int k = 0; k < argc; k++) {
    if (strcmp(argv[k], "--write_faces") == 0) {
      write_faces_to_vtk = 1;
    }
    if (strcmp(argv[k], "--test_surf_mesh") == 0) {
      test_surf_mesh = 1;
    }
    if (sscanf(argv[k], "order=%d", &order) == 1) {
      if (order < 2) {
        order = 2;
      }
      if (order > 5) {
        order = 5;
      }
    }
  }

  // Get the communicator
  MPI_Comm comm = MPI_COMM_WORLD;

  // This is all tied to this STEP file
  const char *filename = "bracket_solid.stp";
  double htarget = 8.0;

  // Load in the geometry file
  TMRModel *geo = TMR_LoadModelFromSTEPFile(filename);
  if (geo) {
    geo->incref();

    // Get the volume
    TMRVolume **volume;
    geo->getVolumes(NULL, &volume);

    // Get the faces from the volume
    int num_faces;
    TMRFace **faces;
    volume[0]->getFaces(&num_faces, &faces);

    // Write the surfaces files out, if needed
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0 && write_faces_to_vtk) {
      for (int i = 0; i < num_faces; i++) {
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

    if (test_surf_mesh) {
      int num_verts, num_edges;
      TMRVertex **verts;
      TMREdge **edges;
      geo->getVertices(&num_verts, &verts);
      geo->getEdges(&num_edges, &edges);
      geo->getFaces(&num_faces, &faces);

      TMRModel *geo_surf =
          new TMRModel(num_verts, verts, num_edges, edges, num_faces, faces);
      geo_surf->incref();

      TMRMesh *mesh_surf = new TMRMesh(comm, geo_surf);
      mesh_surf->incref();

      // Mesh the geometry
      TMRMeshOptions opts;
      opts.num_smoothing_steps = 5;
      mesh_surf->mesh(opts, htarget);

      mesh_surf->writeToVTK("surface-mesh.vtk");
      mesh_surf->decref();
      geo_surf->decref();
    }

    TMRMesh *mesh = new TMRMesh(comm, geo);
    mesh->incref();

    // Mesh the geometry
    TMRMeshOptions options;
    options.num_smoothing_steps = 5;
    mesh->mesh(options, htarget);

    // Write the original volume mesh to the file
    mesh->writeToVTK("volume-mesh.vtk", TMRMesh::TMR_HEX);

    // Construct the model from the mesh
    TMRModel *model = mesh->createModelFromMesh();
    model->incref();

    // Create the topology object from the geo-mesh
    TMRTopology *topo = new TMRTopology(comm, model);
    topo->incref();

    // Set the topology
    TMROctForest *forest = new TMROctForest(comm);
    forest->incref();

    // Set the mesh order
    forest->setMeshOrder(order, TMR_GAUSS_LOBATTO_POINTS);

    // Create the random trees
    forest->setTopology(topo);
    forest->createRandomTrees(10, 0, 5);
    forest->balance();
    forest->createNodes();

    // Write the volume mesh
    mesh->writeToBDF("volume-mesh.bdf", TMRMesh::TMR_HEX);

    topo->decref();
    forest->decref();
    model->decref();
    mesh->decref();
    geo->decref();
  }

  TMRFinalize();
  MPI_Finalize();

  return 0;
}

#endif  // TMR_HAS_OPENCASCADE
