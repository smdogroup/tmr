#ifdef TMR_HAS_OPENCASCADE

#include "MITCShell.h"
#include "TACSMeshLoader.h"
#include "TMRMesh.h"
#include "TMROpenCascade.h"
#include "isoFSDTStiffness.h"

// Include the netgen library

namespace nglib {
#include "nglib.h"
}

using namespace nglib;

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  TMRInitialize();

  char filename[256];
  snprintf(filename, sizeof(filename), "misc1.step");

  double htarget = 4.0;
  int surface_only = 0;
  for (int k = 0; k < argc; k++) {
    if (sscanf(argv[k], "h=%lf", &htarget) == 1) {
      if (htarget < 0.1) {
        htarget = 0.1;
      }
    }
    if (sscanf(argv[k], "file=%s", filename) == 1) {
      printf("file=%s\n", filename);
    }
    if (strcmp(argv[k], "surface") == 0) {
      surface_only = 1;
    }
  }

  // Load in the geometry file
  TMRModel *geo = TMR_LoadModelFromSTEPFile(filename);
  if (geo) {
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

    TMRModel *model =
        new TMRModel(num_verts, verts, num_edges, edges, num_faces, faces);

    // Allocate the new mesh
    TMRMesh *mesh = new TMRMesh(MPI_COMM_WORLD, model);
    mesh->incref();

    // Adjust the quality factor
    TMRMeshOptions options;
    options.mesh_type_default = TMR_TRIANGLE;
    options.frontal_quality_factor = 1.25;
    options.num_smoothing_steps = 10;
    options.write_mesh_quality_histogram = 1;
    options.triangularize_print_level = 1;

    // Create the feature size
    TMRPoint a, b;
    a.x = a.y = a.z = -1e3;
    b.x = b.y = b.z = 1e3;
    TMRBoxFeatureSize *fs =
        new TMRBoxFeatureSize(a, b, 0.01 * htarget, htarget);

    // Set the default box size
    const int NUM_BOXES = 280;
    TMRPoint p1[NUM_BOXES], p2[NUM_BOXES];

    double p[3] = {65.0, -150.0, -175.0};
    double d[3] = {220.0, 150.0, 260.0};

    // Create random boxes
    for (int i = 0; i < NUM_BOXES; i++) {
      p1[i].x = p[0] + (d[0] * rand()) / RAND_MAX;
      p1[i].y = p[1] + (d[1] * rand()) / RAND_MAX;
      p1[i].z = p[2] + (d[2] * rand()) / RAND_MAX;

      p2[i].x = p1[i].x + 25.0;
      p2[i].y = p1[i].y + 25.0;
      p2[i].z = p1[i].z + 25.0;

      fs->addBox(p1[i], p2[i], 0.5 * htarget);
    }

    // Print out the refinement volumes to a VTK file
    FILE *rfp = fopen("volume-refine.vtk", "w");
    if (rfp) {
      fprintf(rfp, "# vtk DataFile Version 3.0\n");
      fprintf(rfp, "vtk output\nASCII\n");
      fprintf(rfp, "DATASET UNSTRUCTURED_GRID\n");

      // Write out the points
      fprintf(rfp, "POINTS %d float\n", 8 * NUM_BOXES);
      for (int kk = 0; kk < NUM_BOXES; kk++) {
        for (int k = 0; k < 2; k++) {
          for (int j = 0; j < 2; j++) {
            for (int i = 0; i < 2; i++) {
              double x = p1[kk].x + i * (p2[kk].x - p1[kk].x);
              double y = p1[kk].y + j * (p2[kk].y - p1[kk].y);
              double z = p1[kk].z + k * (p2[kk].z - p1[kk].z);
              fprintf(rfp, "%e %e %e\n", x, y, z);
            }
          }
        }
      }

      fprintf(rfp, "\nCELLS %d %d\n", NUM_BOXES, 9 * NUM_BOXES);

      // Write out the cell connectivities
      for (int k = 0; k < NUM_BOXES; k++) {
        fprintf(rfp, "8 %d %d %d %d %d %d %d %d\n", 8 * k, 8 * k + 1, 8 * k + 3,
                8 * k + 2, 8 * k + 4, 8 * k + 5, 8 * k + 7, 8 * k + 6);
      }

      // All quadrilaterals
      fprintf(rfp, "\nCELL_TYPES %d\n", NUM_BOXES);
      for (int k = 0; k < NUM_BOXES; k++) {
        fprintf(rfp, "%d\n", 12);
      }

      fclose(rfp);
    }

    // Mesh the object of interest
    mesh->mesh(options, fs);
    mesh->writeToVTK("surface-mesh.vtk");

    if (!surface_only) {
      // Get the triangle connectivity
      int ntris;
      const int *tris;
      mesh->getTriConnectivity(&ntris, &tris);

      // Get the mesh points
      TMRPoint *X;
      int npts = mesh->getMeshPoints(&X);

      // Allocate the new mesh
      Ng_Mesh *m = Ng_NewMesh();

      // Count up the number of boundary nodes for this volume
      for (int i = 0; i < npts; i++) {
        double pt[3];
        pt[0] = X[i].x;
        pt[1] = X[i].y;
        pt[2] = X[i].z;
        Ng_AddPoint(m, pt);
      }

      // Add the surface elements
      for (int i = 0; i < ntris; i++) {
        int tri[3];
        tri[0] = tris[3 * i] + 1;
        tri[1] = tris[3 * i + 1] + 1;
        tri[2] = tris[3 * i + 2] + 1;
        Ng_AddSurfaceElement(m, NG_TRIG, tri);
      }

      // Set the meshing parameters
      Ng_RestrictMeshSizeGlobal(m, htarget);
      for (int i = 0; i < NUM_BOXES; i++) {
        double pt1[3] = {p1[i].x, p1[i].y, p1[i].z};
        double pt2[3] = {p2[i].x, p2[i].y, p2[i].z};
        Ng_RestrictMeshSizeBox(m, pt1, pt2, 0.5 * htarget);
      }

      // generate volume mesh
      Ng_Meshing_Parameters mp;
      mp.maxh = htarget;
      mp.fineness = 1;
      mp.second_order = 0;

      // Generate the volume mesh
      Ng_GenerateVolumeMesh(m, &mp);

      // Get the total number of points and tets in the volume
      int num_points = Ng_GetNP(m);
      int num_tet = Ng_GetNE(m);

      // Allocate space to store everything
      TMRPoint *Xt = new TMRPoint[num_points];
      int *tet = new int[4 * num_tet];

      // Retrieve the points
      for (int i = 0; i < num_points; i++) {
        double x[3];
        Ng_GetPoint(m, i + 1, x);
        Xt[i].x = x[0];
        Xt[i].y = x[1];
        Xt[i].z = x[2];
      }

      // Retrieve the tets
      for (int i = 0; i < num_tet; i++) {
        Ng_GetVolumeElement(m, i + 1, &tet[4 * i]);
        for (int k = 0; k < 4; k++) {
          tet[4 * i + k] -= 1;
        }
      }

      // Free the memory and exit from netgen
      Ng_DeleteMesh(m);
      Ng_Exit();

      // Print out the volume mesh to a VTK file
      FILE *fp = fopen("volume-mesh.vtk", "w");
      if (fp) {
        fprintf(fp, "# vtk DataFile Version 3.0\n");
        fprintf(fp, "vtk output\nASCII\n");
        fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

        // Write out the points
        fprintf(fp, "POINTS %d float\n", num_points);
        for (int k = 0; k < num_points; k++) {
          fprintf(fp, "%e %e %e\n", Xt[k].x, Xt[k].y, Xt[k].z);
        }

        fprintf(fp, "\nCELLS %d %d\n", num_tet, 5 * num_tet);

        // Write out the cell connectivities
        for (int k = 0; k < num_tet; k++) {
          fprintf(fp, "4 %d %d %d %d\n", tet[4 * k], tet[4 * k + 1],
                  tet[4 * k + 2], tet[4 * k + 3]);
        }

        // All quadrilaterals
        fprintf(fp, "\nCELL_TYPES %d\n", num_tet);
        for (int k = 0; k < num_tet; k++) {
          fprintf(fp, "%d\n", 10);
        }
        fclose(fp);
      }
    }

    // Decref the objects
    mesh->decref();
    geo->decref();
  }

  TMRFinalize();
  MPI_Finalize();

  return 0;
}

#endif  // TMR_HAS_OPENCASCADE
