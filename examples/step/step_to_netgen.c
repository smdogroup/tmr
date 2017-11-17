#ifdef TMR_HAS_OPENCASCADE

#include "TMROpenCascade.h"
#include "TMRMesh.h"
#include "TACSMeshLoader.h"
#include "MITCShell.h"
#include "isoFSDTStiffness.h"

// Include the netgen library
namespace nglib {
#include "nglib.h"
}

using namespace nglib;

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

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
    options.mesh_type_default = TMR_TRIANGLE;
    options.frontal_quality_factor = 1.25;
    options.num_smoothing_steps = 10;
    options.write_mesh_quality_histogram = 1;
    options.triangularize_print_level = 1;

    // Mesh the object of interest
    mesh->mesh(options, htarget);
    mesh->writeToVTK("surface-mesh.vtk");
    
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
    for ( int i = 0; i < npts; i++ ){
      double pt[3];
      pt[0] = X[i].x;
      pt[1] = X[i].y;
      pt[2] = X[i].z;
      Ng_AddPoint(m, pt);
    }

    // Add the surface elements
    for ( int i = 0; i < ntris; i++ ){
      int tri[3];
      tri[0] = tris[3*i]+1;
      tri[1] = tris[3*i+1]+1;
      tri[2] = tris[3*i+2]+1;
      Ng_AddSurfaceElement(m, NG_TRIG, tri);
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
    TMRPoint *Xt = new TMRPoint[ num_points ];
    int *tet = new int[ 4*num_tet ];

    // Retrieve the points
    for ( int i = 0; i < num_points; i++ ){
      double x[3];
      Ng_GetPoint(m, i+1, x);
      Xt[i].x = x[0];
      Xt[i].y = x[1];
      Xt[i].z = x[2];
    }

    // Retrieve the tets
    for ( int i = 0; i < num_tet; i++ ){
      Ng_GetVolumeElement(m, i+1, &tet[4*i]);
      for ( int k = 0; k < 4; k++ ){
        tet[4*i+k] -= 1;
      }
    }
  
    // Free the memory and exit from netgen
    Ng_DeleteMesh(m);
    Ng_Exit();

    // Print out the volume mesh to a VTK file
    FILE *fp = fopen("volume-mesh.vtk", "w");
    if (fp){
      fprintf(fp, "# vtk DataFile Version 3.0\n");
      fprintf(fp, "vtk output\nASCII\n");
      fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

      // Write out the points
      fprintf(fp, "POINTS %d float\n", num_points);
      for ( int k = 0; k < num_points; k++ ){
        fprintf(fp, "%e %e %e\n", Xt[k].x, Xt[k].y, Xt[k].z);
      }

      fprintf(fp, "\nCELLS %d %d\n", num_tet, 5*num_tet);
    
      // Write out the cell connectivities
      for ( int k = 0; k < num_tet; k++ ){
        fprintf(fp, "4 %d %d %d %d\n", 
                tet[4*k], tet[4*k+1], tet[4*k+2], tet[4*k+3]);
      }

      // All quadrilaterals
      fprintf(fp, "\nCELL_TYPES %d\n", num_tet);
      for ( int k = 0; k < num_tet; k++ ){
        fprintf(fp, "%d\n", 10);
      }
      fclose(fp);
    }

    // Decref the objects
    mesh->decref();
    geo->decref();
  }
  
  TMRFinalize();
  MPI_Finalize();

  return 0;
}

#endif // TMR_HAS_OPENCASCADE
