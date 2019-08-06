#include "TMRBspline.h"
#include "TMRNativeTopology.h"
#include "TMRMesh.h"
#include "TMRQuadForest.h"
#include "TACSAssembler.h"
#include "TACSToFH5.h"
#include "PlaneStressQuad.h"
#include <stdio.h>
#include <math.h>

const int rae2822_npts = 128;
const double rae2822_pts[] =
  {1.0000, 0.0000,
   0.9994, 0.0003,
   0.99759, 0.00069,
   0.99459, 0.00132,
   0.99039, 0.00218,
   0.98502, 0.00326,
   0.97847, 0.00455,
   0.97077, 0.00606,
   0.96194, 0.00775,
   0.952, 0.00964,
   0.94096, 0.0117,
   0.92886, 0.01393,
   0.91574, 0.01627,
   0.9016, 0.01874,
   0.88651, 0.02131,
   0.87048, 0.02397,
   0.85355, 0.0267,
   0.83578, 0.02948,
   0.8172, 0.03231,
   0.79785, 0.03514,
   0.77778, 0.03795,
   0.75705, 0.04075,
   0.7357, 0.04338,
   0.71378, 0.04612,
   0.69134, 0.04857,
   0.66845, 0.05112,
   0.64514, 0.05339,
   0.62149, 0.05547,
   0.59754, 0.05733,
   0.57336, 0.05895,
   0.54901, 0.0603,
   0.52453, 0.06135,
   0.5, 0.06212,
   0.47547, 0.06261,
   0.45099, 0.06286,
   0.42663, 0.06285,
   0.40245, 0.06263,
   0.37851, 0.0622,
   0.35486, 0.06155,
   0.33156, 0.0607,
   0.30866, 0.05967,
   0.28622, 0.05848,
   0.2643, 0.05713,
   0.24295, 0.05556,
   0.22221, 0.05377,
   0.20215, 0.05187,
   0.1828, 0.04987,
   0.16422, 0.04778,
   0.14645, 0.04558,
   0.12952, 0.04321,
   0.11349, 0.04073,
   0.0984, 0.03817,
   0.08427, 0.03552,
   0.07114, 0.0328,
   0.05904, 0.03004,
   0.04801, 0.02726,
   0.03806, 0.02445,
   0.02923, 0.02163,
   0.02153, 0.01875,
   0.01498, 0.01579,
   0.00961, 0.01269,
   0.00541, 0.00945,
   0.00241, 0.00642,
   0.0006, 0.00323,
   0.0000, 0.0000,
   0.0006, -0.00317,
   0.00241, -0.00658,
   0.00541, -0.00957,
   0.00961, -0.01273,
   0.01498, -0.0158,
   0.02153, -0.0188,
   0.02923, -0.0218,
   0.03806, -0.02472,
   0.04801, -0.02761,
   0.05904, -0.03042,
   0.07114, -0.03315,
   0.08427, -0.03584,
   0.0984, -0.03844,
   0.11349, -0.04094,
   0.12952, -0.04333,
   0.14645, -0.04561,
   0.16422, -0.04775,
   0.1828, -0.04977,
   0.20215, -0.05167,
   0.22221, -0.0534,
   0.24295, -0.05498,
   0.2643, -0.05638,
   0.28622, -0.05753,
   0.30866, -0.05843,
   0.33156, -0.059,
   0.35486, -0.05919,
   0.37851, -0.05893,
   0.40245, -0.05817,
   0.42663, -0.05689,
   0.45099, -0.05515,
   0.47547, -0.05297,
   0.5, -0.05044,
   0.52453, -0.04761,
   0.54901, -0.04452,
   0.57336, -0.04127,
   0.59754, -0.03791,
   0.62149, -0.03463,
   0.64514, -0.0311,
   0.66845, -0.0277,
   0.69134, -0.02438,
   0.71378, -0.02118,
   0.7357, -0.01812,
   0.75705, -0.01524,
   0.77778, -0.01256,
   0.79785, -0.01013,
   0.8172, -0.00792,
   0.83578, -0.00594,
   0.85355, -0.00422,
   0.87048, -0.00273,
   0.88651, -0.00149,
   0.9016, -0.00049,
   0.91574, 0.00027,
   0.92886, 0.00081,
   0.94096, 0.00113,
   0.952, 0.00125,
   0.96194, 0.00125,
   0.97077, 0.00113,
   0.97847, 0.00094,
   0.98502, 0.00071,
   0.99039, 0.00048,
   0.99459, 0.00026,
   0.99759, 0.00009,
   0.9994, -0.00001,
   1.000, 0.0000 };

void test_surface_lofter( double htarget ){
  // Create the control points
  int nctl = rae2822_npts;
  TMRPoint line_pts[rae2822_npts];
  memset(line_pts, 0, nctl*sizeof(TMRPoint));
  
  // Create a series of curves
  int num_curves = 5;
  TMRBsplineCurve *lofts[5];
  double chord[5] = {6.0, 4.0, 3.0,  2.0, 1.0};
  double twist[5] = {0, -1.0, -2.0, -4.0, -10.0};
  for ( int k = 0; k < num_curves; k++ ){
    twist[k] *= M_PI/180.0;
  }

  for ( int k = 0; k < num_curves; k++ ){
    for ( int i = 0; i < nctl; i++ ){
      double c = cos(twist[k]);
      double s = sin(twist[k]);
      double x = rae2822_pts[2*i];
      double y = rae2822_pts[2*i+1];
      line_pts[i].x = chord[k]*(c*x + y*s) + 3*k;
      line_pts[i].y = chord[k]*(-s*x + y*c);
      line_pts[i].z = 5.0*k;
    }

    // Create the interpolation object
    TMRCurveInterpolation *interper = 
      new TMRCurveInterpolation(line_pts, nctl);
    interper->incref();

    interper->setNumControlPoints(13 + 2*k);
    lofts[k] = interper->createCurve(4);
    lofts[k]->incref();

    // Free the interpolation object
    interper->decref();
  }  

  // Create the lofter object
  TMRCurveLofter *lofter = new TMRCurveLofter(lofts, num_curves);
  lofter->incref();

  // Create the surface object from the lofted surface
  TMRBsplineSurface *surface = lofter->createSurface(4);
  surface->incref();
  lofter->decref();
  TMRFace *face = new TMRFaceFromSurface(surface);
  face->incref();

  surface->writeToVTK("bspline_surface.vtk");

  // (0,1) -- (1,1)
  //   |        |
  // (0,0) -- (1,0)
  double pts1[] = {0.1, 0.0, 0.4, 0.0};
  double pts2[] = {0.4, 0.0, 0.4, 1.0};
  double pts3[] = {0.4, 1.0, 0.1, 1.0};
  double pts4[] = {0.1, 1.0, 0.1, 0.0};
  TMRBsplinePcurve *p1 = new TMRBsplinePcurve(2, 2, pts1);
  TMRBsplinePcurve *p2 = new TMRBsplinePcurve(2, 2, pts2);
  TMRBsplinePcurve *p3 = new TMRBsplinePcurve(2, 2, pts3);
  TMRBsplinePcurve *p4 = new TMRBsplinePcurve(2, 2, pts4);

  // Create the curves and add them to the surface
  int ncurves = 4;
  TMREdge *curves[4];
  curves[0] = new TMREdgeFromFace(face, p1);
  curves[1] = new TMREdgeFromFace(face, p2);
  curves[2] = new TMREdgeFromFace(face, p3);
  curves[3] = new TMREdgeFromFace(face, p4);

  // Create the boundary curves for the surface
  TMRVertexFromEdge *v1 = new TMRVertexFromEdge(curves[0], 0.0);
  TMRVertexFromEdge *v2 = new TMRVertexFromEdge(curves[1], 0.0);
  TMRVertexFromEdge *v3 = new TMRVertexFromEdge(curves[2], 0.0);
  TMRVertexFromEdge *v4 = new TMRVertexFromEdge(curves[3], 0.0);

  int num_verts = 4;
  TMRVertex *verts[] = {v1, v2, v3, v4};

  // Set the vertices
  curves[0]->setVertices(v1, v1);
  curves[1]->setVertices(v2, v3);
  curves[2]->setVertices(v3, v4);
  curves[3]->setVertices(v4, v1);
  
  // Set the directions of the curves
  int dir[4];
  dir[0] = 1;
  dir[1] = 1;
  dir[2] = 1;
  dir[3] = 1;
  TMREdgeLoop *loop = new TMREdgeLoop(ncurves, curves, dir);
  face->addEdgeLoop(1, loop);

  // Set the communicator
  MPI_Comm comm = MPI_COMM_WORLD;

  // Create the model
  TMRModel *geo = new TMRModel(num_verts, verts,
                               ncurves, curves,
                               1, &face);

  // Allocate the new mesh
  TMRMesh *mesh = new TMRMesh(comm, geo);

  // Mesh the geometry
  TMRMeshOptions options;
  mesh->mesh(options, htarget);
  mesh->writeToVTK("mesh.vtk");

  // Create the mesh
  TMRModel *geo_mesh = mesh->createModelFromMesh();
  geo_mesh->incref();

  TMRTopology *topo = new TMRTopology(comm, geo_mesh);
  topo->incref();

  TMRQuadForest *forest = new TMRQuadForest(comm);
  forest->incref();

  // Set up the forest 
  forest->setTopology(topo);

  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  // Create random trees
  forest->createRandomTrees(10, 0, 4);

  double tbal = MPI_Wtime();
  forest->balance();
  tbal = MPI_Wtime() - tbal;
  printf("[%d] Balance: %f\n", mpi_rank, tbal);
  double tnodes = MPI_Wtime();
  forest->createNodes();
  tnodes = MPI_Wtime() - tnodes;
  printf("[%d] Nodes: %f\n", mpi_rank, tnodes);

  // Allocate the stiffness object
  TacsScalar rho = 2570.0, E = 70e9, nu = 0.3;
  PlaneStressStiffness *stiff = 
    new PlaneStressStiffness(rho, E, nu);

  // Allocate the solid element class
  TACSElement *elem = new PlaneStressQuad<2>(stiff, LINEAR);
  
  // Create the TACSAssembler objects
  TACSAssembler *tacs = NULL;

  // Find the number of nodes for this processor
  const int *range;
  forest->getOwnedNodeRange(&range);
  int num_nodes = range[mpi_rank+1] - range[mpi_rank];

  // Create the mesh
  const int *elem_conn;
  int num_elements = 0;
  forest->getNodeConn(&elem_conn, &num_elements);

  // Get the dependent node information
  const int *dep_ptr, *dep_conn;
  const double *dep_weights;
  int num_dep_nodes = forest->getDepNodeConn(&dep_ptr, &dep_conn,
                                             &dep_weights);

  // Create the associated TACSAssembler object
  int vars_per_node = 2;
  tacs = new TACSAssembler(comm, vars_per_node,
                           num_nodes, num_elements,
                           num_dep_nodes);
  tacs->incref();

  // Set the element ptr
  int *ptr = new int[ num_elements+1 ];
  for ( int i = 0; i < num_elements+1; i++ ){
    ptr[i] = 4*i;
  }
  
  // Set the element connectivity into TACSAssembler
  tacs->setElementConnectivity(elem_conn, ptr);
  delete [] ptr;
    
  // Set the dependent node information
  tacs->setDependentNodes(dep_ptr, dep_conn, dep_weights);

  // Set the elements
  TACSElement **elems = new TACSElement*[ num_elements ];
  for ( int k = 0; k < num_elements; k++ ){
    elems[k] = elem;
  }
    
  // Set the element array
  tacs->setElements(elems);
  delete [] elems;
    
  // Initialize the TACSAssembler object
  tacs->initialize();

  // Create the node vector
  TacsScalar *Xn;
  TACSBVec *X = tacs->createNodeVec();
  X->getArray(&Xn);

  // Get the points
  TMRPoint *Xp;
  forest->getPoints(&Xp);

  // Get all of the local node numbers
  const int *local_nodes;
  int num_local_nodes = forest->getNodeNumbers(&local_nodes);
  
  // Loop over all the nodes
  for ( int i = 0; i < num_local_nodes; i++ ){
    if (local_nodes[i] >= range[mpi_rank] &&
        local_nodes[i] < range[mpi_rank+1]){
      int loc = local_nodes[i] - range[mpi_rank];
      Xn[3*loc] = Xp[i].x;
      Xn[3*loc+1] = Xp[i].y;
      Xn[3*loc+2] = Xp[i].z;
    }
  }
    
  // Set the node locations into TACSAssembler
  tacs->setNodes(X);

  // Create and write out an fh5 file
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS |
                             TACSElement::OUTPUT_STRESSES |
                             TACSElement::OUTPUT_EXTRAS);
  TACSToFH5 *f5 = new TACSToFH5(tacs, TACS_PLANE_STRESS, write_flag);
  f5->incref();
    
  // Write out the solution
  f5->writeToFile("output.f5");
  f5->decref();

  // Deallocate the mesh
  forest->decref();
  topo->decref();
  geo_mesh->decref();

  // Free the objects
  surface->decref();
  for ( int k = 0; k < num_curves; k++ ){
    lofts[k]->decref();
  }
}

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  double htarget = 0.1;
  for ( int i = 0; i < argc; i++ ){
    if (sscanf(argv[i], "h=%lf", &htarget) == 1){
      if (htarget < 0.01){ htarget = 0.01; }
      if (htarget > 1.0){ htarget = 1.0; }
    }
  }

  test_surface_lofter(htarget);

  TMRFinalize();
  MPI_Finalize();
  return (0);
}
