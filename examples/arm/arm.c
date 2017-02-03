#include "TMRTopology.h"
#include "TMRBspline.h"
#include "TMRQuadForest.h"
#include "TACSAssembler.h"
#include "isoFSDTStiffness.h"
#include "PlaneStressQuad.h"
#include "PlaneStressTraction.h"
#include "TACSToFH5.h"
#include <stdio.h>
#include <math.h>

/*
  Create a transfinite interpolation of the surface
*/
class TFIPlanar : public TMRSurface {
 public:
  TFIPlanar( TMRVertex *v[], TMREdge *e[], const int edge_dir[] ){
    for ( int i = 0; i < 4; i++ ){
      v[i]->evalPoint(&pt[i]);
      e[i]->getCurve(&curve[i]);
      dir[i] = edge_dir[i];
    }
  }
  // Get the parameter range for this surface
  void getRange( double *umin, double *vmin,
                 double *umax, double *vmax ){
    *umin = 0.0; *umax = 1.0;
    *vmin = 0.0; *vmax = 1.0;
  }
 
  // Given the parametric point, compute the x,y,z location
  int evalPoint( double u, double v, TMRPoint *X ){
    // Evaluate the curve locations
    TMRPoint c[4];

    // Evaluate the curves along the v-direction
    int fail = 0;
    if (dir[0] > 0){ fail = fail || curve[0]->evalPoint(v, &c[0]); }
    else { fail = fail || curve[0]->evalPoint(1.0-v, &c[0]); }
    
    if (dir[1] > 0){ fail = fail || curve[1]->evalPoint(v, &c[1]); }
    else { fail = fail || curve[1]->evalPoint(1.0-v, &c[1]); }

    if (dir[2] > 0){ fail = fail || curve[2]->evalPoint(u, &c[2]); }
    else { fail = fail || curve[2]->evalPoint(1.0-u, &c[2]); }

    if (dir[3] > 0){ fail = fail || curve[3]->evalPoint(u, &c[3]); }
    else { fail = fail || curve[3]->evalPoint(1.0-u, &c[3]); }

    if (fail){ return fail; }

    // Set the points
    X->x = (1.0-u)*c[0].x + u*c[1].x + (1.0-v)*c[2].x + v*c[3].x
      - ((1.0-u)*(1.0-v)*pt[0].x + u*(1.0-v)*pt[1].x + v*(1.0-u)*pt[2].x + u*v*pt[3].x);
    X->y = (1.0-u)*c[0].y + u*c[1].y + (1.0-v)*c[2].y + v*c[3].y
      - ((1.0-u)*(1.0-v)*pt[0].y + u*(1.0-v)*pt[1].y + v*(1.0-u)*pt[2].y + u*v*pt[3].y);
    X->z = (1.0-u)*c[0].z + u*c[1].z + (1.0-v)*c[2].z + v*c[3].z
      - ((1.0-u)*(1.0-v)*pt[0].z + u*(1.0-v)*pt[1].z + v*(1.0-u)*pt[2].z + u*v*pt[3].z); 

    return 1;
  }
  
  // Perform the inverse evaluation
  int invEvalPoint( TMRPoint p, double *u, double *v ){
    return 1;
  }

  // Given the parametric point, evaluate the first derivative 
  int evalDeriv( double u, double v, 
                 TMRPoint *Xu, TMRPoint *Xv ){
    return 1;
  }
 private:
  int dir[4];
  TMRPoint pt[4];
  TMRCurve *curve[4];
};

/*
  Create a circle centered at the point p
*/
TMRBsplineCurve* createCircle( TMRPoint c, double r ){
  // Set the points and weights for the B-spline circle
  const int nctl = 9;
  const int ku = 3;
  TMRPoint p[nctl];
  double wts[nctl];
  memset(p, 0, nctl*sizeof(TMRPoint));

  // Set the knot locations
  double Tu[] = {
    0.0, 0.0, 0.0,
    0.25, 0.25,
    0.5, 0.5,
    0.75, 0.75, 
    1.0, 1.0, 1.0};
  
  for ( int k = 0; k < nctl; k++ ){
    p[k] = c;
  }

  // Set the weights
  double sqrt2 = 1.0/sqrt(2.0);

  // c + (r,0)
  p[0].x += r;
  wts[0] = 1.0;

  // c + (r,r)
  p[1].x += r;
  p[1].y += r;
  wts[1] = sqrt2;

  // c + (0,r)
  p[2].y += r;
  wts[2] = 1.0;
  
  // c + (-r,r)
  p[3].x -= r;
  p[3].y += r;
  wts[3] = sqrt2;

  // c + (-r,0)
  p[4].x -= r;
  wts[4] = 1.0;

  // c + (-r,-r)
  p[5].x -= r;
  p[5].y -= r;
  wts[5] = sqrt2;

  // c + (0,-r)
  p[6].y -= r;
  wts[6] = 1.0;
  
  // c + (r,-r)
  p[7].x += r;
  p[7].y -= r;
  wts[7] = sqrt2;

  // c + (0,r)
  p[8].x += r;
  wts[8] = 1.0;
  
  // Create the circle
  TMRBsplineCurve *curve =
    new TMRBsplineCurve(nctl, ku, Tu, wts, p);

  // Return the curve
  return curve;
}

TMRBsplineCurve* createSemiCircle( TMRPoint center, double r,
                                   double theta ){
  // Set the points and weights for the B-spline circle
  const int nctl = 5;
  const int ku = 3;
  TMRPoint p[nctl];
  double wts[nctl];
  memset(p, 0, nctl*sizeof(TMRPoint));

  for ( int k = 0; k < nctl; k++ ){
    p[k] = center;
  }

  double Tu[] = {
    0.0, 0.0, 0.0,
    0.5, 0.5,
    1.0, 1.0, 1.0};

  // Set the weights
  double sqrt2 = 1.0/sqrt(2.0);

  // Compute the sine/cosine
  double c = cos(theta);
  double s = sin(theta);

  // Use the transformation
  // [ c  -s ]
  // [ s   c ]
  
  // c + (r,0)
  p[0].x += r*c;
  p[0].y += r*s;
  wts[0] = 1.0;

  // c + (r,r)
  p[1].x += (c - s)*r;
  p[1].y += (c + s)*r;
  wts[1] = sqrt2;

  // c + (0,r)
  p[2].x -= s*r;
  p[2].y += c*r;
  wts[2] = 1.0;
  
  // c + (-r,r)
  p[3].x -= (c + s)*r;
  p[3].y += (c - s)*r;
  wts[3] = sqrt2;

  // c + (-r,0)
  p[4].x -= c*r;
  p[4].y -= s*r; 
  wts[4] = 1.0;
  
  // Create the circle
  TMRBsplineCurve *curve =
    new TMRBsplineCurve(nctl, ku, Tu, wts, p);

  // Return the curve
  return curve;
}

/*
  Create a line between two points
*/
TMRBsplineCurve* createLine( TMRPoint p1, TMRPoint p2 ){
  TMRPoint p[2];
  p[0] = p1;
  p[1] = p2;
  return new TMRBsplineCurve(2, 2, p);
}

/*
  Create a line between two specified vertices
*/
TMRBsplineCurve* createLine( TMRVertex *v1, TMRVertex *v2 ){
  TMRPoint p[2];
  v1->evalPoint(&p[0]);
  v2->evalPoint(&p[1]);
  TMRBsplineCurve *bspline = new TMRBsplineCurve(2, 2, p);
  bspline->setVertices(v1, v2);
  return bspline;
}

/*
  This arm example is set up as follows:
*/
TMRTopology* setUpTopology( MPI_Comm comm, 
                            double r1, double r2, 
                            double L, double t,
                            int write_vtk_files=0 ){
  TMRPoint p1;
  p1.x = p1.y = p1.z = 0.0;

  TMRPoint p2;
  p2.x = L;
  p2.y = p2.z = 0.0;

  TMRPoint p3;
  p3.x = 0.0;
  p3.y = r1+t;
  p3.z = 0.0;

  TMRPoint p4;
  p4.x = L;
  p4.y = r2+t;
  p4.z = 0.0;

  TMRPoint p5;
  p5.x = 0.0;
  p5.y = -(r1+t);
  p5.z = 0.0;

  TMRPoint p6;
  p6.x = L;
  p6.y = -(r2+t);
  p6.z = 0.0;

  // Set the curves that form the outline of the bracket
  TMRBsplineCurve *inner1 = createCircle(p1, r1);
  TMRBsplineCurve *inner2 = createCircle(p2, r2);
  TMRBsplineCurve *outer1 = createSemiCircle(p1, r1+t, 0.5*M_PI);
  TMRBsplineCurve *outer2 = createSemiCircle(p2, r2+t, 1.5*M_PI);
  TMRBsplineCurve *line1 = createLine(p3, p4);
  TMRBsplineCurve *line2 = createLine(p5, p6);

  // Set the names of the curves
  inner1->setAttribute("inner1");
  inner2->setAttribute("inner2");

  // Write out the geometry
  if (write_vtk_files){
    inner1->writeToVTK("innerA.vtk");
    inner2->writeToVTK("innerB.vtk");
    outer1->writeToVTK("outerA.vtk");
    outer2->writeToVTK("outerB.vtk");
    line1->writeToVTK("lineA.vtk");
    line2->writeToVTK("lineB.vtk");
  }

  // Create the vertices within the model
  const int num_vertices = 18;
  TMRVertex *v[18];
  v[0] = new TMRVertexFromCurve(inner1, 0.0);
  v[1] = new TMRVertexFromCurve(inner1, 0.25);
  v[2] = new TMRVertexFromCurve(inner1, 0.5);
  v[3] = new TMRVertexFromCurve(inner1, 0.75);

  v[4] = new TMRVertexFromCurve(inner2, 0.0);
  v[5] = new TMRVertexFromCurve(inner2, 0.25);
  v[6] = new TMRVertexFromCurve(inner2, 0.5);
  v[7] = new TMRVertexFromCurve(inner2, 0.75);

  v[8] = new TMRVertexFromCurve(outer1, 0.0);
  v[9] = new TMRVertexFromCurve(outer1, 0.5);
  v[10] = new TMRVertexFromCurve(outer1, 1.0);

  v[11] = new TMRVertexFromCurve(outer2, 0.0);
  v[12] = new TMRVertexFromCurve(outer2, 0.5);
  v[13] = new TMRVertexFromCurve(outer2, 1.0);

  v[14] = new TMRVertexFromCurve(line1, 0.45);
  v[15] = new TMRVertexFromCurve(line1, 0.7);

  v[16] = new TMRVertexFromCurve(line2, 0.45);
  v[17] = new TMRVertexFromCurve(line2, 0.7);

  // Create the edges within the model
  const int num_edges = 29;
  TMREdge *e[29];
  // First inner edges
  e[0] = new TMREdge(new TMRSplitCurve(inner1, v[0], v[1]));
  e[1] = new TMREdge(new TMRSplitCurve(inner1, v[1], v[2]));
  e[2] = new TMREdge(new TMRSplitCurve(inner1, v[2], v[3]));
  e[3] = new TMREdge(new TMRSplitCurve(inner1, v[3], v[0]));

  // Second set of inner edges
  e[4] = new TMREdge(new TMRSplitCurve(inner2, v[4], v[5]));
  e[5] = new TMREdge(new TMRSplitCurve(inner2, v[5], v[6]));
  e[6] = new TMREdge(new TMRSplitCurve(inner2, v[6], v[7]));
  e[7] = new TMREdge(new TMRSplitCurve(inner2, v[7], v[4]));

  // Outer edges -- all the way around
  e[8] = new TMREdge(new TMRSplitCurve(outer1, v[8], v[9]));
  e[9] = new TMREdge(new TMRSplitCurve(outer1, v[9], v[10]));

  // Lower line segment edges
  e[10] = new TMREdge(new TMRSplitCurve(line2, v[10], v[16]));
  e[11] = new TMREdge(new TMRSplitCurve(line2, v[16], v[17]));
  e[12] = new TMREdge(new TMRSplitCurve(line2, v[17], v[11]));

  // Outer edges
  e[13] = new TMREdge(new TMRSplitCurve(outer2, v[11], v[12]));
  e[14] = new TMREdge(new TMRSplitCurve(outer2, v[12], v[13]));

  // Upper line segments  
  e[15] = new TMREdge(new TMRSplitCurve(line1, v[15], v[13]));
  e[16] = new TMREdge(new TMRSplitCurve(line1, v[14], v[15]));
  e[17] = new TMREdge(new TMRSplitCurve(line1, v[8], v[14]));

  // The inner line segments - most of these requrie the creation of a new line
  e[18] = new TMREdge(createLine(v[1], v[8]));
  e[19] = new TMREdge(createLine(v[2], v[9]));
  e[20] = new TMREdge(createLine(v[3], v[10]));

  e[21] = new TMREdge(createLine(v[0], v[16]));
  e[22] = new TMREdge(createLine(v[6], v[17]));

  e[23] = new TMREdge(createLine(v[7], v[11]));
  e[24] = new TMREdge(createLine(v[4], v[12]));
  e[25] = new TMREdge(createLine(v[5], v[13]));

  e[26] = new TMREdge(createLine(v[6], v[15]));
  e[27] = new TMREdge(createLine(v[0], v[14]));

  // Create the last edge joining the two segments
  e[28] = new TMREdge(createLine(v[0], v[6]));

  if (write_vtk_files){
    for ( int i = 0; i < num_edges; i++ ){
      TMRCurve *curve;
      e[i]->getCurve(&curve);
      char filename[128];
      sprintf(filename, "edge%d.vtk", i);
      curve->writeToVTK(filename);
    }
  }

  // Now, create all of the faces
  const int num_faces = 10;
  TMRFace *f[10];
  
  // Arrays to store the vertices/edges for each
  int dir0[4] = {1, -1, 1, 1};
  TMRVertex *v0[4] = {v[0], v[14], v[1], v[8]}; 
  TMREdge *e0[4] = {e[0], e[17], e[27], e[18]};
  f[0] = new TMRFace(new TFIPlanar(v0, e0, dir0), e0, dir0);

  int dir1[4] = {1, 1, 1, 1};
  TMRVertex *v1[4] = {v[1], v[8], v[2], v[9]}; 
  TMREdge *e1[4] = {e[1], e[8], e[18], e[19]};
  f[1] = new TMRFace(new TFIPlanar(v1, e1, dir1), e1, dir1);

  int dir2[4] = {1, 1, 1, 1};
  TMRVertex *v2[4] = {v[2], v[9], v[3], v[10]}; 
  TMREdge *e2[4] = {e[2], e[9], e[19], e[20]};
  f[2] = new TMRFace(new TFIPlanar(v2, e2, dir2), e2, dir2);

  int dir3[4] = {1, 1, 1, 1};
  TMRVertex *v3[4] = {v[3], v[10], v[0], v[16]}; 
  TMREdge *e3[4] = {e[3], e[10], e[20], e[21]};
  f[3] = new TMRFace(new TFIPlanar(v3, e3, dir3), e3, dir3);

  int dir4[4] = {1, 1, 1, 1};
  TMRVertex *v4[4] = {v[0], v[16], v[6], v[17]}; 
  TMREdge *e4[4] = {e[28], e[11], e[21], e[22]};
  f[4] = new TMRFace(new TFIPlanar(v4, e4, dir4), e4, dir4);

  int dir5[4] = {1, 1, 1, 1};
  TMRVertex *v5[4] = {v[6], v[17], v[7], v[11]}; 
  TMREdge *e5[4] = {e[6], e[12], e[22], e[23]};
  f[5] = new TMRFace(new TFIPlanar(v5, e5, dir5), e5, dir5);

  int dir6[4] = {1, 1, 1, 1};
  TMRVertex *v6[4] = {v[7], v[11], v[4], v[12]}; 
  TMREdge *e6[4] = {e[7], e[13], e[23], e[24]};
  f[6] = new TMRFace(new TFIPlanar(v6, e6, dir6), e6, dir6);

  int dir7[4] = {1, 1, 1, 1};
  TMRVertex *v7[4] = {v[4], v[12], v[5], v[13]}; 
  TMREdge *e7[4] = {e[4], e[14], e[24], e[25]};
  f[7] = new TMRFace(new TFIPlanar(v7, e7, dir7), e7, dir7);

  int dir8[4] = {1, -1, 1, 1};
  TMRVertex *v8[4] = {v[5], v[13], v[6], v[15]}; 
  TMREdge *e8[4] = {e[5], e[15], e[25], e[26]};
  f[8] = new TMRFace(new TFIPlanar(v8, e8, dir8), e8, dir8);

  int dir9[4] = {1, 1, 1, 1};
  TMRVertex *v9[4] = {v[0], v[6], v[14], v[15]}; 
  TMREdge *e9[4] = {e[27], e[26], e[28], e[16]};
  f[9] = new TMRFace(new TFIPlanar(v9, e9, dir9), e9, dir9);

  if (write_vtk_files){
    for ( int i = 0; i < num_faces; i++ ){
      TMRSurface *surf;
      f[i]->getSurface(&surf);
      char filename[128];
      sprintf(filename, "surf%d.vtk", i);
      surf->writeToVTK(filename);
    }
  }

  // Create the topology object with the vertices, faces and edge objects
  TMRTopology *topo = 
    new TMRTopology(comm, v, num_vertices, 
                    e, num_edges, f, num_faces);

  return topo;
}

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  // Set the communicator
  MPI_Comm comm = MPI_COMM_WORLD;

  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  // Create the dimensions of the part
  double L = 10.0;
  double t = 1.5;
  double r1 = 2.0;
  double r2 = 1.0;

  // Create the topology
  TMRTopology *topo = setUpTopology(comm, r1, r2, L, t);
  topo->incref();

  TMRQuadForest *forest = new TMRQuadForest(comm);
  forest->incref();

  forest->setTopology(topo);
  forest->createTrees(5);
  forest->repartition();
  forest->balance(1);
  forest->repartition();
  forest->createNodes(2);

    // Allocate the stiffness object
  TacsScalar rho = 2570.0, E = 70e9, nu = 0.3;
  PlaneStressStiffness *stiff = 
    new PlaneStressStiffness(rho, E, nu);

  // Allocate the solid element class
  TACSElement *elem = new PlaneStressQuad<2>(stiff, LINEAR, mpi_rank);

  // Find the number of nodes for this processor
  const int *range;
  forest->getOwnedNodeRange(&range);
  int num_nodes = range[mpi_rank+1] - range[mpi_rank];

  // Create the mesh
  double tmesh = MPI_Wtime();
  int *elem_conn, num_elements = 0;
  forest->createMeshConn(&elem_conn, &num_elements);
  tmesh = MPI_Wtime() - tmesh;
  printf("[%d] Mesh: %f\n", mpi_rank, tmesh);

  // Get the dependent node information
  const int *dep_ptr, *dep_conn;
  const double *dep_weights;

  // Create/retrieve the dependent node information
  double tdep = MPI_Wtime();
  int num_dep_nodes = 
    forest->getDepNodeConn(&dep_ptr, &dep_conn,
                           &dep_weights);
  tdep = MPI_Wtime() - tdep;
  printf("[%d] Dependent nodes: %f\n", mpi_rank, tdep);

  // Create the associated TACSAssembler object
  int vars_per_node = 2;
  TACSAssembler *tacs = new TACSAssembler(comm, vars_per_node,
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
  delete [] elem_conn;
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

  // Set the boundary conditions  
  int dim;
  TMRQuadrant *bcs;
  TMRQuadrantArray *inner1 = forest->getNodesWithAttribute("inner1");
  inner1->getArray(&bcs, &dim);
  for ( int i = 0; i < dim; i++ ){
    int node = bcs[i].tag;
    tacs->addBCs(1, &node);
  }
  delete inner1;

  // Initialize the model
  tacs->initialize();

  // Get the quadrants from the forest
  TMRQuadrantArray *quadrants;
  forest->getQuadrants(&quadrants);

  // Add the loads to the
  TMRQuadrantArray *inner2 = forest->getQuadsWithAttribute("inner2");
  inner2->getArray(&bcs, &dim);
  TACSAuxElements *aux = new TACSAuxElements(dim);

  double tx = 0.0;
  double ty = 100.0;
  PSQuadTraction<2> *trac[4];
  for ( int k = 0; k < 4; k++ ){
    trac[k] = new PSQuadTraction<2>(k, tx, ty);
    trac[k]->incref();
  }

  // Set up the boundary tractions 
  for ( int i = 0; i < dim; i++ ){
    int face = bcs[i].tag;
    if (face >= 0 && face < 4){
      int use_node_search = 0;
      TMRQuadrant *t = quadrants->contains(&bcs[i], use_node_search);
      if (t){
        aux->addElement(t->tag, trac[face]);
      }
    }
  }
  tacs->setAuxElements(aux);
  delete inner2;

  // Get the nodes
  TMRQuadrantArray *nodes;
  forest->getNodes(&nodes);

  // Get the quadrants associated with the nodes
  int size;
  TMRQuadrant *array;
  nodes->getArray(&array, &size);
 
  // Create the node vector
  TACSBVec *X = tacs->createNodeVec();
  X->incref();

  // Retrieve the node array
  TacsScalar *Xn;
  X->getArray(&Xn);

  // Get the points
  TMRPoint *Xp;
  forest->getPoints(&Xp);

  // Loop over all the nodes
  for ( int i = 0; i < size; i++ ){
    if (array[i].tag >= range[mpi_rank] &&
        array[i].tag < range[mpi_rank+1]){
      int loc = array[i].tag - range[mpi_rank];
      Xn[3*loc] = Xp[i].x;
      Xn[3*loc+1] = Xp[i].y;
      Xn[3*loc+2] = Xp[i].z;
    }
  }
    
  // Set the node locations into TACSAssembler
  tacs->setNodes(X);

  // Create the preconditioner
  TACSBVec *res = tacs->createVec();
  TACSBVec *ans = tacs->createVec();
  FEMat *mat = tacs->createFEMat();

  // Increment the reference count to the matrix/vectors
  res->incref();
  ans->incref();
  mat->incref();

  // Allocate the factorization
  int lev = 4500;
  double fill = 10.0;
  int reorder_schur = 1;
  PcScMat *pc = new PcScMat(mat, lev, fill, reorder_schur);
  pc->incref();

  // Assemble and factor the stiffness/Jacobian matrix
  double alpha = 1.0, beta = 0.0, gamma = 0.0;
  tacs->assembleJacobian(alpha, beta, gamma, res, mat);
  mat->applyBCs();
  pc->factor();

  int gmres_iters = 10; // Number of GMRES iterations 
  int nrestart = 2; // Number of allowed restarts
  int is_flexible = 1; // Is a flexible preconditioner?
  GMRES *gmres = new GMRES(mat, pc, gmres_iters, 
                           nrestart, is_flexible);

  res->set(1.0);
  res->applyBCs();
  gmres->solve(res, ans);
  tacs->setVariables(ans);

  // Create and write out an fh5 file
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS |
                             TACSElement::OUTPUT_STRESSES |
                             TACSElement::OUTPUT_EXTRAS);
  TACSToFH5 *f5 = new TACSToFH5(tacs, PLANE_STRESS, write_flag);
  f5->incref();
    
  // Write out the solution
  f5->writeToFile("output.f5");
  f5->decref();
  tacs->decref();
  forest->decref();
  topo->decref();

  TMRFinalize();
  MPI_Finalize();
  
  return (0);
}





