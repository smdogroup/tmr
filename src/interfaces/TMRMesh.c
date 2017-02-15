

TMRCurveMesh::TMRCurveMesh( TMRCurve *_curve ){
  curve = _curve;
  curve->incref();

  npts = 0;
  pts = NULL;
}

/*
  Destroy the mesh for this curve, and free the underlying data
*/
TMRCurveMesh::~TMRCurveMesh(){
  curve->decref();
  if (pts){ delete [] pts; }
}

/*
  Create a mesh
*/
void TMRCurveMesh::mesh( double htarget ){
  double tmin, tmax;
  curve->getRange(&tmin, &tmax);

  // Integrate along the curve to obtain the distance function such
  // that dist(tvals[i]) = int_{tmin}^{tvals[i]} ||d{C(t)}dt||_{2} dt
  int nvals;
  double *dist, *tvals;
  curve->integrate(tmin, tmax, integration_eps, &tvals, &dist, &nvals);
  
  // Compute the number of points along this curve
  npts = 1 + (int)(dist[nvals-1]/htarget);
  if (npts < 2){ npts = 2; }

  // The average distance between points
  double d = dist[nvals-1]/(npts-1);

  // Allocate the parametric points that will be used 
  pts = new double[ npts ];

  // Set the starting/end location of the points
  pts[0] = tmin;
  pts[npts-1] = tmax;

  // Perform the integration so that the points are evenly spaced
  // along the curve
  for ( int j = 1, k = 1; (j < nvals && k < npts-1); j++ ){
    while ((k < npts-1) && 
           (dist[j-1] <= d*k && d*k < dist[j])){
      double u = 0.0;
      if (dist[j] > dist[j-1]){ 
        u = (d*k - dist[j-1])/(dist[j] - dist[j-1]);
      }
      tpts[k-1] = tvals[j-1] + (tvals[j] - tvals[j-1])*u;
      k++;
    }
  }
  
  // Free the integration result
  delete [] tvals;
  delete [] dist;
}

/*
  Create the surface mesh object

  This call does not create the underlying surface. You can set
  meshing options into the class after it is created, and then create
  the mesh.

  Note that the curve/edge meshes must be meshed before calling this
  object.
*/
TMRSurfaceMesh::TMRSurfaceMesh( int ncurves, 
                                TMRCurveMesh **cmeshes,
                                TMRSurface *surf ){
  // Set the data
  num_curves = ncurves;
  curve_meshes = new TMRCurveMesh*[ num_curves ];
  for ( int i = 0; i < num_curves; i++ ){
    curve_meshes[i] = cmeshes[i];
    curve_meshes[i]->incref();
  }

  surface = surf;
  surface->incref();
}

/*
  Create the surface mesh
*/
void TMRSurfaceMesh::mesh( double htarget ){
  // Count up the number of points and segments from the curves that
  // bound the surface





  int nseg = 0;
  int npts = num_vertices;
  for ( int i = 0; i < num_curves; i++ ){
    nseg += curve_npts[i]-1;
    npts += curve_npts[i]-2;
  }

  double *params = new double[ 2*npts ];
  int *seg = new int[ 2*(nseg+1) ];

  int pt = 0;
  for ( int i = 0; i < num_curves; i++ ){
    // Retrieve the vertices associated with this curve
    TMRVertex *v1, *v2;
    if (dir[i] > 0){
      curves[i]->getVertices(&v1, &v2);
    }
    else {
      curves[i]->getVertices(&v2, &v1);
    }

    // Perform the inverse evaluation
    TMRPoint p;
    v1->evalPoint(&p);
    surf->invEvalPoint(p, &params[2*pt], &params[2*pt+1]);
    seg[2*pt] = pt;
    seg[2*pt+1] = pt+1;
    pt++;

    // Find the point on the curve
    if (dir[i] > 0){
      for ( int j = 1; j < curve_npts[i]-1; j++ ){
        double t = curve_tpts[i][j-1];
        curves[i]->evalPoint(t, &p);
        surf->invEvalPoint(p, &params[2*pt], &params[2*pt+1]);
        seg[2*pt] = pt;
        seg[2*pt+1] = pt+1;
        pt++;
      }
    }
    else {
      for ( int j = curve_npts[i]-2; j >= 1; j-- ){
        double t = curve_tpts[i][j-1];
        curves[i]->evalPoint(t, &p);
        surf->invEvalPoint(p, &params[2*pt], &params[2*pt+1]);
        seg[2*pt] = pt;
        seg[2*pt+1] = pt+1;
        pt++;
      } 
    }
  }
  seg[2*(pt-1)+1] = 0;



  TMRTriangularize *tri = 
    new TMRTriangularize(npts, params, nseg, seg, surf);
  tri->incref();

  // Create the mesh using the frontal algorithm
  tri->frontal(htarget);

  // Write the result to a file
  tri->writeToVTK("triangle.vtk");

  // Extract the triangularization of the domain
  double *pts;
  int ntris, *tris;
  tri->getMesh(&ntris, &tris, &pts);

  // Free the triangle mesh
  tri->decref();

  // Compute the triangle edges and neighbors in the dual mesh
  int num_tri_edges;
  int *tri_edges, *tri_neighbors, *dual_edges;
  computeTriEdges(nnodes, ntris, tris, 
                  &numtriedges, &triedges,
                  &trineighbors, &dualedges);


  // Recombine the mesh into a quadrilateral mesh


}







/*
  Compute a node to triangle or node to quad data structure
*/
void TMRSurfaceMesh::computeNodeToElems( int nnodes, int nelems, int numelemnodes,
                                         const int conn[], 
                                         int **_ptr, int **_nodetoelems ){
  // Set the pointer
  int *ptr = new int[ nnodes ];
  memset(ptr, 0, (nnodes+1)*sizeof(int));

  // Count up the references
  const int conn_size = nelems*numelemnodes;
  for ( int i = 0; i < conn_size; i++ ){
    if (conn[i] >= 0){
      ptr[conn[i]+1]++;
    }
  }

  // Set the pointer into the quad array
  for ( int i = 0; i < nnodes; i++ ){
    ptr[i+1] += ptr[i];
  }

  // Compute the node to quads
  int *nodetoelems = new int[ ptr[nnodes] ];
  const int *conn_ptr = conn;
  for ( int i = 0; i < nelems; i++ ){
    for ( int j = 0; j < numelemnodes; j++ ){
      int node = conn_ptr[0];
      if (node >= 0){
        nodetoelems[ptr[node]] = i;
        ptr[node]++;
        conn_ptr++;
      }
    }
  }

  // Reset the pointer array
  for ( int i = nnodes-1; i >= 0; i-- ){
    ptr[i+1] = ptr[i];
  }
  ptr[0] = 0;

  // Set the output points
  *_ptr = ptr;
  *_nodetoelems = nodetoelems;
}

/*
  Compute all of the edges within the triangular mesh
*/
void TMRSurfaceMesh::computeTriEdges( int nnodes, int ntris, const int conn[],
                                      int *numtriedges, int **_triedges,
                                      int **_trineighbors, int **_dualedges ){
  // Compute the edges in the triangular mesh
  int *ptr;
  int *nodestotris;
  computeNodeToElems(nnodes, ntris, 3, tris, &ptr, &nodestotris);
  
  // Now compute the neighbors for each triangle
  int *triedgenums = new int[ 3*ntris ];
  for ( int i = 0; i < 3*ntris; i++ ){
    triedgenums[i] = -1;
  }

  // Quck reference from the edge index to the local node numbering
  const int enodes[3][2] = {{1, 2}, {2, 0}, {0, 1}};

  // Allocate the array for the triangle neighbors
  int *trineighbors = new int[ 3*ntris ];
  
  int count = 0;
  int ne = 0;
  for ( int i = 0; i < ntris; i++ ){
    // Search through each edge of the each triangle
    for ( int j = 0; j < 3; j++ ){
      if (triedgenums[3*i+j] < 0){
        triedgenums[3*i+j] = ne;

        // Triangle edges that have no neighbors are labeled with a -1
        trineighbors[3*i+j] = -1;

        // Search for the neighboring quad that shares this edge
        int kp = ptr[conn[3*i + enodes[j][0]]];
        int kpend = ptr[conn[3*i + enodes[j][0]]+1];
        for ( ; kp < kpend; kp++ ){
          // Find the potential quad neighbor
          int n = nodetotris[kp];

          // Don't count the same edge twice
          if (n == i){ continue; }

          // Flag to indicate that we have found the other edge (there
          // will only be at most one other match since this is
          // planar)
          int quit = 0;

          // Search over all the edges on this quad, and see
          // if they match
          for ( int e = 0; e < 3; e++ ){  
            // Check if the adjacent edge matches in either direction
            if ((conn[3*i+enodes[j][0]] == conns[3*n+enodes[e][0]] &&
                 conn[3*i+enodes[j][1]] == conns[3*n+enodes[e][1]]) ||
                (conn[3*i+enodes[j][0]] == conns[3*n+enodes[e][1]] &&
                 conn[3*i+enodes[j][1]] == conns[3*n+enodes[e][0]])){
              // Label the other edge that shares this same node
              triedgenums[3*n+e] = ne;

              // Set the triangle neighbors
              trineighbors[3*n+e] = i;
              trineighbors[3*i+j] = n;

              quit = 1;
            }
          }
          if (quit){ break; }
        }

        // Increment the edge number
        ne++;
      }
    }
  }

  // Free the data that is no longer required
  delete ptr;
  delete [] nodetotris;

  // Now we have a unique list of edge numbers and the total number of
  // edges, we can construct the unique edge list
  int *triedges = new int[ 2*ne ];
  int *dualedges = new int[ 2*ne ];
  for ( int i = 0; i < ntris; i++ ){
    for ( int j = 0; j < 3; j++ ){
      // Get the unique edge number
      int n = triedgenums[3*i+j];
      tridges[2*n] = conn[3*i+enodes[j][0]];
      triedges[2*n+1] = conn[3*i+enodes[j][1]];

      // Set the dual edge numbers - connecting triangles to other
      // triangles. Note that some elements of this array will be -1.
      dualedges[2*n] = i;
      dualedges[2*n+1] = trineighbors[3*i+j];
    }
  } 

  delete [] triedgenums;

  // Set the number of triangle edges and the triangle edges themselves
  *numtriedges = ne;
  *_triedges = triedges; 
  *_trineighbors = trineighbors;
  *_dualedges = dualedges;
}

/*
  Compute the connectivity between the edges
*/
void TMRSurfaceMesh::computeQuadEdges( int nnodes, int nquads, 
                                       const int quads[],
                                       int *_nquadedges,
                                       int **_quadedges ){
  // Compute the edges in the quadrilateral mesh
  int *ptr;
  int *nodestoquads;
  computeNodeToElems(nnodes, nquads, 4, quads, &ptr, &nodetoquads);
  
  // Now compute the neighbors for each quad
  int *quadedgenums = new int[ 4*nquads ];
  for ( int i = 0; i < 4*nquads; i++ ){
    quadedgenums[i] = -1;
  }

  // Quck reference from the quad index to the edge
  const int enodes[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};

  int count = 0;
  int ne = 0;
  for ( int i = 0; i < nquads; i++ ){
    // Search through each edge of the each quadrilateral
    for ( int j = 0; j < 4; j++ ){
      if (quadedgenums[4*i+j] < 0){
        quadedgenums[4*i+j] = ne;

        // Search for the neighboring quad that shares this edge
        int kp = ptr[quads[4*i + enodes[j][0]]];
        int kpend = ptr[quads[4*i + enodes[j][0]]+1];
        for ( ; kp < kpend; kp++ ){
          // Find the potential quad neighbor
          int n = nodetoquads[kp];

          // Don't count the same edge twice
          if (n == i){ continue; }

          // Flag to indicate that we have found the other edge (there
          // will only be at most one other match since this is planar)
          int quit = 0;

          // Search over all the edges on this quad, and see
          // if they match
          for ( int e = 0; e < 4; e++ ){  
            // Check if the adjacent edge matches in either direction
            if ((quads[4*i+enodes[j][0]] == quads[4*n+enodes[e][0]] &&
                 quads[4*i+enodes[j][1]] == quads[4*n+enodes[e][1]]) ||
                (quads[4*i+enodes[j][0]] == quads[4*n+enodes[e][1]] &&
                 quads[4*i+enodes[j][1]] == quads[4*n+enodes[e][0]])){
              // Label the other edge that shares this same node
              quadedgenums[4*n+e] = ne;
              quit = 1;
            }
          }
          if (quit){ break; }
        }

        // Increment the edge number
        ne++;
      }
    }
  }

  // Free the pointers
  delete [] ptr;
  delete [] nodetoquads;

  // Now we have a unique list of edge numbers and the total number of
  // edges, we can construct the unique edge list
  int *quadedges = new int[ 2*ne ];
  for ( int i = 0; i < nquads; i++ ){
    for ( int j = 0; j < 4; j++ ){
      // Get the unique edge number
      int n = quadedgenums[4*i+j];
      quadedges[2*n] = quads[4*i+enodes[j][0]];
      quadedges[2*n+1] = quads[4*i+enodes[j][1]];
    }
  } 

  // Free the data
  delete [] quadedgenums;

  *_nquadedges = ne;
  *_quadedges = quadedges;
}

/*
  Add the motion of the point in parameter space
*/
inline void addParamMovement( double alpha,
                              TMRPoint *Xu, TMRPoint *Xv,
                              TMRPoint *d, double *delta ){
  // Compute the sum - in the parameter space - of the motion 
  // contributed by each node
  double g11 = Xu->dot(Xu);
  double g12 = Xu->dot(Xv);
  double g22 = Xv->dot(Xv);

  // Compute the right hand sides
  double b1 = d->dot(Xu);
  double b2 = d->dot(Xv);

  // Compute the inverse of the determinant
  double invdet = alpha/(g11*g22 - g12*g12);

  // Add the contribution
  delta[0] += invdet*(g22*b1 - g12*b2);
  delta[1] += invdet*(g11*b2 - g12*b1);
}

/*
  Apply Laplacian smoothing
*/
void TMRTriangulation::laplacianSmoothing( int nsmooth,
                                           int num_edges,
                                           const int *edge_list,
                                           int num_pts,
                                           double *prm,
                                           TMRPoint *p ){
  int *count = new int[ num_pts ];
  double *new_params = new double[ 2*num_pts ];
  TMRPoint *Xu = new TMRPoint[ num_pts ];
  TMRPoint *Xv = new TMRPoint[ num_pts ];

  for ( int iter = 0; iter < nsmooth; iter++ ){
    memset(count, 0, num_pts*sizeof(int));
    memset(new_params, 0, 2*num_pts*sizeof(double));

    // Evaluate the derivatives w.r.t. the parameter locations
    for ( int i = 0; i < num_pts; i++ ){
      if (ptmarkers[i] == 0){
        surface->evalDeriv(prm[2*i], prm[2*i+1], &Xu[i], &Xv[i]);
      }
    }    

    // Loop over all the edges
    for ( int i = 0; i < num_edges; i++ ){
      int n1 = edge_list[2*i];
      int n2 = edge_list[2*i+1];

      // Compute the difference between the points along the
      // specified edge
      TMRPoint d;
      d.x = p[n2].x - p[n1].x;
      d.y = p[n2].y - p[n1].y;
      d.z = p[n2].z - p[n1].z;

      // Add the movement of the node in parameter space
      if (ptmarkers[n1] == 0){
        addParamMovement(1.0, &Xu[n1], &Xv[n1], &d, 
                         &new_params[2*n1]);
        count[n1]++;
      }

      // Add the movement of the second node in parameter space
      if (ptmarkers[n2] == 0){
        addParamMovement(-1.0, &Xu[n2], &Xv[n2], &d, 
                         &new_params[2*n2]);
        count[n2]++;
      }
    }

    // Set the locations for the new points, keep in place
    for ( int i = 0; i < npts; i++ ){
      if (ptmarkers[i] == 0 && count[i] > 0){
        prm[2*i] += new_params[2*i]/count[i];
        prm[2*i+1] += new_params[2*i+1]/count[i];
        surface->evalPoint(prm[2*i], prm[2*i+1], &p[i]);
      }
    }
  }

  delete [] count;
  delete [] new_params;
  delete [] Xu;
  delete [] Xv;
}

/*
  Apply the spring smoothing analogy
*/
void TMRTriangulation::springSmoothing( int nsmooth, double alpha,
                                        int num_edges, const int *edge_list,
                                        int num_pts, double *prm,
                                        TMRPoint *p ){
  double *len = new double[ num_edges ];
  double *new_params = new double[ 2*num_pts ];
  TMRPoint *Xu = new TMRPoint[ num_pts ];
  TMRPoint *Xv = new TMRPoint[ num_pts ];

  for ( int iter = 0; iter < nsmooth; iter++ ){
    double sum = 0.0;
    for ( int i = 0; i < num_edges; i++ ){
      int n1 = edge_list[2*i];
      int n2 = edge_list[2*i+1];

      // Compute the difference between the points along the
      // specified edge
      TMRPoint d;
      d.x = p[n2].x - p[n1].x;
      d.y = p[n2].y - p[n1].y;
      d.z = p[n2].z - p[n1].z;

      len[i] = sqrt(d.dot(d));
      sum += len[i];
    }
    double len0 = 0.9*sum/nedges;

    // Evaluate the derivatives w.r.t. the parameter locations
    for ( int i = 0; i < num_pts; i++ ){
      if (ptmarkers[i] == 0){
        surface->evalDeriv(prm[2*i], prm[2*i+1], &Xu[i], &Xv[i]);
      }
    }

    memset(new_params, 0, 2*num_pts*sizeof(double));

    for ( int i = 0; i < num_edges; i++ ){
      int n1 = edge_list[2*i];
      int n2 = edge_list[2*i+1];

      // Compute the difference between the points along the
      // specified edge
      TMRPoint d;
      d.x = p[n2].x - p[n1].x;
      d.y = p[n2].y - p[n1].y;
      d.z = p[n2].z - p[n1].z;

      double scale = (len0 - len[i])/len[i];

      // Add the movement of the node in parameter space
      if (ptmarkers[n1] == 0){
        addParamMovement(-scale, &Xu[n1], &Xv[n1], &d, 
                         &new_params[2*n1]);
      }

      // Add the movement of the second node in parameter space
      if (ptmarkers[n2] == 0){
        addParamMovement(scale, &Xu[n2], &Xv[n2], &d, 
                         &new_params[2*n2]);
      }
    }

    // Set the locations for the new points, keep in place
    for ( int i = 0; i < npts; i++ ){
      if (ptmarkers[i] == 0){
        prm[2*i] += alpha*new_params[2*i];
        prm[2*i+1] += alpha*new_params[2*i+1];
        surface->evalPoint(prm[2*i], prm[2*i+1], &p[i]);
      }
    }
  }

  delete [] new_params;
  delete [] len;
  delete [] Xu;
  delete [] Xv;
}

/*
  Perform spring smoothing specifically for a quadrilateral mesh.

  This code also connects springs across the faces of quadrilateral
  elements in an attempt to achieve higher mesh quality. This is 
  not always successful.
*/
void TMRTriangulation::springQuadSmoothing( int nsmooth,
                                            double alpha,
                                            int num_quads,
                                            const int *quad_list,
                                            int num_edges,
                                            const int *edge_list,
                                            int num_pts,
                                            double *prm,
                                            TMRPoint *p ){
  double *len = new double[ num_edges ];
  double *new_params = new double[ 2*num_pts ];
  TMRPoint *Xu = new TMRPoint[ num_pts ];
  TMRPoint *Xv = new TMRPoint[ num_pts ];

  // Compute the squrare root of 2
  const double sqrt2 = sqrt(2.0);

  for ( int iter = 0; iter < nsmooth; iter++ ){
    double sum = 0.0;
    for ( int i = 0; i < num_edges; i++ ){
      int n1 = edge_list[2*i];
      int n2 = edge_list[2*i+1];

      // Compute the difference between the points along the
      // specified edge
      TMRPoint d;
      d.x = p[n2].x - p[n1].x;
      d.y = p[n2].y - p[n1].y;
      d.z = p[n2].z - p[n1].z;

      len[i] = sqrt(d.dot(d));
      sum += len[i];
    }
    double len0 = 0.8*sum/nedges;

    // Evaluate the derivatives w.r.t. the parameter locations
    for ( int i = 0; i < num_pts; i++ ){
      if (ptmarkers[i] == 0){
        surface->evalDeriv(prm[2*i], prm[2*i+1], &Xu[i], &Xv[i]);
      }
    }

    memset(new_params, 0, 2*num_pts*sizeof(double));

    for ( int i = 0; i < num_edges; i++ ){
      int n1 = edge_list[2*i];
      int n2 = edge_list[2*i+1];

      // Compute the difference between the points along the
      // specified edge
      TMRPoint d;
      d.x = p[n2].x - p[n1].x;
      d.y = p[n2].y - p[n1].y;
      d.z = p[n2].z - p[n1].z;

      double scale = (len0 - len[i])/len[i];

      // Add the movement of the node in parameter space
      if (ptmarkers[n1] == 0){
        addParamMovement(-scale, &Xu[n1], &Xv[n1], &d, 
                         &new_params[2*n1]);
      }

      // Add the movement of the second node in parameter space
      if (ptmarkers[n2] == 0){
        addParamMovement(scale, &Xu[n2], &Xv[n2], &d, 
                         &new_params[2*n2]);
      }
    }

    for ( int i = 0; i < num_quads; i++ ){
      // Add the contribution from the first cross-quad member
      int n1 = quad_list[4*i];
      int n2 = quad_list[4*i+2];

      TMRPoint d;
      d.x = p[n2].x - p[n1].x;
      d.y = p[n2].y - p[n1].y;
      d.z = p[n2].z - p[n1].z;
      double ld = sqrt(d.dot(d));
      double scale = (sqrt2*len0 - ld)/ld;
      if (ptmarkers[n1] == 0){
        addParamMovement(-scale, &Xu[n1], &Xv[n1], &d, 
                         &new_params[2*n1]);
      }
      if (ptmarkers[n2] == 0){
        addParamMovement(scale, &Xu[n2], &Xv[n2], &d, 
                         &new_params[2*n2]);
      }

      // Add the contribution from the second cross-quad member
      n1 = quad_list[4*i+1];
      n2 = quad_list[4*i+3];

      d.x = p[n2].x - p[n1].x;
      d.y = p[n2].y - p[n1].y;
      d.z = p[n2].z - p[n1].z;
      ld = sqrt(d.dot(d));
      scale = (sqrt2*len0 - ld)/ld;
      if (ptmarkers[n1] == 0){
        addParamMovement(-scale, &Xu[n1], &Xv[n1], &d, 
                         &new_params[2*n1]);
      }
      if (ptmarkers[n2] == 0){
        addParamMovement(scale, &Xu[n2], &Xv[n2], &d, 
                         &new_params[2*n2]);
      }
    }

    // Set the locations for the new points, keep in place
    for ( int i = 0; i < npts; i++ ){
      if (ptmarkers[i] == 0){
        prm[2*i] += alpha*new_params[2*i];
        prm[2*i+1] += alpha*new_params[2*i+1];
        surface->evalPoint(prm[2*i], prm[2*i+1], &p[i]);
      }
    }
  }

  delete [] new_params;
  delete [] len;
  delete [] Xu;
  delete [] Xv;
}

/*
  Get the quadrilateral elemnet obtained by combining the triangles
  t1 and t2 together

        . --- .
      / |   /
    /   | /
  . --- . 
*/
int TMRTriangulation::getRecombinedQuad( const int tris[],
                                         const int trineighbors[],
                                         int t1, int t2,
                                         int quad[] ){
  int fail = 0;

  // Find the common edge between the two tirangles
  const int shift[3] = {2, 0, 1};
  int e1 = 0, e2 = 0;
  for ( ; e1 < 3; e1++ ){
    if (trineighbors[3*t1+shift[e1]] == t2){ break; }
  }
  for ( ; e2 < 3; e2++ ){
    if (trineighbors[3*t2+shift[e2]] == t1){ break; }
  }

  if (e1 >= 3 || e2 >= 3){
    fail = 1;
    return fail;
  }

  // Order the triangle, omitting the common edge
  for ( int j = 0, i = 0; i < 3; i++, j++ ){
    quad[j] = tris[3*t1+i];
    if (i == e1){
      j++;
    }
  }

  // Set the node contributed by the second triangle
  e2 += 2;
  if (e2 >= 3){ e2 -= 3; }
  quad[e1+1] = tris[3*t2+e2];

  return fail;
}

/*
  Compute the quality of a quadrilateral element
*/
double TMRTriangulation::computeQuadQuality( const int *quad,
                                             const TMRPoint *p ){
  // Compute the maximum of fabs(0.5*M_PI - alpha)
  double max_val = 0.0;

  for ( int k = 0; k < 4; k++ ){
    int prev = k-1;
    if (prev < 0){ prev = 3; }
    int next = k+1;
    if (next > 3){ next = 0; }

    TMRPoint a;
    a.x = p[quad[k]].x - p[quad[prev]].x;
    a.y = p[quad[k]].y - p[quad[prev]].y;
    a.z = p[quad[k]].z - p[quad[prev]].z;  
    
    TMRPoint b;
    b.x = p[quad[next]].x - p[quad[k]].x;
    b.y = p[quad[next]].y - p[quad[k]].y;
    b.z = p[quad[next]].z - p[quad[k]].z;
   
    // Compute the internal angle between the 
    double alpha = M_PI - acos(a.dot(b)/sqrt(a.dot(a)*b.dot(b)));
    double val = fabs(0.5*M_PI - alpha);
    if (val > max_val){
      max_val = val;
    }
  }

  // Compute the quality
  double eta = 1.0 - (2.0/M_PI)*max_val;
  if (eta < 0.0){ 
    eta = 0.0; 
  }

  return eta;
}

/*
  Compute the recombined quality
*/
double TMRTriangulation::computeRecombinedQuality( const int tris[], 
                                                   const int trineighbors[],
                                                   int t1, int t2,
                                                   const TMRPoint *p ){
  // Find the combined quadrilateral from the two given triangles
  int quad[4];
  int fail = getRecombinedQuad(tris, trineighbors, t1, t2, quad);
  if (fail){
    return 0.0;
  }

  return computeQuadQuality(quad, p);
}

/*
  Compute the quality of a quadrilateral element
*/
double TMRTriangulation::computeTriQuality( const int *tri,
                                            const TMRPoint *p ){
  // Compute the maximum of fabs(M_PI/3 - alpha)
  double max_val = 0.0;

  for ( int k = 0; k < 3; k++ ){
    int prev = k-1;
    if (prev < 0){ prev = 2; }
    int next = k+1;
    if (next > 2){ next = 0; }

    TMRPoint a;
    a.x = p[tri[k]].x - p[tri[prev]].x;
    a.y = p[tri[k]].y - p[tri[prev]].y;
    a.z = p[tri[k]].z - p[tri[prev]].z;  
    
    TMRPoint b;
    b.x = p[tri[next]].x - p[tri[k]].x;
    b.y = p[tri[next]].y - p[tri[k]].y;
    b.z = p[tri[next]].z - p[tri[k]].z;
   
    // Compute the internal angle between the 
    double alpha = M_PI - acos(a.dot(b)/sqrt(a.dot(a)*b.dot(b)));
    double val = fabs(M_PI/3.0 - alpha);
    if (val > max_val){
      max_val = val;
    }
  }

  // Compute the quality
  double eta = 1.0 - (3.0/M_PI)*max_val;
  if (eta < 0.0){ 
    eta = 0.0; 
  }

  return eta;
}
/*
  Recombine the triangulation into a quadrilateral mesh
*/
void TMRTriangulation::recombine(){
  // Count up the number of dual edges
  int ndualedges = 0;
  for ( int i = 0; i < nedges; i++ ){
    if (dualedges[2*i] >= 0 && dualedges[2*i+1] >= 0){
      ndualedges++;
    }
  }

  // Allocate the weights
  double *weights = new double[ ndualedges ];
  int *graphedges = new int[ 2*ndualedges ];

  // Compute the weight associated with each edge by combputing the
  // recombined quality
  for ( int i = 0, e = 0; i < nedges; i++ ){
    int t1 = dualedges[2*i];
    int t2 = dualedges[2*i+1];

    if (t1 >= 0 && t2 >= 0){
      // Compute the weight for this recombination
      double quality = computeRecombinedQuality(t1, t2, pts);
      double weight = (1.0 - quality)*(1.0 + 1.0/(quality + 0.01));
      graphedges[2*e] = t1;
      graphedges[2*e+1] = t2;
      weights[e] = weight;
      e++;
    }
  }

  int *match = new int[ 2*ntris ];
  TMR_PerfectMatchGraph(ntris, ndualedges, graphedges, weights, match);

  delete [] weights;
  delete [] graphedges;

  // Get the required options
  nquads = ntris/2;
  quads = new int[ 4*nquads ];

  // Recombine the quads
  for ( int i = 0; i < nquads; i++ ){
    int fail = getRecombinedQuad(match[2*i], match[2*i+1], &quads[4*i]);
    if (fail){
      printf("Recombined quadrilateral %d between triangles %d and %d failed\n", 
             i, match[2*i], match[2*i+1]);
    }
  }
  delete [] match;
}

/*
  Print the quad quality
*/
void TMRTriangulation::printQuadQuality(){
  const int nbins = 20;
  int total = 0;
  int bins[nbins];
  memset(bins, 0, nbins*sizeof(int));
  for ( int i = 0; i < nquads; i++ ){
    double quality = computeQuadQuality(&quads[4*i], pts);

    int k = 0;
    for ( ; k < nbins; k++ ){
      if (quality < 1.0*(k+1)/nbins){
        break;
      }
    }
    bins[k]++;
  }

  for ( int i = 0; i < nbins; i++ ){
    total += bins[i];
  }

  printf("Quality   # elements   percentage\n");
  for ( int k = 0; k < nbins; k++ ){
    printf("< %.2f    %10d   %10.3f\n",
           1.0*(k+1)/nbins, bins[k], 100.0*bins[k]/total);
  }
}

/*
  Print the triangle quality
*/
void TMRTriangulation::printTriQuality(){
  const int nbins = 20;
  int total = 0;
  int bins[nbins];
  memset(bins, 0, nbins*sizeof(int));
  for ( int i = 0; i < ntris; i++ ){
    double quality = computeTriQuality(&tris[3*i], pts);

    int k = 0;
    for ( ; k < nbins; k++ ){
      if (quality < 1.0*(k+1)/nbins){
        break;
      }
    }
    bins[k]++;
  }

  for ( int i = 0; i < nbins; i++ ){
    total += bins[i];
  }

  printf("Quality   # elements   percentage\n");
  for ( int k = 0; k < nbins; k++ ){
    printf("< %.2f    %10d   %10.3f\n",
           1.0*(k+1)/nbins, bins[k], 100.0*bins[k]/total);
  }
}

/*
  Write the quadrilateral mesh to a VTK file
*/
void TMRTriangulation::writeQuadToVTK( const char *filename ){
  FILE *fp = fopen(filename, "w");
  if (fp){
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "vtk output\nASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
      
    // Write out the points
    fprintf(fp, "POINTS %d float\n", npts);
    for ( int k = 0; k < npts; k++ ){
      fprintf(fp, "%e %e %e\n", pts[k].x, pts[k].y, pts[k].z);
    }
    
    // Write out the cell values
    fprintf(fp, "\nCELLS %d %d\n", nquads, 5*nquads);
    for ( int k = 0; k < nquads; k++ ){
      fprintf(fp, "4 %d %d %d %d\n", quads[4*k], quads[4*k+1], 
              quads[4*k+2], quads[4*k+3]);
    }

    // All quadrilaterals
    fprintf(fp, "\nCELL_TYPES %d\n", nquads);
    for ( int k = 0; k < nquads; k++ ){
      fprintf(fp, "%d\n", 9);
    }

    // Print out the rest as fields one-by-one
    fprintf(fp, "CELL_DATA %d\n", nquads);
    fprintf(fp, "SCALARS quality float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for ( int k = 0; k < nquads; k++ ){
      fprintf(fp, "%e\n", computeQuadQuality(&quads[4*k], pts));
    }

    fclose(fp);
  }
}

/*
  Write the output to a VTK file
*/
void TMRTriangulation::writeToVTK( const char *filename ){
  FILE *fp = fopen(filename, "w");
  if (fp){
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "vtk output\nASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
      
    // Write out the points
    fprintf(fp, "POINTS %d float\n", npts);
    for ( int k = 0; k < npts; k++ ){
      fprintf(fp, "%e %e %e\n", pts[k].x, pts[k].y, pts[k].z);
    }
    
    // Write out the cell values
    fprintf(fp, "\nCELLS %d %d\n", ntris, 4*ntris);
    for ( int k = 0; k < ntris; k++ ){
      fprintf(fp, "3 ");
      for ( int j = 0; j < 3; j++ ){
        int node = tris[3*k+j];
        fprintf(fp, "%d ", node);
      }
      fprintf(fp, "\n");
    }

    // All quadrilaterals
    fprintf(fp, "\nCELL_TYPES %d\n", ntris);
    for ( int k = 0; k < ntris; k++ ){
      fprintf(fp, "%d\n", 5);
    }

        // Print out the rest as fields one-by-one
    fprintf(fp, "CELL_DATA %d\n", ntris);
    fprintf(fp, "SCALARS quality float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for ( int k = 0; k < ntris; k++ ){
      fprintf(fp, "%e\n", computeTriQuality(&tris[3*k], pts));
    }

    fclose(fp);
  }
}
