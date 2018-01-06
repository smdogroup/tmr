#include "TMRMeshSmoothing.h"
#include <math.h>

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
void laplacianSmoothing( int nsmooth, int num_fixed_pts,
                         int num_edges, const int *edge_list,
                         int num_pts, double *prm, TMRPoint *p,
                         TMRFace *face ){
  int *count = new int[ num_pts ];
  double *new_params = new double[ 2*num_pts ];
  TMRPoint *Xu = new TMRPoint[ num_pts ];
  TMRPoint *Xv = new TMRPoint[ num_pts ];

  for ( int iter = 0; iter < nsmooth; iter++ ){
    memset(count, 0, num_pts*sizeof(int));
    memset(new_params, 0, 2*num_pts*sizeof(double));

    // Evaluate the derivatives w.r.t. the parameter locations
    for ( int i = num_fixed_pts; i < num_pts; i++ ){
      face->evalDeriv(prm[2*i], prm[2*i+1], &Xu[i], &Xv[i]);
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
      if (n1 >= num_fixed_pts){
        addParamMovement(1.0, &Xu[n1], &Xv[n1], &d, 
                         &new_params[2*n1]);
        count[n1]++;
      }

      // Add the movement of the second node in parameter space
      if (n2 >= num_fixed_pts){
        addParamMovement(-1.0, &Xu[n2], &Xv[n2], &d, 
                         &new_params[2*n2]);
        count[n2]++;
      }
    }

    // Set the locations for the new points, keep in place
    for ( int i = num_fixed_pts; i < num_pts; i++ ){
      if (count[i] > 0){
        prm[2*i] += new_params[2*i]/count[i];
        prm[2*i+1] += new_params[2*i+1]/count[i];
        face->evalPoint(prm[2*i], prm[2*i+1], &p[i]);
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
void springSmoothing( int nsmooth, double alpha, int num_fixed_pts,
                      int num_edges, const int *edge_list,
                      int num_pts, double *prm, TMRPoint *p,
                      TMRFace *face ){
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
    double len0 = 0.9*sum/num_edges;

    // Evaluate the derivatives w.r.t. the parameter locations
    for ( int i = num_fixed_pts; i < num_pts; i++ ){
      face->evalDeriv(prm[2*i], prm[2*i+1], &Xu[i], &Xv[i]);
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
      if (n1 >= num_fixed_pts){
        addParamMovement(-scale, &Xu[n1], &Xv[n1], &d, 
                         &new_params[2*n1]);
      }

      // Add the movement of the second node in parameter space
      if (n2 >= num_fixed_pts){
        addParamMovement(scale, &Xu[n2], &Xv[n2], &d, 
                         &new_params[2*n2]);
      }
    }

    // Set the locations for the new points, keep in place
    for ( int i = num_fixed_pts; i < num_pts; i++ ){
      prm[2*i] += alpha*new_params[2*i];
      prm[2*i+1] += alpha*new_params[2*i+1];
      face->evalPoint(prm[2*i], prm[2*i+1], &p[i]);
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
void springQuadSmoothing( int nsmooth, double alpha, int num_fixed_pts,
                          int num_quads, const int *quad_list,
                          int num_edges, const int *edge_list,
                          int num_pts, double *prm, TMRPoint *p,
                          TMRFace *face ){
  double *len = new double[ num_edges ];
  double *new_params = new double[ 2*num_pts ];
  TMRPoint *Xu = new TMRPoint[ num_pts ];
  TMRPoint *Xv = new TMRPoint[ num_pts ];

  // Compute the squrare root of 2
  const double sqrt2 = 1.1*sqrt(2.0);

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
    double len0 = sum/num_edges;

    // Evaluate the derivatives w.r.t. the parameter locations
    for ( int i = num_fixed_pts; i < num_pts; i++ ){
      face->evalDeriv(prm[2*i], prm[2*i+1], &Xu[i], &Xv[i]);
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
      if (n1 >= num_fixed_pts){
        addParamMovement(-scale, &Xu[n1], &Xv[n1], &d, 
                         &new_params[2*n1]);
      }

      // Add the movement of the second node in parameter space
      if (n2 >= num_fixed_pts){
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
      if (n1 >= num_fixed_pts){
        addParamMovement(-scale, &Xu[n1], &Xv[n1], &d, 
                         &new_params[2*n1]);
      }
      if (n2 >= num_fixed_pts){
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
      if (n1 >= num_fixed_pts){
        addParamMovement(-scale, &Xu[n1], &Xv[n1], &d, 
                         &new_params[2*n1]);
      }
      if (n2 >= num_fixed_pts){
        addParamMovement(scale, &Xu[n2], &Xv[n2], &d, 
                         &new_params[2*n2]);
      }
    }

    // Set the locations for the new points, keep in place
    for ( int i = num_fixed_pts; i < num_pts; i++ ){
      prm[2*i] += alpha*new_params[2*i];
      prm[2*i+1] += alpha*new_params[2*i+1];
      face->evalPoint(prm[2*i], prm[2*i+1], &p[i]);
    }
  }

  delete [] new_params;
  delete [] len;
  delete [] Xu;
  delete [] Xv;
}

/*
  Smooth the mesh
*/
void quadSmoothing( int nsmooth, int num_fixed_pts,
                    int num_pts, const int *ptr, const int *pts_to_quads,
                    int num_quads, const int *quads,
                    double *prm, TMRPoint *p,
                    TMRFace *face ){
  // Compute the min/max degree
  int min_degree = 0;
  int max_degree = 0;
  for ( int i = num_fixed_pts; i < num_pts; i++ ){
    int N = ptr[i+1] - ptr[i];
    if (i == num_fixed_pts){
      min_degree = N;
      max_degree = N; 
    }
    if (N < min_degree){
      min_degree = N;
    }
    if (N > max_degree){
      max_degree = N;
    }
  }

  for ( int iter = 0; iter < nsmooth; iter++ ){
    for ( int degree = min_degree; degree <= max_degree; degree++ ){
      for ( int i = num_fixed_pts; i < num_pts; i++ ){
        int N = ptr[i+1] - ptr[i];

        if (N == degree){
          // Evaluate the derivatives w.r.t. the parameter locations so
          // that we can take movement in the physical plane and convert
          // it to movement in the parametric coordinates
          TMRPoint Xu, Xv;
          face->evalDeriv(prm[2*i], prm[2*i+1], &Xu, &Xv);

          // Normalize the directions Xu, Xv to form a locally-orthonormal 
          // coordinate frame aligned with the surface
          TMRPoint xdir, ydir;

          // Normalize the x-direction
          double xnorm = sqrt(Xu.dot(Xu));
          xdir.x = Xu.x/xnorm;
          xdir.y = Xu.y/xnorm;
          xdir.z = Xu.z/xnorm;

          // Remove the component of the x-direction from Xv
          double dot = xdir.dot(Xv);
          ydir.x = Xv.x - dot*xdir.x;
          ydir.y = Xv.y - dot*xdir.y;
          ydir.z = Xv.z - dot*xdir.z;

          double ynorm = sqrt(ydir.dot(ydir));
          ydir.x = ydir.x/ynorm;
          ydir.y = ydir.y/ynorm;
          ydir.z = ydir.z/ynorm;

          if (N > 0){
            // Loop over the quadrilaterals that reference this point 
            double A = 0.0, B = 0.0;  
            for ( int qp = ptr[i]; qp < ptr[i+1]; qp++ ){
              const int *quad = &quads[4*pts_to_quads[qp]];

              // Pick out the influence triangle points from the quadrilateral 
              // This consists of the base point i and the following two
              int ijk[3];
              if (quad[0] == i){
                ijk[0] = quad[0];  ijk[1] = quad[1];  ijk[2] = quad[3];
              }
              else if (quad[1] == i){
                ijk[0] = quad[1];  ijk[1] = quad[2];  ijk[2] = quad[0];
              }
              else if (quad[2] == i){
                ijk[0] = quad[2];  ijk[1] = quad[3];  ijk[2] = quad[1];
              }
              else {
                ijk[0] = quad[3];  ijk[1] = quad[0];  ijk[2] = quad[2];
              }

              // Now compute the geometric quantities
              // p = yj - yk, q = xk - xj
              double xi = xdir.dot(p[ijk[0]]);
              double yi = ydir.dot(p[ijk[0]]);
              double xj = xdir.dot(p[ijk[1]]);
              double yj = ydir.dot(p[ijk[1]]);
              double xk = xdir.dot(p[ijk[2]]);
              double yk = ydir.dot(p[ijk[2]]);
              double p = yj - yk;
              double q = xk - xj;
              double r = xj*yk - xk*yj;
              double a = 0.5*(p*xi + q*yi + r);
              double b = sqrt(p*p + q*q); 
              A += a;
              B += b;
            }

            double hbar = 2.0*A/B;
            double bbar = B/N;

            // Set the weights
            double w1 = 1.0/(hbar*hbar);
            double w2 = 4.0/(bbar*bbar);

            // The parameters for the Jacobian/right-hand-side
            double s1 = 0.0, s2 = 0.0, s3 = 0.0, s4 = 0.0, s5 = 0.0;

            for ( int qp = ptr[i]; qp < ptr[i+1]; qp++ ){
              const int *quad = &quads[4*pts_to_quads[qp]];

              // Pick out the influence triangle points from the quadrilateral 
              // This consists of the base point i and the following two
              int ijk[3];
              if (quad[0] == i){
                ijk[0] = quad[0];  ijk[1] = quad[1];  ijk[2] = quad[3];
              }
              else if (quad[1] == i){
                ijk[0] = quad[1];  ijk[1] = quad[2];  ijk[2] = quad[0];
              }
              else if (quad[2] == i){
                ijk[0] = quad[2];  ijk[1] = quad[3];  ijk[2] = quad[1];
              }
              else {
                ijk[0] = quad[3];  ijk[1] = quad[0];  ijk[2] = quad[2];
              }

              // Now compute the geometric quantities
              // p = yj - yk, q = xk - xj
              double xi = xdir.dot(p[ijk[0]]);
              double yi = ydir.dot(p[ijk[0]]);
              double xj = xdir.dot(p[ijk[1]]);
              double yj = ydir.dot(p[ijk[1]]);
              double xk = xdir.dot(p[ijk[2]]);
              double yk = ydir.dot(p[ijk[2]]);
              double p = yj - yk;
              double q = xk - xj;
              double r = xj*yk - xk*yj;
              double a = 0.5*(p*xi + q*yi + r);
              double b = sqrt(p*p + q*q); 

              // Other quantities derived from the in-plane triangle data
              double xm = 0.5*(xj + xk);
              double ym = 0.5*(yj + yk);
              double binv2 = 1.0/(b*b);

              // Sum up the contributions to the s terms
              s1 += binv2*(w1*p*p + w2*q*q);
              s2 += binv2*p*q*(w1 - w2);
              s3 += binv2*(w1*p*(hbar*b - 2*a) - 
                           w2*q*((xi - xm)*q - (yi - ym)*p));
              s4 += binv2*(w1*q*q + w2*p*p);
              s5 += binv2*(w1*q*(hbar*b - 2*a) - 
                           w2*p*((yi - ym)*p - (xi - xm)*q));
            }

            // Compute the updates in the physical plane
            double det = s1*s4 - s2*s2;
            double lx = 0.0, ly = 0.0;
            if (det != 0.0){
              det = 1.0/det;
              lx = det*(s3*s4 - s2*s5);
              ly = det*(s1*s5 - s2*s3);
            }

            // Check that the requested move direction is well-defined
            if (lx == lx && ly == ly){
              // Add up the displacements along the local coordinate directions
              TMRPoint dir;
              dir.x = lx*xdir.x + ly*ydir.x;
              dir.y = lx*xdir.y + ly*ydir.y;
              dir.z = lx*xdir.z + ly*ydir.z;

              // Add the parameter movement along the specified direction
              // and compute the update
              addParamMovement(1.0, &Xu, &Xv, &dir, &prm[2*i]);
              face->evalPoint(prm[2*i], prm[2*i+1], &p[i]);
            }
          }
        }
      }
    }
  }
}
