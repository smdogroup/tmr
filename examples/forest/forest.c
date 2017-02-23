#include "TMRQuadForest.h"
#include "specialFSDTStiffness.h"
#include "TACSAssembler.h"
#include "PlaneStressQuad.h"
#include "TACSToFH5.h"
#include "TACSShellTraction.h"
#include "TensorToolbox.h"
#include "MITCShell.h"
#include "KSFailure.h"
#include "Compliance.h"
#include "TACSMeshLoader.h"
#include "tacslapack.h"

// Set the element order = 3rd order, quadratic elements
const int ELEMENT_ORDER = 3;

class ProblemGeometry {
 public:
  virtual void getConnectivity( int *num_nodes, 
                                const int **conn, int *num_faces ) = 0;
  virtual void getLocation( int face, int x, int y,
                            TacsScalar *X ) = 0; 
  virtual int isBoundary( int face, int x, int y ) = 0;
};

class SquarePlate : public ProblemGeometry {
 public:
  SquarePlate(){}

  void getConnectivity( int *_num_nodes,
                        const int **_conn,
                        int *_num_faces ){
    *_num_nodes = 8;
    *_num_faces = 5;
    *_conn = conn;
  }
  void getLocation( int face, int x, int y,
                    TacsScalar *X ){
    // Compute the u/v locations on the face
    const double dh = 1.0/(1 << TMR_MAX_LEVEL);
    double u = dh*x, v = dh*y;

    // Evaluate the shape functions
    double N[4];
    N[0] = (1.0 - u)*(1.0 - v);
    N[1] = u*(1.0 - v);
    N[2] = (1.0 - u)*v;
    N[3] = u*v;

    // Compute the node location
    X[0] = (N[0]*Xpts[3*conn[4*face]] +
            N[1]*Xpts[3*conn[4*face+1]] +
            N[2]*Xpts[3*conn[4*face+2]] +
            N[3]*Xpts[3*conn[4*face+3]]);
    X[1] = (N[0]*Xpts[3*conn[4*face]+1] +
            N[1]*Xpts[3*conn[4*face+1]+1] +
            N[2]*Xpts[3*conn[4*face+2]+1] +
            N[3]*Xpts[3*conn[4*face+3]+1]);
    X[2] = (N[0]*Xpts[3*conn[4*face]+2] +
            N[1]*Xpts[3*conn[4*face+1]+2] +
            N[2]*Xpts[3*conn[4*face+2]+2] +
            N[3]*Xpts[3*conn[4*face+3]+2]);
  }
  int isBoundary( int face, int x, int y ){
    const int hmax = 1 << TMR_MAX_LEVEL;
    if (face == 0 && y == 0){ return 1; }
    if (face == 2 && y == 0){ return 1; }
    if (face == 3 && y == hmax){ return 1; }
    if (face == 4 && x == hmax){ return 1; }
    return 0;
  }

 private:
  static const double Xpts[24];
  static const int conn[20];
};

const double SquarePlate::Xpts[] = {0.0, 0.0, 0.0,
                                    1.0, 0.0, 0.0,
                                    0.3, 0.7, 0.0,
                                    0.8, 0.25, 0.0,
                                    0.25, 0.2, 0.0,
                                    0.75, 0.6, 0.0,
                                    0.0, 1.0, 0.0,
                                    1.0, 1.0, 0.0};

const int SquarePlate::conn[] = {0, 1, 4, 3,
                                 2, 4, 5, 3, 
                                 6, 0, 2, 4, 
                                 2, 5, 6, 7,
                                 3, 1, 5, 7};


class MeshProblem : public ProblemGeometry {
 public:
  MeshProblem( int _num_nodes, int _num_faces,
               const int *elem_conn, 
               const TacsScalar *_Xpts ){
    num_nodes = _num_nodes;
    num_faces = _num_faces;

    conn = new int[ 4*num_faces ];
    memcpy(conn, elem_conn, 4*num_faces*sizeof(int));

    Xpts = new TacsScalar[ 3*num_nodes ];
    memcpy(Xpts, _Xpts, 3*num_nodes*sizeof(TacsScalar));
  }
  ~MeshProblem(){
    delete [] conn;
    delete [] Xpts;
  }

  void getConnectivity( int *_num_nodes,
                        const int **_conn,
                        int *_num_faces ){
    *_num_nodes = num_nodes;
    *_num_faces = num_faces;
    *_conn = conn;
  }
  void getLocation( int face, int x, int y,
                    TacsScalar *X ){
    // Compute the u/v locations on the face
    const double dh = 1.0/(1 << TMR_MAX_LEVEL);
    double u = dh*x, v = dh*y;

    // Evaluate the shape functions
    double N[4];
    N[0] = (1.0 - u)*(1.0 - v);
    N[1] = u*(1.0 - v);
    N[2] = (1.0 - u)*v;
    N[3] = u*v;

    // Compute the node location
    X[0] = (N[0]*Xpts[3*conn[4*face]] +
            N[1]*Xpts[3*conn[4*face+1]] +
            N[2]*Xpts[3*conn[4*face+2]] +
            N[3]*Xpts[3*conn[4*face+3]]);
    X[1] = (N[0]*Xpts[3*conn[4*face]+1] +
            N[1]*Xpts[3*conn[4*face+1]+1] +
            N[2]*Xpts[3*conn[4*face+2]+1] +
            N[3]*Xpts[3*conn[4*face+3]+1]);
    X[2] = (N[0]*Xpts[3*conn[4*face]+2] +
            N[1]*Xpts[3*conn[4*face+1]+2] +
            N[2]*Xpts[3*conn[4*face+2]+2] +
            N[3]*Xpts[3*conn[4*face+3]+2]);
  }
  int isBoundary( int face, int x, int y ){
    TacsScalar X[3];
    getLocation(face, x, y, X);
    if (X[1] < 0.0015){
      return 1;
    }

    return 0;
  }

 private:
  TacsScalar *Xpts;
  int *conn;
  int num_faces, num_nodes;
};

/*
  Compute the Jacobian transformation at a point within the element
*/
void computeJacobianTrans( const TacsScalar Xpts[],
                           const double Na[], const double Nb[], 
                           TacsScalar Xd[], TacsScalar J[] ){
  memset(Xd, 0, 9*sizeof(TacsScalar));

  // Compute the derivative along the coordinate directions
  const double *na = Na, *nb = Nb;
  const TacsScalar *X = Xpts;
  for ( int i = 0; i < 9; i++ ){
    Xd[0] += X[0]*na[0];
    Xd[1] += X[1]*na[0];
    Xd[2] += X[2]*na[0];

    Xd[3] += X[0]*nb[0];
    Xd[4] += X[1]*nb[0];
    Xd[5] += X[2]*nb[0];

    na++;  nb++;
    X += 3;
  }

  // Compute the cross-product with the normal
  Tensor::crossProduct3D(&Xd[6], &Xd[0], &Xd[3]);
  Tensor::normalize3D(&Xd[6]);

  // Compute the transpose of the Jacobian transformation
  FElibrary::jacobian3d(Xd, J);
}

/*
  Evaluate the enrichment functions
*/
void evalEnrichmentFuncs( const double pt[], double N[] ){
  // Compute the cubic enrichment shape functions along the two
  // coordinate directions
  double ca = (1.0 + pt[0])*pt[0]*(1.0 - pt[0]);
  double cb = (1.0 + pt[1])*pt[1]*(1.0 - pt[1]);

  // Set the shape functions themselves
  N[0] = ca;
  N[1] = pt[1]*ca;
  N[2] = pt[1]*pt[1]*ca;
  N[3] = cb;
  N[4] = pt[0]*cb;
  N[5] = pt[0]*pt[0]*cb;
  N[6] = ca*cb;
}

/*
  The enriched shape functions for the tensor-product quadratic
  Lagrange shape functions.
*/
void evalEnrichmentFuncs( const double pt[],
                          double N[], double Na[], double Nb[] ){
  // Compute the cubic enrichment shape functions along the two
  // coordinate directions
  double ca = (1.0 + pt[0])*pt[0]*(1.0 - pt[0]);
  double cb = (1.0 + pt[1])*pt[1]*(1.0 - pt[1]);

  // Compute the derivative of the enrichment functions w.r.t. 
  // the two coordinate directions
  double da = 1.0 - 3*pt[0]*pt[0];
  double db = 1.0 - 3*pt[1]*pt[1];

  // Set the shape functions themselves
  N[0] = ca;
  N[1] = pt[1]*ca;
  N[2] = pt[1]*pt[1]*ca;
  N[3] = cb;
  N[4] = pt[0]*cb;
  N[5] = pt[0]*pt[0]*cb;
  N[6] = ca*cb;

  // Set the derivatives of the enrichment functions with respect to
  // the first and second coordinate directions
  Na[0] = da;
  Na[1] = pt[1]*da;
  Na[2] = pt[1]*pt[1]*da;
  Na[3] = 0.0;
  Na[4] = cb;
  Na[5] = 2.0*pt[0]*cb;
  Na[6] = da*cb;
  
  Nb[0] = 0.0;
  Nb[1] = ca;
  Nb[2] = 2.0*pt[1]*ca;
  Nb[3] = db;
  Nb[4] = pt[0]*db;
  Nb[5] = pt[0]*pt[0]*db;
  Nb[6] = ca*db;
}

/*
  Given the values of the derivatives of one component of the
  displacement or rotation field at the nodes, compute the
  reconstruction over the element by solving a least-squares problem

  input:
  Xpts:    the element node locations
  uvals:   the solution at the nodes
  uderiv:  the derivative of the solution in x/y/z at the nodes
  
  output:
  ubar:    the values of the coefficients on the enrichment functions
*/
void computeElemRecon( const TacsScalar Xpts[],
                       const TacsScalar uvals[],
                       const TacsScalar uderiv[],
                       TacsScalar ubar[] ){
  // Set up the least squares problem at the nodes
  TacsScalar A[18*7], b[18*6];
  int nrhs = 6;

  memset(b, 0, sizeof(b));
  memset(A, 0, sizeof(A));

  const double wvals[] = {0.5, 1.0, 0.5};

  for ( int c = 0, jj = 0; jj < 3; jj++ ){
    for ( int ii = 0; ii < 3; ii++, c += 2 ){
      // Set the parametric location within the element
      double pt[2];
      pt[0] = -1.0 + 1.0*ii;
      pt[1] = -1.0 + 1.0*jj;

      // Compute the quadratic shape functions at this point
      double N[9], Na[9], Nb[9];
      FElibrary::biLagrangeSF(N, Na, Nb, pt, 3);
      
      // Evaluate the Jacobian transformation at this point
      TacsScalar Xd[9], J[9];
      computeJacobianTrans(Xpts, Na, Nb, Xd, J);

      // Normalize the first direction
      TacsScalar d1[3], d2[3];
      d1[0] = Xd[0];  d1[1] = Xd[1];  d1[2] = Xd[2];
      Tensor::normalize3D(d1);

      // Compute d2 = n x d1
      Tensor::crossProduct3D(d2, &Xd[6], d1);

      // First, compute the contributions to the righ-hand-side. The
      // right vector contains the difference between the prescribed
      // derivative and the contribution to the derivative from the
      // quadratic shape function terms
      const TacsScalar *ud = &uderiv[18*(ii + 3*jj)];
      for ( int k = 0; k < nrhs; k++ ){
        b[18*k+c] = 
          wvals[ii]*wvals[jj]*(d1[0]*ud[0] + d1[1]*ud[1] + d1[2]*ud[2]);
        b[18*k+c+1] = 
          wvals[ii]*wvals[jj]*(d2[0]*ud[0] + d2[1]*ud[1] + d2[2]*ud[2]);
        ud += 3;
      }

      // Compute the derivatives from the interpolated solution
      TacsScalar Ud[6*2];
      memset(Ud, 0, sizeof(Ud));
      for ( int i = 0; i < 9; i++ ){
        for ( int k = 0; k < 6; k++ ){
          Ud[2*k] += uvals[6*i + k]*Na[i];
          Ud[2*k+1] += uvals[6*i + k]*Nb[i];
        }
      }
      
      // Set the right-hand side
      for ( int k = 0; k < nrhs; k++ ){
        TacsScalar d[3];
        d[0] = Ud[2*k]*J[0] + Ud[2*k+1]*J[1];
        d[1] = Ud[2*k]*J[3] + Ud[2*k+1]*J[4];
        d[2] = Ud[2*k]*J[6] + Ud[2*k+1]*J[7];

        b[18*k+c] -= 
          wvals[ii]*wvals[jj]*(d1[0]*d[0] + d1[1]*d[1] + d1[2]*d[2]);
        b[18*k+c+1] -= 
          wvals[ii]*wvals[jj]*(d2[0]*d[0] + d2[1]*d[1] + d2[2]*d[2]);
      }
      
      // Now, evaluate the terms for the left-hand-side
      // that contribute to the
      double Nr[7], Nar[7], Nbr[7];
      evalEnrichmentFuncs(pt, Nr, Nar, Nbr);

      // Add the contributions to the the enricment 
      for ( int i = 0; i < 7; i++ ){
        // Evaluate the
        TacsScalar d[3];
        d[0] = Nar[i]*J[0] + Nbr[i]*J[1];
        d[1] = Nar[i]*J[3] + Nbr[i]*J[4];
        d[2] = Nar[i]*J[6] + Nbr[i]*J[7];

        A[18*i+c] = 
          wvals[ii]*wvals[jj]*(d1[0]*d[0] + d1[1]*d[1] + d1[2]*d[2]);
        A[18*i+c+1] = 
          wvals[ii]*wvals[jj]*(d2[0]*d[0] + d2[1]*d[1] + d2[2]*d[2]);
      }
    }
  }

  // Singular values
  TacsScalar s[7];
  int m = 18, n = 7;
  double rcond = -1.0;
  int rank;
  int lwork = 10*18;
  TacsScalar work[10*18];
  int info;

  // Using LAPACK, compute the least squares solution
  LAPACKdgelss(&m, &n, &nrhs, A, &m, b, &m, s, 
               &rcond, &rank, work, &lwork, &info);

  // Copy over the ubar solution
  for ( int i = 0; i < n; i++ ){
    for ( int j = 0; j < nrhs; j++ ){
      ubar[nrhs*i + j] = b[m*j + i];
    }
  }
}

/*
  Compute the local derivative weights
*/
void computeLocalWeights( TACSAssembler *tacs, 
                          TACSBVec **wlocal ){
  // Create the weight vector - the weights on the averages
  // Set the weights on each of the nodes
  TACSBVec *weights = new TACSBVec(tacs->getVarMap(), 1,
                                   NULL, tacs->getBVecDistribute(),
                                   tacs->getBVecDepNodes());
  weights->incref();

  // Set the local element weights - this assumes that all the
  // elements are bi-quadratic 9-noded elements
  TacsScalar welem[9]; 
  for ( int i = 0; i < 9; i++ ){
    welem[i] = 1.0;
  }

  // Add unit weights to all the elements. This will sum up
  // the number of times each node is referenced
  int nelems = tacs->getNumElements();
  for ( int i = 0; i < nelems; i++ ){
    int len = 0;
    const int *nodes;
    tacs->getElement(i, &nodes, &len);
    weights->setValues(len, nodes, welem, ADD_VALUES);
  }

  // Finish setting all of the values
  weights->beginSetValues(ADD_VALUES);
  weights->endSetValues(ADD_VALUES);

  // Distribute the values
  weights->beginDistributeValues();
  weights->endDistributeValues();

  // Return the local weight values
  *wlocal = weights;
}

/*
  Given the input solution, compute and set the derivatives 
  in the vector ux. Note that this vector must be 3*

  input:
  TACSAssembler:  the TACSAssembler object
  uvec:           the solution vector
  wlocal:         the local element-weight values

  output:
  uderiv:         the approximate derivatives at the nodes
*/
void computeNodeDeriv( TACSAssembler *tacs, TACSBVec *uvec, 
                       TACSBVec *wlocal, TACSBVec **_uderiv ){
  // Allocate a vector for the derivatives
  TACSBVec *uderiv = 
    new TACSBVec(tacs->getVarMap(), 3*tacs->getVarsPerNode(), NULL,
                 tacs->getBVecDistribute(), tacs->getBVecDepNodes());

  // The number of variables at each node
  int vars_per_node = tacs->getVarsPerNode();
  int deriv_per_node = 3*vars_per_node;

  // Get the number of elements
  int nelems = tacs->getNumElements();

  // Allocate space for the element-wise values and derivatives
  TacsScalar *uelem = new TacsScalar[ 9*vars_per_node ];
  TacsScalar *delem = new TacsScalar[ 9*deriv_per_node ];

  // Perform the reconstruction for the local
  for ( int elem = 0; elem < nelems; elem++ ){
    // Get the element nodes
    int len = 0;
    const int *nodes;
    tacs->getElement(elem, &nodes, &len);

    // Get the local weight values
    TacsScalar welem[9];
    wlocal->getValues(len, nodes, welem);

    // Get the local element adjoint variables
    uvec->getValues(len, nodes, uelem);

    // Get the node locations for the element
    TacsScalar Xpts[3*9];
    tacs->getElement(elem, Xpts);
    
    // Compute the derivative of the 6 components of the adjoint
    // variables along the 3-coordinate directions
    TacsScalar *d = delem;

    // Compute the contributions to the derivative from this side of
    // the element    
    for ( int jj = 0; jj < 3; jj++ ){
      for ( int ii = 0; ii < 3; ii++ ){
        double pt[2];
        pt[0] = -1.0 + 1.0*ii;
        pt[1] = -1.0 + 1.0*jj;

        // Evaluate the the quadratic shape functions at this point
        double N[9], Na[9], Nb[9];
        FElibrary::biLagrangeSF(N, Na, Nb, pt, 3);
      
        // Evaluate the Jacobian transformation at this point
        TacsScalar Xd[9], J[9];
        computeJacobianTrans(Xpts, Na, Nb, Xd, J);

        // Compute the derivatives from the interpolated solution
        TacsScalar Ud[6*2];
        memset(Ud, 0, sizeof(Ud));
        for ( int i = 0; i < 9; i++ ){
          for ( int k = 0; k < 6; k++ ){
            Ud[2*k] += uelem[6*i + k]*Na[i];
            Ud[2*k+1] += uelem[6*i + k]*Nb[i];
          }
        }

        // Evaluate the x/y/z derivatives of each value at the
        // coordinate locations
        TacsScalar winv = 1.0/welem[3*jj + ii];
        for ( int k = 0; k < vars_per_node; k++ ){
          d[0] = winv*(Ud[2*k]*J[0] + Ud[2*k+1]*J[1]);
          d[1] = winv*(Ud[2*k]*J[3] + Ud[2*k+1]*J[4]);
          d[2] = winv*(Ud[2*k]*J[6] + Ud[2*k+1]*J[7]);
          d += 3;
        }
      }
    }

    // Add the values
    uderiv->setValues(len, nodes, delem, ADD_VALUES);
  }

  // Free the element values
  delete [] uelem;
  delete [] delem;

  // Add the values in parallel
  uderiv->beginSetValues(ADD_VALUES);
  uderiv->endSetValues(ADD_VALUES);

  // Distribute the values so that we can call getValues
  uderiv->beginDistributeValues();
  uderiv->endDistributeValues();

  *_uderiv = uderiv;
}

/*
  Refine the mesh with the given refinement specified on each
  processor

  This function takes the relative refinement on each local processor,
  and reduces this to the root processor. The code then computes the
  local elements which need to be refined and then applies the
  refinement to the quadtree forest. This quadtree forest can then
  be used for creating the next mesh.

  input:
  forest:      the quadtree forest
  refine:      the refinement flag for each element
  min_refine:  the minimum quadrant refinement level
  max_refine:  the maximum quadrant refinement level
*/
void refineQuadMesh( TMRQuadForest *forest, const int *refine,
                     const int min_refine, const int max_refine ){
  // Set the refinement on each of the quadtrees
  TMRQuadtree **trees;
  int ntrees = forest->getQuadtrees(&trees);
    
  // Loop through all trees and refine them
  for ( int i = 0, offset = 0; i < ntrees; i++ ){
    if (trees[i]){
      int nlocal = trees[i]->getNumElements();
      trees[i]->refine(&refine[offset], min_refine, max_refine);
      offset += nlocal;
    }
  }
}

/*
  The following function performs a mesh refinement based on a strain
  energy criteria. It is based on the following relationship for
  linear finite-element analysis

  a(u-uh,u-uh) = a(u,u) - a(uh,uh) 

  where a(u,u) is the bilinear strain energy functional, u is the
  exact solution, and uh is the discretized solution at any mesh
  level. This relies on the relationship that a(uh, u - uh) = 0 which
  is satisfied due to the method of Galerkin/Ritz.

  The following function computes a localized error indicator using
  the element-wise strain energy. The code computes a higher-order
  reconstructed solution using a cubic enrichment functions. These
  enrichment functions expand the original solution space and are
  computed based on a least-squares approximation with nodal gradient
  values. The localized error indicator is evaluated as follows:

  err = [sum_{i=1}^{4} ae(uCe, uCe) ] - ae(ue, ue)

  where uCe is the element-wise cuibc element reconstruction projected
  onto a uniformly refined mesh.

  input:
  tacs:        the TACSAssembler object
  uvec:        the solution variables in the mesh
  partition:   the element partition
  forest:      the forest of quadtrees 
  target_err:  target absolute error value
  min_refine:  minimum refinement
  max_refnie:  maximum refinement

  returns:     the strain energy error
*/
TacsScalar strainEnergyRefine( TACSAssembler *tacs,
                               TACSBVec *uvec,
                               TMRQuadForest *forest,
                               double target_err,
                               int min_refine, int max_refine ){
  // Get the communicator
  MPI_Comm comm = tacs->getMPIComm();
  
  // Ensure that the local variables have been set
  tacs->setVariables(uvec);

  // Perform a local refinement of the nodes based on the strain energy
  // within each element
  int nelems = tacs->getNumElements();
  TacsScalar *SE_error = new TacsScalar[ nelems ];
  
  // Compute the reconstructed solution
  TACSBVec *wlocal;
  computeLocalWeights(tacs, &wlocal);
  wlocal->incref();

  // Compute the nodal derivatives
  TACSBVec *uderiv;
  computeNodeDeriv(tacs, uvec, wlocal, &uderiv);
  uderiv->incref();

  // Zero the time-derivatives: this assumes a steady-state
  TacsScalar dvars[6*9], ddvars[6*9];
  memset(dvars, 0, sizeof(dvars));
  memset(ddvars, 0, sizeof(ddvars));

  // Keep track of the total error
  TacsScalar SE_total_error = 0.0;

  // For each element in the mesh, compute the original strain energy
  for ( int i = 0; i < nelems; i++ ){
    // Set the simulation time
    double time = 0.0;

    // Get the node locations and variables
    TacsScalar Xpts[3*9], uelem[6*9];
    TACSElement *elem = tacs->getElement(i, Xpts, uelem, NULL, NULL);
    
    // Get the nodes from TACSAssembler
    int len;
    const int *nodes;
    tacs->getElement(i, &nodes, &len);

    // Evaluate the element residual
    TacsScalar res[6*9];
    memset(res, 0, elem->numVariables()*sizeof(TacsScalar));
    elem->addResidual(time, res, Xpts, uelem, dvars, ddvars);
    
    // Take the inner product to get the strain energy
    SE_error[i] = 0.0;
    for ( int j = 0; j < elem->numVariables(); j++ ){
      SE_error[i] += res[j]*uelem[j];
    }

    // Compute the solution on the refined mesh
    TacsScalar delem[18*9];
    uderiv->getValues(len, nodes, delem);

    // 7 enrichment functions for each degree of freedom
    TacsScalar ubar[7*6];
    computeElemRecon(Xpts, uelem, delem, ubar);

    TacsScalar SE_refine = 0.0;
    for ( int ii = 0; ii < 2; ii++ ){
      for ( int jj = 0; jj < 2; jj++ ){
        TacsScalar rXpts[3*9], ruelem[6*9];
        memset(rXpts, 0, 3*9*sizeof(TacsScalar));
        memset(ruelem, 0, 6*9*sizeof(TacsScalar));

        for ( int m = 0; m < 3; m++ ){
          for ( int n = 0; n < 3; n++ ){
            double pt[2];
            pt[0] = -1.0 + 0.5*(2*ii + n);
            pt[1] = -1.0 + 0.5*(2*jj + m);

            // Evaluate the locations of the new nodes
            double N[9], Nr[7];
            FElibrary::biLagrangeSF(N, pt, 3);
            evalEnrichmentFuncs(pt, Nr);
            
            // Set the values of the variables at this point
            for ( int k = 0; k < 9; k++ ){
              rXpts[3*(n + 3*m)] += Xpts[3*k]*N[k];
              rXpts[3*(n + 3*m)+1] += Xpts[3*k+1]*N[k];
              rXpts[3*(n + 3*m)+2] += Xpts[3*k+2]*N[k];

              // Evaluate the interpolation part of the reconstruction
              ruelem[6*(n + 3*m)] += uelem[6*k]*N[k];
              ruelem[6*(n + 3*m)+1] += uelem[6*k+1]*N[k];
              ruelem[6*(n + 3*m)+2] += uelem[6*k+2]*N[k];
              ruelem[6*(n + 3*m)+3] += uelem[6*k+3]*N[k];
              ruelem[6*(n + 3*m)+4] += uelem[6*k+4]*N[k];
              ruelem[6*(n + 3*m)+5] += uelem[6*k+5]*N[k];
            }

            // Add the portion from the enrichment functions
            for ( int k = 0; k < 7; k++ ){
              // Evaluate the interpolation part of the reconstruction
              ruelem[6*(n + 3*m)] += ubar[6*k]*Nr[k];
              ruelem[6*(n + 3*m)+1] += ubar[6*k+1]*Nr[k];
              ruelem[6*(n + 3*m)+2] += ubar[6*k+2]*Nr[k];
              ruelem[6*(n + 3*m)+3] += ubar[6*k+3]*Nr[k];
              ruelem[6*(n + 3*m)+4] += ubar[6*k+4]*Nr[k];
              ruelem[6*(n + 3*m)+5] += ubar[6*k+5]*Nr[k];
            }
          }
        }

        // Compute the element residual
        memset(res, 0, elem->numVariables()*sizeof(TacsScalar));
        elem->addResidual(time, res, rXpts, ruelem, dvars, ddvars);
        
        // Take the inner product to get the strain energy
        for ( int j = 0; j < elem->numVariables(); j++ ){
          SE_refine += res[j]*ruelem[j];
        }
      }
    }
    // SE_refine - SE_error should always be a positive quantity
    SE_error[i] = SE_refine - SE_error[i];

    // Add up th etotal error
    SE_total_error += SE_error[i];
  }

  // Count up the total strain energy 
  TacsScalar SE_temp = 0.0;
  MPI_Allreduce(&SE_total_error, &SE_temp, 1, TACS_MPI_TYPE, MPI_SUM, comm);
  SE_total_error = SE_temp;

  // Go through and flag which element should be refined
  int *refine_local = new int[ nelems ];
  memset(refine_local, 0, nelems*sizeof(int));
  for ( int i = 0; i < nelems; i++ ){
    if (SE_error[i] >= target_err){
      refine_local[i] = 1;
    }
  }

  // Free some of the data that is no longer required
  delete [] SE_error;
  uderiv->decref();
  wlocal->decref();

  // refine the quadrant mesh based on the local refinement values
  refineQuadMesh(forest, refine_local,
                 min_refine, max_refine);
  delete [] refine_local;

  // Return the error
  return SE_total_error;
}

/*
  Refine the mesh using the original solution and the adjoint solution

  input:
  tacs:        the TACSAssembler object
  uvec:        the solution variables
  adj:         the adjoint solution variables
  forest:      the forest of quadtrees 
  target_err:  absolute value of the target error
  min_refine:  minimum refinement
  max_refnie:  maximum refinement
*/
TacsScalar adjointRefine( TACSAssembler *tacs,
                          TACSAssembler *refine,
                          TACSBVec *uvec, TACSBVec *adjvec,
                          TMRQuadForest *forest,
                          double target_err,
                          int min_refine, int max_refine,
                          TacsScalar *adj_corr ){
  // Get the communicator
  MPI_Comm comm = tacs->getMPIComm();
  
  // Ensure that the local variables have been set
  tacs->setVariables(uvec);

  // Perform a local refinement of the nodes based on the strain energy
  // within each element
  int nelems = tacs->getNumElements();

  // Get the number of variables per node (and the number of
  // derivatives per node required for the reconstruction)
  int vars_per_node = tacs->getVarsPerNode();
  int deriv_per_node = 3*vars_per_node;
  
  // Create the refined residual vector
  TACSBVec *residual = refine->createVec();
  residual->incref();

  // Get the auxiliary elements (surface tractions) associated with the
  // element class
  TACSAuxElements *aux_elements = refine->getAuxElements();
  int naux = 0, aux_count = 0;
  TACSAuxElem *aux = NULL;
  if (aux_elements){
    aux_elements->sort();
    naux = aux_elements->getAuxElements(&aux);
  }

  // For each element in the mesh, compute the original strain energy
  for ( int i = 0; i < nelems; i++ ){
    // Set the simulation time
    double time = 0.0;

    // Get the node locations and variables
    TacsScalar vars[6*9], dvars[6*9], ddvars[6*9];
    tacs->getElement(i, NULL, vars, dvars, ddvars);

    // For each element on the refined mesh, retrieve the
    // local element residual and the element
    for ( int ii = 0; ii < 2; ii++ ){
      for ( int jj = 0; jj < 2; jj++ ){
        // Set the element number on the refined mesh
        int elem_num = 4*i + jj + 2*ii;

        // The refined node locations and element variables
        TacsScalar Xpts[3*9];
        TACSElement *elem = refine->getElement(elem_num, Xpts);

        // Set the local element variables
        TacsScalar evars[6*9];
        memset(evars, 0, elem->numVariables()*sizeof(TacsScalar));

        // Perform the interpolation
        for ( int m = 0; m < 3; m++ ){
          for ( int n = 0; n < 3; n++ ){
            double pt[2];
            pt[0] = -1.0 + 0.5*(2*ii + n);
            pt[1] = -1.0 + 0.5*(2*jj + m);

            // Evaluate the locations of the new nodes
            double N[9];
            FElibrary::biLagrangeSF(N, pt, 3);
           
            // Set the values of the variables at this point
            for ( int k = 0; k < 9; k++ ){
              // Evaluate the interpolation part of the reconstruction
              for ( int kk = 0; kk < 6; kk++ ){
                evars[6*(n + 3*m)+kk] += vars[6*k+kk]*N[k];
              }
            }
          }
        }

        // Compute the quadratic element residual
        TacsScalar res[6*9];
        memset(res, 0, elem->numVariables()*sizeof(TacsScalar));
        elem->addResidual(time, res, Xpts, evars, dvars, ddvars);

        while (aux_count < naux && aux[aux_count].num == elem_num){
          aux[aux_count].elem->addResidual(time, res, Xpts, 
                                           evars, dvars, ddvars);
          aux_count++;
        }

        // Get refined element numbers
        int len;
        const int *nodes;
        refine->getElement(elem_num, &nodes, &len);
        
        // Add the solution error to the residual
        residual->setValues(len, nodes, res, ADD_VALUES);
      }
    }
  }

  // Add all the values together
  residual->beginSetValues(ADD_VALUES);
  residual->endSetValues(ADD_VALUES);

  // Distribute the residual entries
  residual->beginDistributeValues();
  residual->endDistributeValues();

  // Compute the interpolated and reconstructed solution
  TACSBVec *wlocal;
  computeLocalWeights(tacs, &wlocal);
  wlocal->incref();

  // Distribute the adjoint vector derivatives
  adjvec->beginDistributeValues();
  adjvec->endDistributeValues();

  // Compute the nodal derivatives of the adjoint vector
  TACSBVec *adjderiv;
  computeNodeDeriv(tacs, adjvec, wlocal, &adjderiv);
  adjderiv->incref();

  // Free the wlocal vector
  wlocal->decref();  wlocal = NULL;

  // Set local values for the quadratic and cubic adjoint
  // contributions
  TACSBVec *qadjvec = refine->createVec();
  qadjvec->incref();

  // Compute the weights on the refined mesh
  computeLocalWeights(refine, &wlocal);
  wlocal->incref();

  // Allocate the refinement array - which elements will be refined
  int *refine_local = new int[ nelems ];
  memset(refine_local, 0, nelems*sizeof(int));

  for ( int i = 0; i < nelems; i++ ){
    // Get the element node numbers
    int len = 0;
    const int *nodes;
    tacs->getElement(i, &nodes, &len);
        
    // Get the values of the adjoint
    TacsScalar aelem[6*9];
    adjvec->getValues(len, nodes, aelem);

    // Get the derivative of the solution and the adjoint
    TacsScalar dadjelem[18*9];
    adjderiv->getValues(len, nodes, dadjelem);

    // Retrieve the element node locations
    TacsScalar Xpts[3*9];
    tacs->getElement(i, Xpts);

    // 7 enrichment functions for each degree of freedom
    TacsScalar adjbar[7*6];
    computeElemRecon(Xpts, aelem, dadjelem, adjbar);

    // Compute the remaining error for this element
    TacsScalar elem_remain = 0.0;

    for ( int ii = 0; ii < 2; ii++ ){
      for ( int jj = 0; jj < 2; jj++ ){
        // Get the values from the residual
        int elem_num = 4*i + jj + 2*ii;

        // Get the element object
        TACSElement *elem = refine->getElement(elem_num);

        // Get the element node numbers
        refine->getElement(elem_num, &nodes, &len);

        // Get the local part of the residual
        TacsScalar res[6*9], wvals[9];
        residual->getValues(len, nodes, res);
        wlocal->getValues(len, nodes, wvals);

        // The quadratic and cubuic reconstruction of the adjoint
        TacsScalar qadjelem[6*9], cadjelem[6*9];
        memset(qadjelem, 0, 6*9*sizeof(TacsScalar));

        for ( int m = 0; m < 3; m++ ){
          for ( int n = 0; n < 3; n++ ){
            double pt[2];
            pt[0] = -1.0 + 0.5*(2*ii + n);
            pt[1] = -1.0 + 0.5*(2*jj + m);

            // Evaluate the locations of the new nodes
            double N[9];
            FElibrary::biLagrangeSF(N, pt, 3);
           
            // Set the values of the variables at this point
            for ( int k = 0; k < 9; k++ ){
              // Evaluate the interpolation part of the reconstruction
              for ( int kk = 0; kk < 6; kk++ ){
                qadjelem[6*(n + 3*m)+kk] += aelem[6*k+kk]*N[k];
              }
            }
          }
        }

        // Copy over the quadratic part of the adjoint solution
        memcpy(cadjelem, qadjelem, 6*9*sizeof(TacsScalar));
        for ( int m = 0; m < 3; m++ ){
          for ( int n = 0; n < 3; n++ ){
            double pt[2];
            pt[0] = -1.0 + 0.5*(2*ii + n);
            pt[1] = -1.0 + 0.5*(2*jj + m);

            // Evaluate the locations of the new nodes
            double Nr[7];
            evalEnrichmentFuncs(pt, Nr);
           
            // Set the values of the variables at this point
            for ( int k = 0; k < 7; k++ ){
              // Evaluate the interpolation part of the reconstruction
              for ( int kk = 0; kk < 6; kk++ ){
                cadjelem[6*(n + 3*m)+kk] += adjbar[6*k+kk]*Nr[k];
              }
            }
          }
        }

        for ( int j = 0; j < elem->numVariables(); j++ ){
          elem_remain += fabs((cadjelem[j] - qadjelem[j])*res[j]);
          qadjelem[j] = wvals[j/6]*(cadjelem[j] - qadjelem[j]);
        }

        // Add the local residual values
        qadjvec->setValues(len, nodes, qadjelem, ADD_VALUES);
      }
    }

    // If the predicted element error exceeds the target element
    // error, then refine this element
    if (elem_remain >= target_err){
      refine_local[i] = 1;
    }
  }

  qadjvec->beginSetValues(ADD_VALUES);
  qadjvec->endSetValues(ADD_VALUES);
  
  qadjvec->beginDistributeValues();
  qadjvec->endDistributeValues();

  // Set the total remaining error and correction
  TacsScalar *rvals, *avals;
  int vsize = qadjvec->getArray(&avals);
  residual->getArray(&rvals);
 
 TacsScalar total_err_remain = 0.0;
  for ( int i = 0; i < vsize; i++ ){
    total_err_remain += fabs(avals[i]*rvals[i]);
  }

  // Sum up the total error contribution
  MPI_Allreduce(MPI_IN_PLACE, &total_err_remain, 1,
                TACS_MPI_TYPE, MPI_SUM, comm);

  // Compute the adjoint-based residual correction
  TacsScalar total_corr = qadjvec->dot(residual);

  // Free some of the data that is no longer required
  wlocal->decref();
  residual->decref();
  qadjvec->decref();
  adjderiv->decref();

  // Compute the total number of elements
  int ntotal;
  MPI_Allreduce(&nelems, &ntotal, 1, MPI_INT, MPI_SUM, comm);

  // refine the quadrant mesh based on the local refinement values
  refineQuadMesh(forest, refine_local,
                 min_refine, max_refine);

  delete [] refine_local;

  // Set the adjoint residual correction
  *adj_corr = total_corr;

  // Return the error
  return total_err_remain;
}

/*
  Create the TACSAssembler object based on the TMR forest that is
  supplied.
*/
TACSAssembler *createTACSAssembler( MPI_Comm comm, 
                                    TACSElement *element,
                                    TACSElement *trac,
                                    TMRQuadForest *forest,
                                    ProblemGeometry *problem ){
  int mpi_size, mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // Create the TACSCreator object
  int vars_per_node = 6;

  // Create the nodes
  forest->createNodes(ELEMENT_ORDER);

  // Create the element mesh
  int *conn, num_elements = 0;
  forest->createMeshConn(&conn, &num_elements);

  // Find the number of nodes for this processor
  const int *range;
  forest->getOwnedNodeRange(&range);
  int num_nodes = range[mpi_rank+1] - range[mpi_rank];

  // Get the dependent node information
  const int *dep_ptr, *dep_conn;
  const double *dep_weights;
  int num_dep_nodes = 
    forest->getDepNodeConn(&dep_ptr, &dep_conn, &dep_weights);

  // Create the associated TACSAssembler object
  TACSAssembler *tacs = 
    new TACSAssembler(comm, vars_per_node, num_nodes, 
                      num_elements, num_dep_nodes);

  // Set the element ptr
  int nodes_per_elem = ELEMENT_ORDER*ELEMENT_ORDER;
  int *ptr = new int[ nodes_per_elem*num_elements ];
  for ( int i = 0; i < num_elements+1; i++ ){
    ptr[i] = nodes_per_elem*i;
  }
    
  // Set the element connectivity into TACSAssembler
  tacs->setElementConnectivity(conn, ptr);
  delete [] conn;
  delete [] ptr;
  
  // Set the dependent node information
  tacs->setDependentNodes(dep_ptr, dep_conn, dep_weights);
  
  // Create the elements
  TACSElement **elements = new TACSElement*[ num_elements ];
  for ( int i = 0; i < num_elements; i++ ){
    elements[i] = element;
  }
  tacs->setElements(elements);
  delete [] elements;

  // Create the node vector
  TACSBVec *X = tacs->createNodeVec();
  X->incref();

  // Get the node array
  TacsScalar *Xpts;
  X->getArray(&Xpts);
    
  // Fill in the positions
  TMRQuadtree **quad;
  int num_faces = forest->getQuadtrees(&quad);
      
  for ( int face = 0; face < num_faces; face++ ){
    if (quad[face]){
      // Retrieve the node quadrants
      TMRQuadrantArray *nodes;
      quad[face]->getNodes(&nodes);
      
      // Get the array of node quadrants from this face
      int size;
      TMRQuadrant *array;
      nodes->getArray(&array, &size);
      
      // Set the node locations
      for ( int i = 0; i < size; i++ ){
        int node = array[i].tag;
        if (node >= range[mpi_rank] && node < range[mpi_rank+1]){
          int index = node - range[mpi_rank];
          problem->getLocation(face, array[i].x, array[i].y, 
                               &Xpts[3*index]);
          
          // Set the boundary conditions based on the spatial location
          if (problem->isBoundary(face, array[i].x, array[i].y)){
            tacs->addBCs(1, &node);
          }
        }
      }
    }
  }

  // Initialize TACS
  tacs->initialize();

  // Set the node locations
  tacs->setNodes(X);
  X->decref();
  
  // Create the auxiliary element class
  int nelems = tacs->getNumElements();
  TACSAuxElements *aux = new TACSAuxElements(nelems);
  for ( int i = 0; i < nelems; i++ ){
    aux->addElement(i, trac);
  }
  tacs->setAuxElements(aux);

  return tacs;
}

/*
  Run an adaptive analysis
*/
int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  // Set the communicator
  MPI_Comm comm = MPI_COMM_WORLD;

  // Get the MPI communicator rank
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Set the thickness to use in the calculation
  double t = 2.25e-3;
  int orthotropic_flag = 0;

  // Set the refinement levels
  int min_refine = 1;
  int max_refine = TMR_MAX_LEVEL;

  // Create the default problem 
  ProblemGeometry *problem = new SquarePlate();
  char problem_name[64];
  sprintf(problem_name, "plate");
  
  // Set the initial function value
  TacsScalar fval_init = 0.0;

  // Set the target relative error for the functional
  int test_element = 0;
  double target_rel_err = 0.01;
  int use_energy_norm = 0;
  double ks_weight = 10.0;
  int max_iters = 8;

  // Set the traction values
  TacsScalar tx = 0.0, ty = 0.0, tz = 100.0e3;

  for ( int k = 0; k < argc; k++ ){
    if (strcmp(argv[k], "use_energy_norm") == 0){
      target_rel_err = 0.001;
      use_energy_norm = 1;
    }
    double rho = 0.0;
    if (sscanf(argv[k], "ks_weight=%lf", &rho) == 1){
      ks_weight = rho;
    }
    if (strcmp(argv[k], "test_element") == 0){
      test_element = 1;
    }
    int niter = 0;
    if (sscanf(argv[k], "max_iters=%d", &niter) == 1){
      if (niter > 0 && niter < 15){
        max_iters = niter;
      }
    }
    if (strcmp(argv[k], "wing") == 0){
      // Set the min refinement and maximum number of iterations
      max_iters = 3;
      min_refine = 0;

      // Set the thickness and traction values
      t = 0.05;
      tx = 0.0, ty = 0.0;
      tz = 1000.0;
      target_rel_err = 0.001;
   
      // Set the new problem name
      sprintf(problem_name, "wing");

      // Create the TACSMeshLoader class
      TACSMeshLoader *mesh = new TACSMeshLoader(comm);
      mesh->incref();
      mesh->scanBDFFile("wing/CRM_box_2nd.bdf");
      
      // Extract the connectivity
      int nnodes, nelems;
      const int *elem_node_conn;
      const double *Xpts;
      mesh->getConnectivity(&nnodes, &nelems, NULL, 
                            &elem_node_conn, &Xpts);
      
      delete problem;
      problem = new MeshProblem(nnodes, nelems, 
                                elem_node_conn, Xpts);
      mesh->decref();
    }
  }

  // Create a prefix for the type of problem we're running
  char prefix[128];
  if (use_energy_norm){
    sprintf(prefix, "%s/energy_norm", problem_name);
  }
  else {
    sprintf(prefix, "%s/ks%.0f", problem_name, ks_weight);
  }

  OrthoPly *ply = NULL;
  if (orthotropic_flag){
    // Set the material properties to use
    double rho = 1.0;
    double E1 = 100.0e9;
    double E2 = 5.0e9;
    double nu12 = 0.25;
    double G12 = 10.0e9;
    double G13 = 10.0e9;
    double G23 = 4.0e9;

    double Xt = 100.0e6;
    double Xc = 50.0e6;
    double Yt = 2.5e6;
    double Yc = 10.0e6;
    double S12 = 8.0e6;

    ply = new OrthoPly(t, rho, E1, E2, nu12, 
                       G12, G23, G13, 
                       Xt, Xc, Yt, Yc, S12);
    printf("Using orthotropic material properties: \n");
  }
  else {
    // Set the material properties to use
    double rho = 2700.0;
    double E = 70e9;
    double nu = 0.3;
    double ys = 350e6;
  
    ply = new OrthoPly(t, rho, E, nu, ys);
    printf("Using isotropic material properties: \n");
  }

  // Create the stiffness relationship
  double kcorr = 5.0/6.0;
  FSDTStiffness *stiff = 
    new specialFSDTStiffness(ply, orthotropic_flag, t, kcorr);
  
  TacsScalar axis[] = {1.0, 0.0, 0.0};
  stiff->setRefAxis(axis);

  if (rank == 0){
    ply->printProperties();
    stiff->printStiffness();
  }

  // Create the traction class 
  TACSElement *trac = new TACSShellTraction<ELEMENT_ORDER>(tx, ty, tz);
  trac->incref();

  // Create the shell element
  TACSElement *element = new MITCShell<ELEMENT_ORDER>(stiff);
  element->incref();

  // Print out the parameters to the screen
  if (rank == 0){
    printf("target_rel_err =  %f\n", target_rel_err);
    printf("use_energy_norm = %d\n", use_energy_norm);
    printf("ks_weight =       %f\n", ks_weight);
  }

  TMRQuadForest *forest = new TMRQuadForest(comm);
    
  // Set the connectivity
  int num_nodes, num_faces;
  const int *conn;
  problem->getConnectivity(&num_nodes, &conn, &num_faces);
  forest->setConnectivity(num_nodes, conn, num_faces);
    
  // Allocate the trees
  forest->createTrees(min_refine);

  // Create the problem file
  FILE *fp = NULL;
  if (rank == 0){
    char filename[256];
    sprintf(filename, "%serror_history.dat", prefix);
    fp = fopen(filename, "w");
    fprintf(fp, "Variables = iter, nelems, nnodes, fval, fcorr, error\n");
  }

  for ( int iter = 0; iter < max_iters; iter++ ){
    // Balance the forest so that we can use it!
    forest->balance(1);
    
    TACSAssembler *tacs = 
      createTACSAssembler(comm, element, trac, forest, problem);
    tacs->incref();

    // Create the KS function
    TACSKSFailure *ks_func = new TACSKSFailure(tacs, ks_weight);
    ks_func->incref();
    ks_func->setKSFailureType(TACSKSFailure::CONTINUOUS);

    // Allocate the compliance functional
    TACSCompliance *comp_func = new TACSCompliance(tacs);
    comp_func->incref();

    // Set the actual function to use
    TACSFunction *functional = ks_func;
    if (use_energy_norm){
      functional = comp_func;
    }
        
    // Create the preconditioner
    TACSBVec *res = tacs->createVec();
    TACSBVec *ans = tacs->createVec();
    TACSBVec *tmp = tacs->createVec();
    TACSBVec *adjoint = tacs->createVec();
    FEMat *mat = tacs->createFEMat();

    // Increment the reference count to the matrix/vectors
    res->incref();
    ans->incref();
    tmp->incref();
    adjoint->incref();
    mat->incref();
    
    // Allocate the factorization
    int lev = 4500;
    double fill = 10.0;
    int reorder_schur = 1;
    PcScMat *pc = new PcScMat(mat, lev, fill, reorder_schur);
    pc->incref();

    // Create the GMRES object
    int gmres_iters = 10; 
    int nrestart = 2; // Number of allowed restarts
    int is_flexible = 0; // Is a flexible preconditioner?
    GMRES *gmres = new GMRES(mat, pc, gmres_iters, nrestart, is_flexible);
    gmres->incref();
    
    // Assemble and factor the stiffness/Jacobian matrix
    double alpha = 1.0, beta = 0.0, gamma = 0.0;
    tacs->assembleJacobian(alpha, beta, gamma, res, mat);

    // Factor the Jacobian matrix
    pc->factor();
    gmres->solve(res, ans);
    ans->scale(-1.0);

    // Set the solution
    tacs->setVariables(ans);

    // Check the residual error for the solution
    mat->mult(ans, tmp);
    tmp->axpy(1.0, res);

    // Print the solution norm
    TacsScalar norm = tmp->norm()/res->norm();
    if (rank == 0){
      printf("Solution residual norm: %15.5e\n", norm);
    }

    // Solve the adjoint equation
    TacsScalar fval;
    tacs->evalFunctions(&functional, 1, &fval);

    if (iter == 0){
      fval_init = fval;
    }

    res->zeroEntries();
    tacs->addSVSens(alpha, beta, gamma, &functional, 1, &res);
    gmres->solve(res, adjoint);
    adjoint->scale(-1.0);
 
    // Create an TACSToFH5 object for writing output to files
    unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                               TACSElement::OUTPUT_DISPLACEMENTS |
                               TACSElement::OUTPUT_STRAINS |
                               TACSElement::OUTPUT_STRESSES |
                               TACSElement::OUTPUT_EXTRAS);
    TACSToFH5 *f5 = new TACSToFH5(tacs, SHELL, write_flag);
    f5->incref();

    // Write the displacements
    char filename[256];
    sprintf(filename, "%s_solution%d.f5", prefix, iter);
    f5->writeToFile(filename);

    // Delete the viewer
    f5->decref();

    // Set the target error using the factor: max(1.0, 16*(2**-iter))
    double factor = 16.0/(1 << iter);
    if (factor < 1.0){ factor = 1.0; }

    // Get the total number of elements
    int nnodes = tacs->getNumNodes();
    int nelems = tacs->getNumElements();
    int totals[2] = {nnodes, nelems};
    MPI_Allreduce(MPI_IN_PLACE, totals, 2, MPI_INT, MPI_SUM, comm);
    nnodes = totals[0];
    nelems = totals[1];
    
    // Set the absolute element-level error based on the relative
    // error that is requested
    double target_abs_err = factor*target_rel_err*fval/nelems;

    // Compute the computable error contribution
    TacsScalar error_total = 0.0;

    if (use_energy_norm){
      // Perform the refinement
      error_total = strainEnergyRefine(tacs, ans, forest,
                                       target_abs_err,
                                       min_refine, max_refine);
    
      if (fp){
        fprintf(fp, "%d %d %d %15.10e %15.10e %15.10e\n",
                iter, nelems, nnodes, 
                fval/fval_init, 
                (fval+error_total)/fval_init, 
                error_total/fval_init);
        fflush(fp);
      }
    }
    else {
      // Duplicate the forest
      TMRQuadForest *dup = forest->duplicate();

      // Refine the quadtrees uniformly
      TMRQuadtree **trees;
      int ntrees = dup->getQuadtrees(&trees);
      for ( int i = 0; i < ntrees; i++ ){
        if (trees[i]){
          trees[i]->refine();
        }
      }

      // Balance the tree
      dup->balance(1);

      // Create the uniformly refined TACSAssembler object
      TACSAssembler *refine = 
        createTACSAssembler(comm, element, trac, dup, problem);
      refine->incref();
      
      // Free the duplicated elements
      delete dup;
      
      // Perform the error assessment
      TacsScalar adjcorr;
      error_total = 
        adjointRefine(tacs, refine, ans, adjoint, forest,
                      target_abs_err, min_refine, max_refine, &adjcorr);
      
      // Deallocate the refined object
      refine->decref();

      if (fp){
        fprintf(fp, "%d %d %d %15.10e %15.10e %15.10e\n",
                iter, nelems, nnodes, 
                fval/fval_init, 
                (fval+adjcorr)/fval_init, 
                error_total/fval_init);
        fflush(fp);
      }
    }
    
    // Free everything
    tacs->decref();
    gmres->decref();
    pc->decref();
    ks_func->decref();
    comp_func->decref();
    res->decref();
    ans->decref();
    tmp->decref();
    adjoint->decref();
    mat->decref();
  }

  if (fp){
    fclose(fp);
  }

  if (rank == 0){
    delete forest;
  }

  element->decref();
  trac->decref();
  delete problem;

  TMRFinalize();
  MPI_Finalize();
  return (0);
}
