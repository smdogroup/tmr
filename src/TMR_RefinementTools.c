#include "TMR_RefinementTools.h"
#include "TensorToolbox.h"
#include "tacslapack.h"

// Include the stdlib set/string classes
#include <string>
#include <set>

/*
  Create a multgrid object for a forest of octrees
*/
void TMR_CreateTACSMg( int num_levels, 
                       TACSAssembler *tacs[],
                       TMROctForest *forest[],
                       TACSMg **_mg, int use_pairs ){
  // Get the communicator
  MPI_Comm comm = tacs[0]->getMPIComm();

  // Create the multigrid object
  int zero_guess = 0;
  double omega = 1.0;
  int mg_sor_iters = 1;
  int mg_sor_symm = 1;
  int mg_iters_per_level = 1;
  TACSMg *mg = new TACSMg(comm, num_levels);
  
  // Create the intepolation/restriction objects between mesh levels
  for ( int level = 0; level < num_levels-1; level++ ){
    // Create the interpolation object
    TACSBVecInterp *interp = 
      new TACSBVecInterp(tacs[level+1]->getVarMap(),
                         tacs[level]->getVarMap(),
                         tacs[level]->getVarsPerNode());

    // Set the interpolation
    forest[level]->createInterpolation(forest[level+1], interp);
    
    // Initialize the interpolation
    interp->initialize();

    // If we're using pairs, then identify the nodes that have
    // multiple dof/node in the octree
    int npairs = 0;
    int *pairs = NULL;
    if (use_pairs){
      TMROctantArray *nodes;
      forest[level]->getNodes(&nodes);
      
      // Get the nodes
      int size;
      TMROctant *array;
      nodes->getArray(&array, &size);
      for ( int i = 0; i < size; i++ ){
        if (array[i].level == 2){
          npairs++;
        }
      }

      // Set the DOF associted with the pairs
      pairs = new int[ npairs ];
      npairs = 0;
      for ( int i = 0; i < size; i++ ){
        if (array[i].level == 2){
          pairs[npairs] = array[i].tag;
          npairs++;
        }
      }
    }

    // Create the matrix
    PMat *mat = tacs[level]->createMat();
    PSOR *pc = new PSOR(mat, zero_guess, omega, 
                        mg_sor_iters, mg_sor_symm, pairs, npairs);

    // Set the interpolation and TACS object within the multigrid object
    mg->setLevel(level, tacs[level], interp, mg_iters_per_level,
                 mat, pc);

    if (pairs){ delete [] pairs; }
  }

  // Set the lowest level - with no interpolation object
  mg->setLevel(num_levels-1, tacs[num_levels-1], NULL);

  // Return the multigrid object
  *_mg = mg;
} 

/*
  Create the TACS multigrid objects for the quadrilateral case
*/
void TMR_CreateTACSMg( int num_levels, 
                       TACSAssembler *tacs[],
                       TMRQuadForest *forest[],
                       TACSMg **_mg ){
  // Get the communicator
  MPI_Comm comm = tacs[0]->getMPIComm();

  // Create the multigrid object
  double omega = 1.0;
  int mg_sor_iters = 1;
  int mg_sor_symm = 1;
  int mg_iters_per_level = 1;
  TACSMg *mg = new TACSMg(comm, num_levels, omega, 
                          mg_sor_iters, mg_sor_symm);
  
  // Create the intepolation/restriction objects between mesh levels
  for ( int level = 0; level < num_levels-1; level++ ){
    // Create the interpolation object
    TACSBVecInterp *interp = 
      new TACSBVecInterp(tacs[level+1]->getVarMap(),
                         tacs[level]->getVarMap(),
                         tacs[level]->getVarsPerNode());

    // Set the interpolation
    forest[level]->createInterpolation(forest[level+1], interp);
    
    // Initialize the interpolation
    interp->initialize();

    // Set the interpolation and TACS object within the multigrid object
    mg->setLevel(level, tacs[level], interp, mg_iters_per_level);
  }

  // Set the lowest level - with no interpolation object
  mg->setLevel(num_levels-1, tacs[num_levels-1], NULL);

  // Return the multigrid object
  *_mg = mg;
} 

/*
  Compute the transpose of the Jacobian transformation at a point
  within the element.
*/
static void computeJacobianTrans2D( const TacsScalar Xpts[],
                                    const double Na[], const double Nb[], 
                                    TacsScalar Xd[], TacsScalar J[],
                                    const int num_nodes ){
  memset(Xd, 0, 9*sizeof(TacsScalar));

  // Compute the derivative along the coordinate directions
  const double *na = Na, *nb = Nb;
  const TacsScalar *X = Xpts;
  for ( int i = 0; i < num_nodes; i++ ){
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
  Compute the transpose of the 3D Jacobian transformation at a point
  within the element
*/
static void computeJacobianTrans3D( const TacsScalar Xpts[],
                                    const double Na[], 
                                    const double Nb[], 
                                    const double Nc[],
                                    TacsScalar Xd[], TacsScalar J[],
                                    const int num_nodes ){
  memset(Xd, 0, 9*sizeof(TacsScalar));

  // Compute the derivative along the coordinate directions
  const double *na = Na, *nb = Nb, *nc = Nc;
  const TacsScalar *X = Xpts;
  for ( int i = 0; i < num_nodes; i++ ){
    Xd[0] += X[0]*na[0];
    Xd[1] += X[1]*na[0];
    Xd[2] += X[2]*na[0];

    Xd[3] += X[0]*nb[0];
    Xd[4] += X[1]*nb[0];
    Xd[5] += X[2]*nb[0];

    Xd[6] += X[0]*nc[0];
    Xd[7] += X[1]*nc[0];
    Xd[8] += X[2]*nc[0];

    na++;  nb++;  nc++;
    X += 3;
  }

  // Compute the transpose of the Jacobian transformation
  FElibrary::jacobian3d(Xd, J);
}

/*
  Evaluate the enrichment function for a second-order shell problem
*/
static void eval2ndEnrichmentFuncs2D( const double pt[], double N[] ){
  // Compute the cubic enrichment shape functions along the two
  // coordinate directions
  double ca = (1.0 + pt[0])*(1.0 - pt[0]);
  double cb = (1.0 + pt[1])*(1.0 - pt[1]);

  N[0] = ca;
  N[1] = pt[1]*ca;
  N[2] = cb;
  N[3] = pt[0]*cb;
  N[4] = ca*cb;
}

/*
  Evaluate the derivative of the enrichment functions for a
  second-order problem.
*/
static void eval2ndEnrichmentFuncs2D( const double pt[], 
                                      double N[], double Na[], double Nb[] ){
  
  // Compute the cubic enrichment shape functions along the two
  // coordinate directions
  double ca = (1.0 + pt[0])*(1.0 - pt[0]);
  double cb = (1.0 + pt[1])*(1.0 - pt[1]);

  // Compute the derivatives
  double da = -2.0*pt[0];
  double db = -2.0*pt[1];

  // Evaluate the shape functions
  N[0] = ca;
  N[1] = pt[1]*ca;
  N[2] = cb;
  N[3] = pt[0]*cb;
  N[4] = ca*cb;

  // Evaluate the derivatives of the shape functions
  Na[0] = da;
  Na[1] = pt[1]*da;
  Na[2] = 0.0;
  Na[3] = cb;
  Na[4] = da*cb;

  Nb[0] = 0.0;
  Nb[1] = ca;
  Nb[2] = db;
  Nb[3] = pt[0]*db;
  Nb[4] = ca*db;
} 

/*
  Evaluate the enrichment functions for a third-order shell problem
*/
static void eval3rdEnrichmentFuncs2D( const double pt[], double N[] ){
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
static void eval3rdEnrichmentFuncs2D( const double pt[],
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
static void computeElemRecon2D( const int order,
                                const int vars_per_node,
                                const TacsScalar Xpts[],
                                const TacsScalar uvals[],
                                const TacsScalar uderiv[],
                                TacsScalar ubar[],
                                TacsScalar *tmp ){
  // The number of enrichment functions: 5 if we have a second-order
  // mesh 7 if we have a third-order mesh
  const int nenrich = (order == 2 ? 5 : 7);

  // The number of equations
  const int neq = 2*order*order;

  // The number of derivatives per node
  const int deriv_per_node = 3*vars_per_node;  

  // Set up the least squares problem at the nodes
  int nrhs = vars_per_node;

  TacsScalar *A = &tmp[0];
  TacsScalar *b = &tmp[nenrich*neq];

  // Set the weighs 
  double wvals[3];
  if (order == 2){
    wvals[0] = wvals[1] = 1.0;
  }
  else if (order == 3){
    wvals[0] = wvals[2] = 0.5;
    wvals[1] = 1.0;
  }

  for ( int c = 0, jj = 0; jj < order; jj++ ){
    for ( int ii = 0; ii < order; ii++, c += 2 ){
      // Set the parametric location within the element
      double pt[2];
      pt[0] = -1.0 + (2.0*ii)/(order-1);
      pt[1] = -1.0 + (2.0*jj)/(order-1);

      // Compute the element shape functions at this point
      double N[9], Na[9], Nb[9];
      FElibrary::biLagrangeSF(N, Na, Nb, pt, order);
      
      // Evaluate the Jacobian transformation at this point
      TacsScalar Xd[9], J[9];
      computeJacobianTrans2D(Xpts, Na, Nb, Xd, J, order*order);

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
      const TacsScalar *ud = &uderiv[deriv_per_node*(ii + order*jj)];

      for ( int k = 0; k < vars_per_node; k++ ){
        b[neq*k+c] = 
          wvals[ii]*wvals[jj]*(d1[0]*ud[0] + d1[1]*ud[1] + d1[2]*ud[2]);
        b[neq*k+c+1] = 
          wvals[ii]*wvals[jj]*(d2[0]*ud[0] + d2[1]*ud[1] + d2[2]*ud[2]);
        ud += 3;
      }
      
      // Add the contribution from the nodes
      for ( int k = 0; k < vars_per_node; k++ ){
        // Evaluate the derivatives along the parametric directions
        TacsScalar Ua = 0.0, Ub = 0.0;
        for ( int i = 0; i < order*order; i++ ){
          Ua += uvals[vars_per_node*i + k]*Na[i];
          Ub += uvals[vars_per_node*i + k]*Nb[i];
        }

        // Compute the derivative along the x,y,z directions
        TacsScalar d[3];
        d[0] = Ua*J[0] + Ub*J[1];
        d[1] = Ua*J[3] + Ub*J[4];
        d[2] = Ua*J[6] + Ub*J[7];

        b[neq*k+c] -= 
          wvals[ii]*wvals[jj]*(d1[0]*d[0] + d1[1]*d[1] + d1[2]*d[2]);
        b[neq*k+c+1] -= 
          wvals[ii]*wvals[jj]*(d2[0]*d[0] + d2[1]*d[1] + d2[2]*d[2]);
      }

      // xi,X = [X,xi]^{-1}
      // U,X = U,xi*xi,X = U,xi*J^{T}
      
      // Now, evaluate the terms for the left-hand-side
      // that contribute to the derivative
      double Nr[7], Nar[7], Nbr[7];
      if (order == 2){
        eval2ndEnrichmentFuncs2D(pt, Nr, Nar, Nbr);
      }
      else {
        eval3rdEnrichmentFuncs2D(pt, Nr, Nar, Nbr);
      }

      // Add the contributions to the the enricment 
      for ( int i = 0; i < nenrich; i++ ){
        TacsScalar d[3];
        d[0] = Nar[i]*J[0] + Nbr[i]*J[1];
        d[1] = Nar[i]*J[3] + Nbr[i]*J[4];
        d[2] = Nar[i]*J[6] + Nbr[i]*J[7];

        A[neq*i+c] = 
          wvals[ii]*wvals[jj]*(d1[0]*d[0] + d1[1]*d[1] + d1[2]*d[2]);
        A[neq*i+c+1] = 
          wvals[ii]*wvals[jj]*(d2[0]*d[0] + d2[1]*d[1] + d2[2]*d[2]);
      }
    }
  }

  // Singular values
  TacsScalar s[7];
  int m = neq, n = nenrich;
  double rcond = -1.0;
  int rank;

  // Length of the work array
  int lwork = 10*18;
  TacsScalar work[10*18];
  
  // LAPACK output status
  int info;

  // Using LAPACK, compute the least squares solution
  LAPACKdgelss(&m, &n, &nrhs, A, &m, b, &m, s, 
               &rcond, &rank, work, &lwork, &info);

  // Copy over the ubar solution
  for ( int i = 0; i < nenrich; i++ ){
    for ( int j = 0; j < vars_per_node; j++ ){
      ubar[vars_per_node*i + j] = b[m*j + i];
    }
  }
}

/*
  Compute the local derivative weights
*/
static void computeLocalWeights( TACSAssembler *tacs, 
                                 TACSBVec *weights,
                                 const int *element_nums=NULL,
                                 int num_elements=-1 ){
  // Zero the entries in the array
  weights->zeroEntries();

  // Set the local element weights
  int max_nodes = tacs->getMaxElementNodes();
  TacsScalar *welem = new TacsScalar[ max_nodes ];
  for ( int i = 0; i < max_nodes; i++ ){
    welem[i] = 1.0;
  }

  if (num_elements < 0){
    // Add unit weights to all the elements. This will sum up the
    // number of times each node is referenced.
    int nelems = tacs->getNumElements();
    for ( int i = 0; i < nelems; i++ ){
      int len = 0;
      const int *nodes;
      tacs->getElement(i, &nodes, &len);   
      
      // Compute the local weights
      for ( int j = 0; j < len; j++ ){
        welem[j] = 1.0;
        if (nodes[j] < 0){
          welem[j] = 0.0;
        }
      }
      
      weights->setValues(len, nodes, welem, TACS_ADD_VALUES);
    }
  }
  else {
    // Add up the weights for only those elements in the list
    for ( int i = 0; i < num_elements; i++ ){
      int elem = element_nums[i];

      // Get the node numbers for this elements
      int len = 0;
      const int *nodes;
      tacs->getElement(elem, &nodes, &len);
      
      // Compute the local weights
      for ( int j = 0; j < len; j++ ){
        welem[j] = 1.0;
        if (nodes[j] < 0){
          welem[j] = 0.0;
        }
      }
      
      weights->setValues(len, nodes, welem, TACS_ADD_VALUES);
    }
  }

  delete [] welem;

  // Finish setting all of the values
  weights->beginSetValues(TACS_ADD_VALUES);
  weights->endSetValues(TACS_ADD_VALUES);

  // Distribute the values
  weights->beginDistributeValues();
  weights->endDistributeValues();
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
static void computeNodeDeriv2D( const int order,
                                TACSAssembler *tacs, TACSBVec *uvec, 
                                TACSBVec *weights, TACSBVec *uderiv,
                                const int *element_nums=NULL,
                                int num_elements=-1 ){
  // Zero the nodal derivatives
  uderiv->zeroEntries();

  // The number of variables at each node
  int vars_per_node = tacs->getVarsPerNode();

  // Number of derivatives per node - x,y,z derivatives
  int deriv_per_node = 3*vars_per_node;

  // Get the number of elements
  int nelems = tacs->getNumElements();

  // Allocate space for the element-wise values and derivatives
  TacsScalar *Ud = new TacsScalar[ 2*vars_per_node ];
  TacsScalar *uelem = new TacsScalar[ order*order*vars_per_node ];
  TacsScalar *delem = new TacsScalar[ order*order*deriv_per_node ];

  // We're only going to iterate over the elements in the list, not
  // the entire array
  if (element_nums){
    nelems = num_elements;
  }

  // Perform the reconstruction for the local
  for ( int index = 0; index < nelems; index++ ){
    // Set the element number, depending on whether we're using the
    // full set of elements, or only a subset of the elements
    int elem = index;
    if (element_nums){
      elem = element_nums[index];      
    }

    // Get the element nodes
    int len = 0;
    const int *nodes;
    tacs->getElement(elem, &nodes, &len);

    // Get the local weight values for this element
    TacsScalar welem[9];
    weights->getValues(len, nodes, welem);

    // Get the local element variables
    uvec->getValues(len, nodes, uelem);

    // Get the node locations for the element
    TacsScalar Xpts[3*9];
    tacs->getElement(elem, Xpts);
    
    // Compute the derivative of the components of the variables along
    // each of the 3-coordinate directions
    TacsScalar *d = delem;

    // Compute the contributions to the derivative from this side of
    // the element    
    for ( int jj = 0; jj < order; jj++ ){
      for ( int ii = 0; ii < order; ii++ ){
        double pt[2];
        pt[0] = -1.0 + 2.0*ii/(order-1);
        pt[1] = -1.0 + 2.0*jj/(order-1);

        // Evaluate the the quadratic shape functions at this point
        double N[9], Na[9], Nb[9];
        FElibrary::biLagrangeSF(N, Na, Nb, pt, order);
      
        // Evaluate the Jacobian transformation at this point
        TacsScalar Xd[9], J[9];
        computeJacobianTrans2D(Xpts, Na, Nb, Xd, J, order*order);

        // Compute the derivatives from the interpolated solution
        memset(Ud, 0, 2*vars_per_node*sizeof(TacsScalar));
        for ( int k = 0; k < vars_per_node; k++ ){
          const TacsScalar *ue = &uelem[k];
          for ( int i = 0; i < order*order; i++ ){
            Ud[2*k] += ue[0]*Na[i];
            Ud[2*k+1] += ue[0]*Nb[i];
            ue += vars_per_node;
          }
        }

        // Evaluate the x/y/z derivatives of each value at the
        // independent nodes
        TacsScalar winv = 1.0/welem[ii + jj*order];
        if (nodes[ii + jj*order] >= 0){
          for ( int k = 0; k < vars_per_node; k++ ){
            d[0] = winv*(Ud[2*k]*J[0] + Ud[2*k+1]*J[1]);
            d[1] = winv*(Ud[2*k]*J[3] + Ud[2*k+1]*J[4]);
            d[2] = winv*(Ud[2*k]*J[6] + Ud[2*k+1]*J[7]);
            d += 3;
          }
        }
        else {
          for ( int k = 0; k < vars_per_node; k++ ){
            d[0] = d[1] = d[2] = 0.0;
            d += 3;
          }
        }
      }
    }

    // Add the values of the derivatives
    uderiv->setValues(len, nodes, delem, TACS_ADD_VALUES);
  }

  // Free the element values
  delete [] Ud;
  delete [] uelem;
  delete [] delem;

  // Add the values in parallel
  uderiv->beginSetValues(TACS_ADD_VALUES);
  uderiv->endSetValues(TACS_ADD_VALUES);

  // Distribute the values so that we can call getValues
  uderiv->beginDistributeValues();
  uderiv->endDistributeValues();
}

/*
  Reconstruct the solution on a more refined mesh
*/
void addRefinedSolution( const int order,
                         TACSAssembler *tacs,
                         TACSAssembler *tacs_refine,
                         TACSBVec *vec,
                         TACSBVec *vecDeriv,
                         TACSBVec *vec_refine,
                         const int *element_nums=NULL,
                         int num_elements=-1  ){  
  // Number of variables/derivatives per node
  const int vars_per_node = tacs->getVarsPerNode();
  const int deriv_per_node = 3*vars_per_node;

  // The number of equations: 2 times the number of nodes for each element
  const int neq = 2*order*order;

  // The number of enrichment functions: 5 if we have a second-order
  // mesh 7 if we have a third-order mesh
  const int nenrich = (order == 2 ? 5 : 7);
  
  // Allocate space for the element reconstruction problem
  TacsScalar *tmp = new TacsScalar[ neq*(nenrich + vars_per_node) ];
 
  // Element solution on the coarse TACS mesh
  TacsScalar *uelem = new TacsScalar[ vars_per_node*order*order ];
  TacsScalar *delem = new TacsScalar[ deriv_per_node*order*order ];
  TacsScalar *ubar = new TacsScalar[ vars_per_node*nenrich ];

  // Refined element solution
  TacsScalar *uref = new TacsScalar[ vars_per_node*order*order ];

  // The maximum number of nodes for any element
  TacsScalar Xpts[3*9];

  // Number of local elements in the coarse version of TACS
  int nelems = tacs->getNumElements();
  if (element_nums){
    nelems = num_elements;
  }

  for ( int index = 0; index < nelems; index++ ){
    // Get the element number
    int elem = index;
    if (element_nums){
      elem = element_nums[index];
    }
    
    // Get the node numbers and node locations for this element
    int len;
    const int *nodes;
    tacs->getElement(elem, &nodes, &len);
    tacs->getElement(elem, Xpts);

    // Get the derivatives at the nodes
    vec->getValues(len, nodes, uelem);
    vecDeriv->getValues(len, nodes, delem);

    // Compute the reconstruction weights for the enrichment functions
    computeElemRecon2D(order, vars_per_node,
                       Xpts, uelem, delem, ubar, tmp);

    // Loop over the refined elements
    for ( int ii = 0; ii < 2; ii++ ){
      for ( int jj = 0; jj < 2; jj++ ){
        // Get the refined element nodes
        const int *refine_nodes;
        tacs_refine->getElement(4*elem + jj + 2*ii,
                                &refine_nodes, NULL);

        // Zero the refined element contribution
        memset(uref, 0, vars_per_node*order*order*sizeof(TacsScalar));

        // Compute the element order
        for ( int m = 0; m < order; m++ ){
          for ( int n = 0; n < order; n++ ){
            // Set the new parameter point in the refined element
            double pt[2];
            pt[0] = ii - 1.0 + 1.0*n/(order-1);
            pt[1] = jj - 1.0 + 1.0*m/(order-1);

            // Evaluate the shape functions and the enrichment
            // functions at the new parametric point
            double N[9], Nr[7];
            FElibrary::biLagrangeSF(N, pt, order);
            if (order == 2){
              eval2ndEnrichmentFuncs2D(pt, Nr);
            }
            else {
              eval3rdEnrichmentFuncs2D(pt, Nr);
            }
            
            // Set the values of the variables at this point
            for ( int i = 0; i < vars_per_node; i++ ){
              const TacsScalar *ue = &uelem[i];
              TacsScalar *u = &uref[vars_per_node*(n + order*m) + i];

              for ( int k = 0; k < order*order; k++ ){
                u[0] += N[k]*ue[vars_per_node*k];
              }
            }

            // Add the portion from the enrichment functions
            for ( int i = 0; i < vars_per_node; i++ ){
              const TacsScalar *ue = &ubar[i];
              TacsScalar *u = &uref[vars_per_node*(n + order*m) + i];

              for ( int k = 0; k < nenrich; k++ ){
                u[0] += Nr[k]*ue[vars_per_node*k];
              }
            }
          }
        }

        // Zero the contribution if it goes to a dependent node
        for ( int i = 0; i < order*order; i++ ){
          if (refine_nodes[i] < 0){
            for ( int j = 0; j < vars_per_node; j++ ){
              uref[vars_per_node*i + j] = 0.0;
            }
          }
        }

        // Add the contributions to the element
        vec_refine->setValues(len, refine_nodes, uref, TACS_ADD_VALUES);
      }
    }
  }

  // Free allocated data
  delete [] tmp;
  delete [] uelem;
  delete [] delem;
  delete [] ubar;
  delete [] uref;
}

/*
  Compute the reconstructed solution on the uniformly refined mesh.
*/
void TMR_ComputeReconSolution( TMRQuadForest *forest,
                               TACSAssembler *tacs,
                               TACSAssembler *tacs_refined,
                               TACSBVec *_uvec,
                               TACSBVec *_uvec_refined ){
  // Get the mesh order
  const int order = forest->getMeshOrder();

  // Retrieve the variables from the TACSAssembler object
  TACSBVec *uvec = _uvec; 
  if (!uvec){
    uvec = tacs->createVec();
    uvec->incref();
    tacs->getVariables(uvec);
  }
  TACSBVec *uvec_refined = _uvec_refined;
  if (!uvec_refined){
    uvec_refined = tacs_refined->createVec();
    uvec_refined->incref();
  }

  // Distribute the solution vector answer
  uvec->beginDistributeValues();
  uvec->endDistributeValues();

  // Allocate a vector for the derivatives
  int vars_per_node = tacs->getVarsPerNode();
  TACSBVec *uderiv = 
    new TACSBVec(tacs->getVarMap(), 3*vars_per_node,
                 tacs->getBVecDistribute(), tacs->getBVecDepNodes());
  uderiv->incref();

  // Create the weight vector - the weights are the number of times
  // each node is referenced by adjacent elements, including
  // inter-process references.
  TACSBVec *weights = new TACSBVec(tacs->getVarMap(), 1,
                                   tacs->getBVecDistribute(),
                                   tacs->getBVecDepNodes());
  weights->incref();

  // Get the underlying topology object
  TMRTopology *topo = forest->getTopology();

  // Compute the max size of the element array
  int nelems = tacs->getNumElements();
  int *face_elem_nums = new int[ nelems ];


  // Loop over all of the faces and uniquely sort the faces
  int num_faces = topo->getNumFaces();
  std::set<std::string> face_attr_set;
  for ( int face_num = 0; face_num < num_faces; face_num++ ){
    TMRFace *face;
    topo->getFace(face_num, &face);
    const char *attr = face->getAttribute();
    if (attr){
      std::string str(attr);
      face_attr_set.insert(str);
    }
    else {
      std::string str("");
      face_attr_set.insert(str);
    }
  }

  // Loop over all of the faces
  std::set<std::string>::iterator it;
  for ( it = face_attr_set.begin(); it != face_attr_set.end(); it++ ){
    // Get the quads with the given face number
    const char *attr = NULL;
    if (!(*it).empty()){
      attr = (*it).c_str();
    }
    TMRQuadrantArray *quad_array = forest->getQuadsWithAttribute(attr);

    // Get the quadrants for this face
    int num_face_elems;
    TMRQuadrant *array;
    quad_array->getArray(&array, &num_face_elems);

    // Create an array of the element numbers for this face
    for ( int i = 0; i < num_face_elems; i++ ){
      face_elem_nums[i] = array[i].tag;
    }
    
    // Free the quadrant array - it is no longer required
    delete quad_array;

    // Compute the nodal weights for the derivatives
    computeLocalWeights(tacs, weights, 
                        face_elem_nums, num_face_elems);

    // Compute the nodal derivatives
    computeNodeDeriv2D(order, tacs, uvec, weights, uderiv,
                       face_elem_nums, num_face_elems);

    // Compute the refined solution
    addRefinedSolution(order, tacs, tacs_refined,
                       uvec, uderiv, uvec_refined,
                       face_elem_nums, num_face_elems);
  }

  // Free the temp array
  delete [] face_elem_nums;

  // Free the weights on the coarse mesh
  weights->decref();

  // Add the values
  uvec_refined->beginSetValues(TACS_ADD_VALUES);
  uvec_refined->endSetValues(TACS_ADD_VALUES);

  // Create a vector for the refined weights
  TACSBVec *weights_refined = new TACSBVec(tacs_refined->getVarMap(), 1,
                                           tacs_refined->getBVecDistribute(),
                                           tacs_refined->getBVecDepNodes());
  weights_refined->incref();

  // Compute the nodal weights for the refined mesh
  computeLocalWeights(tacs_refined, weights_refined);

  TacsScalar *u, *w;
  uvec_refined->getArray(&u);
  int size = weights_refined->getArray(&w);
  
  for ( int i = 0; i < size; i++ ){
    TacsScalar winv = 1.0/w[i];
    for ( int j = 0; j < vars_per_node; j++ ){
      u[j] *= winv;
    }
    u += vars_per_node; 
  }

  // Free the refined weights
  weights_refined->decref();

  // Distribute the values
  uvec_refined->beginDistributeValues();
  uvec_refined->endDistributeValues();

  // The solution was not passed as an argument
  if (!_uvec){
    uvec->decref();
  }

  // The refined solution vector was not passed as an argument
  if (!_uvec_refined){
    tacs_refined->setVariables(uvec_refined);
    uvec_refined->decref();
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
  forest:      the forest of quadtrees 
  target_err:  target absolute error value
  min_level:   minimum refinement level
  max_level:   maximum refinement level

  returns:     predicted strain energy error
*/
TacsScalar TMR_StrainEnergyRefine( TACSAssembler *tacs,
                                   TMRQuadForest *forest,
                                   double target_err,
                                   int min_level, int max_level ){
  // The maximum number of nodes
  const int max_num_nodes = 9;

  // Set the order of the mesh and the number of enrichment shape functions
  const int order = forest->getMeshOrder();
  const int nenrich = (order == 2 ? 5 : 7);

  // Get the communicator
  MPI_Comm comm = tacs->getMPIComm();
  
  // Retrieve the variables from the TACSAssembler object
  TACSBVec *uvec = tacs->createVec();
  tacs->getVariables(uvec);
  uvec->incref();
  uvec->beginDistributeValues();
  uvec->endDistributeValues();

  // Number of local elements
  const int nelems = tacs->getNumElements();

  // Get the number of variables per node
  const int vars_per_node = tacs->getVarsPerNode();
  const int deriv_per_node = 3*vars_per_node;

  // The number of equations: 2 times the number of nodes for each element
  const int neq = 2*order*order;

  // Compute the contribution to the error in the energy norm from
  // each element
  TacsScalar *SE_error = new TacsScalar[ nelems ];
  
  // Create the weight vector - the weights are the number of times
  // each node is referenced by adjacent elements, including
  // inter-process references.
  TACSBVec *weights = new TACSBVec(tacs->getVarMap(), 1,
                                   tacs->getBVecDistribute(),
                                   tacs->getBVecDepNodes());
  weights->incref();
  computeLocalWeights(tacs, weights);

  // Compute the nodal derivatives
  TACSBVec *uderiv = 
    new TACSBVec(tacs->getVarMap(), 3*vars_per_node,
                 tacs->getBVecDistribute(), tacs->getBVecDepNodes());
  uderiv->incref();
  computeNodeDeriv2D(order, tacs, uvec, weights, uderiv);
  weights->decref();

  // Allocate space for the element reconstruction problem
  TacsScalar *tmp = new TacsScalar[ neq*(nenrich + vars_per_node) ];
  TacsScalar *ubar = new TacsScalar[ vars_per_node*nenrich ];
  TacsScalar *delem = new TacsScalar[ deriv_per_node*order*order ];

  // Allocate arrays needed for the reconstruction
  TacsScalar *vars_elem = new TacsScalar[ vars_per_node*max_num_nodes ];

  // The interpolated variables on the refined mesh
  TacsScalar *dvars = new TacsScalar[ vars_per_node*max_num_nodes ];
  TacsScalar *ddvars = new TacsScalar[ vars_per_node*max_num_nodes ];
  TacsScalar *vars_interp = new TacsScalar[ vars_per_node*max_num_nodes ];

  // Allocate the element residual array
  TacsScalar *res = new TacsScalar[ vars_per_node*max_num_nodes ];

  // Keep track of the total error
  TacsScalar SE_total_error = 0.0;

  for ( int i = 0; i < nelems; i++ ){
    // The simulation time -- we assume time-independent analysis
    double time = 0.0;

    // Get the node locations and variables
    TacsScalar Xpts[3*max_num_nodes];
    TACSElement *elem = tacs->getElement(i, Xpts, vars_elem, dvars, ddvars);
    
    // Get the node numbers for this element
    int len;
    const int *nodes;
    tacs->getElement(i, &nodes, &len);

    // Evaluate the element residual
    memset(res, 0, elem->numVariables()*sizeof(TacsScalar));
    elem->addResidual(time, res, Xpts, vars_elem, dvars, ddvars);
    
    // Take the inner product to get the strain energy
    SE_error[i] = 0.0;
    for ( int j = 0; j < elem->numVariables(); j++ ){
      SE_error[i] += res[j]*vars_elem[j];
    }

    // Compute the solution on the refined mesh
    uderiv->getValues(len, nodes, delem);

    // Compute the enrichment functions for each degree of freedom
    computeElemRecon2D(order, vars_per_node,
                       Xpts, vars_elem, delem, ubar, tmp);

    TacsScalar SE_refine = 0.0;
    for ( int jj = 0; jj < 2; jj++ ){
      for ( int ii = 0; ii < 2; ii++ ){
        // Retrieve the refined element
        TacsScalar rXpts[3*max_num_nodes];
        memset(rXpts, 0, 3*order*order*sizeof(TacsScalar));

        // Set the variables to zero
        memset(vars_interp, 0, vars_per_node*order*order*sizeof(TacsScalar));

        for ( int m = 0; m < order; m++ ){
          for ( int n = 0; n < order; n++ ){
            double pt[2];
            pt[0] = ii - 1.0 + 1.0*n/(order-1);
            pt[1] = jj - 1.0 + 1.0*m/(order-1);

            // Evaluate the locations of the new nodes
            double N[9], Nr[7];
            FElibrary::biLagrangeSF(N, pt, order);
            if (order == 2){
              eval2ndEnrichmentFuncs2D(pt, Nr);
            }
            else {
              eval3rdEnrichmentFuncs2D(pt, Nr);
            }
            
            // Set the values of the variables at this point
            for ( int k = 0; k < order*order; k++ ){
              rXpts[3*(n + m*order)] += Xpts[3*k]*N[k];
              rXpts[3*(n + m*order)+1] += Xpts[3*k+1]*N[k];
              rXpts[3*(n + m*order)+2] += Xpts[3*k+2]*N[k];
            }

            // Evaluate the interpolation part of the reconstruction
            for ( int k = 0; k < order*order; k++ ){
              for ( int kk = 0; kk < vars_per_node; kk++ ){
                vars_interp[vars_per_node*(n + m*order)+kk] += 
                  vars_elem[vars_per_node*k+kk]*N[k];
              }
            }

            // Add the portion from the enrichment functions
            for ( int k = 0; k < nenrich; k++ ){
              // Evaluate the interpolation part of the reconstruction
              for ( int kk = 0; kk < vars_per_node; kk++ ){
                vars_interp[vars_per_node*(n + m*order)+kk] += 
                  ubar[vars_per_node*k+kk]*Nr[k];
              }
            }
          }
        }

        // Compute the element residual
        memset(res, 0, elem->numVariables()*sizeof(TacsScalar));
        elem->addResidual(time, res, rXpts, 
                          vars_interp, dvars, ddvars);
        
        // Take the inner product to get the strain energy
        for ( int j = 0; j < elem->numVariables(); j++ ){
          SE_refine += res[j]*vars_interp[j];
        }
      }
    }

    // SE_refine - SE_error should always be a positive quantity
    SE_error[i] = SE_refine - SE_error[i];

    // Add up the total error
    SE_total_error += SE_error[i];
  }
  
  // Free the global vectors
  uvec->decref();
  uderiv->decref();

  // Free the element-related data
  delete [] tmp;
  delete [] ubar;
  delete [] delem;
  delete [] vars_elem;
  delete [] dvars;
  delete [] ddvars;
  delete [] vars_interp;
  delete [] res;

  // Count up the total strain energy 
  TacsScalar SE_temp = 0.0;
  MPI_Allreduce(&SE_total_error, &SE_temp, 1, TACS_MPI_TYPE, MPI_SUM, comm);
  SE_total_error = SE_temp;

  // Go through and flag which element should be refined
  int *refine = new int[ nelems ];
  memset(refine, 0, nelems*sizeof(int));
  for ( int i = 0; i < nelems; i++ ){
    if (SE_error[i] >= target_err){
      refine[i] = 1;
    }
  }

  forest->refine(refine, min_level, max_level);

  // Free the data that is no longer required
  delete [] SE_error;

  // refine the quadrant mesh based on the local refinement values
  delete [] refine;

  // Return the error
  return SE_total_error;
}

/*
  Refine the mesh using the original solution and the adjoint solution

  input:
  tacs:          the TACSAssembler object
  tacs_refine:   the uniformly refined TACSAssembler object
  adjvec:        the adjoint solution variables
  forest:        the forest of quadtrees 
  target_err:    absolute value of the target error
  min_level:     minimum refinement
  max_level:     maximum refinement

  output:
  adj_corr:      adjoint-based functional correction

  returns: 
  absolute functional error estimate
*/
TacsScalar TMR_AdjointRefine( TACSAssembler *tacs,
                              TACSAssembler *tacs_refine,
                              TACSBVec *adjoint,
                              TMRQuadForest *forest,
                              double target_error,
                              int min_level, int max_level,
                              TacsScalar *adj_corr ){
  // The maximum number of nodes
  const int max_num_nodes = 9;

  // Get the order of the mesh and the number of enrichment shape functions
  const int order = forest->getMeshOrder();
  const int nenrich = (order == 2 ? 5 : 7);

  // Get the communicator
  MPI_Comm comm = tacs->getMPIComm();
  
  // Perform a local refinement of the nodes based on the strain energy
  // within each element
  const int nelems = tacs->getNumElements();

  // Get the number of variables per node 
  const int vars_per_node = tacs->getVarsPerNode();
  
  // Create the refined residual vector
  TACSBVec *adjoint_refine = tacs_refine->createVec();
  adjoint_refine->incref();

  // Allocate the refinement array - which elements will be refined
  int *refine = new int[ nelems ];
  memset(refine, 0, nelems*sizeof(int));

  // Keep track of the total error remaining from each element
  // indicator and the adjoint error correction
  TacsScalar total_error_remain = 0.0;
  TacsScalar total_adjoint_corr = 0.0;

  // Allocate arrays needed for the reconstruction
  TacsScalar *vars_elem = new TacsScalar[ vars_per_node*max_num_nodes ];
  TacsScalar *adj_elem = new TacsScalar[ vars_per_node*max_num_nodes ];

  // The interpolated variables on the refined mesh
  TacsScalar *dvars = new TacsScalar[ vars_per_node*max_num_nodes ];
  TacsScalar *ddvars = new TacsScalar[ vars_per_node*max_num_nodes ];
  TacsScalar *vars_interp = new TacsScalar[ vars_per_node*max_num_nodes ];
  TacsScalar *adj_elem_interp = new TacsScalar[ vars_per_node*max_num_nodes ];
  TacsScalar *adj_elem_refine = new TacsScalar[ vars_per_node*max_num_nodes ];

  // Allocate the element residual array
  TacsScalar *res = new TacsScalar[ vars_per_node*max_num_nodes ];

  // Get the auxiliary elements (surface tractions) associated with the
  // element class
  TACSAuxElements *aux_elements = tacs_refine->getAuxElements();
  int num_aux_elems = 0;
  TACSAuxElem *aux = NULL;
  if (aux_elements){
    aux_elements->sort();
    num_aux_elems = aux_elements->getAuxElements(&aux);
  }

  // Reconstruct the adjoint solution on the finer mesh
  TMR_ComputeReconSolution(forest, tacs, tacs_refine,
                           adjoint, adjoint_refine);
  
  // For each element in the mesh, compute the residual on the refined
  // mesh based on its and store its value in the global solution
  int aux_count = 0;
  for ( int elem = 0; elem < nelems; elem++ ){
    // Set the simulation time
    double time = 0.0;

    // Get the element variable values on the coarse mesh
    tacs->getElement(elem, NULL, vars_elem);

    // Get the element node numbers
    int elem_len = 0;
    const int *elem_nodes;
    tacs->getElement(elem, &elem_nodes, &elem_len);
        
    // Get the values of the adjoint on the coarse mesh
    adjoint->getValues(elem_len, elem_nodes, adj_elem);

    // Keep track of the remaining error from this element
    TacsScalar elem_error_remain = 0.0;

    // For each element on the refined mesh, compute the interpolated
    // solution and sum up the local contribution to the adjoint
    for ( int ii = 0; ii < 2; ii++ ){
      for ( int jj = 0; jj < 2; jj++ ){
        // Set the element number on the refined mesh
        int elem_num = 4*elem + jj + 2*ii;

        // Zero the interpolations
        memset(vars_interp, 0, vars_per_node*order*order*sizeof(TacsScalar));
        memset(adj_elem_interp, 0, vars_per_node*order*order*sizeof(TacsScalar));

        // Perform the interpolation
        for ( int m = 0; m < order; m++ ){
          for ( int n = 0; n < order; n++ ){
            double pt[2];
            pt[0] = ii - 1.0 + 1.0*n/(order-1);
            pt[1] = jj - 1.0 + 1.0*m/(order-1);

            // Evaluate the locations of the new nodes
            double N[max_num_nodes];
            FElibrary::biLagrangeSF(N, pt, order);
           
            // Evaluate the interpolation part of the reconstruction
            for ( int k = 0; k < order*order; k++ ){
              for ( int kk = 0; kk < vars_per_node; kk++ ){
                vars_interp[vars_per_node*(n + m*order)+kk] += 
                  vars_elem[vars_per_node*k+kk]*N[k];
              }
              for ( int kk = 0; kk < vars_per_node; kk++ ){
                adj_elem_interp[vars_per_node*(n + m*order)+kk] += 
                  adj_elem[vars_per_node*k+kk]*N[k];
              }
            }
          }
        }

        // Get the node locations and the velocity/acceleration
        // variables (these shuold be zero and are not used...)
        TacsScalar Xpts[3*max_num_nodes];
        TACSElement *elem = tacs_refine->getElement(elem_num, Xpts,
                                                    NULL, dvars, ddvars);

        // Compute the residual for this element on the refined mesh
        memset(res, 0, elem->numVariables()*sizeof(TacsScalar));
        elem->addResidual(time, res, Xpts, 
                          vars_interp, dvars, ddvars);

        while (aux_count < num_aux_elems && aux[aux_count].num == elem_num){
          aux[aux_count].elem->addResidual(time, res, Xpts, 
                                           vars_interp, dvars, ddvars);
          aux_count++;
        }

        // Get the element and node numbers for the refined mesh
        int len = 0;
        const int *nodes;
        tacs_refine->getElement(elem_num, &nodes, &len);

        // Get the adjoint variables for the refined mesh
        adjoint_refine->getValues(len, nodes, adj_elem_refine);
        
        // Add in the contribution to the error from this element
        for ( int j = 0; j < elem->numVariables(); j++ ){
          elem_error_remain += (adj_elem_refine[j] - adj_elem_interp[j])*res[j];
        }
      }
    }

    // Add the contribution to the total remaining error and the
    // adjoint correction.
    total_error_remain += fabs(elem_error_remain);
    total_adjoint_corr += elem_error_remain;

    // Refine the element if the element error exceeds the target error
    if (fabs(elem_error_remain) > target_error){
      refine[elem] = 1;
    }
  }

  // Sum up the contributions across all processors
  TacsScalar temp[2];
  temp[0] = total_error_remain;
  temp[1] = total_adjoint_corr;
  MPI_Allreduce(MPI_IN_PLACE, temp, 2, TACS_MPI_TYPE, MPI_SUM, comm);
  total_error_remain = temp[0];
  total_adjoint_corr = temp[1];

  // Set the adjoint reconstruction as the variables
  tacs_refine->setVariables(adjoint_refine);

  // Free the data that is no longer required
  delete [] vars_elem;
  delete [] adj_elem;
  delete [] dvars;
  delete [] ddvars;
  delete [] vars_interp;
  delete [] adj_elem_interp;
  delete [] adj_elem_refine;
  delete [] res;
  adjoint_refine->decref();

  // refine the quadrant mesh based on the local refinement values
  forest->refine(refine, min_level, max_level);
  delete [] refine;

  // Set the adjoint residual correction
  if (adj_corr){
    *adj_corr = total_adjoint_corr;
  }

  // Return the error
  return total_error_remain;
}
