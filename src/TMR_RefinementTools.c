#include "TMR_RefinementTools.h"
#include "TensorToolbox.h"
#include "tacslapack.h"

/*
  Compute the Jacobian transformation at a point within the element
*/
static void computeJacobianTrans( const TacsScalar Xpts[],
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
static void evalEnrichmentFuncs( const double pt[], double N[] ){
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
static void evalEnrichmentFuncs( const double pt[],
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
static void computeElemRecon( const TacsScalar Xpts[],
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
static void computeLocalWeights( TACSAssembler *tacs, 
                                 TACSBVec **wlocal ){
  // Create the weight vector - the weights on the averages
  // Set the weights on each of the nodes
  TACSBVec *weights = new TACSBVec(tacs->getVarMap(), 1,
                                   NULL, tacs->getBVecDistribute(),
                                   tacs->getBVecDepNodes());
  weights->incref();

  // Set the local element weights
  int max_nodes = tacs->getMaxElementNodes();
  TacsScalar *welem = new TacsScalar[ max_nodes ];
  for ( int i = 0; i < max_nodes; i++ ){
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

  delete [] welem;

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
static void computeNodeDeriv( TACSAssembler *tacs, TACSBVec *uvec, 
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
  // Get the communicator
  MPI_Comm comm = tacs->getMPIComm();
  
  // Retrieve the variables from the TACSAssembler object
  TACSBVec *uvec = tacs->createVec();
  tacs->getVariables(uvec);
  uvec->incref();

  // Perform a local refinement of the nodes based on the strain
  // energy within each element
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
  
  uvec->decref();
  uderiv->decref();
  wlocal->decref();

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

  // Free some of the data that is no longer required
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
                              TACSBVec *adjvec,
                              TMRQuadForest *forest,
                              double target_err,
                              int min_level, int max_level,
                              TacsScalar *adj_corr ){
  // Get the communicator
  MPI_Comm comm = tacs->getMPIComm();
  
  // Perform a local refinement of the nodes based on the strain energy
  // within each element
  int nelems = tacs->getNumElements();

  // Get the number of variables per node (and the number of
  // derivatives per node required for the reconstruction)
  int vars_per_node = tacs->getVarsPerNode();
  int deriv_per_node = 3*vars_per_node;
  
  // Create the refined residual vector
  TACSBVec *residual = tacs_refine->createVec();
  residual->incref();

  // Get the auxiliary elements (surface tractions) associated with the
  // element class
  TACSAuxElements *aux_elements = tacs_refine->getAuxElements();
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
        TACSElement *elem = tacs_refine->getElement(elem_num, Xpts);

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
        tacs_refine->getElement(elem_num, &nodes, &len);
        
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
  TACSBVec *qadjvec = tacs_refine->createVec();
  qadjvec->incref();

  // Compute the weights on the refined mesh
  computeLocalWeights(tacs_refine, &wlocal);
  wlocal->incref();

  // Allocate the refinement array - which elements will be refined
  int *refine = new int[ nelems ];
  memset(refine, 0, nelems*sizeof(int));

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
        TACSElement *elem = tacs_refine->getElement(elem_num);

        // Get the element node numbers
        tacs_refine->getElement(elem_num, &nodes, &len);

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
      refine[i] = 1;
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

  // refine the quadrant mesh based on the local refinement values
  forest->refine(refine, min_level, max_level);
  delete [] refine;

  // Set the adjoint residual correction
  if (adj_corr){
    *adj_corr = total_corr;
  }

  // Return the error
  return total_err_remain;
}

