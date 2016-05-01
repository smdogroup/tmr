#include "TMRForrest.h"
#include "isoFSDTStiffness.h"
#include "TACSCreator.h"
#include "TACSAssembler.h"
#include "PlaneStressQuad.h"
#include "TACSToFH5.h"
#include "TACSShellTraction.h"
#include "TensorToolbox.h"
#include "MITCShell.h"

// Create the nodal locations for the mesh
const double test_Xpts[] = {0.0, 0.0, 0.0,
                            1.0, 0.0, 0.0,
                            0.3, 0.7, 0.0, 
                            0.8, 0.2, 0.0, 
                            0.25, 0.21, 0.0,
                            0.75, 0.6, 0.0,
                            0.0, 1.0, 0.0,
                            1.0, 1.0, 0.0};

const int test_conn[] = {0, 1, 4, 3,
                         2, 4, 5, 3, 
                         6, 0, 2, 4, 
                         2, 5, 6, 7,
                         3, 1, 5, 7};

/*
  Retrieve the x/y location on the face based on the u/v coordinates
*/
void get_location( int face, double u, double v, 
                   double *X ){
  double N[4];
  N[0] = (1.0 - u)*(1.0 - v);
  N[1] = u*(1.0 - v);
  N[2] = (1.0 - u)*v;
  N[3] = u*v;

  X[0] = (N[0]*test_Xpts[3*test_conn[4*face]] +
          N[1]*test_Xpts[3*test_conn[4*face+1]] +
          N[2]*test_Xpts[3*test_conn[4*face+2]] +
          N[3]*test_Xpts[3*test_conn[4*face+3]]);
  X[1] = (N[0]*test_Xpts[3*test_conn[4*face]+1] +
          N[1]*test_Xpts[3*test_conn[4*face+1]+1] +
          N[2]*test_Xpts[3*test_conn[4*face+2]+1] +
          N[3]*test_Xpts[3*test_conn[4*face+3]+1]);
  X[2] = (N[0]*test_Xpts[3*test_conn[4*face]+2] +
          N[1]*test_Xpts[3*test_conn[4*face+1]+2] +
          N[2]*test_Xpts[3*test_conn[4*face+2]+2] +
          N[3]*test_Xpts[3*test_conn[4*face+3]+2]);
}

/*
  The enriched shape functions for the tensor-product quadratic
  Lagrange shape functions.
*/
void evalEnrichmentFuncs( double pt[],
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
  reconstruction over the element by solving a least-squares 
  problem
*/
void computeReconn( const TacsScalar Xpts[],
                    const TacsScalar uvals[],
                    const TacsScalar ux[],
                    const TacsScalar uy[],
                    const TacsScalar uz[],
                    TacsScalar ubar[] ){
  // Set up the least squares problem at the nodes
  TacsScalar A[18*7], b[18*6];
  int nrhs = 6;

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
      TacsScalar Xd[9];
      memset(Xd, 0, sizeof(Xd));

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
      TacsScalar J[9];
      FElibrary::jacobian3d(Xd, J);

      // Normalize the first direction
      TacsScalar d1[3], d2[3];
      d1[0] = Xd[0];  d1[1] = Xd[1];  d2[2] = Xd[2];
      Tensor::normalize3D(d1);

      // Compute d2 = n x d1
      Tensor::crossProduct3D(d2, &Xd[6], d1);

      // First, compute the contributions to the righ-hand-side. The b
      // vector contains the difference between the prescribed
      // derivative and the contribution to the derivative from the
      // quadratic shape function terms
      for ( int k = 0; k < nrhs; k++ ){
        b[18*k+c] = 
          wvals[ii]*wvals[jj]*(d1[0]*ux[0] + d1[1]*uy[0] + d1[2]*uz[0]);
        b[18*k+c+1] = 
          wvals[ii]*wvals[jj]*(d2[0]*ux[0] + d2[1]*uy[0] + d2[2]*uz[0]);
        ux++, uy++, uz++;
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
        d[0] = Nar[i]*J[3] + Nbr[i]*J[4];
        d[0] = Nar[i]*J[6] + Nbr[i]*J[7];

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
  int lwork = 5*18;
  TacsScalar work[5*18];

  // Using LAPACK, compute the least squares solution
  // ubar, res, rank, s = np.linalg.lstsq(A, b)  
  // LAPACKgless(&m, &n, &nrhs, A, &m, b, &m, s, &rcond, &rank, work, &lwork);
}

/*
  Compute reconstruction
*/
void compute( TACSAssembler *tacs, BVec *adj ){
  int size = tacs->getNumNodes() + tacs->getNumDependentNodes();

  
  // Create the weight vector - the weights on the averages
  // Set the weights on each of the nodes
  BVec *weights = new BVec(tacs->getVarMap(), 1);

  // Retrieve the distributed object from TACSAssembler
  BVecDistribute *vecDist = tacs->getBVecDistribute();

  // Allocate space for the local weights
  TacsScalar *wlocal = new TacsScalar[ size ];
  memset(wlocal, 0, size*sizeof(TacsScalar));

  // Set the local element weights - this assumes that all the
  // elements are bi-quadratic 9-noded elements
  TacsScalar w[9]; 
  for ( int i = 0; i < 9; i++ ){
    w[i] = 1.0;
  }

  // Get the number of elements
  int nelems = tacs->getNumElements();
  for ( int i = 0; i < nelems; i++ ){
    // Add all the values
    tacs->addValues(1, i, w, wlocal);
  }

  // Add in the dependent weight values
  tacs->addDependentResidual(1, wlocal);

  // Add all the weigths together
  vecDist->beginReverse(wlocal, weights, BVecDistribute::ADD);
  vecDist->endReverse(wlocal, weights, BVecDistribute::ADD);

  // Send the weights back
  vecDist->beginForward(weights, wlocal);
  vecDist->endForward(weights, wlocal);

  // Set the weights
  tacs->setDependentVariables(1, wlocal);
  
  // Set the local adjoint variables
  TacsScalar *adjlocal = new TacsScalar[ 6*size ];
  vecDist->beginForward(adj, adjlocal);
  vecDist->endForward(adj, adjlocal);
  tacs->setDependentVariables(6, adjlocal);

  // Now allocate the derivative vectors
  const int adj_per_node = 18;
  BVec *deriv = new BVec(tacs->getVarMap(), adj_per_node);

  // Allocate space for the local variables
  TacsScalar *derivlocal = new TacsScalar[ adj_per_node*size ];
  memset(derivlocal, 0, adj_per_node*size*sizeof(TacsScalar));

  // Perform the reconstruction for the local
  for ( int i = 0; i < nelems; i++ ){
    // Get the local element adjoint variables
    TacsScalar Xpts[3*9], advars[6*9];
    tacs->getValues(vars_per_node, i, adjlocal, advars);

    // Get the node locations for the element
    tacs->getElement(i, Xpts, NULL, NULL, NULL);
    
    // Compute the derivative of the 6 components of the adjoint
    // variables along the 3-coordinate directions
    TacsScalar adjderiv[adj_per_node*9];

    // Compute the contributions to the derivative from this side of
    // the element    
    for ( int jj = 0; jj < 3; jj++ ){
      for ( int ii = 0; ii < 3; ii++ ){
        double pt[2];
        pt[0] = -1.0 + 1.0*ii;
        pt[1] = -1.0 + 1.0*jj;

        double N[9], Na[9], Nb[9];

      }
    }

    // Add the values
    tacs->addValues(adj_per_node, adjderiv, derivlocal);
  }


  tacs->addDependentResidual(adj_per_node, derivlocal);
  
  // Add all the weigths together
  vecDist->beginReverse(derivlocal, , BVecDistribute::ADD);
  vecDist->endReverse(wlocal, weights, BVecDistribute::ADD);

  // Send the weights back
  vecDist->beginForward(weights, wlocal);
  vecDist->endForward(weights, wlocal);

  // Now we can reconstruct the adjoint solution on the full set of nodes


}

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);

  srand(time(NULL));

  const int ORDER = 3;

  // Set the communicator
  MPI_Comm comm = MPI_COMM_WORLD;

  // Get the MPI communicator rank
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Create the TACSCreator object
  int vars_per_node = 6;
  TACSCreator *creator = new TACSCreator(comm, vars_per_node);
  creator->incref();

  if (rank == 0){
    TMRQuadForrest *forrest = new TMRQuadForrest(MPI_COMM_SELF);
    
    // Set the connectivity
    int num_nodes = 8;
    int num_faces = 5;
    forrest->setConnectivity(num_nodes, 
                             test_conn, num_faces);
    
    // Allocate the trees (random trees for now)
    int nrand = 50;
    int min_refine = 0;
    int max_refine = 8;
    // forrest->createRandomTrees(nrand, min_refine, max_refine);
    forrest->createTrees(4);

    // Balance the forrest so that we can use it!
    forrest->balance(1);
    
    // Create the nodes
    forrest->createNodes(ORDER);
    
    // Extract the mesh
    int nnodes, nelems;
    int *elem_ptr, *elem_conn;
    forrest->getMesh(&nnodes, &nelems,
                     &elem_ptr, &elem_conn);
    
    // Set the element id numbers
    int *elem_id_nums = new int[ nelems ];
    memset(elem_id_nums, 0, nelems*sizeof(int));
    
    // Create the Xpts array
    double *Xpts = new double[ 3*nnodes ];
    
    // Fill in the positions
    TMRQuadtree **quad;
    forrest->getQuadtrees(&quad);
    
    int num_bcs = 0;
    int max_bcs = 10000;
    int *bc_nodes = new int[ max_bcs ];
    
    int elem = 0;
    for ( int face = 0; face < num_faces; face++ ){
      // Retrieve the node quadrants
      TMRQuadrantArray *nodes;
      quad[face]->getNodes(&nodes);
      
      // Get the array of node quadrants from this face
      int size;
      TMRQuadrant *array;
      nodes->getArray(&array, &size);
      
      // Iterate through and evaluate the x/y/z locations
      const double dh = 1.0/(1 << TMR_MAX_LEVEL);
      for ( int i = 0; i < size; i++ ){
        int node = array[i].tag;
        if (node >= 0){
          double u = dh*array[i].x;
          double v = dh*array[i].y;         
          get_location(face, u, v, &Xpts[3*node]);
          
          // Set the boundary conditions based on the spatial location
          if (((Xpts[3*node] < 1e-6 || Xpts[3*node] > 0.999999) ||
               (Xpts[3*node+1] < 1e-6 || Xpts[3*node+1] > 0.999999))
              && num_bcs < max_bcs){
            bc_nodes[num_bcs] = node;
            num_bcs++;
          }
        }
      }
    }
    
    // Set the connectivity
    creator->setGlobalConnectivity(nnodes, nelems,
                                   elem_ptr, elem_conn,
                                   elem_id_nums);
    
    // Set the boundary conditions
    creator->setBoundaryConditions(num_bcs, bc_nodes, NULL, NULL);
    
    // Set the nodal locations
    creator->setNodes(Xpts);
    
    // Set the dependent nodes
    int num_dep_nodes;
    const int *dep_ptr;
    const int *dep_conn;
    const double *dep_weights;
    forrest->getDependentNodes(&num_dep_nodes, &dep_ptr,
                               &dep_conn, &dep_weights);
    creator->setDependentNodes(num_dep_nodes, dep_ptr,
                               dep_conn, dep_weights);
    
    // Partition the mesh
    creator->partitionMesh();

    delete forrest;
  }

  // Set the reference direction - this defines the local x-axis
  TacsScalar axis[] = {1.0, 0.0, 0.0};

  // Create the material properties
  TacsScalar rho = 2700.0; 
  TacsScalar E = 70e9;
  TacsScalar nu = 0.3;
  TacsScalar kcorr = 5.0/6.0;
  TacsScalar ys = 350e6;
  TacsScalar thickness = 0.0075;
  isoFSDTStiffness *stiff = new isoFSDTStiffness(rho, E, nu, kcorr,
                                                 ys, thickness);
  stiff->setRefAxis(axis);

  // Create the shell element
  TACSElement *element = new MITCShell<ORDER>(stiff);
  element->incref();

  // Test the element implementation on the root processor
  if (rank == 0){
    TestElement *test = new TestElement(element);
    test->incref();
    test->setPrintLevel(2);
    test->testResidual();
    test->testJacobian();
    test->decref();
  }

  // Set the elements internally within TACS
  creator->setElements(&element, 1);
  
  // Create the TACSAssembler object
  TACSAssembler *tacs = creator->createTACS();
  tacs->incref();

  // Dereference the element
  element->decref();

  // Create the traction class 
  TacsScalar tx = 0.0, ty = 0.0, tz = 100.0e3;
  TACSElement *trac = new TACSShellTraction<ORDER>(tx, ty, tz);

  // Create the auxiliary element class
  int nelems = tacs->getNumElements();
  TACSAuxElements *aux = new TACSAuxElements(nelems);
  for ( int i = 0; i < nelems; i++ ){
    aux->addElement(i, trac);
  }
  tacs->setAuxElements(aux);

  // Create the preconditioner
  BVec *res = tacs->createVec();
  BVec *ans = tacs->createVec();
  BVec *tmp = tacs->createVec();
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
  tacs->assembleJacobian(res, mat, alpha, beta, gamma);
  pc->factor();
  pc->applyFactor(res, ans);
  tacs->setVariables(ans);

  mat->mult(ans, tmp);
  tmp->axpy(-1.0, res);

  // Print the solution norm
  TacsScalar norm = tmp->norm()/res->norm();
  if (rank == 0){
    printf("Solution residual norm: %15.5e\n", norm);
  }

  // Create an TACSToFH5 object for writing output to files
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS |
                             TACSElement::OUTPUT_STRAINS |
                             TACSElement::OUTPUT_STRESSES |
                             TACSElement::OUTPUT_EXTRAS);
  // TACSToFH5 * f5 = new TACSToFH5(tacs, PLANE_STRESS, write_flag);
  TACSToFH5 * f5 = new TACSToFH5(tacs, SHELL, write_flag);
  f5->incref();
  f5->writeToFile("forrest.f5");

  // Free everything
  f5->decref();

  tacs->decref();
  creator->decref();

  MPI_Finalize();
  return (0);
}
