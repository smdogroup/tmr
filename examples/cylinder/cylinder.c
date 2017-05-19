#include "TMRGeometry.h"
#include "TMRMesh.h"
#include "TMR_TACSCreator.h"
#include "TMRQuadForest.h"
#include "TACSAssembler.h"
#include "Compliance.h"
#include "specialFSDTStiffness.h"
#include "MITCShell.h"
#include "TACSShellTraction.h"
#include "TACSToFH5.h"
#include "TACSMg.h"
#include "tacslapack.h"
#include "KSFailure.h"
#include <stdio.h>
#include <math.h>
#include "TMR_RefinementTools.h"

#ifdef TMR_HAS_OPENCASCADE

#include "TMROpenCascade.h"
#include <BRepPrimAPI_MakeCylinder.hxx>

/*
  Get the element and the auxiliary elements
*/
class TMRCylinderCreator : public TMRQuadTACSCreator {
 public:
  TMRCylinderCreator( MPI_Comm comm,
                      TMRBoundaryConditions *_bcs,
                      double _alpha, double _beta, double _R,
                      double _load, FSDTStiffness *stiff ):
  TMRQuadTACSCreator(_bcs){
    int rank;
    MPI_Comm_rank(comm, &rank);
    alpha = _alpha;
    beta = _beta;
    R = _R;
    load = _load;
    elem2 = new MITCShell<2>(stiff); elem2->incref();
    elem3 = new MITCShell<3>(stiff); elem3->incref();
    // elem2->setComponentNum(rank);
    // elem3->setComponentNum(rank);
  }
  ~TMRCylinderCreator(){
    elem2->decref();
    elem3->decref();
  }

  /*
    Create/set all of the elements within the new TACSAssembler object
  */
  void createElements( int order,
                       TMRQuadForest *forest,
                       int num_elements,
                       TACSElement **elements ){
    if (order == 2){
      for ( int i = 0; i < num_elements; i++ ){
        elements[i] = elem2;
      }
    }
    else if (order == 3){
      for ( int i = 0; i < num_elements; i++ ){
        elements[i] = elem3;
      }
    }
  }

  /*
    Create all of the auxiliary elements within TACS
  */ 
  TACSAuxElements *createAuxElements( int order,
                                      TMRQuadForest *forest ){
    TACSAuxElements *aux = new TACSAuxElements();
    
    // Get the quadrants
    TMRQuadrantArray *quadrants;
    forest->getQuadrants(&quadrants);

    // get the quadrant array
    int size; 
    TMRQuadrant *quads;
    quadrants->getArray(&quads, &size);

    // Get the points from the forest
    TMRPoint *Xp;
    forest->getPoints(&Xp);

    for ( int i = 0; i < size; i++ ){
      // Get the side-length of the quadrant
      const int32_t h = 1 << (TMR_MAX_LEVEL - quads[i].level - (order-2));

      // Set the tractions based on the node location
      TacsScalar tx[9], ty[9], tz[9];
      for ( int jj = 0, n = 0; jj < order; jj++ ){
        for ( int ii = 0; ii < order; ii++, n++ ){
          // Find the nodal index and determine the (x,y,z) location
          TMRQuadrant node;
          node.face = quads[i].face;
          node.x = quads[i].x + h*ii;
          node.y = quads[i].y + h*jj;
          forest->transformNode(&node);
          int index = forest->getNodeIndex(&node);
        
          // Set the pressure load
          double z = Xp[index].z;
          double y = -R*atan2(Xp[index].y, Xp[index].x);
          double p = -load*sin(beta*z)*sin(alpha*y);

          tx[n] = p*Xp[index].x/R;
          ty[n] = p*Xp[index].y/R;
          tz[n] = 0.0;
        }
      }
      
      // Create the traction class 
      TACSElement *trac = NULL;
      if (order == 2){
        trac = new TACSShellTraction<2>(tx, ty, tz);
      }
      else {
        trac = new TACSShellTraction<3>(tx, ty, tz);
      }

      aux->addElement(i, trac);
    }

    return aux;
  }

 private:
  double alpha, beta;
  double R;
  double load;
  
  TACSElement *elem2, *elem3;
};

/*
  Compute the coefficients of a single term in a Fourier series for a
  specially orthotropic cylinder subjectedt to a sinusoidally varying
  pressure distribution specified as follows:

  p = sin(alpha*y)*cos(beta*x)

  Note that y = r*theta

  The coefficients U, V, W, theta and phi are for the expressions:
  
  u(x,y) = U*sin(alpha*y)*cos(beta*x)
  v(x,y) = V*cos(alpha*y)*sin(beta*x)
  w(x,y) = W*sin(alpha*y)*sin(beta*x)
  psi_x(x,y) = theta*sin(alpha*y)*cos(beta*x)
  psi_y(x,y) = phi*cos(alpha*y)*sin(beta*x)

  Note that u, v and w correspond to the axial, tangential and normal
  displacements in a reference frame attached to the shell surface -
  they are not the displacements expressed in the global reference
  frame. Likewise, psi_x and psi_y are the rotations of the normal
  along the x-direction and the tangential y-direction.
*/
void compute_coefficients( double *U, double *V, double *W, 
                           double *theta, double *phi,
                           double alpha, double beta, double ainv, 
                           double A11, double A12, double A22, double A33,
                           double D11, double D12, double D22, double D33,
                           double bA11, double bA22, double load ){

  double A[5*5], A2[5*5]; // 5 x 5 system of linear equations
  int ipiv[5];
  double rhs[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  rhs[2] = - load;
 
  double B[8*5];
  memset(B, 0, sizeof(B));
  double * b = B;

  // Assign the columns for u
  b[0] = -beta;
  b[2] = alpha;
  b[5] = -alpha*ainv;
  b += 8;

  // Assign columns for v
  b[1] = -alpha;
  b[2] = beta;
  b[4] = alpha*ainv;
  b[6] = -ainv;
  b += 8;

  // Assign the columns for w
  b[1] = ainv;
  b[4] = -ainv*ainv;
  b[6] = alpha;
  b[7] = beta;
  b += 8;

  // Assign the columns for psi_x
  b[3] = -beta;
  b[5] = alpha;
  b[7] = 1.0;
  b += 8;

  // Assign the columns for psi_y
  b[4] = -alpha; 
  b[5] = beta;
  b[6] = 1.0;

  for ( int j = 0; j < 5; j++ ){
    double * bj = &B[8*j];
    for ( int i = 0; i < 5; i++ ){
      double * bi = &B[8*i];

      A2[i + 5*j] = 
        - ((bi[0]*(A11*bj[0] + A12*bj[1]) +  
            bi[1]*(A12*bj[0] + A22*bj[1]) + bi[2]*A33*bj[2]) +
           (bi[3]*(D11*bj[3] + D12*bj[4]) +
            bi[4]*(D12*bj[3] + D22*bj[4]) + bi[5]*D33*bj[5]) +
           bi[6]*bA11*bj[6] + bi[7]*bA22*bj[7]);
    }
  }

  // The first equation for u
  A[0]  = -(A11*beta*beta + A33*alpha*alpha) - D33*ainv*ainv*alpha*alpha;
  A[5]  = -(A33 + A12)*alpha*beta;
  A[10] = A12*beta*ainv;
  A[15] = D33*ainv*alpha*alpha;
  A[20] = D33*ainv*alpha*beta;

  // The second equation for v
  A[1]  = -(A12 + A33)*alpha*beta;
  A[6]  = -(A33*beta*beta + A22*alpha*alpha) - ainv*ainv*bA11 - D22*ainv*ainv*alpha*alpha;
  A[11] = (A22 + bA11)*ainv*alpha + D22*alpha*ainv*ainv*ainv;
  A[16] = D12*ainv*alpha*beta;
  A[21] = bA11*ainv + D22*ainv*alpha*alpha;

  // The third equation for w
  A[2]  = A12*beta*ainv;
  A[7]  = (bA11 + A22)*alpha*ainv + D22*alpha*ainv*ainv*ainv;
  A[12] = -(bA11*alpha*alpha + bA22*beta*beta) - A22*ainv*ainv - D22*ainv*ainv*ainv*ainv;
  A[17] = -bA22*beta - D12*beta*ainv*ainv;
  A[22] = -bA11*alpha - D22*alpha*ainv*ainv;

  // Fourth equation for theta
  A[3]  = D33*ainv*alpha*alpha;
  A[8]  = D12*ainv*alpha*beta;
  A[13] = -bA22*beta - D12*beta*ainv*ainv;
  A[18] = -(D11*beta*beta + D33*alpha*alpha) - bA22;
  A[23] = -(D12 + D33)*alpha*beta;

  // Fifth equation for phi
  A[4]  =  D33*ainv*alpha*beta;
  A[9]  =  bA11*ainv + D22*ainv*alpha*alpha;
  A[14] = -bA11*alpha - D22*alpha*ainv*ainv;
  A[19] = -(D33 + D12)*alpha*beta;
  A[24] = -(D33*beta*beta + D22*alpha*alpha) - bA11;

  int info = 0;
  int n = 5;
  int nrhs = 1;
  LAPACKgesv(&n, &nrhs, A2, &n, ipiv, rhs, &n, &info);

  // Solve for the coefficients 
  *U = rhs[0];
  *V = rhs[1];
  *W = rhs[2];
  *theta = rhs[3];
  *phi = rhs[4];
}

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  MPI_Comm comm = MPI_COMM_WORLD;
  int mpi_size, mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  int order = 2;
  int orthotropic_flag = 0;
  double htarget = 10.0;
  int strain_energy_refine = 0;
  for ( int i = 0; i < argc; i++ ){
    if (sscanf(argv[i], "h=%lf", &htarget) == 1){
      printf("htarget = %f\n", htarget);
    }
    if (sscanf(argv[i], "order=%d", &order) == 1){
      if (order < 2){ order = 2; }
      if (order > 3){ order = 3; }
    }
  }

  // Set the shell geometry parameters
  double t = 1.0;
  double L = 100.0;
  double R = 100.0/M_PI;

  // Set the alpha/beta parameters
  double alpha = 4.0/R;
  double beta = 3*M_PI/L;

  // Set the load parameter
  double load = 3e6;

  // Set the yield stress
  double ys = 350e6;

  OrthoPly *ply = NULL;
  if (orthotropic_flag){
    // Set the material properties to use
    double rho = 1.0;
    double E1 = 100.0e9;
    double E2 =   5.0e9;
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
  }
  else {
    // Set the material properties to use
    double rho = 2700.0;
    double E = 70e9;
    double nu = 0.3;
  
    ply = new OrthoPly(t, rho, E, nu, ys);
  }

  // Create the stiffness relationship
  double kcorr = 5.0/6.0;
  FSDTStiffness *stiff = 
    new specialFSDTStiffness(ply, orthotropic_flag, t, kcorr);
  stiff->incref();
  ply->incref();

  // Set the boundary conditions
  int clamped_bc_nums[] = {0, 1, 2, 5};
  int restrained_bc_nums[] = {0, 1, 5};
  TMRBoundaryConditions *bcs = new TMRBoundaryConditions();
  bcs->addBoundaryCondition("Clamped", 4, clamped_bc_nums);
  bcs->addBoundaryCondition("Restrained", 3, restrained_bc_nums);

  // Set up the creator object - this facilitates creating the
  // TACSAssembler objects for different geometries
  TMRCylinderCreator *creator = new TMRCylinderCreator(comm, bcs, alpha, beta, R,
                                                       load, stiff);

  // Create the geometry
  BRepPrimAPI_MakeCylinder cylinder(R, L);
  TopoDS_Compound compound;
  BRep_Builder builder;
  builder.MakeCompound(compound);
  builder.Add(compound, cylinder.Shape());

  // Load in the geometry
  TMRModel *geo = TMR_LoadModelFromCompound(compound);

  if (geo){
    geo->incref();

    int num_verts;
    TMRVertex **verts;
    geo->getVertices(&num_verts, &verts);
    verts[0]->setAttribute("Clamped");
    verts[1]->setAttribute("Restrained");

    int num_edges;
    TMREdge **edges;
    geo->getEdges(&num_edges, &edges);
    edges[0]->setAttribute("Restrained");
    edges[2]->setAttribute("Restrained");

    // Only include the cylinder face
    int num_faces;
    TMRFace **faces;
    geo->getFaces(&num_faces, &faces);
    num_faces = 1;

    TMRModel *face_geo = new TMRModel(num_verts, verts,
                                      num_edges, edges,
                                      num_faces, &faces[0]);

    // Allocate the new mesh
    TMRMesh *mesh = new TMRMesh(comm, face_geo);
    mesh->incref();
    
    TMRMeshOptions options;
    options.frontal_quality_factor = 1.25;
    options.mesh_type_default = TMR_UNSTRUCTURED;
    mesh->mesh(options, htarget);
    mesh->writeToVTK("cylinder-mesh.vtk");

    TMRModel *model = mesh->createModelFromMesh();
    model->incref();

    const int MAX_REFINE = 4;
    TMRQuadForest *forest[MAX_REFINE+2];
    forest[0] = new TMRQuadForest(comm);
    forest[0]->incref();

    TMRTopology *topo = new TMRTopology(comm, model);
    forest[0]->setTopology(topo);
    forest[0]->createTrees(0);
    forest[0]->repartition();
    
    // The target relative error on the compliance
    double target_rel_err = 1e-4;

    FILE *fp = NULL;
    if (mpi_rank == 0){
      fp = fopen("cylinder_refine.dat", "w");
      fprintf(fp, "Variables = \"%s\" \"%s\" \"%s\" \"%s\" \"%s\" \"%s\"\n",
              "iter", "nelems", "nnodes", "fval", "fval corr", "abs_err");
    }

    for ( int iter = 0; iter < MAX_REFINE; iter++ ){
      // Create the TACSAssembler object associated with this forest
      TACSAssembler *tacs[MAX_REFINE+2];
      forest[0]->balance(1);
      forest[0]->repartition();
      tacs[0] = creator->createTACS(order, forest[0]);
      tacs[0]->incref();

      // Set the number of levels: Cap it with a value of 4
      int num_levels = 2 + iter;
      if (num_levels > 4){ num_levels = 4; }

      if (order == 3){
        // Create the coarser versions of tacs
        forest[1] = forest[0]->duplicate();
        forest[1]->incref();
        forest[1]->balance();
        tacs[1] = creator->createTACS(2, forest[1]);
        tacs[1]->incref();
      }

      // Create all of the coarser levels
      for ( int i = order-1; i < num_levels; i++ ){
        forest[i] = forest[i-1]->coarsen();
        forest[i]->incref();
        forest[i]->balance();
        tacs[i] = creator->createTACS(2, forest[i]);
        tacs[i]->incref();
      }

      // Create the multigrid object for TACS
      TACSMg *mg;
      TMR_CreateTACSMg(num_levels, tacs, forest, &mg);
      mg->incref();
      
      // Create the vectors 
      TACSBVec *ans = tacs[0]->createVec();
      TACSBVec *res = tacs[0]->createVec();
      ans->incref();
      res->incref();

      // Assemble the matrix
      mg->assembleJacobian(1.0, 0.0, 0.0, res);
      mg->factor();
    
      // Set up the solver
      int gmres_iters = 100; 
      int nrestart = 1;
      int is_flexible = 1;
      GMRES *gmres = new GMRES(mg->getMat(0), mg, 
                               gmres_iters, nrestart, is_flexible);
      gmres->incref();
      gmres->setMonitor(new KSMPrintStdout("GMRES", mpi_rank, 10));
      gmres->setTolerances(1e-10, 1e-30);

      // Get solution and store in ans
      gmres->solve(res, ans);
      ans->scale(-1.0);
      tacs[0]->setVariables(ans);

      // Evaluate the function of interest
      TacsScalar fval = -res->dot(ans);

      // Create and write out an fh5 file
      unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                                 TACSElement::OUTPUT_DISPLACEMENTS |
                                 TACSElement::OUTPUT_EXTRAS);
      TACSToFH5 *f5 = new TACSToFH5(tacs[0], TACS_SHELL, write_flag);
      f5->incref();

      // Write out the solution
      char outfile[128];
      sprintf(outfile, "output%02d.f5", iter);
      f5->writeToFile(outfile);
      f5->decref();

      int nelems = tacs[0]->getNumElements();

      // Set the target error using the factor: max(1.0, 16*(2**-iter))
      double factor = 1.0;

      // Determine the total number of elements across all processors
      int nelems_total = 0;
      MPI_Allreduce(&nelems, &nelems_total, 1, MPI_INT, MPI_SUM, comm);

      // The target absolute error per element
      double target_abs_err = 0.0;

      // The absolute error estimate and the corrected function value
      // estimate with the adjoint error correction
      TacsScalar abs_err = 0.0, fval_est = 0.0;

      if (strain_energy_refine){
        // Set the absolute element-level error based on the relative
        // error that is requested
        target_abs_err = factor*target_rel_err*fabs(fval)/nelems_total;

        // Perform the strain energy refinement
        abs_err = TMR_StrainEnergyRefine(tacs[0], forest[0], 
                                         target_abs_err);
        fval_est = fval + abs_err;
      }
      else {
        double ks_weight = 100;
        TACSKSFailure *ks_func = new TACSKSFailure(tacs[0], ks_weight);
        ks_func->setKSFailureType(TACSKSFailure::CONTINUOUS);
        ks_func->incref();

        // Set the pointer to the TACS function
        TACSFunction *func = ks_func;

        // Evaluate the function of interest
        tacs[0]->evalFunctions(&func, 1, &fval);

        // Set the absolute element-level error based on the relative
        // error that is requested
        target_abs_err = factor*target_rel_err*fabs(fval)/nelems_total;

        // Create the adjoint vector
        TACSBVec *adjvec = tacs[0]->createVec();
        adjvec->incref();

        // Compute the derivative df/du
        res->zeroEntries();
        tacs[0]->addSVSens(1.0, 0.0, 0.0, &func, 1, &res);

        // Solve for the adjoint vector
        gmres->solve(res, adjvec);
        adjvec->scale(-1.0);

        // Duplicate the forest
        int min_level = 0, max_level = TMR_MAX_LEVEL;
        TMRQuadForest *forest_refine = forest[0]->duplicate();
        forest_refine->incref();

        // Uniformly refine the mesh everywhere
        forest_refine->refine();
        forest_refine->balance();

        // Create the new, refined version of TACS
        TACSAssembler *tacs_refine = creator->createTACS(order, forest_refine);
        tacs_refine->incref();

        // Perform the adjoint refinement
        TacsScalar adj_corr;
        abs_err = TMR_AdjointRefine(tacs[0], tacs_refine,
                                    adjvec, forest[0],
                                    target_abs_err, min_level, max_level,
                                    &adj_corr);

        // Compute the function estimate
        fval_est = fval + adj_corr;

        // Delete all the info that is no longer needed
        tacs_refine->decref();
        forest_refine->decref();
        adjvec->decref();
        func->decref();
      }

      if (mpi_rank == 0){
        printf("Relative factor         = %e\n", factor);
        printf("Function value          = %e\n", fval);
        printf("Target element error    = %e\n", target_abs_err);
        printf("Absolute error estimate = %e\n", abs_err);
        printf("Relative error estimate = %e\n", fabs(abs_err/fval));

        const int *range;
        forest[0]->getOwnedNodeRange(&range);
        int nnodes = range[mpi_size];
        fprintf(fp, "%d %d %d %15.10e %15.10e %15.10e\n",
                iter, nelems, nnodes, fval, fval_est, abs_err);
      }

      // Decrease the reference count
      tacs[0]->decref();
      for ( int i = 1; i < num_levels; i++ ){
        forest[i]->decref();
        tacs[i]->decref();
      }

      res->decref();
      ans->decref();
      mg->decref();
      gmres->decref();
    }

    if (fp){ fclose(fp); }

    forest[0]->decref();
    model->decref();
    mesh->decref();
    geo->decref();
  }

  TMRFinalize();
  MPI_Finalize();
  return (0);
}

#endif // TMR_HAS_OPENCASCADE
