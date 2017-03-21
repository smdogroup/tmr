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
  TMRCylinderCreator( double _alpha, double _beta, double _R,
                      double _load, FSDTStiffness *stiff ){
    alpha = _alpha;
    beta = _beta;
    R = _R;
    load = _load;
    elem2 = new MITCShell<2>(stiff); elem2->incref();
    elem3 = new MITCShell<3>(stiff); elem3->incref();
  }
  ~TMRCylinderCreator(){
    elem2->decref();
    elem3->decref();
  }

  TACSElement *createElement( int order,
                              TMRQuadForest *forest,
                              TMRQuadrant quad ){
    if (order == 2){ return elem2; }
    if (order == 3){ return elem3; }
  }
  TACSElement *createAuxElement( int order,
                                 TMRQuadForest *forest,
                                 TMRQuadrant quad ){
    // Get the points from the forest
    TMRPoint *Xp;
    forest->getPoints(&Xp);

    // Get the side-length of the quadrant
    const int32_t h = 1 << (TMR_MAX_LEVEL - quad.level - (order-2));

    // Set the tractions based on the node location
    TacsScalar tx[9], ty[9], tz[9];
    for ( int j = 0, n = 0; j < order; j++ ){
      for ( int i = 0; i < order; i++, n++ ){
        // Find the nodal index and determine the (x,y,z) location
        TMRQuadrant node;
        node.x = quad.x + h*i;
        node.y = quad.y + h*j;
        node.face = quad.face;
        forest->transformNode(&node);
        int index = forest->getNodeIndex(&node);
        
        // Set the pressure load
        double z = Xp[index].z;
        double theta =-R*atan2(Xp[index].y, Xp[index].x);
        double p = -load*sin(beta*z)*sin(alpha*theta);

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

    return trac;
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

  int orthotropic_flag = 0;
  double htarget = 10.0;
  for ( int i = 0; i < argc; i++ ){
    if (sscanf(argv[i], "h=%lf", &htarget) == 1){
      printf("htarget = %f\n", htarget);
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
  double load = 1.0e3;

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

  // Set up the creator object - this facilitates creating the
  // TACSAssembler objects for different geometries
  TMRCylinderCreator *creator = new TMRCylinderCreator(alpha, beta, R,
                                                       load, stiff);

  // Create the geometry
  BRepPrimAPI_MakeCylinder cylinder(R, L);
  TopoDS_Compound compound;
  BRep_Builder builder;
  builder.MakeCompound(compound);
  builder.Add(compound, cylinder.Shape());

  // Load in the geometry
  TMRModel *geo = TMR_LoadModelFromCompound(compound);

  int num_edges;
  TMREdge **edges;
  geo->getEdges(&num_edges, &edges);
  edges[0]->setAttribute("Edge1");
  edges[1]->setAttribute("Edge2");
  edges[2]->setAttribute("Edge3");

  // Set the boundary conditions
  int bcs[] = {0, 1, 2, 3, 4, 5};
  creator->addBoundaryCondition("Edge1", 6, bcs);
  creator->addBoundaryCondition("Edge3", 6, bcs);

  if (geo){
    geo->incref();
    // Allocate the new mesh
    TMRMesh *mesh = new TMRMesh(comm, geo);
    mesh->incref();
    
    TMRMeshOptions options;
    options.frontal_quality_factor = 1.25;
    mesh->mesh(options, htarget);
    mesh->writeToVTK("cylinder-mesh.vtk");

    TMRModel *model = mesh->createModelFromMesh();
    model->incref();

    const int MAX_REFINE = 5;
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
    }

    for ( int iter = 0; iter < MAX_REFINE; iter++ ){
      // Create the TACSAssembler object associated with this forest
      TACSAssembler *tacs[MAX_REFINE+2];
      forest[0]->balance();
      forest[0]->repartition();
      tacs[0] = creator->createTACS(3, forest[0]);
      tacs[0]->incref();

      // Create the coarser versions of tacs
      forest[1] = forest[0]->duplicate();
      forest[1]->balance();
      forest[1]->incref();
      tacs[1] = creator->createTACS(2, forest[1]);
      tacs[1]->incref();

      // Set the number of levels: Cap it with a value of 4
      int num_levels = 2 + iter;
      if (num_levels > 4){ num_levels = 4; }

      // Create all of the coarser levels
      for ( int i = 2; i < num_levels; i++ ){
        forest[i] = forest[i-1]->coarsen();
        forest[i]->incref();
        forest[i]->balance();
        tacs[i] = creator->createTACS(2, forest[i]);
        tacs[i]->incref();
      }

      // Create the interpolation
      TACSBVecInterp *interp[MAX_REFINE+1];

      for ( int level = 0; level < num_levels-1; level++ ){
        // Create the interpolation object
        interp[level] = new TACSBVecInterp(tacs[level+1]->getVarMap(),
                                           tacs[level]->getVarMap(),
                                           tacs[level]->getVarsPerNode());
        interp[level]->incref();

        // Set the interpolation
        forest[level]->createInterpolation(forest[level+1], interp[level]);
    
        // Initialize the interpolation
        interp[level]->initialize();
      }

      // Create the multigrid object
      double omega = 1.0;
      int mg_sor_iters = 1;
      int mg_sor_symm = 1;
      int mg_iters_per_level = 1;
      TACSMg *mg = new TACSMg(comm, num_levels, omega, 
                              mg_sor_iters, mg_sor_symm);
      mg->incref();

      for ( int level = 0; level < num_levels; level++ ){
        if (level < num_levels-1){
          mg->setLevel(level, tacs[level], interp[level], mg_iters_per_level);
        }
        else {
          mg->setLevel(level, tacs[level], NULL);
        }
      }
      
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

      // The compliance function
      TACSFunction *func = new TACSCompliance(tacs[0]);
      func->incref();

      // Evaluate the function of interest
      TacsScalar fval = 0.0;
      tacs[0]->evalFunctions(&func, 1, &fval);

      // Delete the compliance function
      func->decref();

      // Create and write out an fh5 file
      unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                                 TACSElement::OUTPUT_DISPLACEMENTS);
      TACSToFH5 *f5 = new TACSToFH5(tacs[0], SHELL, write_flag);
      f5->incref();
      
      // Write out the solution
      char outfile[128];
      sprintf(outfile, "output%02d.f5", iter);
      f5->writeToFile(outfile);
      f5->decref();

      // Duplicate and refine the mesh
      TMRQuadForest *dup = forest[0]->duplicate();
      dup->incref();
      dup->refine(NULL);
      dup->balance();
      TACSAssembler *dup_tacs = creator->createTACS(3, dup);
      dup_tacs->incref();
      TMR_ComputeReconSolution(3, tacs[0], dup_tacs);

      // Create and write out an fh5 file
      TACSToFH5 *f5_dup = new TACSToFH5(dup_tacs, SHELL, write_flag);
      f5_dup->incref();
      
      // Write out the solution
      sprintf(outfile, "refined_output%02d.f5", iter);
      f5_dup->writeToFile(outfile);
      f5_dup->decref();

      dup->decref();
      dup_tacs->decref();

      // Set the target error using the factor: max(1.0, 16*(2**-iter))
      double factor = 1.0; // 16.0/(1 << iter);
      if (factor < 1.0){ factor = 1.0; }

      // Determine the total number of elements across all processors
      int nelems_total = 0;
      int nelems = tacs[0]->getNumElements();
      MPI_Allreduce(&nelems, &nelems_total, 1, MPI_INT, MPI_SUM, comm);

      // Set the absolute element-level error based on the relative
      // error that is requested
      double target_abs_err = factor*target_rel_err*fval/nelems_total;

      // Perform the strain energy refinement
      TacsScalar abs_err = TMR_StrainEnergyRefine(tacs[0], forest[0], 
                                                  target_abs_err);

      if (mpi_rank == 0){
        printf("Absolute error estimate = %e\n", abs_err);

        const int *range;
        forest[0]->getOwnedNodeRange(&range);
        int nnodes = range[mpi_size-1];
        fprintf(fp, "%d %d %d %15.10e %15.10e %15.10e\n",
                iter, nelems, nnodes, fval, (fval + abs_err)/fval, abs_err);
      }

      // Decrease the reference count
      tacs[0]->decref();
      for ( int i = 1; i < num_levels; i++ ){
        forest[i]->decref();
        tacs[i]->decref();
        interp[i-1]->decref();
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
