#include "TMRGeometry.h"
#include "TMRBspline.h"
#include "TMRMesh.h"
#include "TMRQuadForest.h"
#include "TACSAssembler.h"
#include "PlaneStressQuad.h"
#include "TACSToFH5.h"
#include <stdio.h>
#include <math.h>

#ifdef TMR_HAS_OPENCASCADE

#include "TMROpenCascade.h"
#include <BRepPrimAPI_MakeCylinder.hxx>

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  MPI_Comm comm = MPI_COMM_WORLD;
  int mpi_size, mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  double htarget = 10.0;
  for ( int i = 0; i < argc; i++ ){
    if (sscanf(argv[i], "h=%lf", &htarget) == 1){
      printf("htarget = %f\n", htarget);
    }
  }

  // Create the cylinder and load the compund and mesh it
  double R = 100.0;
  double H = 200.0;
  BRepPrimAPI_MakeCylinder cylinder(R, H);
  TopoDS_Compound compound;
  BRep_Builder builder;
  builder.MakeCompound(compound);
  builder.Add(compound, cylinder.Shape());

  // Allocate the stiffness object
  TacsScalar rho = 2570.0, E = 70e9, nu = 0.3;
  PlaneStressStiffness *stiff = 
    new PlaneStressStiffness(rho, E, nu);

  // Allocate the solid element class
  TACSElement *elem = new PlaneStressQuad<3>(stiff, LINEAR, mpi_rank);

  // Load in the geometry file
  TMRModel *geo = TMR_LoadModelFromCompound(compound);
  if (geo){
    geo->incref();

    // Get the faces that have been created - if any
    // and write them all to different VTK files
    int num_faces;
    TMRFace **faces;
    geo->getFaces(&num_faces, &faces);

    // Allocate the new mesh
    TMRMesh *mesh = new TMRMesh(comm, geo);
    mesh->incref();
    mesh->mesh(htarget);
    mesh->writeToVTK("cylinder-mesh.vtk");

    TMRModel *model = mesh->createModelFromMesh();
    model->incref();

    TMRQuadForest *forest = new TMRQuadForest(comm);
    forest->incref();

    TMRTopology *topo = new TMRTopology(model);
    forest->setTopology(topo);
    // forest->createTrees(1);
    forest->createRandomTrees(10, 0, 5);
    forest->balance();
    forest->createNodes(3);

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
    TACSAssembler *tacs = 
      new TACSAssembler(comm, vars_per_node,
                        num_nodes, num_elements, num_dep_nodes);
    tacs->incref();

    // Set the element ptr
    int *ptr = new int[ num_elements+1 ];
    for ( int i = 0; i < num_elements+1; i++ ){
      ptr[i] = 9*i;
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

    // Initialize the TACSAssembler object
    tacs->initialize();

    // Get the nodes
    TMRQuadrantArray *nodes;
    forest->getNodes(&nodes);

    // Get the quadrants associated with the nodes
    int size;
    TMRQuadrant *array;
    nodes->getArray(&array, &size);

    // Get the points
    TMRPoint *Xp;
    forest->getPoints(&Xp);

    TACSBVec *X = tacs->createNodeVec();
    TacsScalar *Xn;
    X->getArray(&Xn);

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

    tacs->setNodes(X);

    // Create and write out an fh5 file
    unsigned int write_flag = TACSElement::OUTPUT_NODES;
    TACSToFH5 *f5 = new TACSToFH5(tacs, PLANE_STRESS, write_flag);
    f5->incref();
      
    // Write out the solution
    f5->writeToFile("output.f5");
    f5->decref();

    forest->decref();
    model->decref();
    mesh->decref();
    geo->decref();
  }

  TMRFinalize();
  MPI_Finalize();
  return (0);
}

#endif // TMR_HAS_OPENCASCADE
