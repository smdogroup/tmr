#ifdef TMR_HAS_OPENCASCADE

#include "TMROpenCascade.h"
#include "TMRMesh.h"
#include "TACSMeshLoader.h"
#include "Solid.h"
#include "TMROctForest.h"
#include "TACS3DTraction.h"
#include "TMR_TACSTopoCreator.h"
#include "TMR_RefinementTools.h"

/*
  Add the 3D traction
*/
void addFaceTractions( int order,
                       TMROctForest *forest,
                       const char *attr, 
                       TACSAssembler *tacs,
                       TacsScalar Tr[3] ){
  // Create the tractions on each surface
  TACSElement *trac[6];
  for ( int face = 0; face < 6; face++ ){
    if (order == 2){
      trac[face] = new TACS3DTraction<2>(face, Tr[0], Tr[1], Tr[2]);
    }
    else {
      trac[face] = new TACS3DTraction<3>(face, Tr[0], Tr[1], Tr[2]);
    }
    trac[face]->incref();
  }

  // Retrieve the array of octants from the array
  TMROctantArray *octants;
  forest->getOctants(&octants);
  TMROctant *first;
  octants->getArray(&first, NULL);
  
  // Get the octants with the specified attribute
  TMROctantArray *octs = forest->getOctsWithAttribute(attr);

  // Get the octant 
  int size;
  TMROctant *array;
  octs->getArray(&array, &size);

  // Get the auxilary elements from TACS
  TACSAuxElements *aux = tacs->getAuxElements();
  if (!aux){
    aux = new TACSAuxElements();
  }

  for ( int i = 0; i < size; i++ ){
    // Get the face index
    int face_index = array[i].tag;

    // Get the local octant index in the array
    int use_node_search = 0;
    TMROctant *me = octants->contains(&array[i], use_node_search);
    int element_num = me - first;

    // Add the element to the auxiliary elements
    aux->addElement(element_num, trac[face_index]);
  }

  // Set the auxiliary elements into TACS
  tacs->setAuxElements(aux);

  delete octs;
  for ( int face = 0; face < 6; face++ ){
    trac[face]->decref();
  }
}

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  // Get the communicator
  MPI_Comm comm = MPI_COMM_WORLD;

  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  // This is all tied to this STEP file
  const char *filename = "bracket_solid.stp";
  double htarget = 5.0;

  // Create the boundary conditions 
  int bc_nums[] = {0, 1, 2};
  TMRBoundaryConditions *bcs = new TMRBoundaryConditions();
  bcs->addBoundaryCondition("Fixed", 3, bc_nums);

  // Material properties
  TMRStiffnessProperties properties;

  // Load in the geometry file
  TMRModel *geo = TMR_LoadModelFromSTEPFile(filename);
  if (geo){
    geo->incref();

    // Get the volume
    TMRVolume **volume;
    geo->getVolumes(NULL, &volume);

    // Get the faces from the volume
    TMRFace **faces;
    volume[0]->getFaces(NULL, &faces, NULL);

    // Set the upper/lower face numbers. These are based on the
    // ordering in the STEP file. This will hopefully be preserved
    // independent of how the STEP file is loaded.
    int lower_face_num = 1;
    int upper_face_num = 4;
    TMRFace *lower = faces[lower_face_num];
    TMRFace *upper = faces[upper_face_num];
    upper->setMaster(lower);

    // Reset the master orientations based on the volume object
    volume[0]->updateOrientation();

    // Set the attributes associated with the boundary conditions and
    // the loading condition
    int fixed_face_one = 15;
    int fixed_face_two = 16;
    faces[fixed_face_one]->setAttribute("Fixed");
    faces[fixed_face_two]->setAttribute("Fixed");

    int load_applied_on_face = 14;
    faces[load_applied_on_face]->setAttribute("Load");

    // Mesh the bracket
    TMRMesh *mesh = new TMRMesh(comm, geo);
    mesh->incref();

    // Mesh the geometry
    TMRMeshOptions options;
    options.num_smoothing_steps = 50;
    mesh->mesh(options, htarget);

    TMRModel *model = mesh->createModelFromMesh();
    model->incref();

    // Create the topology object from the geo-mesh
    TMRTopology *topo = new TMRTopology(comm, model);
    topo->incref();

    // Set the topology
    static const int MAX_NUM_LEVELS = 5;
    TMROctForest *forest[MAX_NUM_LEVELS];
    TMROctForest *filter[MAX_NUM_LEVELS];
    TMROctTACSTopoCreator *creator[MAX_NUM_LEVELS];
    TACSAssembler *tacs[MAX_NUM_LEVELS];

    int order = 2;

    // Create the filter
    filter[0] = new TMROctForest(comm);
    filter[0]->incref();

    filter[0]->setTopology(topo);
    filter[0]->createTrees(2);
    filter[0]->balance();
    filter[0]->repartition();

    // Duplicate and refine the filter so that the quadrants are
    // aligned locally on this processor
    forest[0] = new TMROctForest(comm);
    forest[0]->incref();

    forest[0]->setTopology(topo);
    forest[0]->createRandomTrees(10, 0, 5);
    forest[0]->balance();
    forest[0]->repartition();

    // Create the levels of refinement
    creator[0] = new TMROctTACSTopoCreator(bcs, properties, filter[0]);
    tacs[0] = creator[0]->createTACS(order, forest[0]);   

    int num_levels = 4;
    // Create all of the coarser levels
    for ( int i = 1; i < num_levels; i++ ){
      forest[i] = forest[i-1]->coarsen();
      forest[i]->incref();
      forest[i]->balance();
      forest[i]->repartition();
      tacs[i] = creator[0]->createTACS(order, forest[i]);
      tacs[i]->incref();
    }

    // Create the multigrid object for TACS
    TACSMg *mg;
    TMR_CreateTACSMg(num_levels, tacs, forest, &mg);
    mg->incref();

    // Create the tractions on the face
    TacsScalar Tr[] = {100.0, 0.0, 0.0};
    addFaceTractions(order, forest[0], "Load", tacs[0], Tr);

    /*
      TMRTopoProblem *prob = new TMRTopoProblem( int _nlevels, 
                    TACSAssembler *_tacs[],
                    TMROctForest *_filter[], 
                    TACSVarMap *_filter_maps[],
                    TACSBVecIndices *_filter_indices[],
                    TACSMg *_mg,
                    double _target_mass,
                    const char *_prefix );

    // Create the topology optimization object
    ParOpt *opt = new ParOpt(prob, max_num_bfgs);

    // Set the optimization parameters
    opt->setMaxMajorIterations(2000);
    opt->setOutputFrequency(1);
  
    // Set the Hessian reset frequency
    int max_num_bfgs = 20;
    opt->setBFGSUpdateType(LBFGS::DAMPED_UPDATE);
    opt->setHessianResetFreq(max_num_bfgs);

    // Set the barrier parameter information
    opt->setBarrierPower(1.0);
    opt->setBarrierFraction(0.25);

    // Set the log/output file
    char outfile[256];
    sprintf(outfile, "%s//paropt_output%d.out", prefix, k);
    opt->setOutputFile(outfile);
  
    // Set the history/restart file
    char new_restart_file[256];
    sprintf(new_restart_file, "%s//paropt_restart%d.bin", prefix, k);
    opt->optimize(new_restart_file);

    // Free the optimization problem data
    delete opt;
    delete prob;

    // Free the analysis/mesh data
    forest->decref();
    filter->decref();
    */

    topo->decref();
    model->decref();
    mesh->decref();
    geo->decref();
  }

  TMRFinalize();
  MPI_Finalize();

  return 0;
}

#endif // TMR_HAS_OPENCASCADE
