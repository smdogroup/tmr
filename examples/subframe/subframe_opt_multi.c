#ifdef TMR_HAS_OPENCASCADE

#include "TMROpenCascade.h"
#include "TMRMesh.h"
#include "TACSMeshLoader.h"
#include "Solid.h"
#include "TMROctForest.h"
#include "TACS3DTraction.h"
#include "TMR_TACSTopoCreator.h"
#include "TMR_RefinementTools.h"
#include "TMRMultiLoadTopoProblem.h"

/* /\* */
/*   Add the 3D traction */
/* *\/ */
/* void addFaceTractions( int order, */
/*                        TMROctForest *forest, */
/*                        const char *attr,  */
/*                        TACSAssembler *tacs, */
/*                        TacsScalar Tr[3] ){ */
/*   // Create the tractions on each surface */
/*   TACSElement *trac[6]; */
/*   for ( int face = 0; face < 6; face++ ){ */
/*     if (order == 2){ */
/*       trac[face] = new TACS3DTraction<2>(face, Tr[0], Tr[1], Tr[2]); */
/*     } */
/*     else { */
/*       trac[face] = new TACS3DTraction<3>(face, Tr[0], Tr[1], Tr[2]); */
/*     } */
/*     trac[face]->incref(); */
/*   } */

/*   // Retrieve the array of octants from the array */
/*   TMROctantArray *octants; */
/*   forest->getOctants(&octants); */
/*   TMROctant *first; */
/*   octants->getArray(&first, NULL); */
  
/*   // Get the octants with the specified attribute */
/*   TMROctantArray *octs = forest->getOctsWithAttribute(attr); */

/*   // Get the octant  */
/*   int size; */
/*   TMROctant *array; */
/*   octs->getArray(&array, &size); */

/*   // Get the auxilary elements from TACS */
/*   TACSAuxElements *aux = tacs->getAuxElements(); */
/*   if (!aux){ */
/*     aux = new TACSAuxElements(); */
/*   } */

/*   for ( int i = 0; i < size; i++ ){ */
/*     // Get the face index */
/*     int face_index = array[i].tag; */

/*     // Get the local octant index in the array */
/*     int use_node_search = 0; */
/*     TMROctant *me = octants->contains(&array[i], use_node_search); */
/*     int element_num = me - first; */

/*     // Add the element to the auxiliary elements */
/*     aux->addElement(element_num, trac[face_index]); */
/*   } */

/*   // Set the auxiliary elements into TACS */
/*   tacs->setAuxElements(aux); */

/*   delete octs; */
/*   for ( int face = 0; face < 6; face++ ){ */
/*     trac[face]->decref(); */
/*   } */
/* } */

void addFaceTractions( int order,
                       TMROctForest *forest,
                       const char *attr, 
                       TACSAuxElements **aux,
                       TACSAssembler *tacs,
                       TacsScalar Tr[3],
		       int load_num,
		       int case_num){
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

  
  //-------------------------------------------
  if (case_num == 0){
    aux[0] = new TACSAuxElements();
    aux[0]->incref();
  }
  
  /* // Get the auxilary elements from TACS */
  
  /* TACSAuxElements *_aux = tacs->getAuxElements(); */
  

  /* if (!_aux){ */
  /*   printf("load_n: %d %d\n", load_num,case_num); */
  /*   _aux = new TACSAuxElements(); */
  /* } */

  for ( int i = 0; i < size; i++ ){
    // Get the face index
    int face_index = array[i].tag;

    // Get the local octant index in the array
    int use_node_search = 0;
    TMROctant *me = octants->contains(&array[i], use_node_search);
    int element_num = me - first;

    // Add the element to the auxiliary elements
    aux[0]->addElement(element_num, trac[face_index]);
  }

  /* // Set the auxiliary elements into TACS */
  /* tacs->setAuxElements(aux); */
  //*aux = _aux;
  delete octs;
  for ( int face = 0; face < 6; face++ ){
    trac[face]->decref();
  }
}

/*
  Based on the forest, create a filter and a hierarchy of
  TACSAssembler functions
*/
void createTopoProblem( int num_levels,
                        TMRBoundaryConditions *bcs,
                        TMRStiffnessProperties *props,
                        TMROctForest *forest,
                        TMROctForest *filter,
                        double target_mass,
                        const char *prefix,
                        TACSAssembler **_tacs,
                        TACSVarMap **_filter_map,
                        TMRMultiLoadTopoProblem **_prob,
                        int nloads){
  // Just use 2nd order meshes throughout
  const int order = 2;

  // Allocate the arrays to store the different levels
  TMROctForest **forest_levs = new TMROctForest*[ num_levels ];
  TMROctForest **filter_levs = new TMROctForest*[ num_levels ];
  TACSAssembler **tacs = new TACSAssembler*[ num_levels ];

  // Store the design variable mapping information
  TACSVarMap **filter_maps = new TACSVarMap*[ num_levels ];
  TACSBVecIndices **filter_ext_indices = new TACSBVecIndices*[ num_levels ];

  // Set the forest/filter for the most refined mesh  
  forest_levs[0] = forest;
  filter_levs[0] = filter; 
   
  TACSAuxElements **aux = new TACSAuxElements*[nloads];
  for (int i = 0; i < nloads; i++){
    aux[i] = NULL;
  }

  for ( int level = 0; level < num_levels; level++ ){
    if (level > 0){
      forest_levs[level] = forest_levs[level-1]->coarsen();
      filter_levs[level] = filter_levs[level-1]->coarsen();
    }

    // Set the current forest/filter
    forest = forest_levs[level];
    filter = filter_levs[level];

    forest->incref();
    filter->incref();
     
    // Balance and repartition the forest and filter and create nodes
    forest->balance();
    forest->repartition();    

    filter->balance();
    filter->repartition();

    TMROctTACSTopoCreator *creator = 
      new TMROctTACSTopoCreator(bcs, *props, filter);
    creator->incref();

    // Create tacs
    tacs[level] = creator->createTACS(order, forest);

    // Extract the filter indices/node maps from the creator class
    creator->getMap(&filter_maps[level]);
    creator->getIndices(&filter_ext_indices[level]);
    filter_maps[level]->incref();
    filter_ext_indices[level]->incref();

    // Set the filter map for the most refined filter
    if (_filter_map && level == 0){
      *_filter_map = filter_maps[level];
    }
    if (_tacs && level == 0){
      *_tacs = tacs[level];
    }

    // Delete the creator class
    creator->decref();
    int load_case_number[] = {12,12,12,12,12};
    char load_name[256];
    if (level == 0){
      // Create the aux elements on the face
      TacsScalar Tr[] = { 0.040673,-0.009505,0.015694,
                         -0.034705,-7.636785,-0.337983,
                         -0.040452,0.000000,-0.006853,
                         -0.435466,1.197862,0.040010,
                         -0.117987,-0.852858,0.038409,
                         -0.117987,-0.852858,0.038409,
                         0.396190,-0.876838,-0.137934,
                         0.396190,-0.876838,-0.137934,
                         0.002894,0.008102,-0.076684,
                         0.002894,0.008102,-0.076684,
                         -0.004341,-0.005498,0.094046,
                         -0.004341,-0.005498,0.094046,

                          -0.012379,0.002874,0.007295,
                         0.010610,2.416503,0.087093,
                         0.012821,0.001105,0.001547,
                         0.017463,-0.122240,-0.004642,
                         0.019311,0.197140,-0.009549,
                         0.019311,0.197140,-0.009549,
                         -0.066633,0.199899,0.032680,
                         -0.066633,0.199899,0.032680,
                         -0.000868,-0.003183,0.024018,
                         -0.000868,-0.003183,0.024018,
                         0.000868,0.002604,-0.009549,
                         0.000868,0.002604,-0.009549,
                               
                         -0.023873,-0.028073,0.003537,
                         2.586931,8.661345,0.372025,
                         0.019010,-0.002653,0.006410,
                         2.289842,-6.467217,-0.068525,
                         -0.098676,-1.374674,0.062389,
                         -0.098676,-1.374674,0.062389,
                         0.292421,-1.349846,-0.224302,
                         0.292421,-1.349846,-0.224302,
                         0.000289,-0.025465,0.037908,
                         0.000289,-0.025465,0.037908,
                         0.000000,0.024597,-0.023439,
                         0.000000,0.024597,-0.023439,

                          0.009063,0.011273,0.011273,
                         -0.714208,-2.947903,-0.119366,
                         -0.007295,0.001547,-0.003095,
                         -0.841090,2.104161,0.036915,
                         0.058993,0.449454,-0.014430,
                         0.058993,0.449454,-0.014430,
                         -0.202870,0.445846,0.058781,
                         -0.202870,0.445846,0.058781,
                         0.000289,0.007234,0.010417,
                         0.000289,0.007234,0.010417,
                         -0.000579,-0.007524,0.004341,
                         -0.000579,-0.007524,0.004341,

                         -0.036694,-0.185902,0.042662,
                         -0.271227,1.323418,0.405845,
                         -0.035810,-0.074935,0.007958,
                         -0.314773,0.579368,0.142576,
                         0.017613,0.197989,-0.037136,
                         0.017613,0.197989,-0.037136,
                         -0.089551,0.198201,0.177405,
                         -0.089551,0.198201,0.177405,
                         0.038776,-0.003762,-0.684945,
                         0.038776,-0.003762,-0.684945,
                         -0.050061,-0.000289,0.769152,
                         -0.050061,-0.000289,0.769152};
      
      for (int i = 0; i < nloads; i++){
        for (int j = 0; j < load_case_number[i]; j++){
          sprintf(load_name,"Load%d",j);
          addFaceTractions(order, forest, load_name, &aux[i],
                           tacs[0],&Tr[3*j+3*12*(i)],
    			   i,j);
        }
      }
      /* //TacsScalar Tr[] = {100.0, 0.0, 0.0}; */
      
      /* TacsScalar Tr[] = {-.294, .440, .015, */
      /*                         .198, .094, .042, */
      /*                    -.020, 0.0, -.003, */
      /*                    -.218, .599, .020, */
      /*                    -.059, -.426, .019, */
      /*                         .198, -.438, -.069, */
      /*                         .001, .004, -.038, */
      /*                         -.002, -.003, .047}; */
      /* for (int p = 0; p < num_load_case; p++){ */
      /*        for (int j = 0; j < load_case_number[p]; j++){ */
      /*          sprintf(load_name,"Load%d",j); */
      /*          addFaceTractions(order, forest, load_name, tacs[0], &Tr[3*j]); */
      /*        } */
      /* } */
        //addFaceTractions(order, forest, "Load", tacs[0], Tr);
    }
  }

  // Create the multigrid object for TACS
  TACSMg *mg;
  TMR_CreateTACSMg(num_levels, tacs, forest_levs, &mg);
  
  // Create the topology optimization problem
  TMRMultiLoadTopoProblem *prob = 
    new TMRMultiLoadTopoProblem(num_levels, tacs, filter_levs,
                                filter_maps, filter_ext_indices, mg,
                                target_mass, prefix, aux, nloads);
  *_prob = prob;

  // Free the filter maps
  for ( int level = 0; level < num_levels; level++ ){
    forest_levs[level]->decref();
    filter_levs[level]->decref();
    filter_maps[level]->decref();
    filter_ext_indices[level]->decref(); 
  }
}

/*
  Main function for the bracket optimization problem
*/
int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  // Get the communicator
  MPI_Comm comm = MPI_COMM_WORLD;

  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  // Compute the approximate volume of the part
  
  //double r = 7.5;
  //double a = 50.0;
  //double t = 25.0;
  //double volume = 3*a*a*t - 3*M_PI*r*r*t;

  double r = 100.;
  double a = 910.;
  double t = 350.;
  double volume = a*r*t;

  // Set the target mass fraction
  double target_mass_fraction = 0.2;

  char prefix[256];
  sprintf(prefix, "results_subframe/multi");

  for ( int k = 0; k < argc; k++ ){
    if (sscanf(argv[k], "prefix=%s", prefix) == 0){
      if (mpi_rank == 0){
        printf("Using prefix = %s\n", prefix);
      }
    }
  }

  // Compute the target mass
  double target_mass = target_mass_fraction*volume;

  if (mpi_rank == 0){
    printf("Approximate volume = %f\n", volume);
    printf("Target mass        = %f\n", target_mass);
  }

  // This is all tied to this STEP file
  const char *filename = "subframe1.stp";
  double htarget = 8.5;

  // Create the boundary conditions 
  int bc_nums[] = {0, 1, 2};
  TMRBoundaryConditions *bcs = new TMRBoundaryConditions();
  bcs->incref();
  bcs->addBoundaryCondition("Fixed", 3, bc_nums);

  // Material properties
  TMRStiffnessProperties properties;
  properties.rho = 1.0;
  properties.E = 70e9;
  properties.nu = 0.3;
  properties.q = 5.0;
  
  int nloads = 2;
    
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
    int lower_face_num = 31;
    int upper_face_num = 32;
    TMRFace *lower = faces[lower_face_num];
    TMRFace *upper = faces[upper_face_num];
    upper->setMaster(lower);

    // Reset the master orientations based on the volume object
    volume[0]->updateOrientation();

    // Set the attributes associated with the boundary conditions and
    // the loading condition
    int face_fixed[] = {30,29,16,15};
    for (int p = 0; p < 4; p++){
      faces[face_fixed[p]]->setAttribute("Fixed");
    }
    
    int face_load[] = {28,27,18,17,
		       26,24,25,23,22,21,20,19};
  
    for (int p = 0; p < 12; p++){
      char load_name[256];
      sprintf(load_name, "Load%d",p);
      faces[face_load[p]]->setAttribute(load_name);
    }
    // Mesh the bracket
    TMRMesh *mesh = new TMRMesh(comm, geo);
    mesh->incref();
    
    // Mesh the geometry
    TMRMeshOptions options;
    options.num_smoothing_steps = 50;
    options.frontal_quality_factor = 1.5;
    mesh->mesh(options, htarget);

    mesh->writeToVTK("subframe_mesh.vtk", TMRMesh::TMR_HEXES);

    TMRModel *model = mesh->createModelFromMesh();
    model->incref();

    // Create the topology object from the geo-mesh
    TMRTopology *topo = new TMRTopology(comm, model);
    topo->incref();

    // The forest used to represent the mesh and the filter
    // respectively
    TMROctForest *forest;

    // Create the filter
    forest = new TMROctForest(comm);
    forest->incref();

    // Set the geometry 
    forest->setTopology(topo);

    // create the trees for the mesh
    int min_level = 1; // The minimum level of refinement
    int max_level = 4; // The maximum level of refinement
    forest->createTrees(1);
    forest->balance();

    ParOptBVecWrap *old_design_vars = NULL;
    TACSVarMap *old_filter_map = NULL;
    TMROctForest *old_filter = NULL;

    int max_iterations = 1;
    for ( int iter = 0; iter < max_iterations; iter++ ){
      // Create the new filter and possibly interpolate from the old to
      // the new filter
      TACSVarMap *filter_map = NULL;
      TMROctForest *filter = forest->coarsen();
      filter->incref();

      // Set the number of levels of multigrid to use
      int num_levels = 2 + iter;
      if (num_levels > max_level+1){
        num_levels = max_level+1;
      }

      // Create the topology optimization object
      TACSAssembler *tacs;
      TMRMultiLoadTopoProblem *prob;
      
      createTopoProblem(num_levels, bcs, &properties,
                        forest, filter, target_mass,
                        prefix, &tacs, &filter_map, &prob,
                        nloads);
      tacs->incref();
      filter_map->incref();

      // Create the topology optimization object
      int max_num_bfgs = 10;
      ParOpt *opt = new ParOpt(prob, max_num_bfgs);

      // THe new design variable values
      ParOptBVecWrap *new_design_vars = 
        dynamic_cast<ParOptBVecWrap*>(prob->createDesignVec());

      // Set the design variables using the old design
      // variable values interpolated from the old filter
      if (old_filter){
        // Create the interpolation object
        TACSBVecInterp *interp = new TACSBVecInterp(old_filter_map,
                                                    filter_map, 1);
        interp->incref();

        // Create the interpolation between the two filter
        filter->createInterpolation(old_filter, interp);
        interp->initialize();

        // Decrease the reference counts to the data that is
        // no longer required
        old_filter_map->decref();
        old_filter->decref();

        // Perform the interpolation from the old set of values to the
        // new set of values
        TACSBVec *old_vec = old_design_vars->vec;
        TACSBVec *new_vec = new_design_vars->vec;
        interp->mult(old_vec, new_vec);
        prob->setInitDesignVars(new_design_vars);
        interp->decref();

        // Free the old design variable values and reset the pointer
        delete old_design_vars;
      }

      // Set the old design variable vector to what was once
      // the new design variable values
      old_design_vars = new_design_vars;
      
      // check the gradients
      opt->checkGradients(1e-6);
      exit(0);
      //printf("Here\n");
      
      //opt->setGradientCheckFrequency(100,1e-6);
      // Set the optimization parameters
      int max_opt_iters = 250;
      opt->setMaxMajorIterations(max_opt_iters);
      prob->setIterationCounter(max_opt_iters*iter);
      opt->setOutputFrequency(1);
        
      // Set the Hessian reset frequency
      opt->setBFGSUpdateType(LBFGS::DAMPED_UPDATE);
      opt->setHessianResetFreq(max_num_bfgs);
      
      // Set the barrier parameter information
      opt->setBarrierPower(1.0);
      opt->setBarrierFraction(0.25);
      
      // Set the log/output file
      char outfile[256];
      sprintf(outfile, "%s//paropt_output%d.out", prefix, iter);
      opt->setOutputFile(outfile);
      
      // Set the history/restart file
      char new_restart_file[256];
      sprintf(new_restart_file, "%s//paropt_restart%d.bin", prefix, iter);
      opt->optimize(new_restart_file);
      //exit(0);
      // Get the final values of the design variables
      ParOptVec *old;
      opt->getOptimizedPoint(&old, NULL, NULL, NULL, NULL);
      old_design_vars->copyValues(old);

      // Set the old filter/variable map for the next time through the
      // loop so that we can interpolate design variable values
      old_filter_map = filter_map;
      old_filter = filter;
      
      // Free the optimization problem data
      delete opt;
      delete prob;

      // Create the visualization for the object
      unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                                 TACSElement::OUTPUT_DISPLACEMENTS |
                                 TACSElement::OUTPUT_STRAINS |
                                 TACSElement::OUTPUT_STRESSES |
                                 TACSElement::OUTPUT_EXTRAS); 
      TACSToFH5 *f5 = new TACSToFH5(tacs, TACS_SOLID, write_flag);
      f5->incref();
      sprintf(outfile, "%s//tacs_output%d.f5", prefix, iter);
      f5->writeToFile(outfile);
      f5->decref();

      // Create the refinement array
      int num_elements = tacs->getNumElements();
      int *refine = new int[ num_elements ];
      memset(refine, 0, num_elements*sizeof(int));
      //exit(0);
      // Refine based solely on the value of the density variable
      TACSElement **elements = tacs->getElements();
      for ( int i = 0; i < num_elements; i++ ){
        TACSConstitutive *c = elements[i]->getConstitutive();

        // In the case of the oct stiffness objects, they always return
        // the density as the first design variable, regardless of the
        // parameter point used in the element
        double pt[3] = {0.0, 0.0, 0.0};
        TacsScalar rho = c->getDVOutputValue(0, pt);

        // Refine things differently depending on whether
        // the density is above or below a threshold
        if (rho > 0.5){
          refine[i] = 1;
        }
        else if (rho < 0.05){
          refine[i] = -1;
        }
      }
      
      // Refine the forest
      forest->refine(refine, min_level, max_level);
      delete [] refine;

      // Free tacs and the filter map
      filter_map->decref();
      tacs->decref();
    }
 
    // Free the analysis/mesh data
    forest->decref();
    topo->decref();
    model->decref();
    mesh->decref();
    geo->decref();
  }

  bcs->decref();

  TMRFinalize();
  MPI_Finalize();

  return 0;
}

#endif // TMR_HAS_OPENCASCADE
