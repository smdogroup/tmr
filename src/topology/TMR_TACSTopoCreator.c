#include "TMR_TACSTopoCreator.h"
#include "TMROctStiffness.h"
#include "FElibrary.h"
#include "Solid.h"

/*
  Set up a creator class for the given filter problem
*/
TMROctTACSTopoCreator::TMROctTACSTopoCreator( TMRBoundaryConditions *_bcs,
                                              TMRStiffnessProperties _properties,
                                              TMROctForest *_filter ):
TMROctTACSCreator(_bcs){
  // Reference the filter
  filter = _filter;
  filter->incref();
  
  MPI_Comm comm = filter->getMPIComm();
  MPI_Comm_rank(comm, &mpi_rank);

  // Set the material properties
  properties = _properties;

  // Balance the filter and create nodes. Create the dependent node
  // connectivity
  filter->balance();
  filter->createNodes(2);
  filter->createDepNodeConn();

  // Get the node range for the filter design variables
  const int *filter_range;
  filter->getOwnedNodeRange(&filter_range);

  // Get a sorted list of the external node numbers
  int *ext_node_nums;
  int num_ext_nodes = filter->getExtNodeNums(&ext_node_nums);

  // Set up the variable map for the design variable numbers
  int num_filter_local = filter_range[mpi_rank+1] - filter_range[mpi_rank];
  filter_map = new TACSVarMap(comm, num_filter_local);
  filter_map->incref();

  // Set up the external filter indices for this filter.  The indices
  // objects steals the array for the external nodes.
  filter_indices = new TACSBVecIndices(&ext_node_nums, num_ext_nodes);
  filter_indices->incref();
  filter_indices->setUpInverse();
}

/*
  Free the creator object
*/
TMROctTACSTopoCreator::~TMROctTACSTopoCreator(){
  filter->decref();
  filter_map->decref();
  filter_indices->decref();
}

  // Get the underlying information about the
void TMROctTACSTopoCreator::getForest( TMROctForest **_filter ){
  *_filter = filter;
}

void TMROctTACSTopoCreator::getMap( TACSVarMap **_map ){
  *_map = filter_map;
}

void TMROctTACSTopoCreator::getIndices( TACSBVecIndices **_indices ){
  *_indices = filter_indices;
}

TACSElement* TMROctTACSTopoCreator::createElement( int order, 
                                                   TMROctForest *_forest,
                                                   TMROctant octant ){
  /*
  // Find the enclosing octant
  TMROctant *oct = filter->findEnclosing(&octant);

  // Get the side-length of the element
  const int32_t h = 1 << (TMR_MAX_LEVEL - octant.level);

  // Find the side length of the octant in the filter that contains
  // the element octant
  const int32_t hoct = 1 << (TMR_MAX_LEVEL - oct->level);
    
  // Get the u/v/w values within the filter octant
  double pt[3];
  pt[0] = -1.0 + 2.0*(octant.x + 0.5*h - oct->x)/hoct;
  pt[1] = -1.0 + 2.0*(octant.y + 0.5*h - oct->y)/hoct;
  pt[2] = -1.0 + 2.0*(octant.z + 0.5*h - oct->z)/hoct;
        
  // Get the Lagrange shape functions
  double N[8];
  FElibrary::triLagrangeSF(N, pt, 2);
  
  // Keep track of the weights for each node
  int nweights = 0;
  TMRIndexWeight weights[32];

  // Get the dependent node information for this mesh
  const int *dep_ptr, *dep_conn;
  const double *dep_weights;
  int num_dep_nodes = filter->getDepNodeConn(&dep_ptr, &dep_conn,
                                             &dep_weights);

  // Get the node range for the filter design variables
  const int *filter_range;
  filter->getOwnedNodeRange(&filter_range);
  int num_filter_local = filter_range[mpi_rank+1] - filter_range[mpi_rank];

  // Get the octant array for the nodes
  TMROctantArray *filter_nodes;
  filter->getNodes(&filter_nodes);

  // Loop over the adjacent nodes within the filter
  for ( int kk = 0; kk < 2; kk++ ){
    for ( int jj = 0; jj < 2; jj++ ){
      for ( int ii = 0; ii < 2; ii++ ){
        // Set the weights
        double wval = N[ii + 2*jj + 4*kk];

        // Compute the location of the node
        TMROctant p;
        p.block = octant.block;
        p.x = oct->x + hoct*ii;
        p.y = oct->y + hoct*jj;
        p.z = oct->z + hoct*kk;

        // Search for the node p within the octant
        const int use_node_search = 1;
        TMROctant *t = filter_nodes->contains(&p, use_node_search);

        // Get the node number
        int node = t->tag;
        if (node >= 0){
          if (node >= filter_range[mpi_rank] && 
              node < filter_range[mpi_rank+1]){
            node = node - filter_range[mpi_rank];
          }
          else {
            node = num_filter_local + filter_indices->findIndex(node);
          }
          weights[nweights].index = node;
          weights[nweights].weight = wval;
          nweights++;
        }
        else {
          node = -node-1;
          for ( int jp = dep_ptr[node]; jp < dep_ptr[node+1]; jp++ ){
            int dep_node = dep_conn[jp];
            if (dep_node >= filter_range[mpi_rank] && 
                dep_node < filter_range[mpi_rank+1]){
              dep_node = dep_node - filter_range[mpi_rank];
            }
            else { 
              dep_node = num_filter_local + filter_indices->findIndex(dep_node);
            }
            weights[nweights].index = dep_node;
            weights[nweights].weight = wval*dep_weights[jp];
            nweights++;
          }
        }
      }
    }
  }

  // Sort and sum the array of weights
  nweights = TMRIndexWeight::uniqueSort(weights, nweights);
        
  // Allocate the stiffness object
  SolidStiffness *stiff = 
    new TMROctStiffness(weights, nweights, 
                        properties.rho, properties.E, properties.nu, 
                        properties.q);
  */

  // Create the solid stiffness object
  SolidStiffness *stiff = new SolidStiffness(properties.rho, properties.E,
                                             properties.nu);

  TACSElement *solid = NULL;
  if (order == 2){
    solid = new Solid<2>(stiff);
  }
  else if (order == 3){
    solid = new Solid<3>(stiff);
  }

  return solid;
}
