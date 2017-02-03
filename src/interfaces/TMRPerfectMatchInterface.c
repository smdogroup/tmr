#include "TMRPerfectMatchInterface.h"
#include <stdio.h>

// Include the perfect matching code
#define PERFECT_MATCHING_DOUBLE
#include "PerfectMatching.h"

int TMR_PerfectMatchGraph( int nnodes, int nedges, 
                           const int *edges, const double *weights,
                           int *match ){
  if (nnodes % 2 == 1){
  	printf("Perfect matching does not exist\n");
  }
  
  // Allocate the perfect matching code
  PerfectMatching *pm = new PerfectMatching(nnodes, nedges);

  // Set the edges in the graph and the weights
  for ( int i = 0; i < nedges; i++ ){
    // Add the dual edge to the perfect matching code
    pm->AddEdge(edges[2*i], edges[2*i+1], weights[i]);
  }

  // Set the options for the perfect matching algorithm
  struct PerfectMatching::Options options;
  options.verbose = false;
  pm->options = options;

  pm->Solve();

  // Allocate the match
  int nmatch = 0;
  for ( int i = 0; i < nnodes; i++ ){
    int j = pm->GetMatch(i);
    if (i < j){
      match[2*nmatch] = i;
      match[2*nmatch+1] = j;
      nmatch++;
    }
  }
  
  delete pm;

  return 0;
}
