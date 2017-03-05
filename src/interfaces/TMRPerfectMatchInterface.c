#include "TMRPerfectMatchInterface.h"
#include <string.h>
#include <stdio.h>

// Include the perfect matching code
#define PERFECT_MATCHING_DOUBLE
#include "PerfectMatching.h"

int TMR_PerfectMatchGraph( int nnodes, int nedges, 
                           const int *edges, const double *weights,
                           int *match ){
  if (nnodes % 2 == 1){
    fprintf(stderr,
            "TMR_PerfectMatchGraph error: Perfect matching does not exist\n");
  }
  
  // Check that all nodes are referenced
  int *count = new int[ nnodes ];
  memset(count, 0, nnodes*sizeof(int));
  for ( int i = 0; i < nedges; i++ ){
    if (edges[2*i] >= 0 && edges[2*i] < nnodes){
      count[edges[2*i]]++;
    }
    else {
      fprintf(stderr, 
              "TMR_PerfectMatchGraph error: Edge %d, node %d out of range\n",
              i, edges[2*i]);
    }
    if (edges[2*i+1] >= 0 && edges[2*i+1] < nnodes){
      count[edges[2*i+1]]++;
    }
    else {
      fprintf(stderr, 
              "TMR_PerfectMatchGraph error: Edge %d, node %d out of range\n",
              i, edges[2*i+1]);
    }
  }

  for ( int i = 0; i < nnodes; i++ ){
    if (count[i] == 0){
      fprintf(stderr,
              "TMR_PerfectMatchGraph error: Node %d not referenced\n", i);
    }
  }

  delete [] count;

  // Allocate the perfect matching code
  PerfectMatching *pm = new PerfectMatching(nnodes, nedges);

  // Set the edges in the graph and the weights
  int max_int = 1 << 10;
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
