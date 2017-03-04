#ifndef TMR_PERFECT_MATCH_INTERFACE
#define TMR_PERFECT_MATCH_INTERFACE

int TMR_PerfectMatchGraph( int nnodes, int nedges, 
                           const int *edges, const double *weights,
                           int *match );

#endif // TMR_PERFECT_MATCH_INTERFACE
