#ifndef TMR_GENERATE_STL_H
#define TMR_GENERATE_STL_H

#include "TMROctree.h"

extern void generate_stl_file( const char *filename,
                               double *x,
                               TMROctree *filter,
                               double cutoff );

#endif
