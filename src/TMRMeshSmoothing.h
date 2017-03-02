#ifndef TMR_MESH_SMOOTHING_H
#define TMR_MESH_SMOOTHING_H

#include "TMRTopology.h"

/*
  Smoothing methods for different meshes
*/
void laplacianSmoothing( int nsmooth, int num_fixed_pts,
                         int num_edges, const int *edge_list,
                         int num_pts, double *prm, TMRPoint *p,
                         TMRFace *face );
void springSmoothing( int nsmooth, double alpha, int num_fixed_pts,
                      int num_edges, const int *edge_list,
                      int num_pts, double *prm, TMRPoint *p,
                      TMRFace *face );
void springQuadSmoothing( int nsmooth, double alpha, int num_fixed_pts,
                          int num_quads, const int *quad_list,
                          int num_edges, const int *edge_list,
                          int num_pts, double *prm, TMRPoint *p,
                          TMRFace *face );
void quadSmoothing( int nsmooth, int num_fixed_pts,
                    int num_pts, const int *ptr, const int *pts_to_quads,
                    int num_quads, const int *quad_list,
                    double *prm, TMRPoint *p,
                    TMRFace *face );

#endif // TMR_MESH_SMOOTHING_H
