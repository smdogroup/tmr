/*
  This file is part of the package TMR for adaptive mesh refinement.

  Copyright (C) 2015 Georgia Tech Research Corporation.
  Additional copyright (C) 2015 Graeme Kennedy.
  All rights reserved.

  TMR is licensed under the Apache License, Version 2.0 (the "License");
  you may not use this software except in compliance with the License.
  You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0
  
  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
*/

#ifndef TMR_MESH_SMOOTHING_H
#define TMR_MESH_SMOOTHING_H

#include "TMRTopology.h"

/*
  Smoothing methods for different meshes
*/
void TMR_LaplacianSmoothing( int nsmooth, int num_fixed_pts,
                             int num_edges, const int *edge_list,
                             int num_pts, double *prm, TMRPoint *p,
                             TMRFace *face );
void TMR_SpringSmoothing( int nsmooth, double alpha, int num_fixed_pts,
                          int num_edges, const int *edge_list,
                          int num_pts, double *prm, TMRPoint *p,
                          TMRFace *face );
void TMR_SpringQuadSmoothing( int nsmooth, double alpha, int num_fixed_pts,
                              int num_quads, const int *quad_list,
                              int num_edges, const int *edge_list,
                              int num_pts, double *prm, TMRPoint *p,
                              TMRFace *face );
void TMR_QuadSmoothing( int nsmooth, int num_fixed_pts,
                        int num_pts, const int *ptr, const int *pts_to_quads,
                        int num_quads, const int *quad_list,
                        double *prm, TMRPoint *p,
                        TMRFace *face );

#endif // TMR_MESH_SMOOTHING_H
