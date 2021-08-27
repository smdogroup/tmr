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

#ifndef TMR_EDGE_MESH_H
#define TMR_EDGE_MESH_H

#include "TMRBase.h"
#include "TMRTopology.h"
#include "TMRMesh.h"

/*
  The mesh for a geometric curve
*/
class TMREdgeMesh : public TMREntity {
 public:
  TMREdgeMesh( MPI_Comm _comm, TMREdge *edge,
               TMRPoint *_X=NULL, int _npts=0 );
  ~TMREdgeMesh();

  // Is this edge mesh degenerate
  int isDegenerate(){ return edge->isDegenerate(); }

  // Retrieve the underlying curve
  void getEdge( TMREdge **_edge );

  // Mesh the geometric object
  void mesh( TMRMeshOptions options,
             TMRElementFeatureSize *fs );

  // Order the mesh points uniquely
  int setNodeNums( int *num );
  int getNodeNums( const int **_vars );

  // Retrieve the mesh points
  void getMeshPoints( int *_npts, const double **_pts, TMRPoint **X );

  // Write the edge points to the file
  void writeToVTK( const char *filename );

  // Get the relative orientation of the edge
  static int getEdgeCopyOrient( TMREdge *e );

 private:
  MPI_Comm comm;
  TMREdge *edge;

  // Keep track if this mesh is prescribed
  int prescribed_mesh;

  // The parametric locations of the points that are obtained from
  // meshing the curve
  int npts; // number of points along the curve
  double *pts; // Parametric node locations
  TMRPoint *X; // Physical node locations
  int *vars; // Global node variable numbers
};

#endif // TMR_EDGE_MESH_H
