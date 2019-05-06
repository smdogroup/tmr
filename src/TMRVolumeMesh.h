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

#ifndef TMR_VOLUME_MESH_H
#define TMR_VOLUME_MESH_H

#include "TMRMesh.h"
#include "TMRFaceMesh.h"

/*
  TMRVolumeMesh class

  This is the class that contains the volume mesh. This mesh
  takes in the arguments needed to form a volume mesh
*/
class TMRVolumeMesh : public TMREntity {
 public:
  TMRVolumeMesh( MPI_Comm _comm, TMRVolume *volume );
  ~TMRVolumeMesh();

  // Create the volume mesh
  int mesh( TMRMeshOptions options );

  // Retrieve the mesh points
  void getMeshPoints( int *_npts, TMRPoint **X );

  // Retrieve the local connectivity from this volume mesh
  int getHexConnectivity( const int **hex );
  int getTetConnectivity( const int **tets );

  // Order the mesh points uniquely
  int setNodeNums( int *num );
  int getNodeNums( const int **_vars );

  // Write the volume mesh to a VTK file
  void writeToVTK( const char *filename );

 private:
  // Create a tetrahedral mesh (if possible)
  int tetMesh( TMRMeshOptions options );

  // Set the node locations based on the surface node locations
  int setNodeLocations( TMRMeshOptions options );

  // Get the mapping for the swept mesh
  int getSweptMapping( int source_plane, int dest_plane,
                       int num_fixed_pts,  int num_quad_pts,
                       double *A, double *b, double *c,
                       double regfactor );

  // The underlying volume
  MPI_Comm comm;
  TMRVolume *volume;

  // Set the number of face loops
  int num_swept_faces;
  TMREdge **swept_edges;
  int *swept_edges_orient;
  TMRFace **swept_faces;

  // Number of points swept through-thickness
  int num_swept_pts;

  // Keep the bottom/top surfaces (master/target) in the mesh for
  // future reference
  TMRFace *target, *source;
  int target_orient, source_orient;

  int num_points; // The number of points
  TMRPoint *X; // The physical node locations
  int *vars; // The global variable numbers

  // Hexahedral mesh information
  int num_hex;
  int *hex;

  // Tetrahedral mesh information
  int num_tet;
  int *tet;
};

#endif // TMR_VOLUME_MESH_H
