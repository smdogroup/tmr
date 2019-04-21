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

#include <math.h>
#include <stdio.h>
#include "TMRMesh.h"
#include "TMRNativeTopology.h"
#include "TMRVolumeMesh.h"

#ifdef TMR_USE_NETGEN
// The namespace is required because of the way nglib is compiled by
// default. This is a funny way to do it, but who am I to complain.
namespace nglib {
#include "nglib.h"
}
#endif // TMR_USE_NETGEN

/*
  Try to create a volume mesh based on the structured/unstructured
  surface meshes
*/
TMRVolumeMesh::TMRVolumeMesh( MPI_Comm _comm,
                              TMRVolume *_volume ){
  comm = _comm;
  volume = _volume;
  volume->incref();

  // Set the number of loops on the source face
  num_swept_faces = 0;
  swept_edges = NULL;
  swept_edges_orient = NULL;
  swept_faces = NULL;
  num_swept_pts = -1;

  // Set the points in the mesh
  num_points = 0;
  num_hex = 0;
  num_tet = 0;

  // Set all the pointer to null
  X = NULL;
  hex = NULL;
  tet = NULL;
  vars = NULL;

  // Set the additional connectivity information that is required for
  // the volume mesh.
  source = target = NULL;
  source_orient = target_orient = 0;
}

TMRVolumeMesh::~TMRVolumeMesh(){
  volume->decref();
  if (X){ delete [] X; }
  if (hex){ delete [] hex; }
  if (tet){ delete [] tet; }
  if (vars){ delete [] vars; }
  if (swept_faces){
    for ( int k = 0; k < num_swept_faces; k++ ){
      swept_faces[k]->decref();
      swept_edges[k]->decref();
    }
    delete [] swept_edges;
    delete [] swept_edges_orient;
    delete [] swept_faces;
  }
}

/*
  Retrieve the mesh points
*/
void TMRVolumeMesh::getMeshPoints( int *_npts, TMRPoint **_X ){
  if (_npts){ *_npts = num_points; }
  if (_X){ *_X = X; }
}

/*
  Retrieve the local connectivity for this volume mesh
*/
int TMRVolumeMesh::getHexConnectivity( const int **_hex ){
  if (_hex){ *_hex = hex; }
  return num_hex;
}

/*
  Retrieve the local connectivity for this volume mesh
*/
int TMRVolumeMesh::getTetConnectivity( const int **_tet ){
  if (_tet){ *_tet = tet; }
  return num_tet;
}

/*
  Write the volume mesh to a VTK file
*/
void TMRVolumeMesh::writeToVTK( const char *filename ){
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  if (mpi_rank == 0 && hex){
    // Write out the connectivity of the mesh to a temp vtk file
    FILE *fp = fopen("volume-mesh.vtk", "w");
    if (fp){
      fprintf(fp, "# vtk DataFile Version 3.0\n");
      fprintf(fp, "vtk output\nASCII\n");
      fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

      // Write out the points
      fprintf(fp, "POINTS %d float\n", num_points);
      for ( int k = 0; k < num_points; k++ ){
        fprintf(fp, "%e %e %e\n", X[k].x, X[k].y, X[k].z);
      }

      // Write out the cell values
      fprintf(fp, "\nCELLS %d %d\n", num_hex, 9*num_hex);
      for ( int k = 0; k < num_hex; k++ ){
        fprintf(fp, "8 %d %d %d %d %d %d %d %d\n",
                hex[8*k], hex[8*k+1], hex[8*k+2], hex[8*k+3],
                hex[8*k+4], hex[8*k+5], hex[8*k+6], hex[8*k+7]);
      }

      // All hex
      fprintf(fp, "\nCELL_TYPES %d\n", num_hex);
      for ( int k = 0; k < num_hex; k++ ){
        fprintf(fp, "%d\n", 12);
      }
      fclose(fp);
    }
  }
}

/*
  Mesh the volume
*/
int TMRVolumeMesh::mesh( TMRMeshOptions options ){
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  if (options.mesh_type_default == TMR_TRIANGLE){
    return tetMesh(options);
  }

  // Keep track of whether the mesh has failed at any time. Try and
  // print out a helpful message.
  int mesh_fail = 0;

  // Get the faces and orientations associated with the volume
  int num_faces;
  TMRFace **faces;
  volume->getFaces(&num_faces, &faces);

  // Set integers for each face to determine whether it is a source or
  // a target face or it is connected to one of them or whether it is
  // structured.
  int *f = new int[ num_faces ];
  memset(f, 0, num_faces*sizeof(int));

  // Set the pointers for the target/source faces and their directions
  // relative to the interior of the hexahedral volume
  target = NULL;
  target_orient = 1;
  source = NULL;
  source_orient = 1;

  // Loop over all the faces to find the sources/targets
  for ( int i = 0; i < num_faces; i++ ){
    TMRFace *src;
    faces[i]->getSource(NULL, &src);
    if (src){
      // To be copied directly, the target orientation should point in
      // to the volume, otherwise its orientation must be flipped.
      target = faces[i];
      target_orient = -target->getOrientation();

      // Find the source face
      source = src;
      f[i] = 1;
      for ( int j = 0; j < num_faces; j++ ){
        if (faces[j] == source){
          f[j] = 1;

          // The natural orientation for the source face should point
          // out of the volume, otherwise its orientation must be
          // flipped.
          source_orient = source->getOrientation();
        }
      }

      // Break here. The meshing algorithm only works with one
      // source/target face pair per volume.
      break;
    }
  }

  // Check if the remaining surface meshes are structured
  for ( int i = 0; i < num_faces; i++ ){
    if (!f[i]){
      TMRFaceMesh *mesh;
      faces[i]->getMesh(&mesh);
      if (!mesh){
        fprintf(stderr,
                "TMRVolumeMesh error: No mesh associated with face %d\n",
                faces[i]->getEntityId());
        mesh_fail = 1;
      }
      else if (mesh->getMeshType() != TMR_STRUCTURED){
        fprintf(stderr,
                "TMRVolumeMesh error: \
Through-thickness meshes must be structured\n");
        fprintf(stderr,
                "Try setting source-target relations on edges and surfaces\n");
        mesh_fail = 1;
      }
    }
  }

  // Check that the source mesh is not a triangular mesh
  if (source){
    TMRFaceMesh *source_mesh;
    source->getMesh(&source_mesh);
    if (source_mesh->getMeshType() == TMR_TRIANGLE){
      fprintf(stderr,
              "TMRVolumeMesh error: Cannot extrude a triangluar mesh\n");
      mesh_fail = 1;
    }
  }

  // Free the f pointer
  delete [] f;

  if (mesh_fail){
    return mesh_fail;
  }

  // Each face that is not either the target or the source, must have
  // only four edges and must be structured. Two of the parallel edges
  // associated with these faces should touch the target and source
  // face, while remaining edges must be parallel and have the same
  // number of nodes. Furthermore, the one parallel edge will be
  // shared by the next face. We loop over the source edge loops and
  // find the connecting faces.

  // Get the number of edge loops
  int num_face_loops = source->getNumEdgeLoops();

  // Count up the number of swept faces
  num_swept_faces = 0;
  for ( int k = 0; k < num_face_loops; k++ ){
    TMREdgeLoop *source_loop;
    source->getEdgeLoop(k, &source_loop);

    // Get the number of edges for this loop
    int nedges;
    source_loop->getEdgeLoop(&nedges, NULL, NULL);
    num_swept_faces += nedges;
  }

  // Allocate space for the swept faces and edges
  swept_faces = new TMRFace*[ num_swept_faces ];
  swept_edges = new TMREdge*[ num_swept_faces ];
  swept_edges_orient = new int[ num_swept_faces ];

  // Set the swept faces/edges data
  for ( int count = 0, k = 0; k < num_face_loops; k++ ){
    // The target and source edge loops
    TMREdgeLoop *source_loop;
    source->getEdgeLoop(k, &source_loop);

    // Get the edges associated with the source loop
    int nedges;
    TMREdge **edges;
    source_loop->getEdgeLoop(&nedges, &edges, NULL);

    for ( int j = 0; j < nedges; j++, count++ ){
      // Search for the other face object that shares the edge object
      // edges[j] with the source face object
      TMRFace *face = NULL;
      for ( int i = 0; i < num_faces; i++ ){
        if (!(faces[i] == target || faces[i] == source)){
          // Get the edge loop associated with the face
          TMREdgeLoop *loop;
          faces[i]->getEdgeLoop(0, &loop);

          // Get the edge loops associated with face i
          int n;
          TMREdge **e;
          const int *edge_orient;
          loop->getEdgeLoop(&n, &e, &edge_orient);

          // Does this edge loop contain edges[j]
          for ( int ii = 0; ii < n; ii++ ){
            if (e[ii] == edges[j]){
              face = faces[i];
              swept_faces[count] = face;

              if (face->getOrientation() > 0){
                int next = ii+1;
                if (next >= n){
                  next = 0;
                }
                swept_edges[count] = e[next];
                swept_edges_orient[count] = edge_orient[next];
              }
              else {
                int next = ii-1;
                if (next < 0){
                  next = n-1;
                }
                swept_edges[count] = e[next];
                swept_edges_orient[count] = -edge_orient[next];
              }

              swept_faces[count]->incref();
              swept_edges[count]->incref();
              break;
            }
          }

          if (face){
            break;
          }
        }
      }

      // Get the number of mesh points along the swept direction
      TMREdgeMesh *mesh = NULL;
      swept_edges[count]->getMesh(&mesh);
      int npts = -1;
      mesh->getMeshPoints(&npts, NULL, NULL);

      // If this is the first time finding the number of points along
      // the depth of this volume, record the number of points.
      // Otherwise, verify that the number of points through the depth
      // of the volume is consistent.
      if (num_swept_pts < 0){
        num_swept_pts = npts;
      }
      else if (num_swept_pts != npts){
        fprintf(stderr,
                "TMRVolumeMesh error: \
Inconsistent number of edge points through-thickness %d != %d\n",
                num_swept_pts, npts);
        mesh_fail = 1;
      }
    }
  }

  // Get the information for the target surface
  TMRFaceMesh *mesh;
  target->getMesh(&mesh);

  // Number of points in the quadrilateral mesh on the surface
  int num_quad_pts = 0;
  TMRPoint *Xtarget;
  mesh->getMeshPoints(&num_quad_pts, NULL, &Xtarget);

  // Get information for the source surface
  source->getMesh(&mesh);

  // Points on the target surface
  TMRPoint *Xsource;
  mesh->getMeshPoints(NULL, NULL, &Xsource);

  // Get the local connectivity on the source surface
  const int *quads;
  int num_quads = mesh->getQuadConnectivity(&quads);

  // The number of hexahedral elements in the mesh
  num_hex = (num_swept_pts-1)*num_quads;
  num_points = num_swept_pts*num_quad_pts;
  hex = new int[ 8*num_hex ];
  X = new TMRPoint[ num_points ];

  // Flip the ordering in the hexahedral elements
  const int flip[] = {0, 3, 2, 1};

  int *h = hex;
  for ( int j = 0; j < num_swept_pts-1; j++ ){
    for ( int i = 0; i < num_quads; i++ ){
      // Set the quadrilateral points in the base layer. Note that the
      // orientation is fliped for the source when the face normal and
      // mesh are pointed outwards so that the hexahedral elements are
      // oriented correctly.
      if (source_orient > 0){
        for ( int k = 0; k < 4; k++ ){
          h[k] = j*num_quad_pts + quads[4*i+flip[k]];
        }
        for ( int k = 0; k < 4; k++ ){
          h[4+k] = (j+1)*num_quad_pts + quads[4*i+flip[k]];
        }
      }
      else {
        for ( int k = 0; k < 4; k++ ){
          h[k] = j*num_quad_pts + quads[4*i+k];
        }
        for ( int k = 0; k < 4; k++ ){
          h[4+k] = (j+1)*num_quad_pts + quads[4*i+k];
        }
      }
      h += 8;
    }
  }

  // Set the new coordinates within the hexahedral mesh
  TMRFaceMesh *target_mesh;
  target->getMesh(&target_mesh);

  // Get the source to target indices
  const int *source_to_target;
  target_mesh->getSourceToTargetMapping(&source_to_target);

  TMRPoint *x = X;
  for ( int j = 0; j < num_swept_pts; j++ ){
    double u = 1.0*j/(num_swept_pts-1);
    for ( int i = 0; i < num_quad_pts; i++ ){
      x[0].x = (1.0 - u)*Xsource[i].x + u*Xtarget[source_to_target[i]].x;
      x[0].y = (1.0 - u)*Xsource[i].y + u*Xtarget[source_to_target[i]].y;
      x[0].z = (1.0 - u)*Xsource[i].z + u*Xtarget[source_to_target[i]].z;
      x += 1;
    }
  }

  return mesh_fail;
}

/*
  Create a tetrahedral mesh
*/
int TMRVolumeMesh::tetMesh( TMRMeshOptions options ){
#ifdef TMR_USE_NETGEN
  nglib::Ng_Init();
  nglib::Ng_Mesh *m = nglib::Ng_NewMesh();

  // First count up the total number of unique points that bound this
  // volume that is about to be meshed
  int num_boundary_pts = 0;

  // Get the faces associated with the volume
  int num_faces;
  TMRFace **faces;
  const int *face_dir;
  volume->getFaces(&num_faces, &faces, &face_dir);

  // Count up the number of boundary nodes for this volume
  for ( int i = 0; i < ; i++ ){
    double pt[3] = {0.0, 0.0, 0.0};
    nglib::Ng_AddPoint(m, pt);
  }

  // Add the surface elements
  for ( int i = 0; i < ntris; i++ ){
    int tri[3];
    tri[0] = tris[3*i]+1;
    tri[1] = tris[3*i+1]+1;
    tri[2] = tris[3*i+2]+1;
    nglib::Ng_AddSurfaceElement(m, NG_TRIG, tri);
  }

  // Set the mesh parameters
  nglib::Ng_Meshing_Parameters mp;
  mp.maxh = htarget;
  mp.fineness = 1;
  mp.second_order = 0;

  // Generate the volume mesh
  nglib::Ng_GenerateVolumeMesh(m, &mp);

  // Get the total number of points and tets in the volume
  num_points = nglib::Ng_GetNP(m);
  num_tet = nglib::Ng_GetNE(m);

  // Allocate space to store everything
  X = new TMRPoint[ num_points ];
  tet = new int[ 4*num_tet ];

  // Retrieve the points
  for ( int i = 0; i < num_points; i++ ){
    double x[3];
    nglib::Ng_GetPoint(m, i+1, x);
    X[i].x = x[0];
    X[i].y = x[1];
    X[i].z = x[2];
  }

  // Retrieve the tets
  for ( int i = 0; i < num_tet; i++ ){
    nglib::Ng_GetVolumeElement(m, i+1, &tet[4*i]);
    for ( int k = 0; k < 4; k++ ){
      tet[4*i+k] -= 1;
    }
  }

  // Free the memory and exit from netgen
  nglib::Ng_DeleteMesh(m);
  nglib::Ng_Exit();
#endif // TMR_USE_NETGEN
  return 0;
}

/*
  Order the mesh points uniquely

  This code first copies the mesh locations from the target and source
  surfaces and then the surrounding surfaces that are structured. The
  code takes into account the orientations of the surrounding surfaces
  by computing the absolute orientations of the surface meshes. The
  orientation of the source/target surfaces is already accounted for
  within the connectivity.
*/
int TMRVolumeMesh::setNodeNums( int *num ){
  if (!vars){
    // Initially, set all of the nodes to zero
    vars = new int[ num_points ];
    for ( int k = 0; k < num_points; k++ ){
      vars[k] = -1;
    }

    // Get the source and target surface mesh and target surface meshes
    TMRFaceMesh *target_mesh;
    target->getMesh(&target_mesh);

    TMRFaceMesh *source_mesh;
    source->getMesh(&source_mesh);

    // Number of points in the quadrilateral mesh on the surface
    int num_quad_pts = 0;
    target_mesh->getMeshPoints(&num_quad_pts, NULL, NULL);

    // Get the nodes on the target face
    const int *target_face_vars;
    target_mesh->getNodeNums(&target_face_vars);

    // Use the source to target index mapping so that the variables will
    // be consistent with the target
    const int *source_to_target;
    target_mesh->getSourceToTargetMapping(&source_to_target);
    for ( int i = 0; i < num_quad_pts; i++ ){
      int target_var = target_face_vars[source_to_target[i]];
      vars[i + num_quad_pts*(num_swept_pts-1)] = target_var;
    }

    // Set the source face variable numbers
    const int *source_face_vars;
    source_mesh->getNodeNums(&source_face_vars);
    for ( int i = 0; i < num_quad_pts; i++ ){
      vars[i] = source_face_vars[i];
    }

    // Now the target and source surfaces of the volume have the
    // correct ordering, but the structured sides are not ordered
    // correctly. Loop through the structured sides (with the same
    // ordering as defined by the edge loops on the source surface)
    // and determine the node ordering.
    int num_source_loops = source->getNumEdgeLoops();

    for ( int k = 0, count = 0; k < num_source_loops; k++ ){
      TMREdgeLoop *source_loop;
      source->getEdgeLoop(k, &source_loop);

      // Get the number of edges for this loop
      int num_loop_edges;
      TMREdge **loop_edges;
      source_loop->getEdgeLoop(&num_loop_edges, &loop_edges, NULL);

      for ( int i = 0; i < num_loop_edges; i++, count++ ){
        // Get the loop edge mesh
        TMREdgeMesh *edge_mesh;
        loop_edges[i]->getMesh(&edge_mesh);

        // Get the swept face mesh (that is structured)
        TMRFaceMesh *swept_mesh;
        swept_faces[count]->getMesh(&swept_mesh);

        const int *swept_face_vars;
        int nswept = swept_mesh->getNodeNums(&swept_face_vars);

        // Get the number of points along this edge
        int npts;
        edge_mesh->getMeshPoints(&npts, NULL, NULL);

        // Loop over the edges on this face
        for ( int ix = 0; ix < npts; ix++ ){
          // Get the index on the source face mesh of the edge at the
          // specified edge index
          int src_face_index =
            source_mesh->getFaceIndexFromEdge(loop_edges[i], ix);

          if (src_face_index < 0 ||
              src_face_index >= num_quad_pts){
            fprintf(stderr,
                    "TMRVolume error: Source surface index %d out of range\n",
                    src_face_index);
          }
          else {
            // Loop over the through-swept nodes
            int iy_index = num_swept_pts-1;
            int orient = swept_edges_orient[count];
            if (orient > 0){
              iy_index = 0;
            }

            for ( int iy = 0; iy < num_swept_pts; iy++, iy_index += orient ){
              int swept_face_index =
                swept_mesh->getStructuredFaceIndex(loop_edges[i], ix,
                                                   swept_edges[count], iy_index);

              if (swept_face_index < 0 || swept_face_index >= nswept){
                fprintf(stderr,
                        "TMRVolume error: Swept face index %d out of range\n",
                        swept_face_index);
              }
              else {
                vars[src_face_index + iy*num_quad_pts] =
                  swept_face_vars[swept_face_index];
              }
            }
          }
        }
      }
    }

    // Now order the variables as they arrive
    int start = *num;
    for ( int k = 0; k < num_points; k++ ){
      if (vars[k] < 0){
        vars[k] = *num;
        (*num)++;
      }
    }

    // Return the number of points that have been allocated
    return *num - start;
  }

  return 0;
}

/*
  Retrieve the global node numbers for the mesh
*/
int TMRVolumeMesh::getNodeNums( const int **_vars ){
  if (_vars){
    *_vars = vars;
  }
  return num_points;
}
