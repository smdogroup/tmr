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

#include "TMRFaceMesh.h"

#include <math.h>
#include <stdio.h>

#include "TMRMesh.h"
#include "TMRMeshSmoothing.h"
#include "TMRNativeTopology.h"
#include "TMRPerfectMatchInterface.h"
#include "TMRTriangularize.h"
#include "tmrlapack.h"

/*
  The triangle nodes and edges are ordered locally as follows. Note
  that the edges are ordered based on the node across the triangle.

       2
      / \
   e1/   \e0
    /     \
   /       \
  0 ------- 1
       e2

  The edge nodes for a given edge index in a triangle
*/
const int tri_edge_nodes[][2] = {{1, 2}, {2, 0}, {0, 1}};

/*
  The edges that are connected to a given node
*/
const int tri_node_edges[][2] = {{1, 2}, {0, 2}, {0, 1}};

/*
  Returns the index number for the (i, j) node location along
  the structured edge.

  This the structured face nodes are ordered by first counting around
  the edges of the face counter-clockwise. After these nodes are
  counted, the interior face nodes are ordered in a cartesian order.
  For nx = 5, and ny = 4, this produces the following face numbering:

  11-- 10-- 9 -- 8 -- 7
  |    |    |    |    |
  12-- 17-- 18-- 19-- 6
  |    |    |    |    |
  13-- 14-- 15-- 16-- 5
  |    |    |    |    |
  0 -- 1 -- 2 -- 3 -- 4
*/
inline int get_structured_index(const int nx, const int ny, const int i,
                                const int j) {
  if (j == 0) {
    return i;
  } else if (i == nx - 1) {
    return nx - 1 + j;
  } else if (j == ny - 1) {
    return 2 * nx + ny - 3 - i;
  } else if (i == 0) {
    return 2 * nx + 2 * ny - 4 - j;
  }

  return 2 * nx + 2 * ny - 4 + (i - 1) + (j - 1) * (nx - 2);
}

/*
  Create the surface mesh object

  This call does not create the underlying surface. You can set
  meshing options into the class after it is created, and then create
  the mesh.

  Note that the curve/edge meshes must be meshed before calling this
  object.
*/
TMRFaceMesh::TMRFaceMesh(MPI_Comm _comm, TMRFace *_face, TMRPoint *_X,
                         int _npts, const int *_quads, int _nquads) {
  comm = _comm;
  face = _face;
  face->incref();
  mesh_type = TMR_NO_MESH;

  num_fixed_pts = 0;
  num_points = 0;
  num_quads = 0;
  num_tris = 0;

  // NULL things that will be used later
  pts = NULL;
  X = NULL;
  vars = NULL;
  quads = NULL;
  tris = NULL;
  source_to_target = NULL;
  copy_to_target = NULL;

  // This is not a prescribed mesh
  prescribed_mesh = 0;

  // Check if the input mesh is consistent
  if (_X && _quads && _npts > 0 && _nquads > 0) {
    prescribed_mesh = 1;
    setPrescribedMesh(_X, _npts, _quads, _nquads);
  }
}

/*
  Free the data associated with the surface mesh
*/
TMRFaceMesh::~TMRFaceMesh() {
  face->decref();
  if (pts) {
    delete[] pts;
  }
  if (X) {
    delete[] X;
  }
  if (vars) {
    delete[] vars;
  }
  if (quads) {
    delete[] quads;
  }
  if (tris) {
    delete[] tris;
  }
  if (source_to_target) {
    delete[] source_to_target;
  }
  if (copy_to_target) {
    delete[] copy_to_target;
  }
}

/*
  Set the prescribed mesh points/locations
*/
void TMRFaceMesh::setPrescribedMesh(const TMRPoint *_X, int _npts,
                                    const int *_quads, int _nquads) {
  // Check that all the edge meshes are already defined...
  for (int k = 0; k < face->getNumEdgeLoops(); k++) {
    // Get the curve information for this loop segment
    TMREdgeLoop *loop;
    face->getEdgeLoop(k, &loop);

    int nedges;
    TMREdge **edges;
    const int *dir;
    loop->getEdgeLoop(&nedges, &edges, &dir);

    for (int i = 0; i < nedges; i++) {
      // Retrieve the underlying curve mesh
      TMREdgeMesh *mesh = NULL;
      edges[i]->getMesh(&mesh);

      if (!mesh) {
        prescribed_mesh = 0;
      }
    }
  }

  // The edge meshes are not defined: cannot prescribe the mesh
  if (!prescribed_mesh) {
    fprintf(stderr,
            "TMRFaceMesh Error: Must prescribe all edge meshes "
            "for a prescribed face mesh\n");
    return;
  }

  // Otherwise, we can now prescribe the mesh
  mesh_type = TMR_UNSTRUCTURED;
  num_points = _npts;
  num_quads = _nquads;

  // Compute a mapping that places the edge nodes first in the correct
  // order
  int *node_mapping = new int[num_points];

  // Set all the variable values to negative
  for (int i = 0; i < num_points; i++) {
    node_mapping[i] = -1;
  }

  // Keep track of whether we encounter any problems
  int fail = 0;

  // Get the face orientation
  int face_orient = face->getOrientation();

  // Retrieve the boundary node numbers from the surface loops
  int pt = 0;
  for (int k = 0; k < face->getNumEdgeLoops(); k++) {
    // Get the curve information for this loop segment
    TMREdgeLoop *loop;
    face->getEdgeLoop(k, &loop);
    int nedges;
    TMREdge **edges;
    const int *edge_orient;
    loop->getEdgeLoop(&nedges, &edges, &edge_orient);

    int edge_index = nedges - 1;
    if (face_orient > 0) {
      edge_index = 0;
    }

    for (int i = 0; i < nedges; i++, edge_index += face_orient) {
      // Retrieve the underlying curve mesh
      TMREdge *edge = edges[edge_index];
      TMREdgeMesh *mesh = NULL;
      edge->getMesh(&mesh);

      // Get the points along the edge
      int npts;
      TMRPoint *Xedge;
      mesh->getMeshPoints(&npts, NULL, &Xedge);

      // Find the orientation on the edge
      int orientation = face_orient * edge_orient[edge_index];

      int index = npts - 1;
      if (orientation > 0) {
        index = 0;
      }

      for (int j = 0; j < npts - 1; j++, index += orientation) {
        int min_index = -1;
        double min_dist = 1e20;
        for (int point = 0; point < num_points; point++) {
          double dist =
              (Xedge[index].x - _X[point].x) * (Xedge[index].x - _X[point].x) +
              (Xedge[index].y - _X[point].y) * (Xedge[index].y - _X[point].y) +
              (Xedge[index].z - _X[point].z) * (Xedge[index].z - _X[point].z);
          if (dist < min_dist) {
            min_dist = dist;
            min_index = point;
          }
        }
        if (node_mapping[min_index] == -1) {
          node_mapping[min_index] = pt;
          pt++;
        } else {
          fail = 1;
        }
      }
    }
  }

  // Abandon the prescribed mesh
  if (fail) {
    prescribed_mesh = 0;
    mesh_type = TMR_NO_MESH;
    num_points = 0;
    num_quads = 0;
    delete[] node_mapping;
    fprintf(stderr, "TMRFaceMesh Error: Failed to map edge nodes\n");
    return;
  }

  // Set the number of fixed points
  num_fixed_pts = pt;

  // Find the un-ordered nodes
  for (int i = 0; i < num_points; i++) {
    if (node_mapping[i] < 0) {
      node_mapping[i] = pt;
      pt++;
    }
  }

  // Copy over the quadrilaterals
  quads = new int[4 * num_quads];
  memcpy(quads, _quads, 4 * num_quads * sizeof(int));
  for (int i = 0; i < 4 * num_quads; i++) {
    quads[i] = node_mapping[quads[i]];
  }

  // Now copy over the node locations
  X = new TMRPoint[num_points];
  for (int i = 0; i < num_points; i++) {
    X[node_mapping[i]] = _X[i];
  }

  // Free the mapping
  delete[] node_mapping;

  // Set the point locations
  pts = new double[2 * num_points];
  for (int i = 0; i < num_points; i++) {
    face->invEvalPoint(X[i], &pts[2 * i], &pts[2 * i + 1]);
  }
}

/*
  Create the surface mesh
*/
void TMRFaceMesh::mesh(TMRMeshOptions options, TMRElementFeatureSize *fs) {
  // Check if the mesh has already been allocated
  if (prescribed_mesh) {
    return;
  }

  int mpi_rank, mpi_size;
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);

  // Set the default mesh type
  TMRFaceMeshType _mesh_type = options.mesh_type_default;
  if (_mesh_type == TMR_NO_MESH) {
    _mesh_type = TMR_STRUCTURED;
  }

  // Get the source face and its orientation relative to this
  // face. Note that the source face may be NULL in which case the
  // source orientation is meaningless.
  TMRFace *source, *copy;
  face->getSource(NULL, &source);
  face->getCopySource(NULL, &copy);

  // If the face mesh for the source or copy source does not yet
  // exist, create it...
  if (source) {
    TMRFaceMesh *face_mesh;
    source->getMesh(&face_mesh);
    if (!face_mesh) {
      face_mesh = new TMRFaceMesh(comm, source);
      face_mesh->mesh(options, fs);
      source->setMesh(face_mesh);
    }
  } else if (copy) {
    TMRFaceMesh *face_mesh;
    copy->getMesh(&face_mesh);
    if (!face_mesh) {
      face_mesh = new TMRFaceMesh(comm, copy);
      face_mesh->mesh(options, fs);
      copy->setMesh(face_mesh);
    }
  }

  // First check if the conditions for a structured mesh are satisfied
  if (_mesh_type == TMR_STRUCTURED) {
    int nloops = face->getNumEdgeLoops();
    if (nloops != 1) {
      _mesh_type = TMR_UNSTRUCTURED;
    }

    // Get the first edge loop and the edges in the loop
    TMREdgeLoop *loop;
    face->getEdgeLoop(0, &loop);
    int nedges;
    TMREdge **edges;
    loop->getEdgeLoop(&nedges, &edges, NULL);
    if (nedges != 4) {
      _mesh_type = TMR_UNSTRUCTURED;
    }
    for (int k = 0; k < nedges; k++) {
      if (edges[k]->isDegenerate()) {
        _mesh_type = TMR_UNSTRUCTURED;
      }
    }

    if (nedges == 4) {
      // Check that parallel edges have the same number of nodes
      int nx1 = 0, nx2 = 0, ny1 = 0, ny2 = 0;
      TMREdgeMesh *mesh;
      for (int k = 0; k < 4; k++) {
        edges[k]->getMesh(&mesh);
        if (!mesh) {
          fprintf(stderr, "TMRFaceMesh Error: Edge mesh does not exist\n");
        }
      }

      edges[0]->getMesh(&mesh);
      mesh->getMeshPoints(&nx1, NULL, NULL);
      edges[2]->getMesh(&mesh);
      mesh->getMeshPoints(&nx2, NULL, NULL);
      if (nx1 != nx2) {
        _mesh_type = TMR_UNSTRUCTURED;
      }

      edges[1]->getMesh(&mesh);
      mesh->getMeshPoints(&ny1, NULL, NULL);
      edges[3]->getMesh(&mesh);
      mesh->getMeshPoints(&ny2, NULL, NULL);
      if (ny1 != ny2) {
        _mesh_type = TMR_UNSTRUCTURED;
      }
    }
  }

  // Record the mesh type
  mesh_type = _mesh_type;

  if (mpi_rank == 0) {
    // Count up the number of points and segments from the curves that
    // bound the surface. Keep track of the number of points = the
    // number of segments.
    int total_num_pts = 0;

    // Get the face orientation
    int face_orient = face->getOrientation();

    // Keep track of the number of closed loop cycles in the domain
    int nloops = face->getNumEdgeLoops();

    // The number of degenerate edges
    int num_degen = 0;

    // Get all of the edges and count up the mesh points
    for (int k = 0; k < nloops; k++) {
      TMREdgeLoop *loop;
      face->getEdgeLoop(k, &loop);
      int nedges;
      TMREdge **edges;
      loop->getEdgeLoop(&nedges, &edges, NULL);

      for (int i = 0; i < nedges; i++) {
        // Count whether this edge is degenerate
        if (edges[i]->isDegenerate()) {
          num_degen++;
        }

        // Check whether the edge mesh exists - it has to!
        TMREdgeMesh *mesh = NULL;
        edges[i]->getMesh(&mesh);
        if (!mesh) {
          fprintf(stderr, "TMRFaceMesh Error: Edge mesh does not exist\n");
        }

        // Get the number of points associated with the curve
        int npts;
        mesh->getMeshPoints(&npts, NULL, NULL);

        // Update the total number of points
        total_num_pts += npts - 1;
      }
    }

    // The number of holes is equal to the number of loops-1. One loop
    // bounds the domain, the other loops cut out holes in the domain.
    // Note that the domain must be contiguous.
    int nholes = nloops - 1;

    // All the boundary loops are closed, therefore, the total number
    // of segments is equal to the total number of points
    int nsegs = total_num_pts;

    // Set the maximum number of extra segments that will be added to
    // handle problematic corners
    const int max_extra_segs = 128;
    const int max_extra_pts = 128;

    // Keep track of the beginning/end of each llop
    int *loop_pt_offset = new int[nloops + 1];

    // Allocate the points and the number of segments based on the
    // number of holes
    double *params = new double[2 * (total_num_pts + nholes + max_extra_pts)];
    int *segments = new int[2 * (nsegs + max_extra_segs)];

    // Start entering the points from the end of the last hole entry in
    // the parameter points array.
    int pt = 0;

    // Set up the degenerate edges
    int *degen = NULL;
    if (num_degen > 0) {
      degen = new int[2 * num_degen];
    }
    num_degen = 0;

    for (int k = 0; k < nloops; k++) {
      // Set the offset to the initial point/segment on this loop
      loop_pt_offset[k] = pt;

      // Get the curve information for this loop segment
      TMREdgeLoop *loop;
      face->getEdgeLoop(k, &loop);
      int nedges;
      TMREdge **edges;
      const int *edge_orient;
      loop->getEdgeLoop(&nedges, &edges, &edge_orient);

      int edge_index = nedges - 1;
      if (face_orient > 0) {
        edge_index = 0;
      }

      for (int i = 0; i < nedges; i++, edge_index += face_orient) {
        // Retrieve the underlying curve mesh
        TMREdge *edge = edges[edge_index];
        TMREdgeMesh *mesh = NULL;
        edge->getMesh(&mesh);

        // Get the mesh points corresponding to this curve
        int npts;
        const double *tpts;
        mesh->getMeshPoints(&npts, &tpts, NULL);

        // Get the orientation of the edge
        int orientation = face_orient * edge_orient[edge_index];

        int index = npts - 1;
        if (orientation > 0) {
          index = 0;
        }

        for (int j = 0; j < npts - 1; j++, index += orientation) {
          int info =
              edge->getParamsOnFace(face, tpts[index], edge_orient[edge_index],
                                    &params[2 * pt], &params[2 * pt + 1]);
          if (info != 0) {
            fprintf(stderr,
                    "TMRFaceMesh Error: getParamsOnFace "
                    "failed with error code %d\n",
                    info);
          } else {
            segments[2 * pt] = pt;
            segments[2 * pt + 1] = pt + 1;
            if (edge->isDegenerate()) {
              degen[2 * num_degen] = pt;
              degen[2 * num_degen + 1] = pt + 1;
              num_degen++;
            }
            pt++;
          }
        }
      }

      // Close off the loop by connecting the segment back to the
      // initial loop point
      segments[2 * (pt - 1) + 1] = loop_pt_offset[k];
    }

    // Set the last loop
    loop_pt_offset[nloops] = pt;

    // Set the total number of fixed points. These are the points that
    // will not be smoothed and constitute the boundary nodes. Note
    // that the Triangularize class removes the holes from the domain
    // automatically.  The boundary points are guaranteed to be
    // ordered first.
    num_fixed_pts = total_num_pts - num_degen;

    if (source) {
      mapSourceToTarget(options, params);
    } else if (copy) {
      mapCopyToTarget(options, params);
    } else if (mesh_type == TMR_STRUCTURED) {
      createStructuredMesh(options, params);
    } else if (mesh_type == TMR_TRIANGLE) {
      // Compute hole points (inside the holes)
      computeHolePts(nloops, total_num_pts, loop_pt_offset, segments, params);

      // Create an unstructured triangular mesh
      createUnstructuredMesh(options, fs, mesh_type, total_num_pts, nholes,
                             params, nsegs, segments, num_degen, degen,
                             &num_points, &pts, &X, &num_quads, &quads,
                             &num_tris, &tris);
    } else if (mesh_type == TMR_UNSTRUCTURED) {
      // Loop over the segments in the mesh to find corners
      // with angles less than 60 degrees. These corners will be
      // cut and replaced with a specified quadrilateral corner pattern

      // Evaluate all of the points around the edge
      TMRPoint *Xparam = new TMRPoint[total_num_pts];
      for (int i = 0; i < total_num_pts; i++) {
        face->evalPoint(params[2 * i], params[2 * i + 1], &Xparam[i]);
      }

      // Go through the edge loops and find corners that will be problematic
      // for the quadrilateral mesh generator. Add extra segments
      // to alleviate the meshing issues in these corners.
      for (int loop = 0; loop < nloops; loop++) {
        for (int p = loop_pt_offset[loop]; p < loop_pt_offset[loop + 1];) {
          int incr = 1;
          int next = p + 1;
          int prev = p - 1;
          if (next >= loop_pt_offset[loop + 1]) {
            next = loop_pt_offset[loop];
          }
          if (prev < loop_pt_offset[loop]) {
            prev = loop_pt_offset[loop + 1] - 1;
          }

          // Compute the difference
          TMRPoint d1, d2;
          d1.x = Xparam[p].x - Xparam[prev].x;
          d1.y = Xparam[p].y - Xparam[prev].y;
          d1.z = Xparam[p].z - Xparam[prev].z;
          d2.x = Xparam[next].x - Xparam[p].x;
          d2.y = Xparam[next].y - Xparam[p].y;
          d2.z = Xparam[next].z - Xparam[p].z;

          // Compute the dot product of the two vectors
          double d1dist = sqrt(d1.dot(d1));
          double d2dist = sqrt(d2.dot(d2));
          double dot = -d1.dot(d2) / (d1dist * d2dist);

          // If the dot product is such that the angle is
          // less than about 75 degrees, add segments to ensure
          // that elements are created on either side of the segment
          if (dot > 0.25) {
            // Set the first point in the new segment list
            segments[2 * nsegs] = pt;
            segments[2 * nsegs + 1] = pt + 1;
            nsegs++;

            // Insert the new point
            TMRPoint Xmid;
            Xmid.x = 0.5 * (Xparam[next].x + Xparam[prev].x);
            Xmid.y = 0.5 * (Xparam[next].y + Xparam[prev].y);
            Xmid.z = 0.5 * (Xparam[next].z + Xparam[prev].z);

            face->invEvalPoint(Xmid, &params[2 * pt], &params[2 * pt + 1]);
            pt++;

            const int max_new_corner_segments = 4;
            for (int i = 0; i < max_new_corner_segments; i++) {
              // Increment the pointers to the next/previous index
              next++;
              prev--;
              if (next >= loop_pt_offset[loop + 1]) {
                next = loop_pt_offset[loop];
              }
              if (prev < loop_pt_offset[loop]) {
                prev = loop_pt_offset[loop + 1] - 1;
              }

              // Compute the mid-point
              Xmid.x = 0.5 * (Xparam[next].x + Xparam[prev].x);
              Xmid.y = 0.5 * (Xparam[next].y + Xparam[prev].y);
              Xmid.z = 0.5 * (Xparam[next].z + Xparam[prev].z);

              // Find the mid-point in the parametric space
              face->invEvalPoint(Xmid, &params[2 * pt], &params[2 * pt + 1]);
              pt++;

              // Find the vector between the two points on the boundary
              TMRPoint d3;
              d3.x = Xparam[next].x - Xparam[prev].x;
              d3.y = Xparam[next].y - Xparam[prev].y;
              d3.z = Xparam[next].z - Xparam[prev].z;

              // If the distance between the next/prev values
              // is less than
              double d3dist = sqrt(d3.dot(d3));
              if (d3dist > 0.75 * (d1dist + d2dist)) {
                break;
              } else if (i + 1 < max_new_corner_segments) {
                // Add the next segment
                segments[2 * nsegs] = pt - 1;
                segments[2 * nsegs + 1] = pt;
                nsegs++;
              }
            }
          }

          p += incr;
        }
      }

      // Reset the total number of points
      total_num_pts = pt;

      // Compute hole points (inside the holes)
      computeHolePts(nloops, total_num_pts, loop_pt_offset, segments, params);

      // Create the unstructured mesh
      createUnstructuredMesh(options, fs, mesh_type, total_num_pts, nholes,
                             params, nsegs, segments, num_degen, degen,
                             &num_points, &pts, &X, &num_quads, &quads,
                             &num_tris, &tris);

      // Free the triangles - these are not needed for this type of mesh
      delete[] tris;
      num_tris = 0;
      tris = NULL;

      // Build connectivity to smooth the quad mesh
      int *pts_to_quad_ptr;
      int *pts_to_quads;
      TMR_ComputeNodeToElems(num_points, num_quads, 4, quads, &pts_to_quad_ptr,
                             &pts_to_quads);

      // Smooth the mesh using a local optimization of node locations
      TMR_QuadSmoothing(options.num_smoothing_steps, num_fixed_pts, num_points,
                        pts_to_quad_ptr, pts_to_quads, num_quads, quads, pts, X,
                        face);

      // Free the connectivity information
      delete[] pts_to_quad_ptr;
      delete[] pts_to_quads;

      if (options.write_post_smooth_quad) {
        char filename[256];
        sprintf(filename, "post_smooth_quad%d.vtk", face->getEntityId());
        writeToVTK(filename);
      }
    }

    if (num_degen > 0) {
      delete[] degen;
    }

    // Free the parameter/segment information
    delete[] loop_pt_offset;
    delete[] params;
    delete[] segments;
  }

  if (mpi_size > 1) {
    // Broadcast the number of points to all the processors
    int temp[3];
    temp[0] = num_points;
    temp[1] = num_quads;
    temp[2] = num_fixed_pts;
    MPI_Bcast(temp, 3, MPI_INT, 0, comm);
    num_points = temp[0];
    num_quads = temp[1];
    num_fixed_pts = temp[2];

    if (mpi_rank != 0) {
      pts = new double[2 * num_points];
      X = new TMRPoint[num_points];
      quads = new int[4 * num_quads];
    }

    // Broadcast the parametric locations and points
    MPI_Bcast(pts, 2 * num_points, MPI_DOUBLE, 0, comm);
    MPI_Bcast(X, num_points, TMRPoint_MPI_type, 0, comm);
    MPI_Bcast(quads, 4 * num_quads, MPI_INT, 0, comm);

    // Broadcast the source to target information
    if (source) {
      if (mpi_rank != 0) {
        source_to_target = new int[num_points];
      }
      MPI_Bcast(source_to_target, num_points, MPI_INT, 0, comm);
    } else if (copy) {
      if (mpi_rank != 0) {
        copy_to_target = new int[num_points];
      }
      MPI_Bcast(copy_to_target, num_points, MPI_INT, 0, comm);
    }
  }
}

/*
  Compute the areas and set the hole point locations
*/
void TMRFaceMesh::computeHolePts(const int nloops, int hole_pt,
                                 const int *loop_pt_offset, const int *segments,
                                 double *params) {
  for (int loop = 0; loop < nloops; loop++) {
    // Compute the area enclosed by the loop. If the area is
    // positive, it is the domain boundary. If the area is negative,
    // we have a hole!  Note that this assumes that the polygon
    // creating the hole is not self-intersecting. (In reality we
    // compute twice the area since we omit the 1/2 factor.)
    double Area = 0.0;
    for (int i = loop_pt_offset[loop]; i < loop_pt_offset[loop + 1]; i++) {
      int s1 = segments[2 * i];
      int s2 = segments[2 * i + 1];
      const double x1 = params[2 * s1];
      const double y1 = params[2 * s1 + 1];
      const double x2 = params[2 * s2];
      const double y2 = params[2 * s2 + 1];
      Area += (x1 * y2 - x2 * y1);
    }

    // Check the area constraint
    if (Area < 0.0) {
      // This is a hole! Compute an approximate position for the hole.
      // Note that this may not work in all cases so beware.
      int s1 = segments[2 * loop_pt_offset[loop]];
      int s2 = segments[2 * loop_pt_offset[loop] + 1];
      const double x1 = params[2 * s1];
      const double y1 = params[2 * s1 + 1];
      const double x2 = params[2 * s2];
      const double y2 = params[2 * s2 + 1];
      const double dx = x2 - x1;
      const double dy = y2 - y1;

      // This is arbitrary and won't work in general if we have a very
      // thin sliver for a hole...
      double frac = 0.01;

      // Set the average location for the hole
      params[2 * hole_pt] = 0.5 * (x1 + x2) + frac * dy;
      params[2 * hole_pt + 1] = 0.5 * (y1 + y2) - frac * dx;

      // Increment the hole pointer
      hole_pt++;
    }
  }
}

/*
  Map the source mesh to the target mesh (this is the target mesh)

  Note that the source and target surface must be contained within the
  same volume. Both faces have outward-facing normals relative to the
  volume. However, the surface meshes are stored in an orientation
  consistent with the underlying surface (which may be reversed
  relative to the face). Therefore, care must be taken to properly
  re-orient the surface mesh.
*/
void TMRFaceMesh::mapSourceToTarget(TMRMeshOptions options,
                                    const double *params) {
  TMRVolume *source_volume;
  TMRFace *source;
  face->getSource(&source_volume, &source);

  // Create the source map of edges and keep track of their local
  // directions relative to the source surface
  std::map<TMREdge *, int> source_edge_orient;
  int source_orient = source->getOrientation();
  for (int k = 0; k < source->getNumEdgeLoops(); k++) {
    TMREdgeLoop *loop;
    source->getEdgeLoop(k, &loop);

    // Get the number of edges from the source loop
    int nedges;
    TMREdge **edges;
    const int *edge_orient;
    loop->getEdgeLoop(&nedges, &edges, &edge_orient);
    for (int j = 0; j < nedges; j++) {
      source_edge_orient[edges[j]] = source_orient * edge_orient[j];
    }
  }

  // Create the target map of edges and keep track of their local
  // directions relative to the target surface
  std::map<TMREdge *, int> target_edge_orient;
  int target_orient = face->getOrientation();
  for (int k = 0; k < face->getNumEdgeLoops(); k++) {
    TMREdgeLoop *loop;
    face->getEdgeLoop(k, &loop);

    // Get the number of edges/edges from the source loop
    int nedges;
    TMREdge **edges;
    const int *edge_orient;
    loop->getEdgeLoop(&nedges, &edges, &edge_orient);
    for (int j = 0; j < nedges; j++) {
      target_edge_orient[edges[j]] = target_orient * edge_orient[j];
    }
  }

  // Compute the relative orientation of the target and face surfaces
  // (and their meshes). The negative sign is due to the fact that the
  // faces are all oriented outwards within the same volume.
  int rel_orient = -target_orient * source_orient;

  // Keep track of the source-to-target edge and target-to-source
  // edge mappings as well as their relative orientations
  std::map<TMREdge *, TMREdge *> target_to_source_edge;

  // Loop over the faces that are within the source volume
  int num_faces;
  TMRFace **faces;
  source_volume->getFaces(&num_faces, &faces);

  for (int i = 0; i < num_faces; i++) {
    // Check that this is not a target or source face
    if (faces[i] != source && faces[i] != face) {
      // Find the source and target edge that are both contained in
      // one of the edge loops for this face
      TMREdge *src_edge = NULL, *tar_edge = NULL;
      for (int k = 0; k < faces[i]->getNumEdgeLoops(); k++) {
        src_edge = tar_edge = NULL;

        // Get the edge loop
        TMREdgeLoop *loop;
        faces[i]->getEdgeLoop(k, &loop);

        // Get the number of edges from the source loop
        int nedges;
        TMREdge **edges;
        const int *edge_orient;
        loop->getEdgeLoop(&nedges, &edges, &edge_orient);

        // Determine which edge is shared
        for (int j = 0; j < nedges; j++) {
          if (target_edge_orient.count(edges[j]) > 0) {
            tar_edge = edges[j];
          }
          if (source_edge_orient.count(edges[j]) > 0) {
            src_edge = edges[j];
          }
        }

        if (src_edge && tar_edge) {
          break;
        }
      }

      // Set the target to source relationship
      if (src_edge && tar_edge) {
        target_to_source_edge[tar_edge] = src_edge;
      }
    }
  }

  // Get the source mesh
  TMRFaceMesh *source_face_mesh;
  source->getMesh(&source_face_mesh);

  // Set the total number of points
  mesh_type = source_face_mesh->mesh_type;
  num_points = source_face_mesh->num_points;
  num_fixed_pts = source_face_mesh->num_fixed_pts;
  num_quads = source_face_mesh->num_quads;
  num_tris = source_face_mesh->num_tris;

  // Allocate the source to target mapping array
  source_to_target = new int[num_points];

  for (int k = 0; k < face->getNumEdgeLoops(); k++) {
    // Get the curve information for this loop segment
    TMREdgeLoop *loop;
    face->getEdgeLoop(k, &loop);

    // Extract the edges and orientations from the loop
    int nedges;
    TMREdge **edges;
    const int *edge_orient;
    loop->getEdgeLoop(&nedges, &edges, &edge_orient);

    for (int i = 0; i < nedges; i++) {
      // Get the edge mesh and the number of points for this edge
      int npts = 0;
      TMREdgeMesh *mesh = NULL;
      edges[i]->getMesh(&mesh);
      mesh->getMeshPoints(&npts, NULL, NULL);

      // Get the source edge
      int source_npts = 0;
      TMREdgeMesh *source_mesh = NULL;
      TMREdge *source_edge = target_to_source_edge[edges[i]];
      source_edge->getMesh(&source_mesh);
      source_mesh->getMeshPoints(&source_npts, NULL, NULL);

      if (npts != source_npts) {
        fprintf(stderr,
                "TMRFaceMesh Error: Inconsistent number of points "
                "on source and target edges\n");
      }

      // Get the orientation of the target edge
      int orient_target = edge_orient[i];
      orient_target *= face->getOrientation();
      int j_target = npts - 1;
      if (orient_target > 0) {
        j_target = 0;
      }

      // Get the orientation of the source edge
      int orient_source = source_edge_orient[source_edge];
      orient_source *= rel_orient;
      int j_source = npts - 1;
      if (orient_source > 0) {
        j_source = 0;
      }

      for (int j = 0; j < npts - 1;
           j++, j_target += orient_target, j_source += orient_source) {
        int source_index =
            source_face_mesh->getFaceIndexFromEdge(source_edge, j_source);
        int target_index = getFaceIndexFromEdge(edges[i], j_target);
        source_to_target[source_index] = target_index;
      }
    }
  }

  // Set the remaining source-to-target mapping
  for (int i = num_fixed_pts; i < num_points; i++) {
    source_to_target[i] = i;
  }

  // Allocate the array for the parametric locations
  pts = new double[2 * num_points];

  // Copy the points from around the boundaries
  for (int i = 0; i < num_fixed_pts; i++) {
    pts[2 * i] = params[2 * i];
    pts[2 * i + 1] = params[2 * i + 1];
  }

  // Compute a least squares transformation between the two
  // surfaces
  double N[16], A[4];
  double sc[2], tc[3];
  memset(N, 0, 16 * sizeof(double));
  memset(A, 0, 4 * sizeof(double));
  sc[0] = sc[1] = 0.0;
  tc[0] = tc[1] = 0.0;

  for (int k = 0; k < num_fixed_pts; k++) {
    sc[0] += source_face_mesh->pts[2 * k];
    sc[1] += source_face_mesh->pts[2 * k + 1];
    tc[0] += pts[2 * k];
    tc[1] += pts[2 * k + 1];
  }
  sc[0] = sc[0] / num_fixed_pts;
  sc[1] = sc[1] / num_fixed_pts;
  tc[0] = tc[0] / num_fixed_pts;
  tc[1] = tc[1] / num_fixed_pts;

  for (int k = 0; k < num_fixed_pts; k++) {
    double uS[2], uT[2];
    uS[0] = source_face_mesh->pts[2 * k] - sc[0];
    uS[1] = source_face_mesh->pts[2 * k + 1] - sc[1];

    // Compute the source->target index number
    int kt = source_to_target[k];
    uT[0] = pts[2 * kt] - tc[0];
    uT[1] = pts[2 * kt + 1] - tc[1];

    // Add the terms to the matrix/right-hand-side
    for (int i = 0; i < 4; i++) {
      int start = i / 2, end = i / 2 + 1;
      for (int j = 2 * start; j < 2 * end; j++) {
        N[i + 4 * j] += uS[i % 2] * uS[j % 2];
      }
      A[i] += uS[i % 2] * uT[i / 2];
    }
  }

  // Factor the least-squares matrix and perform the transformation
  int ipiv[4];
  int n = 4, one = 1, info;
  TmrLAPACKdgetrf(&n, &n, N, &n, ipiv, &info);
  TmrLAPACKdgetrs("N", &n, &one, N, &n, ipiv, A, &n, &info);

  // Set the interior points based on the linear transformation
  for (int k = num_fixed_pts; k < num_points; k++) {
    double uS = source_face_mesh->pts[2 * k] - sc[0];
    double vS = source_face_mesh->pts[2 * k + 1] - sc[1];
    pts[2 * k] = A[0] * uS + A[1] * vS + tc[0];
    pts[2 * k + 1] = A[2] * uS + A[3] * vS + tc[1];
  }

  // Copy the quadrilateral mesh (if any)
  if (num_quads > 0) {
    quads = new int[4 * num_quads];
    memcpy(quads, source_face_mesh->quads, 4 * num_quads * sizeof(int));

    // Adjust the quadrilateral ordering at the boundary
    for (int i = 0; i < 4 * num_quads; i++) {
      if (quads[i] < num_fixed_pts) {
        quads[i] = source_to_target[quads[i]];
      }
    }

    // Flip the orientation of the quads to match the orientation of
    // the face if needed
    if (rel_orient < 0) {
      for (int i = 0; i < num_quads; i++) {
        int tmp = quads[4 * i + 1];
        quads[4 * i + 1] = quads[4 * i + 3];
        quads[4 * i + 3] = tmp;
      }
    }
  }

  // Copy the triangular mesh (if any)
  if (num_tris > 0) {
    tris = new int[3 * num_tris];
    memcpy(tris, source_face_mesh->tris, 3 * num_tris * sizeof(int));

    // Adjust the triangle ordering at the boundary
    for (int i = 0; i < 3 * num_tris; i++) {
      if (tris[i] < num_fixed_pts) {
        tris[i] = source_to_target[tris[i]];
      }
    }

    // Flip the orientation of the triangles to match the
    // orientation of the face
    int orient = -target_orient * source_orient;
    if (orient < 0) {
      for (int i = 0; i < num_tris; i++) {
        int tmp = tris[3 * i + 1];
        tris[4 * i + 1] = tris[4 * i + 2];
        tris[4 * i + 2] = tmp;
      }
    }
  }

  // Evaluate the points
  X = new TMRPoint[num_points];
  for (int i = 0; i < num_points; i++) {
    face->evalPoint(pts[2 * i], pts[2 * i + 1], &X[i]);
  }

  if (num_quads > 0) {
    // Smooth the copied mesh on the new surface
    int *pts_to_quad_ptr;
    int *pts_to_quads;
    TMR_ComputeNodeToElems(num_points, num_quads, 4, quads, &pts_to_quad_ptr,
                           &pts_to_quads);

    // Smooth the mesh using a local optimization of node locations
    TMR_QuadSmoothing(options.num_smoothing_steps, num_fixed_pts, num_points,
                      pts_to_quad_ptr, pts_to_quads, num_quads, quads, pts, X,
                      face);

    // Free the connectivity information
    delete[] pts_to_quad_ptr;
    delete[] pts_to_quads;
  } else if (num_tris > 0) {
    // Compute the triangle edges and neighbors in the dual mesh
    int num_tri_edges;
    int *tri_edges, *tri_neighbors, *dual_edges;
    TMR_ComputePlanarTriEdges(num_points, num_tris, tris, &num_tri_edges,
                              &tri_edges, &tri_neighbors, &dual_edges);

    // Smooth the resulting triangular mesh
    if (options.tri_smoothing_type == TMRMeshOptions::TMR_LAPLACIAN) {
      TMR_LaplacianSmoothing(options.num_smoothing_steps, num_fixed_pts,
                             num_tri_edges, tri_edges, num_points, pts, X,
                             face);
    } else {
      double alpha = 0.1;
      TMR_SpringSmoothing(options.num_smoothing_steps, alpha, num_fixed_pts,
                          num_tri_edges, tri_edges, num_points, pts, X, face);
    }

    delete[] tri_edges;
    delete[] tri_neighbors;
    delete[] dual_edges;
  }
}

/*
  Copy the mesh from the source copy mesh to this mesh. Use inverse
  evaluations to determine the parametric locations of the new mesh
  points.
*/
int TMRFaceMesh::mapCopyToTarget(TMRMeshOptions options, const double *params) {
  // Get the face orientation
  int face_orient = face->getOrientation();

  // Get the relative orientation
  int rel_orient;
  TMRFace *copy;
  face->getCopySource(&rel_orient, &copy);

  // Get the source mesh
  TMRFaceMesh *copy_mesh;
  copy->getMesh(&copy_mesh);

  // Set the total number of points
  mesh_type = copy_mesh->mesh_type;
  num_points = copy_mesh->num_points;
  num_fixed_pts = copy_mesh->num_fixed_pts;
  num_quads = copy_mesh->num_quads;
  num_tris = copy_mesh->num_tris;

  // Set the copy indices
  copy_to_target = new int[num_points];
  for (int i = 0; i < num_points; i++) {
    copy_to_target[i] = -1;
  }

  // Loop over all of the edges, and find their copy edge sources and
  // their orientations on the copied edges
  std::map<TMREdge *, int> copy_orient;
  for (int k = 0; k < copy->getNumEdgeLoops(); k++) {
    TMREdgeLoop *loop;
    copy->getEdgeLoop(k, &loop);

    int nedges;
    TMREdge **edges;
    const int *edge_orient;
    loop->getEdgeLoop(&nedges, &edges, &edge_orient);

    // Set the orientation on each edge
    for (int i = 0; i < nedges; i++) {
      // Index the edge orientation using the source edge
      TMREdge *t = NULL;
      edges[i]->getCopySource(&t);
      if (!t) {
        t = edges[i];
      }

      copy_orient[t] = copy->getOrientation() * edge_orient[i];
    }
  }

  for (int k = 0; k < face->getNumEdgeLoops(); k++) {
    // Get the curve information for this loop segment
    TMREdgeLoop *loop;
    face->getEdgeLoop(k, &loop);

    // Extract the edges and orientations from the loop
    int nedges;
    TMREdge **edges;
    const int *edge_orient;
    loop->getEdgeLoop(&nedges, &edges, &edge_orient);

    for (int i = 0; i < nedges; i++) {
      // Get the edge mesh
      TMREdgeMesh *mesh = NULL;
      edges[i]->getMesh(&mesh);

      int npts;
      mesh->getMeshPoints(&npts, NULL, NULL);

      // Get the edge that we're going to copy. This could be the same
      // edge in this loop (if it is the source)
      TMREdge *copy_edge;
      edges[i]->getCopySource(&copy_edge);
      if (!copy_edge) {
        copy_edge = edges[i];
      }

      if (copy_orient.count(copy_edge) == 0) {
        fprintf(stderr,
                "TMRFaceMesh Error: Copy edge not found in copy source face\n");
      } else {
        // Get the mesh for the x-direction and check its orientation
        int orient_target = edge_orient[i];
        orient_target *= face_orient;
        int j_target = npts - 1;
        if (orient_target > 0) {
          j_target = 0;
        }

        int orient_copy = TMREdgeMesh::getEdgeCopyOrient(edges[i]);
        orient_copy *= copy_orient[copy_edge];
        int j_copy = npts - 1;
        if (orient_copy > 0) {
          j_copy = 0;
        }

        for (int j = 0; j < npts - 1;
             j++, j_target += orient_target, j_copy += orient_copy) {
          int copy_index = copy_mesh->getFaceIndexFromEdge(copy_edge, j_copy);
          int target_index = getFaceIndexFromEdge(edges[i], j_target);
          copy_to_target[copy_index] = target_index;
        }
      }
    }
  }

  int count = 0;
  for (int i = 0; i < num_fixed_pts; i++) {
    if (copy_to_target[i] < 0) {
      count++;
      copy_to_target[i] = 0;
    }
  }
  if (count > 0) {
    fprintf(stderr, "TMRFaceMesh Error: %d errors in mapping copy to target\n",
            count);
  }

  // Allocate the array for the parametric locations
  pts = new double[2 * num_points];
  X = new TMRPoint[num_points];
  quads = new int[4 * num_quads];

  // Copy the points from around the boundaries
  for (int i = 0; i < num_fixed_pts; i++) {
    pts[2 * i] = params[2 * i];
    pts[2 * i + 1] = params[2 * i + 1];
  }

  // Set the edge locations using the locations on the other surface
  for (int i = 0; i < num_fixed_pts; i++) {
    X[copy_to_target[i]] = copy_mesh->X[i];
  }

  // Compute the relative orientations of the faces
  int orient = rel_orient * face->getOrientation() * copy->getOrientation();

  if (mesh_type == TMR_STRUCTURED) {
    // Get the first edge loop and the edges in the loop. There is
    // only one loop since this is a structured mesh
    TMREdgeLoop *loop;
    face->getEdgeLoop(0, &loop);

    // Get the edges associated with the edge loop
    TMREdge **edges;
    const int *edge_orient;
    loop->getEdgeLoop(NULL, &edges, &edge_orient);

    // Get the number of nodes for the x/y edges
    int nx = 0, ny = 0;
    TMREdgeMesh *mesh;

    // If the mesh is structured, copy its indices
    // int xorient = 1, yorient = 1;
    // TMREdge *copy_edge_x, *copy_edge_y;
    if (face_orient > 0) {
      // Get the mesh for the x-direction and check its orientation
      edges[0]->getMesh(&mesh);
      // edges[0]->getCopySource(&copy_edge_x);
      mesh->getMeshPoints(&nx, NULL, NULL);
      // xorient = TMREdgeMesh::getEdgeCopyOrient(edges[0])*edge_orient[0];

      // Get the mesh for the y-direction and check its orientation
      edges[1]->getMesh(&mesh);
      // edges[1]->getCopySource(&copy_edge_y);
      mesh->getMeshPoints(&ny, NULL, NULL);
      // yorient = TMREdgeMesh::getEdgeCopyOrient(edges[1])*edge_orient[1];
    } else {
      // Get the mesh for the y-direction and check its orientation
      edges[0]->getMesh(&mesh);
      // edges[0]->getCopySource(&copy_edge_y);
      mesh->getMeshPoints(&ny, NULL, NULL);
      // yorient = TMREdgeMesh::getEdgeCopyOrient(edges[0])*edge_orient[0];

      // Get the mesh for the x-direction and check its orientation
      edges[1]->getMesh(&mesh);
      // edges[1]->getCopySource(&copy_edge_x);
      mesh->getMeshPoints(&nx, NULL, NULL);
      // xorient = TMREdgeMesh::getEdgeCopyOrient(edges[1])*edge_orient[1];
    }

    // Create the structured quadrilaterals
    int *q = quads;
    for (int j = 0; j < ny - 1; j++) {
      for (int i = 0; i < nx - 1; i++) {
        // Compute the connectivity of the mesh
        q[0] = get_structured_index(nx, ny, i, j);
        q[1] = get_structured_index(nx, ny, i + 1, j);
        q[2] = get_structured_index(nx, ny, i + 1, j + 1);
        q[3] = get_structured_index(nx, ny, i, j + 1);
        q += 4;
      }
    }

    // Use a transfinite interpolation to determine the parametric
    // points where the interior nodes should be placed.
    for (int j = 1; j < ny - 1; j++) {
      for (int i = 1; i < nx - 1; i++) {
        double u = 1.0 * i / (nx - 1);
        double v = 1.0 * j / (ny - 1);

        // Compute the weights on the corners
        double c1 = (1.0 - u) * (1.0 - v);
        double c2 = u * (1.0 - v);
        double c3 = u * v;
        double c4 = (1.0 - u) * v;

        // Compute the weights on the curves
        double w1 = (1.0 - v);
        double w2 = u;
        double w3 = v;
        double w4 = (1.0 - u);

        // New parametric point
        int p = get_structured_index(nx, ny, i, j);

        // Boundary points that we're interpolating from
        int p1 = get_structured_index(nx, ny, i, 0);
        int p2 = get_structured_index(nx, ny, nx - 1, j);
        int p3 = get_structured_index(nx, ny, i, ny - 1);
        int p4 = get_structured_index(nx, ny, 0, j);

        // Evaluate the parametric points based on the transfinite
        // interpolation
        pts[2 * p] =
            ((w1 * pts[2 * p1] + w2 * pts[2 * p2] + w3 * pts[2 * p3] +
              w4 * pts[2 * p4]) -
             (c1 * pts[0] + c2 * pts[2 * (nx - 1)] +
              c3 * pts[2 * (nx + ny - 2)] + c4 * pts[2 * (2 * nx + ny - 3)]));

        pts[2 * p + 1] = ((w1 * pts[2 * p1 + 1] + w2 * pts[2 * p2 + 1] +
                           w3 * pts[2 * p3 + 1] + w4 * pts[2 * p4 + 1]) -
                          (c1 * pts[1] + c2 * pts[2 * (nx - 1) + 1] +
                           c3 * pts[2 * (nx + ny - 2) + 1] +
                           c4 * pts[2 * (2 * nx + ny - 3) + 1]));
      }
    }

    // Allocate and evaluate the new physical point locations
    for (int i = 0; i < num_points; i++) {
      face->evalPoint(pts[2 * i], pts[2 * i + 1], &X[i]);
    }

    double atol = 1e-6;
    for (int i = 0; i < num_points; i++) {
      TMRPoint a = X[i];
      for (int j = 0; j < num_points; j++) {
        TMRPoint b = copy_mesh->X[j];

        if (((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) +
             (a.z - b.z) * (a.z - b.z)) < atol * atol) {
          copy_to_target[j] = i;
          break;
        }
      }
    }

    /*

    // Set the indices depending on the orientation of the edges
    int iy = ny-2;
    if (yorient > 0){
      iy = 1;
    }
    for ( int j = 1; j < ny-1; j++, iy += yorient ){

      int ix = nx-2;
      if (xorient > 0){
        ix = 1;
      }
      for ( int i = 1; i < nx-1; i++, ix += xorient ){
        // Get the index for this (the target face)
        int target_index = get_structured_index(nx, ny, i, j);

        // Get the index on the source face
        int copy_index = copy_mesh->getStructuredFaceIndex(copy_edge_x, ix,
                                                           copy_edge_y, iy);

        // Set the copy-to-target index
        copy_to_target[copy_index] = target_index;

        // Copy over the X position
        X[target_index] = copy_mesh->X[copy_index];

        double u, v;
        int icode = face->invEvalPoint(X[target_index], &u, &v);
        if (icode){
          fprintf(stderr, "TMRFaceMesh Error: Inverse point evaluation "
                  "failed with code %d\n", icode);
        }

        // Set the target index
        pts[2*target_index] = u;
        pts[2*target_index+1] = v;
        face->evalPoint(u, v, &X[target_index]);
      }
    }
    */
  } else {
    for (int i = num_fixed_pts; i < num_points; i++) {
      copy_to_target[i] = i;
      X[i] = copy_mesh->X[i];
      double u, v;
      int icode = face->invEvalPoint(X[i], &u, &v);
      if (icode) {
        fprintf(stderr,
                "TMRFaceMesh Error: Inverse point evaluation "
                "failed with code %d\n",
                icode);
      }

      // Set the target index
      pts[2 * i] = u;
      pts[2 * i + 1] = v;
      face->evalPoint(u, v, &X[i]);
    }

    // Copy over the quadrilaterals
    quads = new int[4 * num_quads];
    for (int i = 0; i < 4 * num_quads; i++) {
      quads[i] = copy_to_target[copy_mesh->quads[i]];
    }

    // Flip the quadrilateral surface orientation, depending on the
    // relative orientations of the copied face
    if (orient < 0) {
      for (int i = 0; i < num_quads; i++) {
        int tmp = quads[4 * i + 1];
        quads[4 * i + 1] = quads[4 * i + 3];
        quads[4 * i + 3] = tmp;
      }
    }
  }

  return 0;
}

/*
  Create a structured mesh
*/
void TMRFaceMesh::createStructuredMesh(TMRMeshOptions options,
                                       const double *params) {
  // Use a straightforward interpolation technique to obtain the
  // structured parametric locations in terms of the boundary
  // point parametric locations. We do not perform checks here
  // since we already know that the surface has four edges and the
  // nodes on those edges can be used for a structured mesh

  // Get the first edge loop and the edges in the loop
  TMREdgeLoop *loop;
  face->getEdgeLoop(0, &loop);

  // Get the edges associated with the edge loop
  TMREdge **edges;
  loop->getEdgeLoop(NULL, &edges, NULL);

  // Get the number of nodes for the x/y edges
  int nx = 0, ny = 0;
  TMREdgeMesh *mesh;

  // Get the face orientation
  int face_orient = face->getOrientation();
  if (face_orient > 0) {
    edges[0]->getMesh(&mesh);
    mesh->getMeshPoints(&nx, NULL, NULL);
    edges[1]->getMesh(&mesh);
    mesh->getMeshPoints(&ny, NULL, NULL);
  } else {
    edges[0]->getMesh(&mesh);
    mesh->getMeshPoints(&ny, NULL, NULL);
    edges[1]->getMesh(&mesh);
    mesh->getMeshPoints(&nx, NULL, NULL);
  }

  // Compute the total number of points
  num_points = nx * ny;
  num_quads = (nx - 1) * (ny - 1);

  // Create the connectivity information
  quads = new int[4 * num_quads];

  int *q = quads;
  for (int j = 0; j < ny - 1; j++) {
    for (int i = 0; i < nx - 1; i++) {
      // Compute the connectivity of the mesh
      q[0] = get_structured_index(nx, ny, i, j);
      q[1] = get_structured_index(nx, ny, i + 1, j);
      q[2] = get_structured_index(nx, ny, i + 1, j + 1);
      q[3] = get_structured_index(nx, ny, i, j + 1);
      q += 4;
    }
  }

  // Now set the parametric locations on the interior
  pts = new double[2 * num_points];

  // Copy the points from around the boundaries
  for (int i = 0; i < num_fixed_pts; i++) {
    pts[2 * i] = params[2 * i];
    pts[2 * i + 1] = params[2 * i + 1];
  }

  // Use a transfinite interpolation to determine the parametric
  // points where the interior nodes should be placed.
  const int c1idx = get_structured_index(nx, ny, 0, 0);
  const int c2idx = get_structured_index(nx, ny, nx - 1, 0);
  const int c3idx = get_structured_index(nx, ny, nx - 1, ny - 1);
  const int c4idx = get_structured_index(nx, ny, 0, ny - 1);

  for (int j = 1; j < ny - 1; j++) {
    double v = 1.0 * j / (ny - 1);
    for (int i = 1; i < nx - 1; i++) {
      double u = 1.0 * i / (nx - 1);

      // Compute the weights on the corners
      double c1 = (1.0 - u) * (1.0 - v);
      double c2 = u * (1.0 - v);
      double c3 = u * v;
      double c4 = (1.0 - u) * v;

      // Compute the weights on the curves
      double w1 = (1.0 - v);
      double w2 = u;
      double w3 = v;
      double w4 = (1.0 - u);

      // New parametric point
      int p = get_structured_index(nx, ny, i, j);

      // Boundary points that we're interpolating from
      int p1 = get_structured_index(nx, ny, i, 0);
      int p2 = get_structured_index(nx, ny, nx - 1, j);
      int p3 = get_structured_index(nx, ny, i, ny - 1);
      int p4 = get_structured_index(nx, ny, 0, j);

      // Evaluate the parametric points based on the transfinite
      // interpolation
      pts[2 * p] = ((w1 * pts[2 * p1] + w2 * pts[2 * p2] + w3 * pts[2 * p3] +
                     w4 * pts[2 * p4]) -
                    (c1 * pts[2 * c1idx] + c2 * pts[2 * c2idx] +
                     c3 * pts[2 * c3idx] + c4 * pts[2 * c4idx]));

      pts[2 * p + 1] = ((w1 * pts[2 * p1 + 1] + w2 * pts[2 * p2 + 1] +
                         w3 * pts[2 * p3 + 1] + w4 * pts[2 * p4 + 1]) -
                        (c1 * pts[2 * c1idx + 1] + c2 * pts[2 * c2idx + 1] +
                         c3 * pts[2 * c3idx + 1] + c4 * pts[2 * c4idx + 1]));
    }
  }

  // Allocate and evaluate the new physical point locations
  X = new TMRPoint[num_points];
  for (int i = 0; i < num_points; i++) {
    face->evalPoint(pts[2 * i], pts[2 * i + 1], &X[i]);
  }

  if (num_quads > 0) {
    // Smooth the copied mesh on the new surface
    int *pts_to_quad_ptr;
    int *pts_to_quads;
    TMR_ComputeNodeToElems(num_points, num_quads, 4, quads, &pts_to_quad_ptr,
                           &pts_to_quads);

    // Smooth the mesh using a local optimization of node locations
    TMR_QuadSmoothing(options.num_smoothing_steps, num_fixed_pts, num_points,
                      pts_to_quad_ptr, pts_to_quads, num_quads, quads, pts, X,
                      face);

    // Free the connectivity information
    delete[] pts_to_quad_ptr;
    delete[] pts_to_quads;
  }
}

/*
  Create the unstructured mesh using the Quad-Blossom algorithm
*/
void TMRFaceMesh::createUnstructuredMesh(
    TMRMeshOptions options, TMRElementFeatureSize *fs,
    TMRFaceMeshType mesh_type, const int total_num_pts, const int nholes,
    const double *params, const int nsegs, const int *segments,
    const int num_degen, const int *degen, int *npts, double **param_pts,
    TMRPoint **Xpts, int *nquads, int **mesh_quads, int *ntris,
    int **mesh_tris) {
  // Set defaults
  *npts = 0;
  *param_pts = NULL;
  *Xpts = NULL;
  *nquads = 0;
  *mesh_quads = NULL;
  *ntris = 0;
  *mesh_tris = NULL;

  // Create the triangularization class
  TMRTriangularize *tri = new TMRTriangularize(total_num_pts + nholes, params,
                                               nholes, nsegs, segments, face);
  tri->incref();

  if (options.write_init_domain_triangle) {
    char filename[256];
    sprintf(filename, "init_domain_triangle%d.vtk", face->getEntityId());
    tri->writeToVTK(filename);
  }

  // Create the mesh using the frontal algorithm
  tri->frontal(options, fs);

  // Free the degenerate triangles and reorder the mesh
  if (num_degen > 0) {
    tri->removeDegenerateEdges(num_degen, degen);
  }

  if (options.write_pre_smooth_triangle) {
    char filename[256];
    sprintf(filename, "pre_smooth_triangle%d.vtk", face->getEntityId());
    tri->writeToVTK(filename);
  }

  // Extract the triangularization
  tri->getMesh(npts, ntris, mesh_tris, param_pts, Xpts);
  tri->decref();

  if (*ntris == 0) {
    fprintf(stderr, "TMRFaceMesh Warning: No triangles for mesh id %d\n",
            face->getEntityId());
  } else {  // *ntris > 0
    // Compute the triangle edges and neighbors in the dual mesh
    int num_tri_edges;
    int *tri_edges, *tri_neighbors, *dual_edges;
    int *node_to_tri_ptr, *node_to_tris;
    TMR_ComputePlanarTriEdges(*npts, *ntris, *mesh_tris, &num_tri_edges,
                              &tri_edges, &tri_neighbors, &dual_edges,
                              &node_to_tri_ptr, &node_to_tris);

    // Smooth the resulting triangular mesh
    if (options.tri_smoothing_type == TMRMeshOptions::TMR_LAPLACIAN) {
      TMR_LaplacianSmoothing(options.num_smoothing_steps, num_fixed_pts,
                             num_tri_edges, tri_edges, *npts, *param_pts, *Xpts,
                             face);
    } else {
      double alpha = 0.1;
      TMR_SpringSmoothing(options.num_smoothing_steps, alpha, num_fixed_pts,
                          num_tri_edges, tri_edges, *npts, *param_pts, *Xpts,
                          face);
    }

    if (options.write_post_smooth_triangle) {
      char filename[256];
      sprintf(filename, "post_smooth_triangle%d.vtk", face->getEntityId());
      writeTrisToVTK(filename, *ntris, *mesh_tris);
    }

    if (mesh_type == TMR_UNSTRUCTURED) {
      // Recombine the mesh into a quadrilateral mesh
      if (*ntris % 2 == 0) {
        recombine(*ntris, *mesh_tris, tri_neighbors, node_to_tri_ptr,
                  node_to_tris, num_tri_edges, dual_edges, nquads, mesh_quads,
                  options);
      } else {
        fprintf(stderr,
                "TMRFaceMesh Error: Odd number of triangles, "
                "cannot perform recombination\n");
      }

      // Simplify the new quadrilateral mesh by removing
      // points/quads with poor quality/connectivity
      for (int k = 0; k < 5; k++) {
        simplifyQuads(0);
      }

      // Free the triangular mesh data
      delete[] tri_edges;
      delete[] tri_neighbors;
      delete[] dual_edges;
      delete[] node_to_tri_ptr;
      delete[] node_to_tris;

      if (options.write_pre_smooth_quad) {
        char filename[256];
        sprintf(filename, "pre_smooth_quad%d.vtk", face->getEntityId());
        writeToVTK(filename);
      }

      // Write out the dual of the final quadrilateral mesh
      if (options.write_quad_dual) {
        int num_quad_edges;
        int *quad_edges;
        int *quad_neighbors, *quad_dual;
        TMR_ComputePlanarQuadEdges(*npts, *nquads, *mesh_quads, &num_quad_edges,
                                   &quad_edges, &quad_neighbors, &quad_dual);

        char filename[256];
        sprintf(filename, "quad_dual%d.vtk", face->getEntityId());
        writeDualToVTK(filename, 4, *nquads, *mesh_quads, num_quad_edges,
                       quad_dual, *Xpts);

        delete[] quad_edges;
        delete[] quad_neighbors;
        delete[] quad_dual;
      }
    } else {
      // Free the triangular mesh data
      delete[] tri_edges;
      delete[] tri_neighbors;
      delete[] dual_edges;
      delete[] node_to_tri_ptr;
      delete[] node_to_tris;
    }
  }
}

/*
  Retrieve the mesh points and parametric locations
*/
void TMRFaceMesh::getMeshPoints(int *_npts, const double **_pts,
                                TMRPoint **_X) {
  if (_npts) {
    *_npts = num_points;
  }
  if (_pts) {
    *_pts = pts;
  }
  if (_X) {
    *_X = X;
  }
}

/*
  Get the source to target mapping
*/
void TMRFaceMesh::getSourceToTargetMapping(const int **_source_to_target) {
  *_source_to_target = source_to_target;
}

/*
  Get the face index along the absolute orientation of the edge
*/
int TMRFaceMesh::getFaceIndexFromEdge(TMREdge *e, int idx) {
  int index = -1;

  if (e && idx >= 0 && mesh_type != TMR_NO_MESH) {
    // Get the face orientation
    int face_orient = face->getOrientation();

    // Keep track of the number of closed loop cycles in the domain
    int nloops = face->getNumEdgeLoops();

    // Keep track of the index
    int pt = 0;

    for (int k = 0; k < nloops; k++) {
      // Get the curve information for this loop segment
      TMREdgeLoop *loop;
      face->getEdgeLoop(k, &loop);

      // Extract the edges and orientations from the loop
      int nedges;
      TMREdge **edges;
      const int *edge_orient;
      loop->getEdgeLoop(&nedges, &edges, &edge_orient);

      // Set the starting edge based on the face orientation
      int edge_index = nedges - 1;
      if (face_orient > 0) {
        edge_index = 0;
      }

      // Keep track of the starting loop point
      int start_loop_pt = pt;

      for (int i = 0; i < nedges; i++, edge_index += face_orient) {
        // Retrieve the underlying edge mesh
        TMREdgeMesh *mesh = NULL;
        TMREdge *edge = edges[edge_index], *copy = NULL;
        edge->getCopySource(&copy);
        edge->getMesh(&mesh);

        // Get the mesh points corresponding to this curve
        int npts;
        mesh->getMeshPoints(&npts, NULL, NULL);

        // Is this the last edge in this loop?
        int last_edge = (i == nedges - 1);

        if (e == edge || e == copy) {
          if (idx < npts) {
            // Check the orientation of the edge. Note that idx is the
            // absolute index along the edge.
            if (face_orient * edge_orient[edge_index] > 0) {
              if (last_edge && idx == npts - 1) {
                index = start_loop_pt;
              } else {
                index = pt + idx;
              }
            } else {
              if (last_edge && idx == 0) {
                index = start_loop_pt;
              } else {
                index = pt + npts - 1 - idx;
              }
            }
          }
          break;
        } else {
          if (edge->isDegenerate()) {
            pt++;
          } else {
            pt += npts - 1;
          }
        }
      }

      if (index != -1) {
        break;
      }
    }
  }

  return index;
}

/*
  Get the index of the structured face node. Note that this only works
  if the face was meshed using the structured mesh code, otherwise
  it returns a negative index.

  input:
  e1, e2:      perpendicular edges
  idx1, idx2:  the indices along each edge
*/
int TMRFaceMesh::getStructuredFaceIndex(TMREdge *e1, int idx1, TMREdge *e2,
                                        int idx2) {
  int index = -1;

  // First ensure that the mesh is actually structured
  if (mesh_type == TMR_STRUCTURED) {
    // Get the face orientation
    int face_orient = face->getOrientation();

    // Get the first edge loop and the edges in the loop. There
    // is only one loop since this is a structured mesh
    TMREdgeLoop *loop;
    face->getEdgeLoop(0, &loop);

    // Get the edges associated with the edge loop
    TMREdge **edges;
    const int *edge_orient;
    loop->getEdgeLoop(NULL, &edges, &edge_orient);

    // Get the number of nodes for the x/y edges
    int nx = 0, ny = 0;
    TMREdgeMesh *mesh;

    // Check for the first edge
    int i = -1, j = -1;

    if (face_orient > 0) {
      edges[0]->getMesh(&mesh);
      mesh->getMeshPoints(&nx, NULL, NULL);
      edges[1]->getMesh(&mesh);
      mesh->getMeshPoints(&ny, NULL, NULL);

      if ((e1 == edges[0] || e1 == edges[2]) &&
          (e2 == edges[1] || e2 == edges[3])) {
        // Check for the index along e1
        if (e1 == edges[0]) {
          if (edge_orient[0] > 0) {
            i = idx1;
          } else {
            i = nx - 1 - idx1;
          }
        } else {  // e1 == edges[2] -- reverse the index
          if (edge_orient[2] > 0) {
            i = nx - 1 - idx1;
          } else {
            i = idx1;
          }
        }

        // Check for the orientation along e2
        if (e2 == edges[1]) {
          if (edge_orient[1] > 0) {
            j = idx2;
          } else {
            j = ny - 1 - idx2;
          }
        } else {  // e2 == edges[3]
          if (edge_orient[3] > 0) {
            j = ny - 1 - idx2;
          } else {
            j = idx2;
          }
        }
      } else if ((e2 == edges[0] || e2 == edges[2]) &&
                 (e1 == edges[1] || e1 == edges[3])) {
        // Check for the index along e1
        if (e2 == edges[0]) {
          if (edge_orient[0] > 0) {
            i = idx2;
          } else {
            i = nx - 1 - idx2;
          }
        } else {  // e2 == edges[2] -- reverse the index
          if (edge_orient[2] > 0) {
            i = nx - 1 - idx2;
          } else {
            i = idx2;
          }
        }

        // Check for the orientation along e2
        if (e1 == edges[1]) {
          if (edge_orient[1] > 0) {
            j = idx1;
          } else {
            j = ny - 1 - idx1;
          }
        } else {  // e1 == edges[3]
          if (edge_orient[3] > 0) {
            j = ny - 1 - idx1;
          } else {
            j = idx1;
          }
        }
      }
    } else {
      edges[0]->getMesh(&mesh);
      mesh->getMeshPoints(&ny, NULL, NULL);
      edges[1]->getMesh(&mesh);
      mesh->getMeshPoints(&nx, NULL, NULL);

      if ((e1 == edges[0] || e1 == edges[2]) &&
          (e2 == edges[1] || e2 == edges[3])) {
        // Check for the index along e1
        if (e1 == edges[0]) {
          if (edge_orient[0] > 0) {
            j = idx1;
          } else {
            j = ny - 1 - idx1;
          }
        } else {  // e1 == edges[2] -- reverse the index
          if (edge_orient[2] > 0) {
            j = ny - 1 - idx1;
          } else {
            j = idx1;
          }
        }

        // Check for the orientation along e2
        if (e2 == edges[1]) {
          if (edge_orient[1] > 0) {
            i = idx2;
          } else {
            i = nx - 1 - idx2;
          }
        } else {  // e2 == edges[3]
          if (edge_orient[3] > 0) {
            i = nx - 1 - idx2;
          } else {
            i = idx2;
          }
        }
      } else if ((e2 == edges[0] || e2 == edges[2]) &&
                 (e1 == edges[1] || e1 == edges[3])) {
        // Check for the index along e1
        if (e2 == edges[0]) {
          if (edge_orient[0] > 0) {
            j = idx2;
          } else {
            j = ny - 1 - idx2;
          }
        } else {  // e2 == edges[2] -- reverse the index
          if (edge_orient[2] > 0) {
            j = ny - 1 - idx2;
          } else {
            j = idx2;
          }
        }

        // Check for the orientation along e2
        if (e1 == edges[1]) {
          if (edge_orient[1] > 0) {
            i = idx1;
          } else {
            i = nx - 1 - idx1;
          }
        } else {  // e1 == edges[3]
          if (edge_orient[3] > 0) {
            i = nx - 1 - idx1;
          } else {
            i = idx1;
          }
        }
      }
    }

    if ((i >= 0 && i < nx) && (j >= 0 && j < ny)) {
      index = get_structured_index(nx, ny, i, j);
    }
  }

  return index;
}

/*
  Set the node numbers internally
*/
int TMRFaceMesh::setNodeNums(int *num) {
  if (!vars) {
    vars = new int[num_points];

    // Set the starting index
    int start = *num;

    // Check if the copy-to-target index has been set
    if (copy_to_target) {
      TMRFace *copy;
      face->getCopySource(NULL, &copy);

      TMRFaceMesh *copy_mesh;
      copy->getMesh(&copy_mesh);

      // Ensure that the variables are set on the copied mesh
      copy_mesh->setNodeNums(num);

      // Set the copy variable nubmers from the target
      for (int i = 0; i < num_points; i++) {
        vars[copy_to_target[i]] = copy_mesh->vars[i];
      }
    } else {
      // Get the face orientation
      int face_orient = face->getOrientation();

      // Retrieve the boundary node numbers from the surface loops
      int pt = 0;
      for (int k = 0; k < face->getNumEdgeLoops(); k++) {
        // Get the curve information for this loop segment
        TMREdgeLoop *loop;
        face->getEdgeLoop(k, &loop);

        int nedges;
        TMREdge **edges;
        const int *edge_orient;
        loop->getEdgeLoop(&nedges, &edges, &edge_orient);

        int edge_index = nedges - 1;
        if (face_orient > 0) {
          edge_index = 0;
        }

        for (int i = 0; i < nedges; i++, edge_index += face_orient) {
          // Retrieve the underlying curve mesh
          TMREdgeMesh *mesh = NULL;
          edges[edge_index]->getMesh(&mesh);

          // Retrieve the variable numbers for this loop
          const int *edge_vars;
          int npts = mesh->getNodeNums(&edge_vars);

          if (edges[edge_index]->isDegenerate()) {
            vars[pt] = edge_vars[0];
          } else {
            // Get the orientation of the edge
            int orientation = face_orient * edge_orient[edge_index];

            int index = npts - 1;
            if (orientation > 0) {
              index = 0;
            }

            for (int j = 0; j < npts - 1; j++, index += orientation, pt++) {
              // Find the point on the curve
              vars[pt] = edge_vars[index];
            }
          }
        }
      }

      // Now order the variables as they arrive
      for (; pt < num_points; pt++) {
        vars[pt] = *num;
        (*num)++;
      }
    }

    // Return the number of points that have been allocated
    return *num - start;
  }

  return 0;
}

/*
  Retrieve the mapping between the local connectivity and the global
  node numbers
*/
int TMRFaceMesh::getNodeNums(const int **_vars) {
  if (_vars) {
    *_vars = vars;
  }
  return num_points;
}

/*
  Get the number of fixed points that are not ordered by this surface
  mesh
*/
int TMRFaceMesh::getNumFixedPoints() { return num_fixed_pts; }

/*
  Get the local quad connectivity
*/
int TMRFaceMesh::getQuadConnectivity(const int **_quads) {
  if (_quads) {
    *_quads = quads;
  }
  return num_quads;
}

/*
  Get the local triangle connectivity
*/
int TMRFaceMesh::getTriConnectivity(const int **_tris) {
  if (_tris) {
    *_tris = tris;
  }
  return num_tris;
}

/*
  Get the quadrilateral elemnet obtained by combining the triangles
  t1 and t2 together

        . --- .
      / |   /
    /   | /
  . --- .
*/
int TMRFaceMesh::getRecombinedQuad(const int triangles[],
                                   const int trineighbors[], int t1, int t2,
                                   int quad[]) {
  int fail = 0;

  // Find the common edge between the two tirangles
  const int shift[3] = {2, 0, 1};
  int e1 = 0, e2 = 0;
  for (; e1 < 3; e1++) {
    if (trineighbors[3 * t1 + shift[e1]] == t2) {
      break;
    }
  }
  for (; e2 < 3; e2++) {
    if (trineighbors[3 * t2 + shift[e2]] == t1) {
      break;
    }
  }

  if (e1 >= 3 || e2 >= 3) {
    fail = 1;
    return fail;
  }

  // Order the triangle, omitting the common edge
  for (int j = 0, i = 0; i < 3; i++, j++) {
    quad[j] = triangles[3 * t1 + i];
    if (i == e1) {
      j++;
    }
  }

  // Set the node contributed by the second triangle
  e2 += 2;
  if (e2 >= 3) {
    e2 -= 3;
  }
  quad[e1 + 1] = triangles[3 * t2 + e2];

  return fail;
}

/*
  Compute the quality of a quadrilateral element
*/
double TMRFaceMesh::computeQuadQuality(const int *quad, const TMRPoint *p) {
  // Compute the maximum of fabs(0.5*M_PI - alpha)
  double max_val = 0.0;

  for (int k = 0; k < 4; k++) {
    int prev = k - 1;
    if (k == 0) {
      prev = 3;
    }
    int next = k + 1;
    if (k == 3) {
      next = 0;
    }

    TMRPoint a;
    a.x = p[quad[k]].x - p[quad[prev]].x;
    a.y = p[quad[k]].y - p[quad[prev]].y;
    a.z = p[quad[k]].z - p[quad[prev]].z;

    TMRPoint b;
    b.x = p[quad[next]].x - p[quad[k]].x;
    b.y = p[quad[next]].y - p[quad[k]].y;
    b.z = p[quad[next]].z - p[quad[k]].z;

    // Compute the internal angle
    double beta = a.dot(b) / sqrt(a.dot(a) * b.dot(b));
    if (beta < -1.0) {
      beta = -1.0;
    }
    if (beta > 1.0) {
      beta = 1.0;
    }
    double alpha = M_PI - acos(beta);
    double val = fabs(0.5 * M_PI - alpha);
    if (val > max_val) {
      max_val = val;
    }
  }

  // Compute the quality
  double eta = 1.0 - (2.0 / M_PI) * max_val;
  if (eta < 0.0) {
    eta = 0.0;
  }

  return eta;
}

/*
  Compute the recombined quality
*/
double TMRFaceMesh::computeRecombinedQuality(const int triangles[],
                                             const int tri_neighbors[], int t1,
                                             int t2, const TMRPoint *p) {
  // Find the combined quadrilateral from the two given triangles
  int quad[4];
  int fail = getRecombinedQuad(triangles, tri_neighbors, t1, t2, quad);
  if (fail) {
    return 0.0;
  }

  return computeQuadQuality(quad, p);
}

/*
  Compute the quality of a quadrilateral element
*/
double TMRFaceMesh::computeTriQuality(const int *tri, const TMRPoint *p) {
  // Compute the maximum of fabs(M_PI/3 - alpha)
  double max_val = 0.0;

  for (int k = 0; k < 3; k++) {
    int prev = k - 1;
    if (prev < 0) {
      prev = 2;
    }
    int next = k + 1;
    if (next > 2) {
      next = 0;
    }

    TMRPoint a;
    a.x = p[tri[k]].x - p[tri[prev]].x;
    a.y = p[tri[k]].y - p[tri[prev]].y;
    a.z = p[tri[k]].z - p[tri[prev]].z;

    TMRPoint b;
    b.x = p[tri[next]].x - p[tri[k]].x;
    b.y = p[tri[next]].y - p[tri[k]].y;
    b.z = p[tri[next]].z - p[tri[k]].z;

    // Compute the internal angle
    double beta = a.dot(b) / sqrt(a.dot(a) * b.dot(b));
    if (beta < -1.0) {
      beta = -1.0;
    }
    if (beta > 1.0) {
      beta = 1.0;
    }
    double alpha = M_PI - acos(beta);
    double val = fabs(M_PI / 3.0 - alpha);
    if (val > max_val) {
      max_val = val;
    }
  }

  // Compute the quality
  double eta = 1.0 - (3.0 / M_PI) * max_val;
  if (eta < 0.0) {
    eta = 0.0;
  }

  return eta;
}

/*
  Recombine the triangulation into a quadrilateral mesh
*/
void TMRFaceMesh::recombine(int ntris, const int triangles[],
                            const int tri_neighbors[],
                            const int node_to_tri_ptr[],
                            const int node_to_tris[], int num_edges,
                            const int dual_edges[], int *_num_quads,
                            int **_new_quads, TMRMeshOptions options) {
  // Allocate the reduced graph weights
  double *weights = new double[num_edges];
  int *graph_edges = new int[2 * num_edges];

  // Compute the weight associated with each edge by combputing the
  // recombined quality
  const double eps = 0.01;

  int edge_num = 0;
  for (int i = 0; i < num_edges; i++) {
    int t1 = dual_edges[2 * i];
    int t2 = dual_edges[2 * i + 1];

    if (t1 >= 0 && t2 >= 0) {
      // Compute the weight for this recombination
      double quality =
          computeRecombinedQuality(triangles, tri_neighbors, t1, t2, X);

      double weight = (1.0 - quality) * (1.0 + 1.0 / (quality + eps));
      graph_edges[2 * edge_num] = t1;
      graph_edges[2 * edge_num + 1] = t2;
      weights[edge_num] = weight;
      edge_num++;
    }
  }

  // Keep track of the number of real edges in the graph
  int num_real_edges = edge_num;

  // Determine the extra edges that connect along the boundaries.
  // between unique triangles that are not already adjacent to one
  // another.
  for (int i = 0; i < ntris; i++) {
    if (tri_neighbors[3 * i] < 0 || tri_neighbors[3 * i + 1] < 0 ||
        tri_neighbors[3 * i + 2] < 0) {
      // Loop over the edges of the triangle
      for (int j = 0; j < 3; j++) {
        // We have a boundary triangle
        if (tri_neighbors[3 * i + j] < 0) {
          // The leading node that we are searching for
          int ij = tri_edge_nodes[j][1];
          int node = triangles[3 * i + ij];

          // Search the triangles adjacent to this node
          const int kpend = node_to_tri_ptr[node + 1];
          for (int kp = node_to_tri_ptr[node]; kp < kpend; kp++) {
            int k = node_to_tris[kp];

            // If this is the same triangle, continue
            if (i == k) {
              continue;
            }

            // If the triangles are actually also neighbors, continue
            if (tri_neighbors[3 * k] == i || tri_neighbors[3 * k + 1] == i ||
                tri_neighbors[3 * k + 2] == i) {
              continue;
            }

            // Find the local node number shared between triangles k
            // and triangle i
            int kj = 0;
            if (triangles[3 * k + 1] == node) {
              kj = 1;
            } else if (triangles[3 * k + 2] == node) {
              kj = 2;
            }

            // Track from the node to the possible edges
            if (tri_neighbors[3 * k + tri_node_edges[kj][0]] < 0 ||
                tri_neighbors[3 * k + tri_node_edges[kj][1]] < 0) {
              graph_edges[2 * edge_num] = i;
              graph_edges[2 * edge_num + 1] = k;
              weights[edge_num] = 0.5 / eps;
              edge_num++;
            }
          }
        }
      }
    }
  }

  // Set the number of edges within the modified dual mesh
  int num_dual_edges = edge_num;

  // Write the dual mesh to a file
  if (options.write_dual_recombine) {
    char filename[256];
    sprintf(filename, "dual_recombine%d.vtk", face->getEntityId());
    writeDualToVTK(filename, 3, ntris, triangles, num_dual_edges, graph_edges,
                   X);
  }

  // Perform the perfect matching
  int *match = new int[ntris / 2];
  int num_match =
      TMR_PerfectMatchGraph(ntris, num_dual_edges, graph_edges, weights, match);
  delete[] weights;

  // The quads formed from the original triangles
  int num_quads_from_tris = 0;

  // The new number of quadrilaterals - after every new quad is added
  int num_new_quads = 0;

  // New points array
  int num_new_points = 0;
  double *new_pts = NULL;
  TMRPoint *new_X = NULL;

  for (int i = 0; i < num_match; i++) {
    if (match[i] < num_real_edges) {
      num_quads_from_tris++;
      num_new_quads++;
    } else {
      // We'll add two extra quads for each extra edge
      num_new_quads += 2;
      num_new_points++;
    }
  }

  // We'll be adding new points so allocate new arrays to handle the
  // new points that will be computed...
  if (num_new_points > 0) {
    num_new_points += num_points;
    new_pts = new double[2 * num_new_points];
    new_X = new TMRPoint[num_new_points];
    memcpy(new_pts, pts, 2 * num_points * sizeof(double));
    memcpy(new_X, X, num_points * sizeof(TMRPoint));
  }

  // Set the number of quadrilateral elements created in the mesh and
  // record the new quadrilateral element connectivity
  int *new_quads = new int[4 * num_new_quads];

  // Recombine the triangles into quadrilateral elements
  num_new_quads = 0;
  for (int i = 0; i < num_match; i++) {
    if (match[i] < num_real_edges) {
      int t1 = graph_edges[2 * match[i]];
      int t2 = graph_edges[2 * match[i] + 1];

      int fail = getRecombinedQuad(triangles, tri_neighbors, t1, t2,
                                   &new_quads[4 * num_new_quads]);
      num_new_quads++;
      if (fail) {
        fprintf(stderr,
                "TMRFaceMesh Error: Quad %d from triangles "
                "%d and %d failed\n",
                i, t1, t2);
      }
    }
  }

  // Set the number of new points
  num_new_points = num_points;

  std::map<int, int> new_point_nums;

  // Add the triangles from edges along the boundary. There should
  // only be a handful of these guys...
  for (int i = 0; i < num_match; i++) {
    if (match[i] >= num_real_edges) {
      // These triangles can only share one common node - the node
      // that we'll now duplicate and move to the interior of the
      // domain. Find the local node index for this shared node in
      // each triangle.
      int t1 = graph_edges[2 * match[i]];
      int t2 = graph_edges[2 * match[i] + 1];

      // j1, j2 are the local nodal indices of the shared node between
      // triangles t1 and t2
      int j1 = 0, j2 = 0;
      for (int flag = 0; !flag && j1 < 3; j1++) {
        for (j2 = 0; j2 < 3; j2++) {
          if (triangles[3 * t1 + j1] == triangles[3 * t2 + j2]) {
            flag = 1;
            break;
          }
        }
        if (flag) {
          break;
        }
      }

      int boundary_pt = triangles[3 * t1 + j1];

      // Go through all previous quadrilaterals and adjust the
      // ordering to reflect the duplicated node location
      for (int k = 0; k < 4 * num_new_quads; k++) {
        if (new_quads[k] == boundary_pt) {
          new_quads[k] = num_new_points;
        }
      }

      // Add the first triangle t1 - this triangle is gauranteed to
      // come first when circling the boundary in the CCW direction.
      int n1 = triangles[3 * t1 + ((j1 + 1) % 3)];
      int n2 = triangles[3 * t1 + ((j1 + 2) % 3)];
      if (new_point_nums.count(n1)) {
        n1 = new_point_nums[n1];
      }
      if (new_point_nums.count(n2)) {
        n2 = new_point_nums[n2];
      }
      new_quads[4 * num_new_quads] = n1;
      new_quads[4 * num_new_quads + 1] = n2;
      new_quads[4 * num_new_quads + 2] = boundary_pt;
      new_quads[4 * num_new_quads + 3] = num_new_points;
      num_new_quads++;

      // Add the connectivity from the second triangle t2. This
      // triangle will always come second when circling the boundary.
      n1 = triangles[3 * t2 + ((j2 + 1) % 3)];
      n2 = triangles[3 * t2 + ((j2 + 2) % 3)];
      if (new_point_nums.count(n1)) {
        n1 = new_point_nums[n1];
      }
      if (new_point_nums.count(n2)) {
        n2 = new_point_nums[n2];
      }
      new_quads[4 * num_new_quads] = n1;
      new_quads[4 * num_new_quads + 1] = n2;
      new_quads[4 * num_new_quads + 2] = num_new_points;
      new_quads[4 * num_new_quads + 3] = boundary_pt;
      num_new_quads++;

      // Add the new boundary points to the map
      new_point_nums[boundary_pt] = num_new_points;

      // Compute the new parameter location by taking the average of
      // the centroid locations for each triangle
      int count = 1;
      new_pts[2 * num_new_points] = pts[2 * boundary_pt];
      new_pts[2 * num_new_points + 1] = pts[2 * boundary_pt + 1];

      int kpend = node_to_tri_ptr[boundary_pt + 1];
      for (int kp = node_to_tri_ptr[boundary_pt]; kp < kpend; kp++) {
        int k = node_to_tris[kp];
        if (k != t1 && k != t2) {
          new_pts[2 * num_new_points] +=
              (pts[2 * triangles[3 * k]] + pts[2 * triangles[3 * k + 1]] +
               pts[2 * triangles[3 * k + 2]]) /
              3.0;
          new_pts[2 * num_new_points + 1] +=
              (pts[2 * triangles[3 * k] + 1] +
               pts[2 * triangles[3 * k + 1] + 1] +
               pts[2 * triangles[3 * k + 2] + 1]) /
              3.0;
          count++;
        }
      }
      if (count > 1) {
        new_pts[2 * num_new_points] = new_pts[2 * num_new_points] / count;
        new_pts[2 * num_new_points + 1] =
            new_pts[2 * num_new_points + 1] / count;
      }
      face->evalPoint(new_pts[2 * num_new_points],
                      new_pts[2 * num_new_points + 1], &new_X[num_new_points]);

      // Increment the number of new points
      num_new_points++;
    }
  }

  delete[] graph_edges;
  delete[] match;

  if (new_pts) {
    num_points = num_new_points;
    delete[] pts;
    delete[] X;
    pts = new_pts;
    X = new_X;
  }

  // Set the quads/output
  *_num_quads = num_new_quads;
  *_new_quads = new_quads;
}

/*
  This code performs several topological improvements to try and avoid
  poor mesh quality.

  First, the code attempts to remove topological triangles on
  boundaries where adjacent boundary edges are included in the same
  quad. Adjacent boundary edges in the same quadrilaterla are not
  always bad (for instance when a quad is at a corner), so we check
  whether to adjust the connectivity based on the angle between the
  two edges. Note that this will only work if the second quadrilateral
  does not have an edge on the boundary and requires that the mesh be
  smoothed afterwards.

           x                            x
         /   \                        / / \
        /     \                     /  /   \
       x       x                   x  /     x
     /   \     /                  /   /     /
    /     \   /         ===>     /   /     /
   /       \ /                  /   /     /
  x -- x -- x                  x -- x -- x

  Next, the code identifies and removes adjacent quadrilaterals that
  look like this:

  x --- x             x ---- x
  | \   |             |      |
  |  x  |     ===>    |      |
  |    \|             |      |
  x --- x             x ---- x

  Finally, the code finds and removes quads that look like this:

  x ---- x ---- x         x ---- x ---- x
  |     / \     |         |      |      |
  |    /   \    |         |      |      |
  x - x     x - x   ===>  x ---- x ---- x
  |    \   /    |         |      |      |
  |     \ /     |         |      |      |
  x ---- x ---- x         x ---- x ---- x
*/
void TMRFaceMesh::simplifyQuads(int dummy_flag) {
  // Compute the node -> quad information
  int *ptr, *pts_to_quads;
  TMR_ComputeNodeToElems(num_points, num_quads, 4, quads, &ptr, &pts_to_quads);

  // First, remove the nodes that are only referred to twice, not on
  // the boundary
  int *new_pt_nums = new int[num_points];
  memset(new_pt_nums, 0, num_points * sizeof(int));

  for (int i = num_fixed_pts; i < num_points; i++) {
    if (ptr[i + 1] - ptr[i] == 2) {
      // Retrieve the pointers to the latest quadrilatral
      int q1 = pts_to_quads[ptr[i]];
      int q2 = pts_to_quads[ptr[i] + 1];
      int *quad1 = &quads[4 * q1];
      int *quad2 = &quads[4 * q2];

      // If any of the other points have been deleted, then
      // skip combining this quadrilateral. This call could
      // probably be repeated to remove this element, but this
      // algorithm cannot handle it.
      int skip_me = 0;
      for (int k = 0; k < 4; k++) {
        if (quad1[k] < 0 || quad2[k] < 0 || new_pt_nums[quad1[k]] < 0 ||
            new_pt_nums[quad2[k]] < 0) {
          skip_me = 1;
          break;
        }
      }

      if (!skip_me) {
        // Find the common node between quads
        int k1 = 0, k2 = 0;
        while (quad1[k1] != i) k1++;
        while (quad2[k2] != i) k2++;

        // Adjust k2 so that it points to the node that
        // will be inserted into the connectivity for the
        // first quadrilateral
        k2 += 2;
        if (k2 >= 4) {
          k2 -= 4;
        }
        quad1[k1] = quad2[k2];

        // Set the point
        quad2[0] = quad2[1] = quad2[2] = quad2[3] = -1;

        // Set the points that we're going to eliminate
        new_pt_nums[i] = -1;
      }
    }
  }

  for (int i = 0; i < num_quads; i++) {
    for (int j1 = 0; j1 < 2; j1++) {
      // Find the local node numbers that are on opposite
      // sides of the quadrilateral
      int j2 = j1 + 2;

      // Compute the points that are on opposite sides of the
      // quadrilateral
      int p1 = quads[4 * i + j1];
      int p2 = quads[4 * i + j2];

      int collapse = 0;
      if (p1 >= 0 && p2 >= 0 && p1 >= num_fixed_pts && p2 >= num_fixed_pts) {
        // If the quad is in the interior
        if ((ptr[p1 + 1] - ptr[p1] == 3 && ptr[p2 + 1] - ptr[p2] == 3)) {
          collapse = 1;
        }
        // If the quad is on the boundary, collapse it if the degree of the
        // nodes are one larger
        if ((quads[4 * i + ((j1 + 1) % 4)] < num_fixed_pts ||
             quads[4 * i + ((j2 + 1) % 4)] < num_fixed_pts) &&
            ((ptr[p1 + 1] - ptr[p1] == 3 && ptr[p2 + 1] - ptr[p2] <= 4) ||
             (ptr[p2 + 1] - ptr[p2] == 3 && ptr[p1 + 1] - ptr[p1] <= 4))) {
          collapse = 1;
        }
      }

      if (collapse) {
        // Check whether any of the quadrilaterals which touch either
        // p1 or p2 have negative indices or are on a boundary.
        int flag = 0;
        for (int kp = ptr[p1]; kp < ptr[p1 + 1]; kp++) {
          int q = pts_to_quads[kp];
          if (quads[4 * q] < 0 || new_pt_nums[quads[4 * q]] < 0 ||
              new_pt_nums[quads[4 * q + 1]] < 0 ||
              new_pt_nums[quads[4 * q + 2]] < 0 ||
              new_pt_nums[quads[4 * q + 3]] < 0) {
            flag = 1;
            break;
          }
        }
        if (flag) {
          break;
        }
        for (int kp = ptr[p2]; kp < ptr[p2 + 1]; kp++) {
          int q = pts_to_quads[kp];
          if (quads[4 * q] < 0 || new_pt_nums[quads[4 * q]] < 0 ||
              new_pt_nums[quads[4 * q + 1]] < 0 ||
              new_pt_nums[quads[4 * q + 2]] < 0 ||
              new_pt_nums[quads[4 * q + 3]] < 0) {
            flag = 1;
            break;
          }
        }
        if (flag) {
          break;
        }

        // Figure out which node to eliminate
        new_pt_nums[p1] = -1;

        // Set the quadrilaterals that refer to the node p1
        // now instead to refer to the node p2
        for (int kp = ptr[p1]; kp < ptr[p1 + 1]; kp++) {
          int q = pts_to_quads[kp];
          for (int k = 0; k < 4; k++) {
            if (quads[4 * q + k] == p1) {
              quads[4 * q + k] = p2;
            }
          }
        }

        // Now eliminate the quad
        quads[4 * i] = quads[4 * i + 1] = quads[4 * i + 2] = quads[4 * i + 3] =
            -1;

        // Quit this quadrilateral
        break;
      }
    }
  }

  // Free the pointer/quad pointer data
  delete[] ptr;
  delete[] pts_to_quads;

  // Remove the points/parameters that have been eliminated
  int pt_num = 0;
  for (; pt_num < num_fixed_pts; pt_num++) {
    new_pt_nums[pt_num] = pt_num;
  }
  for (int i = num_fixed_pts; i < num_points; i++) {
    if (new_pt_nums[i] >= 0) {
      if (i != pt_num) {
        // Copy over the points
        pts[2 * pt_num] = pts[2 * i];
        pts[2 * pt_num + 1] = pts[2 * i + 1];
        X[pt_num] = X[i];
      }
      new_pt_nums[i] = pt_num;
      pt_num++;
    }
  }

  // Set the new number of points
  num_points = pt_num;

  // Set the new connectivity, overwriting the old connectivity
  int quad_num = 0;
  for (int i = 0; i < num_quads; i++) {
    if (quads[4 * i] >= 0) {
      quads[4 * quad_num] = new_pt_nums[quads[4 * i]];
      quads[4 * quad_num + 1] = new_pt_nums[quads[4 * i + 1]];
      quads[4 * quad_num + 2] = new_pt_nums[quads[4 * i + 2]];
      quads[4 * quad_num + 3] = new_pt_nums[quads[4 * i + 3]];
      quad_num++;
    }
  }

  // Set the new number of quadrilaterals
  num_quads = quad_num;

  // Free the new point numbers
  delete[] new_pt_nums;
}

/*
  Add the quad quality
*/
void TMRFaceMesh::addMeshQuality(int nbins, int bins[]) {
  for (int i = 0; i < num_quads; i++) {
    double quality = computeQuadQuality(&quads[4 * i], X);

    int k = 0;
    for (; k < nbins; k++) {
      if (quality < 1.0 * (k + 1) / nbins) {
        break;
      }
    }
    if (k == nbins) {
      k = nbins - 1;
    }
    bins[k]++;
  }

  for (int i = 0; i < num_tris; i++) {
    double quality = computeTriQuality(&tris[3 * i], X);

    int k = 0;
    for (; k < nbins; k++) {
      if (quality < 1.0 * (k + 1) / nbins) {
        break;
      }
    }
    if (k == nbins) {
      k = nbins - 1;
    }
    bins[k]++;
  }
}

/*
  Print the quadrilateral quality
*/
void TMRFaceMesh::printMeshQuality() {
  const int nbins = 20;
  int total = 0;
  int bins[nbins];
  memset(bins, 0, nbins * sizeof(int));
  addMeshQuality(nbins, bins);

  for (int i = 0; i < nbins; i++) {
    total += bins[i];
  }

  printf("Quality   # elements   percentage\n");
  for (int k = 0; k < nbins; k++) {
    printf("< %.2f    %10d   %10.3f\n", 1.0 * (k + 1) / nbins, bins[k],
           100.0 * bins[k] / total);
  }
}

/*
  Print the triangle quality
*/
void TMRFaceMesh::printTriQuality(int ntris, const int triangles[]) {
  const int nbins = 20;
  int total = 0;
  int bins[nbins];
  memset(bins, 0, nbins * sizeof(int));
  for (int i = 0; i < ntris; i++) {
    double quality = computeTriQuality(&triangles[3 * i], X);

    int k = 0;
    for (; k < nbins; k++) {
      if (quality < 1.0 * (k + 1) / nbins) {
        break;
      }
    }
    if (k == nbins) {
      k = nbins - 1;
    }
    bins[k]++;
  }

  for (int i = 0; i < nbins; i++) {
    total += bins[i];
  }

  printf("Quality   # elements   percentage\n");
  for (int k = 0; k < nbins; k++) {
    printf("< %.2f    %10d   %10.3f\n", 1.0 * (k + 1) / nbins, bins[k],
           100.0 * bins[k] / total);
  }
}

/*
  Write the quadrilateral mesh to a VTK file
*/
void TMRFaceMesh::writeToVTK(const char *filename) {
  if (num_quads > 0) {
    FILE *fp = fopen(filename, "w");
    if (fp) {
      fprintf(fp, "# vtk DataFile Version 3.0\n");
      fprintf(fp, "vtk output\nASCII\n");
      fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

      // Write out the points
      fprintf(fp, "POINTS %d float\n", num_points);
      for (int k = 0; k < num_points; k++) {
        fprintf(fp, "%e %e %e\n", X[k].x, X[k].y, X[k].z);
      }

      // Write out the cell values
      fprintf(fp, "\nCELLS %d %d\n", num_quads, 5 * num_quads);
      for (int k = 0; k < num_quads; k++) {
        fprintf(fp, "4 %d %d %d %d\n", quads[4 * k], quads[4 * k + 1],
                quads[4 * k + 2], quads[4 * k + 3]);
      }

      // All quadrilaterals
      fprintf(fp, "\nCELL_TYPES %d\n", num_quads);
      for (int k = 0; k < num_quads; k++) {
        fprintf(fp, "%d\n", 9);
      }

      // Print out the rest as fields one-by-one
      fprintf(fp, "CELL_DATA %d\n", num_quads);
      fprintf(fp, "SCALARS quality float 1\n");
      fprintf(fp, "LOOKUP_TABLE default\n");
      for (int i = 0; i < num_quads; i++) {
        fprintf(fp, "%e\n", computeQuadQuality(&quads[4 * i], X));
      }

      fclose(fp);
    }
  } else if (num_tris > 0) {
    writeTrisToVTK(filename, num_tris, tris);
  }
}

/*
  Write out the segments to a VTK file
*/
void TMRFaceMesh::writeSegmentsToVTK(const char *filename, int npts,
                                     const double *params, int nsegs,
                                     const int segs[]) {
  FILE *fp = fopen(filename, "w");
  if (fp) {
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "vtk output\nASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    // Write out the points
    fprintf(fp, "POINTS %d float\n", npts);
    for (int k = 0; k < npts; k++) {
      fprintf(fp, "%e %e 0.0\n", params[2 * k], params[2 * k + 1]);
    }

    // Write out the cell connectivity
    fprintf(fp, "\nCELLS %d %d\n", nsegs, 3 * nsegs);
    for (int k = 0; k < nsegs; k++) {
      fprintf(fp, "2 %d %d\n", segs[2 * k], segs[2 * k + 1]);
    }

    // Write out the cell types
    fprintf(fp, "\nCELL_TYPES %d\n", nsegs);
    for (int k = 0; k < nsegs; k++) {
      fprintf(fp, "%d\n", 3);
    }
    fclose(fp);
  }
}

/*
  Write the output to a VTK file
*/
void TMRFaceMesh::writeTrisToVTK(const char *filename, int ntris,
                                 const int triangles[]) {
  FILE *fp = fopen(filename, "w");
  if (fp) {
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "vtk output\nASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    // Write out the points
    fprintf(fp, "POINTS %d float\n", num_points);
    for (int k = 0; k < num_points; k++) {
      fprintf(fp, "%e %e %e\n", X[k].x, X[k].y, X[k].z);
    }

    // Write out the cell values
    fprintf(fp, "\nCELLS %d %d\n", ntris, 4 * ntris);
    for (int k = 0; k < ntris; k++) {
      fprintf(fp, "3 %d %d %d\n", triangles[3 * k], triangles[3 * k + 1],
              triangles[3 * k + 2]);
    }

    // All quadrilaterals
    fprintf(fp, "\nCELL_TYPES %d\n", ntris);
    for (int k = 0; k < ntris; k++) {
      fprintf(fp, "%d\n", 5);
    }

    // Print out the rest as fields one-by-one
    fprintf(fp, "CELL_DATA %d\n", ntris);
    fprintf(fp, "SCALARS quality float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int k = 0; k < ntris; k++) {
      fprintf(fp, "%e\n", computeTriQuality(&triangles[3 * k], X));
    }

    fclose(fp);
  }
}

/*
  Write out the dual mesh for visualization
*/
void TMRFaceMesh::writeDualToVTK(const char *filename, int nodes_per_elem,
                                 int nelems, const int elems[],
                                 int num_dual_edges, const int dual_edges[],
                                 const TMRPoint *p) {
  FILE *fp = fopen(filename, "w");
  if (fp) {
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "vtk output\nASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    // Write out the points
    fprintf(fp, "POINTS %d float\n", nelems);
    if (nodes_per_elem == 3) {
      for (int k = 0; k < nelems; k++) {
        fprintf(fp, "%e %e %e\n",
                1.0 / 3.0 *
                    (p[elems[3 * k]].x + p[elems[3 * k + 1]].x +
                     p[elems[3 * k + 2]].x),
                1.0 / 3.0 *
                    (p[elems[3 * k]].y + p[elems[3 * k + 1]].y +
                     p[elems[3 * k + 2]].y),
                1.0 / 3.0 *
                    (p[elems[3 * k]].z + p[elems[3 * k + 1]].z +
                     p[elems[3 * k + 2]].z));
      }
    } else {  // nodes_per_elem == 4
      for (int k = 0; k < nelems; k++) {
        fprintf(fp, "%e %e %e\n",
                0.25 * (p[elems[4 * k]].x + p[elems[4 * k + 1]].x +
                        p[elems[4 * k + 2]].x + p[elems[4 * k + 3]].x),
                0.25 * (p[elems[4 * k]].y + p[elems[4 * k + 1]].y +
                        p[elems[4 * k + 2]].y + p[elems[4 * k + 3]].y),
                0.25 * (p[elems[4 * k]].z + p[elems[4 * k + 1]].z +
                        p[elems[4 * k + 2]].z + p[elems[4 * k + 3]].z));
      }
    }

    // Count up the number of non-degenerate dual edges
    int n = 0;
    for (int k = 0; k < num_dual_edges; k++) {
      if (dual_edges[2 * k] >= 0 && dual_edges[2 * k + 1] >= 0) {
        n++;
      }
    }

    // Write out the cell connectivity
    fprintf(fp, "\nCELLS %d %d\n", n, 3 * n);
    for (int k = 0; k < num_dual_edges; k++) {
      if (dual_edges[2 * k] >= 0 && dual_edges[2 * k + 1] >= 0) {
        fprintf(fp, "2 %d %d\n", dual_edges[2 * k], dual_edges[2 * k + 1]);
      }
    }

    // Write out the cell types
    fprintf(fp, "\nCELL_TYPES %d\n", n);
    for (int k = 0; k < n; k++) {
      fprintf(fp, "%d\n", 3);
    }
    fclose(fp);
  }
}
