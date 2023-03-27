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

#include "TMREgads.h"

#include <map>

#ifdef TMR_HAS_EGADS

using namespace TMR_EgadsInterface;

TMR_EgadsContext::TMR_EgadsContext() {
  EG_open(&ctx);
  ismine = 1;
}

TMR_EgadsContext::TMR_EgadsContext(ego _ctx) {
  ctx = _ctx;
  ismine = 0;
}

TMR_EgadsContext::~TMR_EgadsContext() {
  if (ismine) {
    EG_close(ctx);
  }
}

ego TMR_EgadsContext::getContext() { return ctx; }

TMR_EgadsNode::TMR_EgadsNode(TMR_EgadsContext *_ctx, ego _node) {
  ctx = _ctx;
  ctx->incref();
  node = _node;
}

TMR_EgadsNode::~TMR_EgadsNode() { ctx->decref(); }

int TMR_EgadsNode::evalPoint(TMRPoint *X) {
  double eval[3];
  double params[2] = {0.0, 0.0};
  int icode = EG_evaluate(node, params, eval);

  X->x = eval[0];
  X->y = eval[1];
  X->z = eval[2];

  return icode;
}

// Not great, but the best we can do for now...
int TMR_EgadsNode::getParamOnEdge(TMREdge *edge, double *t) {
  return TMRVertex::getParamOnEdge(edge, t);
}

int TMR_EgadsNode::getParamsOnFace(TMRFace *surface, double *u, double *v) {
  int icode = 0;
  TMR_EgadsFace *f = dynamic_cast<TMR_EgadsFace *>(surface);
  if (f) {
    // Get the face topology
    ego face;
    f->getFaceObject(&face);

    int sense = 0;
    double t = 0.0;
    double params[2];
    icode = EG_getEdgeUV(face, node, sense, t, params);

    *u = params[0];
    *v = params[1];
  }

  if (icode != EGADS_SUCCESS) {
    return TMRVertex::getParamsOnFace(surface, u, v);
  }
  return 0;
}

int TMR_EgadsNode::isSame(TMRVertex *vt) {
  TMR_EgadsNode *v = dynamic_cast<TMR_EgadsNode *>(vt);
  if (v) {
    return (EG_isSame(node, v->node) == EGADS_SUCCESS);
  }
  return 0;
}

void TMR_EgadsNode::getNodeObject(ego *n) { *n = node; }

TMR_EgadsEdge::TMR_EgadsEdge(TMR_EgadsContext *_ctx, ego _edge,
                             int _is_degenerate) {
  ctx = _ctx;
  ctx->incref();
  edge = _edge;
  is_degenerate = _is_degenerate;
}

TMR_EgadsEdge::~TMR_EgadsEdge() { ctx->decref(); }

void TMR_EgadsEdge::getRange(double *tmin, double *tmax) {
  double range[4];

  int periodic;
  EG_getRange(edge, range, &periodic);

  *tmin = range[0];
  *tmax = range[1];
}

int TMR_EgadsEdge::getParamsOnFace(TMRFace *surface, double t, int dir,
                                   double *u, double *v) {
  int icode = 0;
  TMR_EgadsFace *f = dynamic_cast<TMR_EgadsFace *>(surface);
  if (f) {
    // Get the face topology
    ego face;
    f->getFaceObject(&face);

    int sense = -1;
    if (dir > 0) {
      sense = 1;
    } else if (dir == 0) {
      sense = 0;
    }

    double params[2];
    icode = EG_getEdgeUV(face, edge, sense, t, params);

    *u = params[0];
    *v = params[1];
  }

  // Fall through to the underlying implementation
  if (icode != EGADS_SUCCESS) {
    return TMREdge::getParamsOnFace(surface, t, dir, u, v);
  }
  return 0;
}

int TMR_EgadsEdge::evalPoint(double t, TMRPoint *X) {
  double eval[18];

  // Evaluate the position of the object
  double params[2];
  params[0] = t;
  params[1] = 0.0;
  int icode = EG_evaluate(edge, params, eval);

  X->x = eval[0];
  X->y = eval[1];
  X->z = eval[2];
  return icode;
}

int TMR_EgadsEdge::invEvalPoint(TMRPoint X, double *t) {
  // Set the position
  double pos[3];
  pos[0] = X.x;
  pos[1] = X.y;
  pos[2] = X.z;

  // Perform the inverse evaluation
  double params[2];
  double result[3];
  int icode = EG_invEvaluate(edge, pos, params, result);

  // Set the parameters
  *t = params[0];

  return icode;
}

int TMR_EgadsEdge::evalDeriv(double t, TMRPoint *X, TMRPoint *Xt) {
  double eval[18];

  // Evaluate the position of the object
  double params[2];
  params[0] = t;
  params[1] = 0.0;
  int icode = EG_evaluate(edge, params, eval);

  X->x = eval[0];
  X->y = eval[1];
  X->z = eval[2];

  Xt->x = eval[3];
  Xt->y = eval[4];
  Xt->z = eval[5];

  return icode;
}

int TMR_EgadsEdge::eval2ndDeriv(double t, TMRPoint *X, TMRPoint *Xt,
                                TMRPoint *Xtt) {
  double eval[18];

  // Evaluate the position of the object
  double params[2];
  params[0] = t;
  params[1] = 0.0;
  int icode = EG_evaluate(edge, params, eval);

  X->x = eval[0];
  X->y = eval[1];
  X->z = eval[2];

  Xt->x = eval[3];
  Xt->y = eval[4];
  Xt->z = eval[5];

  Xtt->x = eval[6];
  Xtt->y = eval[7];
  Xtt->z = eval[8];

  return icode;
}

int TMR_EgadsEdge::isSame(TMREdge *et) {
  TMR_EgadsEdge *e = dynamic_cast<TMR_EgadsEdge *>(et);
  if (e) {
    return (EG_isSame(edge, e->edge) == EGADS_SUCCESS);
  }
  return 0;
}

void TMR_EgadsEdge::getEdgeObject(ego *e) { *e = edge; }

int TMR_EgadsEdge::isDegenerate() { return is_degenerate; }

/*
  TMR interface to the underlying egads face object
*/
TMR_EgadsFace::TMR_EgadsFace(TMR_EgadsContext *_ctx, int _normal_dir, ego _face)
    : TMRFace(_normal_dir) {
  ctx = _ctx;
  ctx->incref();
  face = _face;
}

TMR_EgadsFace::~TMR_EgadsFace() { ctx->decref(); }

void TMR_EgadsFace::getRange(double *umin, double *vmin, double *umax,
                             double *vmax) {
  double range[4];

  int periodic;
  EG_getRange(face, range, &periodic);

  *umin = range[0];
  *umax = range[1];

  *vmin = range[2];
  *vmax = range[3];
}

int TMR_EgadsFace::evalPoint(double u, double v, TMRPoint *X) {
  double eval[18];

  // Evaluate the position of the object
  double params[2];
  params[0] = u;
  params[1] = v;
  int icode = EG_evaluate(face, params, eval);

  X->x = eval[0];
  X->y = eval[1];
  X->z = eval[2];
  return icode;
}

int TMR_EgadsFace::invEvalPoint(TMRPoint X, double *u, double *v) {
  // Set the position
  double pos[3];
  pos[0] = X.x;
  pos[1] = X.y;
  pos[2] = X.z;

  // Perform the inverse evaluation
  double params[2];
  double result[3];
  int icode = EG_invEvaluate(face, pos, params, result);

  // Set the parameters
  *u = params[0];
  *v = params[1];

  return icode;
}

int TMR_EgadsFace::evalDeriv(double u, double v, TMRPoint *X, TMRPoint *Xu,
                             TMRPoint *Xv) {
  double eval[18];

  // Evaluate the position of the object
  double params[2];
  params[0] = u;
  params[1] = v;
  int icode = EG_evaluate(face, params, eval);

  X->x = eval[0];
  X->y = eval[1];
  X->z = eval[2];

  Xu->x = eval[3];
  Xu->y = eval[4];
  Xu->z = eval[5];

  Xv->x = eval[6];
  Xv->y = eval[7];
  Xv->z = eval[8];

  return icode;
}

int TMR_EgadsFace::eval2ndDeriv(double u, double v, TMRPoint *X, TMRPoint *Xu,
                                TMRPoint *Xv, TMRPoint *Xuu, TMRPoint *Xuv,
                                TMRPoint *Xvv) {
  double eval[18];

  // Evaluate the position of the object
  double params[2];
  params[0] = u;
  params[1] = v;
  int icode = EG_evaluate(face, params, eval);

  X->x = eval[0];
  X->y = eval[1];
  X->z = eval[2];

  Xu->x = eval[3];
  Xu->y = eval[4];
  Xu->z = eval[5];

  Xv->x = eval[6];
  Xv->y = eval[7];
  Xv->z = eval[8];

  Xuu->x = eval[9];
  Xuu->y = eval[10];
  Xuu->z = eval[11];

  Xvv->x = eval[12];
  Xvv->y = eval[13];
  Xvv->z = eval[14];

  Xuv->x = eval[15];
  Xuv->y = eval[16];
  Xuv->z = eval[17];

  return icode;
}

int TMR_EgadsFace::isSame(TMRFace *ft) {
  TMR_EgadsFace *f = dynamic_cast<TMR_EgadsFace *>(ft);
  if (f) {
    return (EG_isSame(face, f->face) == EGADS_SUCCESS);
  }
  return 0;
}

void TMR_EgadsFace::getFaceObject(ego *f) { *f = face; }

/*
  Create the TMRModel by loading in an EGADS model file
*/
TMRModel *TMR_EgadsInterface::TMR_LoadModelFromEGADSFile(const char *filename,
                                                         int print_level) {
  // Create the common context for all egads objects
  TMR_EgadsContext *ctx = new TMR_EgadsContext();

  // Load in the egads model from the file
  int flags = 0;
  ego model;
  int icode = EG_loadModel(ctx->getContext(), flags, filename, &model);
  if (icode != 0) {
    fprintf(stderr, "TMR: Error reading egads file with error code %d\n",
            icode);
    return NULL;
  }

  return TMR_EgadsInterface::TMR_ConvertEGADSModel(model, print_level);
}

/*
  Create a TMR model from an EGADS model object
*/
TMRModel *TMR_EgadsInterface::TMR_ConvertEGADSModel(ego model,
                                                    int print_level) {
  ego context;
  EG_getContext(model, &context);

  // Create the common context for all egads objects
  TMR_EgadsContext *ctx = new TMR_EgadsContext(context);

  // Count the number of objects
  int nverts = 0;
  int nedges = 0;
  int nfaces = 0;
  int nsolids = 0;

  // Forward index to ego map
  std::map<int, ego> verts;
  std::map<int, ego> edges;
  std::map<int, ego> faces;
  std::map<int, ego> solids;

  // Reverse mapping
  std::map<ego, int> vert_map;
  std::map<ego, int> edge_map;
  std::map<ego, int> face_map;

  // Get the values
  ego ref;
  int oclass, mtype;
  double data[4];
  int nbodies;
  ego *bodies;
  int *bsense;

  // Get all the bodies in the object
  int icode = EG_getTopology(model, &ref, &oclass, &mtype, data, &nbodies,
                             &bodies, &bsense);

  if (icode) {
    fprintf(stderr, "TMR: EG_getTopology returned error code %d\n", icode);
    return NULL;
  }
  if (oclass != MODEL) {
    fprintf(stderr, "TMR: Must provide model to TMR_ConvertEGADSModel\n");
    return NULL;
  }

  // Iterate over the bodies
  for (int i = 0; i < nbodies; i++) {
    int nchildren;
    ego *body_children;
    int *body_sense;
    icode = EG_getTopology(bodies[i], &ref, &oclass, &mtype, data, &nchildren,
                           &body_children, &body_sense);

    // Check that we have a body object
    if (oclass != BODY) {
      // This should not occur
      fprintf(stderr, "TMR: EGADS model does not define a body\n");
      return NULL;
    } else if (mtype == WIREBODY) {
      // Cannot convert wirebodies to anything useful
      fprintf(stderr, "TMR: Cannot convert wirebody to TMR object\n");
      return NULL;
    } else if (mtype == SOLIDBODY || mtype == SHEETBODY) {
      solids[nsolids] = bodies[i];
      nsolids++;
    }
  }

  // Iterate over the solid bodies or sheetbodies and extract the
  // faces, edges and nodes (vertices)
  for (int index = 0; index < nbodies; index++) {
    int ntopos;
    ego *topos;

    // Add the faces, edges and nodes
    EG_getBodyTopos(solids[index], NULL, FACE, &ntopos, &topos);
    for (int i = 0; i < ntopos; i++) {
      faces[nfaces] = topos[i];
      face_map[topos[i]] = nfaces;
      nfaces++;
    }
    EG_free(topos);

    EG_getBodyTopos(solids[index], NULL, EDGE, &ntopos, &topos);
    for (int i = 0; i < ntopos; i++) {
      edges[nedges] = topos[i];
      edge_map[topos[i]] = nedges;
      nedges++;
    }
    EG_free(topos);

    EG_getBodyTopos(solids[index], NULL, NODE, &ntopos, &topos);
    for (int i = 0; i < ntopos; i++) {
      verts[nverts] = topos[i];
      vert_map[topos[i]] = nverts;
      nverts++;
    }
    EG_free(topos);
  }

  if (print_level > 0) {
    printf(
        "EGADS model loaded with:\nnverts = %d nedges = %d nfaces = %d "
        "nsolids = %d\n",
        nverts, nedges, nfaces, nsolids);
  }

  // Re-iterate through the list and create the objects needed to
  // define the geometry in TMR
  TMRVertex **all_vertices = new TMRVertex *[nverts];
  memset(all_vertices, 0, nverts * sizeof(TMRVertex *));
  for (int index = 0; index < nverts; index++) {
    all_vertices[index] = new TMR_EgadsNode(ctx, verts[index]);

    // Set the "name" attribute
    int atype, len;
    const int *ints = NULL;
    const double *reals = NULL;
    const char *str = NULL;
    if ((EG_attributeRet(verts[index], "name", &atype, &len, &ints, &reals,
                         &str) == EGADS_SUCCESS) &&
        str != NULL) {
      all_vertices[index]->setName(str);
    }
  }

  TMREdge **all_edges = new TMREdge *[nedges];
  memset(all_edges, 0, nedges * sizeof(TMREdge *));
  for (int index = 0; index < nedges; index++) {
    ego ref;            // reference geometry
    int oclass, mtype;  // object class and type
    int nchildren;
    ego *children;
    int *sense;
    EG_getTopology(edges[index], &ref, &oclass, &mtype, data, &nchildren,
                   &children, &sense);

    int isdegenerate = 0;
    int idx1 = vert_map[children[0]];
    int idx2 = -1;
    if (mtype == ONENODE) {
      idx2 = idx1;
    } else if (mtype == DEGENERATE) {
      idx2 = idx1;
      isdegenerate = 1;
    } else if (nchildren == 2) {
      idx2 = vert_map[children[1]];
    }

    // Set a flag to indicate that the edge is degenerate
    all_edges[index] = new TMR_EgadsEdge(ctx, edges[index], isdegenerate);

    if ((idx1 >= 0 && idx1 < nverts) && (idx2 >= 0 && idx2 < nverts)) {
      all_edges[index]->setVertices(all_vertices[idx1], all_vertices[idx2]);
    } else {
      fprintf(stderr,
              "TMR Error: Unable to set vertices for edge %d "
              "from EGADS file\n",
              index);
    }

    // Set the "name" attribute
    int atype, len;
    const int *ints = NULL;
    const double *reals = NULL;
    const char *str = NULL;
    if ((EG_attributeRet(edges[index], "name", &atype, &len, &ints, &reals,
                         &str) == EGADS_SUCCESS) &&
        str != NULL) {
      all_edges[index]->setName(str);
    }
  }

  TMRFace **all_faces = new TMRFace *[nfaces];
  memset(all_faces, 0, nfaces * sizeof(TMRFace *));
  for (int index = 0; index < nfaces; index++) {
    ego ref;            // reference geometry
    int oclass, mtype;  // object class and type
    int num_loops;
    ego *face_loops;
    int *loop_sense;
    EG_getTopology(faces[index], &ref, &oclass, &mtype, data, &num_loops,
                   &face_loops, &loop_sense);

    // Check if the orientation of the face is flipped relative to the
    // natural orientation of the surface
    int orientation = 1;
    if (mtype == SREVERSE) {
      orientation = -1;
    }

    all_faces[index] = new TMR_EgadsFace(ctx, orientation, faces[index]);

    // Set the "name" attribute
    int atype, len;
    const int *ints = NULL;
    const double *reals = NULL;
    const char *str = NULL;
    if ((EG_attributeRet(faces[index], "name", &atype, &len, &ints, &reals,
                         &str) == EGADS_SUCCESS) &&
        str != NULL) {
      all_faces[index]->setName(str);
    }

    // Create the corresponding loops
    for (int i = 0; i < num_loops; i++) {
      ego loop_ref;                 // reference geometry
      int loop_oclass, loop_mtype;  // object class and type
      int num_edges;
      ego *loop_edges;
      int *edge_sense;
      EG_getTopology(face_loops[i], &loop_ref, &loop_oclass, &loop_mtype, data,
                     &num_edges, &loop_edges, &edge_sense);

      // Create the edge list to go into the EdgeLoop object
      TMREdge **edgs = new TMREdge *[num_edges];
      int *dir = new int[num_edges];

      for (int k = 0; k < num_edges; k++) {
        dir[k] = edge_sense[k];
        int edge_index = edge_map[loop_edges[k]];
        edgs[k] = all_edges[edge_index];
      }

      // Allocate the loop with the given edges/directions
      all_faces[index]->addEdgeLoop(loop_sense[i],
                                    new TMREdgeLoop(num_edges, edgs, dir));

      // Free the allocated data
      delete[] edgs;
      delete[] dir;
    }
  }

  // Create the volumes
  TMRVolume **all_vols = new TMRVolume *[nsolids];
  memset(all_vols, 0, nsolids * sizeof(TMRVolume *));
  for (int index = 0; index < nsolids; index++) {
    ego ref;            // reference geometry
    int oclass, mtype;  // object class and type
    int nchild;
    ego *children;
    int *children_sense;
    EG_getTopology(solids[index], &ref, &oclass, &mtype, data, &nchild,
                   &children, &children_sense);

    if (mtype == SHEETBODY) {
      // Retrieve the children of the shell
      int nshell_faces;
      ego *shell_faces;
      int *face_sense;
      EG_getTopology(children[0], &ref, &oclass, &mtype, data, &nshell_faces,
                     &shell_faces, &face_sense);

      // Allocate the faces arrays/directions
      TMRFace **vol_faces = new TMRFace *[nshell_faces];

      // Now, extract the faces from the underlying shell(s)
      for (int i = 0; i < nshell_faces; i++) {
        int face_index = face_map[shell_faces[i]];

        // Assign the face pointer
        vol_faces[i] = all_faces[face_index];
      }

      // Create the volume object
      all_vols[index] = new TMRVolume(nshell_faces, vol_faces);
      delete[] vol_faces;
    } else if (mtype == SOLIDBODY) {
      // Loop over the closed shells
      int ntotal = 0;
      for (int k = 0; k < nchild; k++) {
        // Retrieve the children of the shell
        int nshell_faces;
        ego *shell_faces;
        int *face_sense;
        EG_getTopology(children[k], &ref, &oclass, &mtype, data, &nshell_faces,
                       &shell_faces, &face_sense);
        ntotal += nshell_faces;
      }

      // Allocate the faces arrays/directions
      TMRFace **vol_faces = new TMRFace *[ntotal];

      for (int k = 0, vol_index = 0; k < nchild; k++) {
        // Retrieve the children of the shell
        int nshell_faces;
        ego *shell_faces;
        int *face_sense;
        EG_getTopology(children[k], &ref, &oclass, &mtype, data, &nshell_faces,
                       &shell_faces, &face_sense);

        // Now, extract the faces from the underlying shell(s)
        for (int i = 0; i < nshell_faces; i++, vol_index++) {
          int face_index = face_map[shell_faces[i]];

          // Assign the face pointer
          vol_faces[vol_index] = all_faces[face_index];
        }
      }

      // Create the volume object
      all_vols[index] = new TMRVolume(ntotal, vol_faces);
      delete[] vol_faces;
    }

    // Set the "name" attribute
    int atype, len;
    const int *ints = NULL;
    const double *reals = NULL;
    const char *str = NULL;
    if ((EG_attributeRet(solids[index], "name", &atype, &len, &ints, &reals,
                         &str) == EGADS_SUCCESS) &&
        str != NULL) {
      all_vols[index]->setName(str);
    }
  }

  TMRModel *geo = new TMRModel(nverts, all_vertices, nedges, all_edges, nfaces,
                               all_faces, nsolids, all_vols);

  // Free the arrays
  delete[] all_vertices;
  delete[] all_edges;
  delete[] all_faces;
  delete[] all_vols;

  return geo;
}

#endif  // TMR_HAS_EGADS
