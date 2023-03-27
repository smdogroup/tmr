#include "TMRBoundaryConditions.h"

/*
  Allocate a new boundary condition set
*/
TMRBoundaryConditions::TMRBoundaryConditions() {
  num_bcs = 0;
  bc_root = NULL;
  bc_current = NULL;
}

/*
  Free the data - including the boundary condition information
*/
TMRBoundaryConditions::~TMRBoundaryConditions() {
  // Free the boundary condition information
  while (bc_root) {
    bc_current = bc_root->next;
    delete bc_root;
    bc_root = bc_current;
  }
}

/*
  Constructor for the BCNode sub-class
*/
TMRBoundaryConditions::BCNode::BCNode(const char *_name, int _num_bcs,
                                      const int *_bc_nums,
                                      const double *_bc_vals) {
  next = NULL;
  name = new char[strlen(_name) + 1];
  strcpy(name, _name);
  num_bcs = _num_bcs;

  if (num_bcs > 0) {
    bc_nums = new int[num_bcs];
    bc_vals = new double[num_bcs];
    memcpy(bc_nums, _bc_nums, num_bcs * sizeof(int));
    if (_bc_vals) {
      memcpy(bc_vals, _bc_vals, num_bcs * sizeof(double));
    } else {
      memset(bc_vals, 0.0, num_bcs * sizeof(double));
    }
  } else {
    num_bcs = -1;
    bc_nums = NULL;
    bc_vals = NULL;
  }
}

/*
  Destructor for the BCNode sub-class
*/
TMRBoundaryConditions::BCNode::~BCNode() {
  delete[] name;
  if (bc_nums) {
    delete[] bc_nums;
  }
  if (bc_vals) {
    delete[] bc_vals;
  }
}

/*
  Add the boundary conditions that will be associated with the
  specified name to the boundary condition linked list.
*/
void TMRBoundaryConditions::addBoundaryCondition(const char *name,
                                                 int num_bc_nums,
                                                 const int bc_nums[],
                                                 const double *bc_vals) {
  num_bcs++;
  BCNode *node = new BCNode(name, num_bc_nums, bc_nums, bc_vals);
  if (!bc_root) {
    bc_root = node;
    bc_current = node;
  } else {
    bc_current->next = node;
    bc_current = bc_current->next;
  }
}

/*
  Get the number of boundary conditions
*/
int TMRBoundaryConditions::getNumBoundaryConditions() { return num_bcs; }

/*
  Retrieve the boundary conditions
*/
void TMRBoundaryConditions::getBoundaryCondition(int bc, const char **_name,
                                                 int *_num_bcs,
                                                 const int **_bc_nums,
                                                 const double **_bc_vals) {
  *_name = NULL;
  *_num_bcs = -1;
  *_bc_nums = NULL;
  *_bc_vals = NULL;

  // Increment until we've found the boundary condition or the end of
  // the linked list
  int count = 0;
  BCNode *node = bc_root;
  while (node && count < bc) {
    node = node->next;
    count++;
  }

  // If the node exists, write out the node
  if (node) {
    *_name = node->name;
    *_num_bcs = node->num_bcs;
    *_bc_nums = node->bc_nums;
    *_bc_vals = node->bc_vals;
  }
}
