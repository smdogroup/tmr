#ifndef TMR_BOUNDARY_CONDITIONS_H
#define TMR_BOUNDARY_CONDITIONS_H

#include "TMRBase.h"

/*
  Specify a list of boundary conditions through a list of names
  and boundary condition information for each node
*/
class TMRBoundaryConditions : public TMREntity {
 public:
  TMRBoundaryConditions();
  ~TMRBoundaryConditions();

  // Add a boundary condition associated with the specified name
  void addBoundaryCondition(const char *name, int num_bcs, const int bc_nums[],
                            const double *_bc_vals = NULL);

  // Get the number of boundary conditions
  int getNumBoundaryConditions();
  void getBoundaryCondition(int bc, const char **_name, int *_num_bcs,
                            const int **_bcs_nums, const double **_bc_vals);

 public:
  // The number of boundary conditions
  int num_bcs;

  // Linked list sub-clas for the boundary conditions -- there are
  // typically only a handful of these objects
  class BCNode {
   public:
    BCNode(const char *_name, int _num_bcs, const int *_bc_nums,
           const double *_bc_vals);
    ~BCNode();
    BCNode *next;
    char *name;
    int num_bcs;
    int *bc_nums;
    double *bc_vals;
  } * bc_root, *bc_current;
};

#endif  // TMR_BOUNDARY_CONDITIONS_H