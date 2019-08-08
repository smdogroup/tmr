Topology Optimization Utilities
===============================

TMR implements the class :class:`~tmr.TMR.TopoProblem` which is designed for topology optimization with ParOpt. 
This class can be set up to contain information about how to evaluate the objective function and constraints using a TACS.Assembler object.
The adjoint derivative calculations for all members are implemented within the code.

* .. autoclass:: tmr.TMR.TopoProblem
      :members:

While you can create the :class:`~tmr.TMR.TopoProblem` from scratch, there are several utilities that can assist you with the creation of this object.
:meth:`~tmr.TopOptUtils.createTopoProblem` generates a problem class with a specified hierarchy of octree or quadtree meshes.
The construction of the finite-element discretization and the topology parametrization are problem-specific, so these are generated through a callback.
Several additional helper functions are implemented to enable the creation of load vectors: :meth:`~tmr.TopOptUtils.computeVertexLoad` creates a load vector from forces at named vertices, 
:meth:`~tmr.TopOptUtils.computeTractionLoad` creates a load vector based on traction on named edges or faces.
:meth:`~tmr.TopOptUtils.interpolateDesignVec` interpolates between design vectors after a refinement step.

* .. automodule:: tmr.TopOptUtils
      :members:

* .. autofunction:: tmr.TMR.convertPVecToVec
