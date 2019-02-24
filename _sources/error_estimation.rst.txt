Error Estimation
****************
There are two mesh refinement strategies that are implemented in :class:`TMR`:
the first is based on the strain energy norm, while the second is an
adjoint-based refinement technique. The strain energy based refinement method is
designed to reduce the solution error in the natural strain energy norm.  This
call takes the form:

.. automodule:: TMR

.. autofunction:: strainEnergyError

The forest can be refined based on the requested target error and the solution
set in :class:`TACS`. 

The adjoint-based refinement method requires an :class:`~TACS.Assembler` object
for both the current level of refinement and a uniformly refined version of the
model. The adjoint-based mesh refinement takes the form:

.. automodule:: TMR

.. autofunction:: adjointError
