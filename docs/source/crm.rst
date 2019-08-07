Common Research Model (CRM) wingbox
===================================

``egads4py`` can be used to create wingbox models that can then be meshed with TMR.
An example of this application can be found under ``examples/egads/crm`` which generates a STEP model of the undeformed CRM wing.
In this example, there are two steps to the creation of the STEP file:

1) The code ``oml.py`` reads in the IGES surface file and extracts the lofted airfoil curves.
   These curves are not all planar, so they are modified slightly and are projected so that they are parallel with the symmetry plane.
   This modified OML geometry is then written to the files ``ucrm_9_oml.egads`` and ``ucrm_9_oml.step``.
2) The code ``crm.py`` reads in the OML geometry.
   This code uses a doubly-connected edge loop (DCEL) data structure to create a planar layout of ribs and spars.
   This planar layout is then extruded in the vertical direction.
   The intersections between the ribs and spars and OML are computed and imprinted and the leading edge, trailing edge, and upper/lower spars and ribs.

.. image:: crm.png