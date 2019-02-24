Installing TMR
**************

Prerequisites
=============
The following items are needed to use TMR:

* Blossom V
* METIS
* numpy
* MPI compiler and mpi4py
* Cython

The following items are optional but needed to use some of the advance features of TMR:

* OpenCASCADE
* Netgen


Steps to compile
================
#. Clone the TMR git repository
#. In the base 'tmr' directory, copy the Makefile.in.info to Makefile.in. Edit
   the Makefile.in to match your compilers, metis and blossom V location,
   etc. If using OpenCASCADE, additional flags are required
#. To compile, from the base directory, run *make* then *make interface*
