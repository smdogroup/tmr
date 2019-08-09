Installing TMR
==============

Prerequisites
-------------

TMR has several dependencies. The current version of TMR requires the following external libraries:

* MPI 2.0 compatible compiler 
* BLAS/LAPACK
* EGADS/egads4py (for geometry and topology information). egads4py is available here: `https://github.com/gjkennedy/egads4py <https://github.com/gjkennedy/egads4py>`_
* TACS (for finite-element analysis). TACS is available here: `https://github.com/gjkennedy/paropt <https://github.com/gjkennedy/paropt>`_
* ParOpt (for topology optimization). ParOpt is available here: `https://github.com/gjkennedy/tacs <https://github.com/gjkennedy/tacs>`_
* Blossom V (for perfect matching in the Quad-Blossom algorithm): `Blossom-V <https://pub.ist.ac.at/~vnk/software/blossom5-v2.05.src.tar.gz>`_
* METIS (for partitioning the coarse hexahedral/quadrilateral meshes): `Metis <http://glaros.dtc.umn.edu/gkhome/metis/metis/overview>`_
* OpenCASCADE

In addition, the python interface for TMR requires the following:

* numpy
* mpi4py
* Cython

Steps to compile
----------------
#. Clone the TMR git repository
#. Make sure your MPI complier works. It is best to use a package manager such as apt-get (linux) or brew (OS X).
#. Follow the instructions to install ParOpt and TACS
#. Install the externals METIS and Blossom V. Note that the Makefile for Blossom V in the extern package creates a library libblossom5.a.
#. Follow the instructions below to install OpenCASCADE
#. In the base 'tmr' directory, copy the Makefile.in.info to Makefile.in. Edit
   the Makefile.in to match your compilers, metis and blossom V location,
   etc. If using OpenCASCADE, additional flags are required
#. To compile, from the base directory, run *make* then *make interface* for the python interface

OpenCASCADE
-----------

There are different ways to install the OpenCASCADE libraries. The instructions below should work for linux and OS X.

1. I downloaded the pre-compiled 64 bit version of OpenCASCADE for linux or OS X from Bob Haimes. For linux I got the following file: OCC680lin64.tgz, and for OS X I got the file: OCC680osx64.tgz. The website is here:

`https://acdl.mit.edu/ESP/ <https://acdl.mit.edu/ESP/>`_


2. I untared/ziped with tar -xzf OCC680lin64.tgz in $HOME/packages/

3. I move OpenCASCADE-6.8.0 to OpenCASCADE

4. I added the following lines to my .bashrc file:

.. code-block::

    # Export the OpenCASCADE root directory
    export CASROOT=$HOME/packages/OpenCASCADE
    export PATH=$CASROOT/bin:$PATH
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CASROOT/lib

    For OSX I added the following lines
    export CASROOT=$HOME/packages/OpenCASCADE
    export CASARCH=i686
    DYLD_LIBRARY_PATH=${CASROOT}/${CASARCH}/lib
    PATH=$CASROOT/bin:$PATH

5. I tried to use the draw command by typing: 

.. code-block::

    sh ./draw.sh

Once the command prompt appeared I typed:

.. code-block::

    Draw[3]> source samples/tcl/bottle.tcl

This showed an image of the bottle example from OpenCASCADE

(This step did not work under OSX, but the final step and compilation
of the TMR library with OpenCASCADE enabled worked anyway.)

6. I added the following lines to the Makefile.in code:

.. code-block::

    # Flags for debugging and regular compilation versions
    TMR_DEBUG_FLAGS = [...] -DTMR_HAS_OPENCASCADE
    TMR_FLAGS = [...] -DTMR_HAS_OPENCASCADE

    OPENCASCADE_INCLUDE = -I${CASROOT}/${CASARCH}/inc
    OPENCASCADE_LIB_PATH = -L${CASROOT}/${CASARCH}/lib -Wl,-rpath,${CASROOT}/${CASARCH}/lib
    OPENCASCADE_LIBS = -lTKBool -lTKernel -lTKFeat -lTKBO -lTKGeomAlgo -lTKMath -lTKOffset -lTKPrim -lTKPShape -lTKTopAlgo -lTKBRep -lTKG2d -lTKG3d -lTKGeomBase -lTKShHealing -lTKSTEP -lTKSTEP209 -lTKSTEPBase -lTKSTEPAttr -lTKXSBase -lTKIGES -lTKFillet -lPTKernel -ldl

If using OpenCASCADE Community Edition (OCE) (not using compiled version):

1. Clone from git repository https://github.com/tpaviot/oce to $HOME/packages/
2. Create an empty directory $HOME/packages/OpenCASCADE
3. Get to $HOME/packages/oce
4. Use the following commands in the directory $HOME/packages/oce

.. code-block::
   
   cmake -DOCE_INSTALL_PREFIX:PATH=$HOME/packages/OpenCASCADE
   
This will build the libraries in $HOME/packages/OpenCASCADE

5. Once cmake is done, a Makefile should be created in $HOME/packages/oce Run the commands in $HOME/packages/oce: make, make install
6. Once completed, the build files should show up in $HOME/packages/OpenCASCADE
7. Add the following lines to my .bashrc file:

.. code-block::

   # Export the OpenCASCADE root directory
   export CASROOT=$HOME/packages/OpenCASCADE
   export PATH=$CASROOT/bin:$PATH
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CASROOT/lib

8. Added the following lines to the tmr/Makefile.in code:

.. code-block::

   # Flags for debugging and regular compilation versions
   TMR_DEBUG_FLAGS = [...] -DTMR_HAS_OPENCASCADE
   TMR_FLAGS = [...] -DTMR_HAS_OPENCASCADE

   OPENCASCADE_INCLUDE = -I${CASROOT}/${CASARCH}/include/oce
   OPENCASCADE_LIB_PATH = -L${CASROOT}/${CASARCH}/lib -Wl,-rpath,${CASROOT}/${CASARCH}/lib
   OPENCASCADE_LIBS = -lTKBool -lTKernel -lTKFeat -lTKBO -lTKGeomAlgo -lTKMath -lTKOffset -lTKPrim -lTKPShape -lTKTopAlgo -lTKBRep -lTKG2d -lTKG3d -lTKGeomBase -lTKShHealing -lTKSTEP -lTKSTEP209 -lTKSTEPBase -lTKSTEPAttr -lTKXSBase -lTKIGES -lTKFillet -lPTKernel -ldl
