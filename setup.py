import os
from subprocess import check_output

# Numpy/mpi4py must be installed prior to installing TACS
import numpy
import mpi4py
import tacs

# Import distutils
from setuptools import setup
from distutils.core import Extension as Ext
from Cython.Build import cythonize

# Convert from local to absolute directories
def get_global_dir(files):
    tmr_root = os.path.abspath(os.path.dirname(__file__))
    new = []
    for f in files:
        new.append(os.path.join(tmr_root, f))
    return new

def get_mpi_flags():
    # Split the output from the mpicxx command
    args = check_output(['mpicxx', '-show']).split()

    # Determine whether the output is an include/link/lib command
    inc_dirs, lib_dirs, libs = [], [], []
    for flag in args:
        if flag[:2] == '-I':
            inc_dirs.append(flag[2:])
        elif flag[:2] == '-L':
            lib_dirs.append(flag[2:])
        elif flag[:2] == '-l':
            libs.append(flag[2:])

    return inc_dirs, lib_dirs, libs

inc_dirs, lib_dirs, libs = get_mpi_flags()

# Relative paths for the include/library directories
rel_inc_dirs = ['src', 'src/interfaces', 'src/topology']
rel_lib_dirs = ['lib']
libs.extend(['tmr'])

# Convert from relative to absolute directories
inc_dirs.extend(get_global_dir(rel_inc_dirs))
lib_dirs.extend(get_global_dir(rel_lib_dirs))

# Add the include directories from OpenCascade
if 'CASARCH' in os.environ:
    inc_dirs.append(os.path.join(os.environ['CASROOT'], 
                                 os.environ['CASARCH'], 'include/oce'))
else:
    inc_dirs.append(os.path.join(os.environ['CASROOT'], 'include/oce'))
    
# This should be made more general so that you can specify
# alternate locations for the installation of AMD/METIS
default_ext_inc = ['extern/blossom5-v2.05.src']
inc_dirs.extend(get_global_dir(default_ext_inc))

# Add the numpy/mpi4py directories
inc_dirs.extend([numpy.get_include(), mpi4py.get_include()])
inc_dirs.extend(tacs.get_include())
inc_dirs.extend(tacs.get_cython_include())

# Add the TACS libraries
tacs_lib_dirs, tacs_libs = tacs.get_libraries()
lib_dirs.extend(tacs_lib_dirs)
libs.extend(tacs_libs)

# Add tmr/lib as a runtime directory
runtime_lib_dirs = get_global_dir(['lib'])
runtime_lib_dirs.extend(tacs_lib_dirs)

exts = []
mod = 'TMR'
exts.append(Ext('tmr.%s'%(mod), sources=['tmr/%s.pyx'%(mod)],
                language='c++',
                include_dirs=inc_dirs, libraries=libs, 
                library_dirs=lib_dirs, runtime_library_dirs=runtime_lib_dirs))

setup(name='tmr',
      version=0.1,
      description='Parallel mesh generation utilities',
      author='Graeme J. Kennedy',
      author_email='graeme.kennedy@ae.gatech.edu',
      ext_modules=cythonize(exts, language='c++', 
                            include_path=inc_dirs))
