import os
from subprocess import check_output
import sys

# Numpy/mpi4py must be installed prior to installing TACS
import numpy
import mpi4py

# Import distutils
from setuptools import setup
from distutils.core import Extension as Ext
from Cython.Build import cythonize
from Cython.Compiler import Options

# Convert from local to absolute directories
def get_global_dir(files):
    tmr_root = os.path.abspath(os.path.dirname(__file__))
    new = []
    for f in files:
        new.append(os.path.join(tmr_root, f))
    return new

def get_mpi_flags():
    # Split the output from the mpicxx command
    args = check_output(['mpicxx', '-show']).decode('utf-8').split()

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
for sufix in ['include/oce', 'inc', 'include']:
    cas_inc = os.path.join(os.environ['CASROOT'], sufix)
    if os.path.isdir(cas_inc):
        inc_dirs.append(cas_inc)
        break

# This should be made more general so that you can specify
# alternate locations for the installation of AMD/METIS
default_ext_inc = ['extern/blossom5-v2.05.src']
inc_dirs.extend(get_global_dir(default_ext_inc))

# Add the numpy/mpi4py/tacs/paropt include directories
inc_dirs.extend([numpy.get_include(), mpi4py.get_include()])

# Add tmr/lib as a runtime directory
runtime_lib_dirs = get_global_dir(['lib'])

# Add the TACS libraries
import tacs
if 'tacs' in sys.modules:
    inc_dirs.extend(tacs.get_include())
    inc_dirs.extend(tacs.get_cython_include())
    tacs_lib_dirs, tacs_libs = tacs.get_libraries()
    lib_dirs.extend(tacs_lib_dirs)
    libs.extend(tacs_libs)
    runtime_lib_dirs.extend(tacs_lib_dirs)

# Add the ParOpt libraries
import paropt
if 'paropt' in sys.modules:
    inc_dirs.extend(paropt.get_include())
    inc_dirs.extend(paropt.get_cython_include())
    paropt_lib_dirs, paropt_libs = paropt.get_libraries()
    lib_dirs.extend(paropt_lib_dirs)
    libs.extend(paropt_libs)
    runtime_lib_dirs.extend(paropt_lib_dirs)

# Add the egads4py libraries
import egads4py
if 'egads4py' in sys.modules:
    inc_dirs.extend(egads4py.get_include())
    inc_dirs.extend(egads4py.get_cython_include())
    egads4py_lib_dirs, egads4py_libs = egads4py.get_libraries()
    lib_dirs.extend(egads4py_lib_dirs)
    libs.extend(egads4py_libs)
    runtime_lib_dirs.extend(egads4py_lib_dirs)

exts = []
mod = 'TMR'
exts.append(Ext('tmr.%s'%(mod), sources=['tmr/%s.pyx'%(mod)],
                include_dirs=inc_dirs, libraries=libs,
                library_dirs=lib_dirs, runtime_library_dirs=runtime_lib_dirs,
                define_macros=[('math_Memory_HeaderFile', '1')]))
for e in exts:
    e.cython_directives = {'embedsignature': True,
                           'binding': True}

setup(name='tmr',
      version=0.1,
      description='Parallel mesh generation utilities',
      author='Graeme J. Kennedy',
      author_email='graeme.kennedy@ae.gatech.edu',
      ext_modules=cythonize(exts, language='c++',
                            include_path=inc_dirs))
