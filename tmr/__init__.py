"""
TMR is a parallel meshing/geometry package for use with TACS
"""

import os


def get_cython_include():
    """
    Get the include directory for the Cython .pxd files in TACS
    """
    return [os.path.abspath(os.path.dirname(__file__))]


def get_include():
    """
    Get the include directory for the Cython .pxd files in TACS
    """
    root_path, tail = os.path.split(os.path.abspath(os.path.dirname(__file__)))

    rel_inc_dirs = ["src", "src/interface", "src/topology", "extern/blossom5-v2.05.src"]

    inc_dirs = []
    for path in rel_inc_dirs:
        inc_dirs.append(os.path.join(root_path, path))

    return inc_dirs


def get_libraries():
    """
    Get the library directories
    """
    root_path, tail = os.path.split(os.path.abspath(os.path.dirname(__file__)))

    rel_lib_dirs = ["lib"]
    libs = ["tmr"]
    lib_dirs = []
    for path in rel_lib_dirs:
        lib_dirs.append(os.path.join(root_path, path))

    return lib_dirs, libs
