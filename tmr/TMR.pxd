# distutils: language=c++

# This file is part of the package TMR for adaptive mesh refinement.

# Copyright (C) 2015 Georgia Tech Research Corporation.
# Additional copyright (C) 2015 Graeme Kennedy.
# All rights reserved.

# TMR is licensed under the Apache License, Version 2.0 (the "License");
# you may not use this software except in compliance with the License.
# You may obtain a copy of the License at

#   http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Import TACS
from tacs import TACS
from tacs.TACS import Vec, VecInterp, Assembler, Mg, Element, Function, Pc
from tacs.TACS cimport (
    Vec,
    VecInterp,
    Assembler,
    Mg,
    Element,
    Function,
    Pc,
    TACSFunction,
    TACSBVec,
    TACSBVecInterp,
    TACSAssembler,
    TACSMg,
    _init_Vec,
    _init_VecInterp,
    _init_Assembler,
    _init_Mg,
    _dynamicTACSMg,
)
from tacs.constitutive import (
    PlaneStressConstitutive,
    SolidConstitutive,
    MaterialProperties,
)
from tacs.constitutive cimport (
    PlaneStressConstitutive,
    SolidConstitutive,
    MaterialProperties,
    TACSMaterialProperties,
)

# Import ParOpt
from paropt.ParOpt import PVec, ProblemBase
from paropt.ParOpt cimport PVec, ProblemBase, ParOptVec, _init_PVec

# Import EGADS
from egads4py.egads import pyego
from egads4py.egads cimport pyego

cdef inline char* tmr_convert_str_to_chars(s):
   if isinstance(s, unicode):
      s = (<unicode>s).encode('utf-8')
   return s

cdef inline str tmr_convert_char_to_str(const char* s):
    if s == NULL:
        return None
    else:
        return s.decode('utf-8')
