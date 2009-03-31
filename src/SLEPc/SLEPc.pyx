# --------------------------------------------------------------------

from petsc4py.PETSc import Error

from petsc4py.PETSc import COMM_NULL
from petsc4py.PETSc import COMM_SELF
from petsc4py.PETSc import COMM_WORLD

# --------------------------------------------------------------------

from petsc4py.PETSc cimport MPI_Comm
from petsc4py.PETSc cimport PetscObject, PetscViewer
from petsc4py.PETSc cimport PetscVec, PetscMat
from petsc4py.PETSc cimport PetscKSP, PetscPC

from petsc4py.PETSc cimport Comm
from petsc4py.PETSc cimport Object, Viewer
from petsc4py.PETSc cimport Vec, Mat
from petsc4py.PETSc cimport KSP, PC

# --------------------------------------------------------------------

cdef extern from *:
   ctypedef unsigned long int size_t

# --------------------------------------------------------------------

cdef extern from *:
    ctypedef char* char_p       "char*"
    ctypedef char* const_char_p "const char*"

cdef inline object cp2str(const_char_p p):
    if p == NULL: return None
    else:         return p

cdef inline char_p str2cp(object s) except ? NULL:
    if s is None: return NULL
    else:         return s

include "allocate.pxi"

# --------------------------------------------------------------------

# Vile hack for raising a exception and not contaminating traceback

cdef extern from *:
    enum: PETSC_ERR_PYTHON "(-1)"

cdef extern from *:
    void pyx_raise"__Pyx_Raise"(object, object, void*)

cdef extern from *:
    void *PyExc_RuntimeError
cdef object PetscError = Error

cdef inline int SETERR(int ierr):
    if (<void*>PetscError):
        pyx_raise(PetscError, ierr, NULL)
    else:
        pyx_raise(<object>PyExc_RuntimeError, ierr, NULL)
    return ierr

cdef inline int CHKERR(int ierr) except -1:
    if ierr == 0: return 0                 # no error
    if ierr == PETSC_ERR_PYTHON: return -1 # error in Python call
    if (<void*>PetscError):
        pyx_raise(PetscError, ierr, NULL)
    else:
        pyx_raise(<object>PyExc_RuntimeError, ierr, NULL)
    return -1

# --------------------------------------------------------------------

cdef extern from "petsc.h":
    ctypedef long   PetscInt
    ctypedef double PetscReal
    ctypedef double PetscScalar
    ctypedef PetscInt    const_PetscInt    "const PetscInt"
    ctypedef PetscReal   const_PetscReal   "const PetscReal"
    ctypedef PetscScalar const_PetscScalar "const PetscScalar"

cdef extern from "compat.h":
    pass

cdef extern from "scalar.h":
    object      PyPetscScalar_FromPetscScalar(PetscScalar)
    PetscScalar PyPetscScalar_AsPetscScalar(object) except*

# ---
cdef inline object toInt(PetscInt value):
    return value
cdef inline PetscInt asInt(object value) except? -1:
    return value

# ---
cdef inline object toReal(PetscReal value):
    return value
cdef inline PetscReal asReal(object value) except? -1:
    return value

# ---
cdef inline object toScalar(PetscScalar value):
    return PyPetscScalar_FromPetscScalar(value)
cdef inline PetscScalar asScalar(object value) except*:
    return PyPetscScalar_AsPetscScalar(value)

# --------------------------------------------------------------------

include "slepcmpi.pxi"
include "slepcsys.pxi"
include "slepcst.pxi"
include "slepcip.pxi"
include "slepceps.pxi"
include "slepcsvd.pxi"

# --------------------------------------------------------------------

__doc__ = \
"""
Scalable Library for Eigenvalue Problem Computations.
"""

DECIDE    = PETSC_DECIDE
IGNORE    = PETSC_IGNORE
DEFAULT   = PETSC_DEFAULT
DETERMINE = PETSC_DETERMINE

include "ST.pyx"
include "IP.pyx"
include "EPS.pyx"
include "SVD.pyx"

# --------------------------------------------------------------------

# --------------------------------------------------------------------

cdef extern from "Python.h":
    int Py_AtExit(void (*)())
    void PySys_WriteStderr(char*,...)

cdef extern from "stdio.h" nogil:
    ctypedef struct FILE
    FILE *stderr
    int fprintf(FILE *, char *, ...)

cdef int initialize(object args) except -1:
    if (<int>SlepcInitializeCalled): return 1
    # initialize SLEPC
    CHKERR( SlepcInitialize(NULL, NULL, NULL, NULL) )
    # register finalization function
    if Py_AtExit(finalize) < 0:
        PySys_WriteStderr("warning: could not register"
                          "SlepcFinalize() with Py_AtExit()")
    return 1 # and we are done, enjoy !!

from petsc4py.PETSc cimport TypeRegistryAdd

cdef extern from *:
    PetscCookie SLEPC_ST_COOKIE  "ST_COOKIE"
    PetscCookie SLEPC_IP_COOKIE  "IP_COOKIE"
    PetscCookie SLEPC_EPS_COOKIE "EPS_COOKIE"
    PetscCookie SLEPC_SVD_COOKIE "SVD_COOKIE"

cdef int register(char path[]) except -1:
    # register Python types
    TypeRegistryAdd(SLEPC_ST_COOKIE,  ST)
    TypeRegistryAdd(SLEPC_IP_COOKIE,  IP)
    TypeRegistryAdd(SLEPC_EPS_COOKIE, EPS)
    TypeRegistryAdd(SLEPC_SVD_COOKIE, SVD)
    return 0

cdef void finalize() nogil:
    # finalize SLEPc
    cdef int ierr = 0
    ierr = SlepcFinalize()
    if ierr != 0:
        fprintf(stderr, "SlepcFinalize() failed "
                "[error code: %d]\n", ierr)
    # and we are done, see you later !!

# --------------------------------------------------------------------

def _initialize(args=None):
    cdef int ready = initialize(args)
    if ready: register(NULL)

def _finalize():
    finalize()

# --------------------------------------------------------------------
if 0: raise RuntimeError # Do not remove this line !!!
