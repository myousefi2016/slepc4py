# --------------------------------------------------------------------

from petsc4py.PETSc cimport Object

# --------------------------------------------------------------------

cdef extern from "slepc.h":
    pass

cdef extern from "slepcst.h":
    struct _p_ST
    ctypedef _p_ST* SlepcST "ST"

cdef extern from "slepcip.h":
    struct _p_IP
    ctypedef _p_IP* SlepcIP "IP"

cdef extern from "slepceps.h":
    struct _p_EPS
    ctypedef _p_EPS* SlepcEPS "EPS"

cdef extern from "slepcsvd.h":
    struct _p_SVD
    ctypedef _p_SVD* SlepcSVD "SVD"

# --------------------------------------------------------------------

ctypedef public api class ST(Object) [type PySlepcST_Type, object PySlepcSTObject]:
    cdef SlepcST st

ctypedef public api class IP(Object) [type PySlepcIP_Type, object PySlepcIPObject]:
    cdef SlepcIP ip

ctypedef public api class EPS(Object) [type PySlepcEPS_Type, object PySlepcEPSObject]:
    cdef SlepcEPS eps

ctypedef public api class SVD(Object) [type PySlepcSVD_Type, object PySlepcSVDObject]:
    cdef SlepcSVD svd

# --------------------------------------------------------------------
