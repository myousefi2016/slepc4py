# Author:  Lisandro Dalcin
# Contact: dalcinl@gmail.com

# -----------------------------------------------------------------------------

cdef extern from "slepc.h":

    struct _p_ST
    ctypedef _p_ST* SlepcST "ST"

    struct _p_IP
    ctypedef _p_IP* SlepcIP "IP"

    struct _p_EPS
    ctypedef _p_EPS* SlepcEPS "EPS"

    struct _p_SVD
    ctypedef _p_SVD* SlepcSVD "SVD"

    struct _p_QEP
    ctypedef _p_QEP* SlepcQEP "QEP"

# -----------------------------------------------------------------------------

from petsc4py.PETSc cimport Object

ctypedef public api class ST(Object) [
    type   PySlepcST_Type,
    object PySlepcSTObject,
    ]:
    cdef SlepcST st

ctypedef public api class IP(Object) [
    type   PySlepcIP_Type,
    object PySlepcIPObject,
    ]:
    cdef SlepcIP ip

ctypedef public api class EPS(Object) [
    type PySlepcEPS_Type,
    object PySlepcEPSObject,
    ]:
    cdef SlepcEPS eps

ctypedef public api class SVD(Object) [
    type   PySlepcSVD_Type,
    object PySlepcSVDObject,
    ]:
    cdef SlepcSVD svd

ctypedef public api class QEP(Object) [
    type   PySlepcQEP_Type,
    object PySlepcQEPObject,
    ]:
    cdef SlepcQEP qep

# -----------------------------------------------------------------------------
