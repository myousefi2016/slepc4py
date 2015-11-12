cdef extern from * :
    enum: PETSC_DECIDE
    enum: PETSC_DEFAULT
    enum: PETSC_DETERMINE

    ctypedef enum PetscBool:
        PETSC_TRUE,  PETSC_YES,
        PETSC_FALSE, PETSC_NO,

    ctypedef enum  PetscNormType "NormType":
        PETSC_NORM_1          "NORM_1"
        PETSC_NORM_2          "NORM_2"
        PETSC_NORM_1_AND_2    "NORM_1_AND_2"
        PETSC_NORM_FROBENIUS  "NORM_FROBENIUS"
        PETSC_NORM_INFINITY   "NORM_INFINITY"
        PETSC_NORM_MAX        "NORM_MAX"

cdef extern from * nogil:
    int PetscMalloc(size_t,void*)
    int PetscFree(void*)
    int PetscMemcpy(void*,void*,size_t)
    int PetscMemzero(void*,size_t)

cdef extern from * nogil:
    MPI_Comm PetscObjectComm(PetscObject)
    int PetscObjectReference(PetscObject)
    int PetscObjectDestroy(PetscObject*)

cdef extern from * nogil:
    int MatGetSize(PetscMat,PetscInt*,PetscInt*)
    int MatGetLocalSize(PetscMat,PetscInt*,PetscInt*)

cdef extern from * nogil:
    enum: SLEPC_VERSION_MAJOR
    enum: SLEPC_VERSION_MINOR
    enum: SLEPC_VERSION_SUBMINOR
    enum: SLEPC_VERSION_PATCH
    enum: SLEPC_VERSION_RELEASE
    char* SLEPC_VERSION_DATE
    char* SLEPC_AUTHOR_INFO
    int SlepcInitialize(int*,char***,char[],char[])
    int SlepcFinalize()
    int SlepcInitializeCalled

cdef inline int PetscINCREF(PetscObject *obj):
    if obj    == NULL: return 0
    if obj[0] == NULL: return 0
    return PetscObjectReference(obj[0])

cdef inline int SlepcCLEAR(PetscObject* obj):
    if obj    == NULL: return 0
    if obj[0] == NULL: return 0
    cdef PetscObject tmp
    tmp = obj[0]; obj[0] = NULL
    return PetscObjectDestroy(&tmp)

cdef extern from * nogil:
    ctypedef enum SlepcFunction "SlepcFunction":
        SLEPC_FUNCTION_NONE
        SLEPC_FUNCTION_EXP
        SLEPC_FUNCTION_LAST
