cdef extern from "slepc.h":

    enum: SLEPC_VERSION_MAJOR
    enum: SLEPC_VERSION_MINOR
    enum: SLEPC_VERSION_SUBMINOR
    enum: SLEPC_VERSION_PATCH
    enum: SLEPC_VERSION_RELEASE
    char* SLEPC_VERSION_DATE
    char* SLEPC_VERSION_PATCH_DATE
    char* SLEPC_AUTHOR_INFO

    int SlepcInitialize(int*,char***,char[],char[])
    int SlepcFinalize()
    int SlepcInitializeCalled

cdef extern from "petsc.h":
    enum: PETSC_DECIDE
    enum: PETSC_DEFAULT
    enum: PETSC_DETERMINE
    enum: PETSC_IGNORE

cdef extern from "petsc.h":
    int PetscObjectReference(PetscObject)
    int PetscObjectDereference(PetscObject)
    int PetscObjectDestroy(PetscObject)

cdef inline int SlepcIncref(PetscObject obj):
    if obj != NULL:
        return PetscObjectReference(obj)
    return 0

cdef inline int SlepcDecref(PetscObject obj):
    if obj != NULL:
        return PetscObjectDereference(obj)
    return 0

cdef inline int SlepcCLEAR(PetscObject* obj):
    if obj == NULL: return 0
    cdef PetscObject tmp = obj[0]
    if tmp == NULL: return 0
    obj[0] = NULL
    return PetscObjectDestroy(tmp)
