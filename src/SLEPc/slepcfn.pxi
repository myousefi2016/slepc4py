cdef extern from * nogil:

    ctypedef char* SlepcFNType "const char*"
    SlepcFNType FNRATIONAL
    SlepcFNType FNEXP
    SlepcFNType FNLOG
    SlepcFNType FNPHI

    int FNCreate(MPI_Comm,SlepcFN*)
    int FNView(SlepcFN,PetscViewer)
    int FNDestroy(SlepcFN*)
    int FNReset(SlepcFN)
    int FNSetType(SlepcFN,SlepcFNType)
    int FNGetType(SlepcFN,SlepcFNType*)

    int FNSetOptionsPrefix(SlepcFN,char[])
    int FNGetOptionsPrefix(SlepcFN,char*[])
    int FNAppendOptionsPrefix(SlepcFN,char[])
    int FNSetFromOptions(SlepcFN)

    int FNSetParameters(SlepcFN,PetscInt,PetscScalar[],PetscInt,PetscScalar[])
    int FNGetParameters(SlepcFN,PetscInt*,PetscScalar*[],PetscInt*,PetscScalar*[])

    int FNEvaluateFunction(SlepcFN,PetscScalar,PetscScalar*)
    int FNEvaluateDerivative(SlepcFN,PetscScalar,PetscScalar*)


cdef object iarray_s(object array, PetscInt* size, PetscScalar** data):
    cdef Py_ssize_t i = 0, n = len(array)
    cdef PetscScalar *a = NULL
    cdef object mem = allocate(n*sizeof(PetscScalar),<void**>&a)
    for i from 0 <= i < n: a[i] = asScalar(array[i])
    if size != NULL: size[0] = <PetscInt> n
    if data != NULL: data[0] = a
    return mem
