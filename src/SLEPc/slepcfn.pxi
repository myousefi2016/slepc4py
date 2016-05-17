cdef extern from * nogil:

    ctypedef char* SlepcFNType "const char*"
    SlepcFNType FNCOMBINE
    SlepcFNType FNRATIONAL
    SlepcFNType FNEXP
    SlepcFNType FNLOG
    SlepcFNType FNPHI
    SlepcFNType FNSQRT
    SlepcFNType FNINVSQRT

    ctypedef enum SlepcFNCombineType "FNCombineType":
        FN_COMBINE_ADD
        FN_COMBINE_MULTIPLY
        FN_COMBINE_DIVIDE
        FN_COMBINE_COMPOSE

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

    int FNSetScale(SlepcFN,PetscScalar,PetscScalar)
    int FNGetScale(SlepcFN,PetscScalar*,PetscScalar*)
    int FNEvaluateFunction(SlepcFN,PetscScalar,PetscScalar*)
    int FNEvaluateDerivative(SlepcFN,PetscScalar,PetscScalar*)
    int FNEvaluateFunctionMat(SlepcFN,PetscMat,PetscMat*)

    int FNRationalSetNumerator(SlepcFN,PetscInt,PetscScalar[])
    int FNRationalGetNumerator(SlepcFN,PetscInt*,PetscScalar*[])
    int FNRationalSetDenominator(SlepcFN,PetscInt,PetscScalar[])
    int FNRationalGetDenominator(SlepcFN,PetscInt*,PetscScalar*[])

    int FNCombineSetChildren(SlepcFN,SlepcFNCombineType,SlepcFN,SlepcFN)
    int FNCombineGetChildren(SlepcFN,SlepcFNCombineType*,SlepcFN*,SlepcFN*)

    int FNPhiSetIndex(SlepcFN,PetscInt)
    int FNPhiGetIndex(SlepcFN,PetscInt*)

cdef object iarray_s(object array, PetscInt* size, PetscScalar** data):
    cdef Py_ssize_t i = 0, n = len(array)
    cdef PetscScalar *a = NULL
    cdef object mem = allocate(n*sizeof(PetscScalar),<void**>&a)
    for i from 0 <= i < n: a[i] = asScalar(array[i])
    if size != NULL: size[0] = <PetscInt> n
    if data != NULL: data[0] = a
    return mem
