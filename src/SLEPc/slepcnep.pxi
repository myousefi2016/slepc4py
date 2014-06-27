cdef extern from * nogil:

    ctypedef char* SlepcNEPType "const char*"
    SlepcNEPType NEPRII
    SlepcNEPType NEPSLP
    SlepcNEPType NEPNARNOLDI

    ctypedef enum SlepcNEPWhich "NEPWhich":
        NEP_LARGEST_MAGNITUDE
        NEP_SMALLEST_MAGNITUDE
        NEP_LARGEST_REAL
        NEP_SMALLEST_REAL
        NEP_LARGEST_IMAGINARY
        NEP_SMALLEST_IMAGINARY
        NEP_TARGET_MAGNITUDE
        NEP_TARGET_REAL
        NEP_TARGET_IMAGINARY

    ctypedef enum SlepcNEPConvergedReason "NEPConvergedReason":
        NEP_CONVERGED_FNORM_ABS
        NEP_CONVERGED_FNORM_RELATIVE
        NEP_CONVERGED_SNORM_RELATIVE
        NEP_DIVERGED_LINEAR_SOLVE
        NEP_DIVERGED_FUNCTION_COUNT
        NEP_DIVERGED_MAX_IT
        NEP_DIVERGED_BREAKDOWN
        NEP_DIVERGED_FNORM_NAN
        NEP_CONVERGED_ITERATING

    ctypedef int (*SlepcNEPFunction)(SlepcNEP,
                                     PetscScalar,
                                     PetscMat*,
                                     PetscMat*,
                                     PetscMatStructure*,
                                     void*) except PETSC_ERR_PYTHON

    ctypedef int (*SlepcNEPJacobian)(SlepcNEP,
                                     PetscScalar,
                                     PetscMat*,
                                     PetscMatStructure*,
                                     void*) except PETSC_ERR_PYTHON

    int NEPCreate(MPI_Comm,SlepcNEP*)
    int NEPDestroy(SlepcNEP*)
    int NEPReset(SlepcNEP)
    int NEPView(SlepcNEP,PetscViewer)

    int NEPSetType(SlepcNEP,SlepcNEPType)
    int NEPGetType(SlepcNEP,SlepcNEPType*)
    int NEPSetTarget(SlepcNEP,PetscScalar)
    int NEPGetTarget(SlepcNEP,PetscScalar*)
    int NEPSetOptionsPrefix(SlepcNEP,char*)
    int NEPGetOptionsPrefix(SlepcNEP,char*[])
    int NEPSetFromOptions(SlepcNEP)
    int NEPAppendOptionsPrefix(SlepcNEP,char*)
    int NEPSetUp(SlepcNEP)
    int NEPSolve(SlepcNEP)

    int NEPSetFunction(SlepcNEP,PetscMat,PetscMat,SlepcNEPFunction,void*)
    int NEPGetFunction(SlepcNEP,PetscMat*,PetscMat*,SlepcNEPFunction*,void**)
    int NEPSetJacobian(SlepcNEP,PetscMat,SlepcNEPJacobian,void*)
    int NEPGetJacobian(SlepcNEP,PetscMat*,SlepcNEPJacobian*,void**)
    int NEPSetSplitOperator(SlepcNEP,PetscInt,PetscMat[],SlepcFN[],PetscMatStructure)
    int NEPGetSplitOperatorTerm(SlepcNEP,PetscInt,PetscMat*,SlepcFN*)
    int NEPGetSplitOperatorInfo(SlepcNEP,PetscInt*,PetscMatStructure*)

    int NEPSetBV(SlepcNEP,SlepcBV)
    int NEPGetBV(SlepcNEP,SlepcBV*)
    int NEPSetTolerances(SlepcNEP,PetscReal,PetscReal,PetscReal,PetscInt,PetscInt)
    int NEPGetTolerances(SlepcNEP,PetscReal*,PetscReal*,PetscReal*,PetscInt*,PetscInt*)

    int NEPSetTrackAll(SlepcNEP,PetscBool)
    int NEPGetTrackAll(SlepcNEP,PetscBool*)

    int NEPSetDimensions(SlepcNEP,PetscInt,PetscInt,PetscInt)
    int NEPGetDimensions(SlepcNEP,PetscInt*,PetscInt*,PetscInt*)
    int NEPSetLagPreconditioner(SlepcNEP,PetscInt)
    int NEPGetLagPreconditioner(SlepcNEP,PetscInt*)
    int NEPSetConstCorrectionTol(SlepcNEP,PetscBool)
    int NEPGetConstCorrectionTol(SlepcNEP,PetscBool*)

    int NEPGetConverged(SlepcNEP,PetscInt*)
    int NEPGetEigenpair(SlepcNEP,PetscInt,PetscScalar*,PetscVec)
    int NEPComputeRelativeError(SlepcNEP,PetscInt,PetscReal*)
    int NEPComputeResidualNorm(SlepcNEP,PetscInt,PetscReal*)
    int NEPGetErrorEstimate(SlepcNEP,PetscInt,PetscReal*)

    int NEPMonitorCancel(SlepcNEP)
    int NEPGetIterationNumber(SlepcNEP,PetscInt*)
    int NEPGetOperationCounters(SlepcNEP,PetscInt*,PetscInt*,PetscInt*)

    int NEPSetInitialSpace(SlepcNEP,PetscInt,PetscVec*)
    int NEPSetWhichEigenpairs(SlepcNEP,SlepcNEPWhich)
    int NEPGetWhichEigenpairs(SlepcNEP,SlepcNEPWhich*)

    int NEPGetConvergedReason(SlepcNEP,SlepcNEPConvergedReason*)

# -----------------------------------------------------------------------------

cdef inline PetscMatStructure matstructure(object structure) \
    except <PetscMatStructure>(-1):
    if   structure is None:  return MAT_DIFFERENT_NONZERO_PATTERN
    elif structure is False: return MAT_DIFFERENT_NONZERO_PATTERN
    elif structure is True:  return MAT_SAME_NONZERO_PATTERN
    else:                    return structure

cdef inline Mat ref_Mat(PetscMat mat):
    cdef Mat ob = <Mat> Mat()
    ob.mat = mat
    PetscINCREF(ob.obj)
    return ob

# -----------------------------------------------------------------------------

cdef inline NEP ref_NEP(SlepcNEP nep):
    cdef NEP ob = <NEP> NEP()
    ob.nep = nep
    PetscINCREF(ob.obj)
    return ob

# -----------------------------------------------------------------------------

cdef int NEP_Function(
    SlepcNEP    nep,
    PetscScalar mu,
    PetscMat*   A,
    PetscMat*   B,
    PetscMatStructure* s,
    void*       ctx,
    ) except PETSC_ERR_PYTHON with gil:
    cdef NEP Nep  = ref_NEP(nep)
    cdef Mat Amat = ref_Mat(A[0])
    cdef Mat Bmat = ref_Mat(B[0])
    (function, args, kargs) = Nep.get_attr('__function__')
    retv = function(Nep, mu, Amat, Bmat, *args, **kargs)
    s[0] = matstructure(retv)
    cdef PetscMat Atmp = NULL, Btmp = NULL
    Atmp = A[0]; A[0] = Amat.mat; Amat.mat = Atmp
    Btmp = B[0]; B[0] = Bmat.mat; Bmat.mat = Btmp
    return 0

# -----------------------------------------------------------------------------

cdef int NEP_Jacobian(
    SlepcNEP    nep,
    PetscScalar mu,
    PetscMat*   J,
    PetscMatStructure* s,
    void*       ctx,
    ) except PETSC_ERR_PYTHON with gil:
    cdef NEP Nep  = ref_NEP(nep)
    cdef Mat Jmat = ref_Mat(J[0])
    (jacobian, args, kargs) = Nep.get_attr('__jacobian__')
    retv = jacobian(Nep, mu, Jmat, *args, **kargs)
    s[0] = matstructure(retv)
    cdef PetscMat Jtmp = NULL
    Jtmp = J[0]; J[0] = Jmat.mat; Jmat.mat = Jtmp
    return 0

