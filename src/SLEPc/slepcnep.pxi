cdef extern from * nogil:

    ctypedef char* SlepcNEPType "const char*"
    SlepcNEPType NEPRII
    SlepcNEPType NEPSLP
    SlepcNEPType NEPNARNOLDI
    SlepcNEPType NEPCISS
    SlepcNEPType NEPINTERPOL
    SlepcNEPType NEPNLEIGS

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
        NEP_ALL
        NEP_WHICH_USER

    ctypedef enum SlepcNEPErrorType "NEPErrorType":
        NEP_ERROR_ABSOLUTE
        NEP_ERROR_RELATIVE
        NEP_ERROR_BACKWARD

    ctypedef enum SlepcNEPRefine "NEPRefine":
        NEP_REFINE_NONE
        NEP_REFINE_SIMPLE
        NEP_REFINE_MULTIPLE

    ctypedef enum SlepcNEPRefineScheme "NEPRefineScheme":
        NEP_REFINE_SCHEME_SCHUR
        NEP_REFINE_SCHEME_MBE
        NEP_REFINE_SCHEME_EXPLICIT

    ctypedef enum SlepcNEPConvergedReason "NEPConvergedReason":
        NEP_CONVERGED_TOL
        NEP_CONVERGED_USER
        NEP_DIVERGED_ITS
        NEP_DIVERGED_BREAKDOWN
        NEP_DIVERGED_LINEAR_SOLVE
        NEP_CONVERGED_ITERATING

    ctypedef int (*SlepcNEPFunction)(SlepcNEP,
                                     PetscScalar,
                                     PetscMat,
                                     PetscMat,
                                     void*) except PETSC_ERR_PYTHON

    ctypedef int (*SlepcNEPJacobian)(SlepcNEP,
                                     PetscScalar,
                                     PetscMat,
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
    int NEPSetRG(SlepcNEP,SlepcRG)
    int NEPGetRG(SlepcNEP,SlepcRG*)
    int NEPSetTolerances(SlepcNEP,PetscReal,PetscInt)
    int NEPGetTolerances(SlepcNEP,PetscReal*,PetscInt*)

    int NEPSetTrackAll(SlepcNEP,PetscBool)
    int NEPGetTrackAll(SlepcNEP,PetscBool*)

    int NEPSetDimensions(SlepcNEP,PetscInt,PetscInt,PetscInt)
    int NEPGetDimensions(SlepcNEP,PetscInt*,PetscInt*,PetscInt*)
    int NEPRIISetLagPreconditioner(SlepcNEP,PetscInt)
    int NEPRIIGetLagPreconditioner(SlepcNEP,PetscInt*)
    int NEPRIISetConstCorrectionTol(SlepcNEP,PetscBool)
    int NEPRIIGetConstCorrectionTol(SlepcNEP,PetscBool*)

    int NEPGetConverged(SlepcNEP,PetscInt*)
    int NEPGetEigenpair(SlepcNEP,PetscInt,PetscScalar*,PetscScalar*,PetscVec,PetscVec)
    int NEPComputeError(SlepcNEP,PetscInt,SlepcNEPErrorType,PetscReal*)
    int NEPErrorView(SlepcNEP,SlepcNEPErrorType,PetscViewer)
    int NEPGetErrorEstimate(SlepcNEP,PetscInt,PetscReal*)

    int NEPMonitorCancel(SlepcNEP)
    int NEPGetIterationNumber(SlepcNEP,PetscInt*)

    int NEPSetInitialSpace(SlepcNEP,PetscInt,PetscVec*)
    int NEPSetWhichEigenpairs(SlepcNEP,SlepcNEPWhich)
    int NEPGetWhichEigenpairs(SlepcNEP,SlepcNEPWhich*)

    int NEPGetConvergedReason(SlepcNEP,SlepcNEPConvergedReason*)

# -----------------------------------------------------------------------------

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
    PetscMat    A,
    PetscMat    B,
    void*       ctx,
    ) except PETSC_ERR_PYTHON with gil:
    cdef NEP Nep  = ref_NEP(nep)
    cdef Mat Amat = ref_Mat(A)
    cdef Mat Bmat = ref_Mat(B)
    (function, args, kargs) = Nep.get_attr('__function__')
    retv = function(Nep, toScalar(mu), Amat, Bmat, *args, **kargs)
    cdef PetscMat Atmp = NULL, Btmp = NULL
    Atmp = A; A = Amat.mat; Amat.mat = Atmp
    Btmp = B; B = Bmat.mat; Bmat.mat = Btmp
    return 0

# -----------------------------------------------------------------------------

cdef int NEP_Jacobian(
    SlepcNEP    nep,
    PetscScalar mu,
    PetscMat    J,
    void*       ctx,
    ) except PETSC_ERR_PYTHON with gil:
    cdef NEP Nep  = ref_NEP(nep)
    cdef Mat Jmat = ref_Mat(J)
    (jacobian, args, kargs) = Nep.get_attr('__jacobian__')
    retv = jacobian(Nep, toScalar(mu), Jmat, *args, **kargs)
    cdef PetscMat Jtmp = NULL
    Jtmp = J; J = Jmat.mat; Jmat.mat = Jtmp
    return 0

