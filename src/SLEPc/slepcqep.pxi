cdef extern from * nogil:

    ctypedef char* SlepcQEPType "const char*"
    SlepcQEPType QEPLINEAR
    SlepcQEPType QEPQARNOLDI

    ctypedef enum SlepcQEPProblemType "QEPProblemType":
        QEP_GENERAL
        QEP_HERMITIAN
        QEP_GYROSCOPIC

    ctypedef enum SlepcQEPWhich "QEPWhich":
        QEP_LARGEST_MAGNITUDE
        QEP_SMALLEST_MAGNITUDE
        QEP_LARGEST_REAL
        QEP_SMALLEST_REAL
        QEP_LARGEST_IMAGINARY
        QEP_SMALLEST_IMAGINARY

    ctypedef enum SlepcQEPConvergedReason "QEPConvergedReason":
        QEP_CONVERGED_TOL
        QEP_DIVERGED_ITS
        QEP_DIVERGED_BREAKDOWN
        QEP_CONVERGED_ITERATING

    int QEPCreate(MPI_Comm,SlepcQEP*)
    int QEPDestroy(SlepcQEP)
    int QEPView(SlepcQEP,PetscViewer)

    int QEPSetType(SlepcQEP,SlepcQEPType)
    int QEPGetType(SlepcQEP,SlepcQEPType*)
    int QEPSetProblemType(SlepcQEP,SlepcQEPProblemType)
    int QEPGetProblemType(SlepcQEP,SlepcQEPProblemType*)
    int QEPSetOperators(SlepcQEP,PetscMat,PetscMat,PetscMat)
    int QEPGetOperators(SlepcQEP,PetscMat*,PetscMat*,PetscMat*)
    int QEPSetOptionsPrefix(SlepcQEP,char*)
    int QEPGetOptionsPrefix(SlepcQEP,char*[])
    int QEPSetFromOptions(SlepcQEP)
    int QEPAppendOptionsPrefix(SlepcQEP,char*)
    int QEPSetUp(SlepcQEP)
    int QEPSolve(SlepcQEP)

    int QEPSetIP(SlepcQEP,SlepcIP)
    int QEPGetIP(SlepcQEP,SlepcIP*)
    int QEPSetTolerances(SlepcQEP,PetscReal,PetscInt)
    int QEPGetTolerances(SlepcQEP,PetscReal*,PetscInt*)

    int QEPSetTrackAll(SlepcQEP,PetscTruth)
    int QEPGetTrackAll(SlepcQEP,PetscTruth*)

    int QEPSetDimensions(SlepcQEP,PetscInt,PetscInt,PetscInt)
    int QEPGetDimensions(SlepcQEP,PetscInt*,PetscInt*,PetscInt*)
    int QEPSetScaleFactor(SlepcQEP,PetscReal)
    int QEPGetScaleFactor(SlepcQEP,PetscReal*)

    int QEPGetConverged(SlepcQEP,PetscInt*)
    int QEPGetEigenpair(SlepcQEP,PetscInt,PetscScalar*,PetscScalar*,PetscVec,PetscVec)
    int QEPComputeRelativeError(SlepcQEP,PetscInt,PetscReal*)
    int QEPComputeResidualNorm(SlepcQEP,PetscInt,PetscReal*)
    int QEPGetErrorEstimate(SlepcQEP,PetscInt,PetscReal*)

    int QEPMonitorCancel(SlepcQEP)
    int QEPGetIterationNumber(SlepcQEP,PetscInt*)
    int QEPGetOperationCounters(SlepcQEP,PetscInt*,PetscInt*,PetscInt*)

    int QEPSetInitialSpace(SlepcQEP,PetscInt,PetscVec*)
    int QEPSetInitialSpaceLeft(SlepcQEP,PetscInt,PetscVec*)
    int QEPSetWhichEigenpairs(SlepcQEP,SlepcQEPWhich)
    int QEPGetWhichEigenpairs(SlepcQEP,SlepcQEPWhich*)
    int QEPSetLeftVectorsWanted(SlepcQEP,PetscTruth)
    int QEPGetLeftVectorsWanted(SlepcQEP,PetscTruth*)

    int QEPGetConvergedReason(SlepcQEP,SlepcQEPConvergedReason*)

    int QEPLinearSetCompanionForm(SlepcQEP,PetscInt)
    int QEPLinearGetCompanionForm(SlepcQEP,PetscInt*)
    int QEPLinearSetExplicitMatrix(SlepcQEP,PetscTruth)
    int QEPLinearGetExplicitMatrix(SlepcQEP,PetscTruth*)
    int QEPLinearSetEPS(SlepcQEP,SlepcEPS)
    int QEPLinearGetEPS(SlepcQEP,SlepcEPS*)

