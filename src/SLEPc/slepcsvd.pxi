cdef extern from * nogil:

    ctypedef char* SlepcSVDType "const char*"
    SlepcSVDType SVDCROSS
    SlepcSVDType SVDCYCLIC
    SlepcSVDType SVDLAPACK
    SlepcSVDType SVDLANCZOS
    SlepcSVDType SVDTRLANCZOS

    ctypedef enum SlepcSVDWhich "SVDWhich":
        SVD_LARGEST
        SVD_SMALLEST

    ctypedef enum SlepcSVDErrorType "SVDErrorType":
        SVD_ERROR_ABSOLUTE
        SVD_ERROR_RELATIVE

    ctypedef enum SlepcSVDConvergedReason "SVDConvergedReason":
        SVD_CONVERGED_TOL
        SVD_CONVERGED_USER
        SVD_DIVERGED_ITS
        SVD_DIVERGED_BREAKDOWN
        SVD_CONVERGED_ITERATING

    int SVDCreate(MPI_Comm,SlepcSVD*)
    int SVDView(SlepcSVD,PetscViewer)
    int SVDDestroy(SlepcSVD*)
    int SVDReset(SlepcSVD)
    int SVDSetType(SlepcSVD,SlepcSVDType)
    int SVDGetType(SlepcSVD,SlepcSVDType*)
    int SVDSetOptionsPrefix(SlepcSVD,char[])
    int SVDAppendOptionsPrefix(SlepcSVD,char[])
    int SVDGetOptionsPrefix(SlepcSVD,char*[])
    int SVDSetFromOptions(SlepcSVD)

    int SVDSetBV(SlepcSVD,SlepcBV,SlepcBV)
    int SVDGetBV(SlepcSVD,SlepcBV*,SlepcBV*)

    int SVDSetOperator(SlepcSVD,PetscMat)
    int SVDGetOperator(SlepcSVD,PetscMat*)

    int SVDSetInitialSpace(SlepcSVD,PetscInt,PetscVec*)

    int SVDSetImplicitTranspose(SlepcSVD,PetscBool)
    int SVDGetImplicitTranspose(SlepcSVD,PetscBool*)
    int SVDSetDimensions(SlepcSVD,PetscInt,PetscInt,PetscInt)
    int SVDGetDimensions(SlepcSVD,PetscInt*,PetscInt*,PetscInt*)
    int SVDSetTolerances(SlepcSVD,PetscReal,PetscInt)
    int SVDGetTolerances(SlepcSVD,PetscReal*,PetscInt*)
    int SVDSetWhichSingularTriplets(SlepcSVD,SlepcSVDWhich)
    int SVDGetWhichSingularTriplets(SlepcSVD,SlepcSVDWhich*)

    int SVDMonitorCancel(SlepcSVD)

    int SVDSetUp(SlepcSVD)
    int SVDSolve(SlepcSVD)
    int SVDGetIterationNumber(SlepcSVD,PetscInt*)
    int SVDGetConvergedReason(SlepcSVD,SlepcSVDConvergedReason*)
    int SVDGetConverged(SlepcSVD,PetscInt*)
    int SVDGetSingularTriplet(SlepcSVD,PetscInt,PetscReal*,PetscVec,PetscVec)
    int SVDComputeError(SlepcSVD,PetscInt,SlepcSVDErrorType,PetscReal*)
    int SVDErrorView(SlepcSVD,SlepcSVDErrorType,PetscViewer)

    int SVDCrossSetEPS(SlepcSVD,SlepcEPS)
    int SVDCrossGetEPS(SlepcSVD,SlepcEPS*)

    int SVDCyclicSetExplicitMatrix(SlepcSVD,PetscBool)
    int SVDCyclicGetExplicitMatrix(SlepcSVD,PetscBool*)
    int SVDCyclicSetEPS(SlepcSVD,SlepcEPS)
    int SVDCyclicGetEPS(SlepcSVD,SlepcEPS*)

    int SVDLanczosSetOneSide(SlepcSVD,PetscBool)

    int SVDTRLanczosSetOneSide(SlepcSVD,PetscBool)
