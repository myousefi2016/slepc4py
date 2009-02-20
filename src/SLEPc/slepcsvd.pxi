cdef extern from "slepcsvd.h" nogil:

    ctypedef char* SlepcSVDType "const char*"
    SlepcSVDType SVDCROSS
    SlepcSVDType SVDCYCLIC
    SlepcSVDType SVDLAPACK
    SlepcSVDType SVDLANCZOS
    SlepcSVDType SVDTRLANCZOS

    ctypedef enum SlepcSVDTransposeMode "SVDTransposeMode":
        SVD_TRANSPOSE_EXPLICIT
        SVD_TRANSPOSE_IMPLICIT

    ctypedef enum SlepcSVDWhich "SVDWhich":
        SVD_LARGEST
        SVD_SMALLEST


    ctypedef enum SlepcSVDConvergedReason "SVDConvergedReason":
        SVD_CONVERGED_ITERATING
        SVD_CONVERGED_TOL
        SVD_DIVERGED_ITS
        SVD_DIVERGED_BREAKDOWN

    int SVDCreate(MPI_Comm,SlepcSVD*)
    int SVDView(SlepcSVD,PetscViewer)
    int SVDDestroy(SlepcSVD)
    int SVDSetType(SlepcSVD,SlepcSVDType)
    int SVDGetType(SlepcSVD,SlepcSVDType*)
    int SVDSetOptionsPrefix(SlepcSVD,char[])
    int SVDAppendOptionsPrefix(SlepcSVD,char[])
    int SVDGetOptionsPrefix(SlepcSVD,char*[])
    int SVDSetFromOptions(SlepcSVD)

    int SVDSetIP(SlepcSVD,SlepcIP)
    int SVDGetIP(SlepcSVD,SlepcIP*)

    int SVDSetOperator(SlepcSVD,PetscMat)
    int SVDGetOperator(SlepcSVD,PetscMat*)

    int SVDSetInitialVector(SlepcSVD,PetscVec)
    int SVDGetInitialVector(SlepcSVD,PetscVec*)

    int SVDSetTransposeMode(SlepcSVD,SlepcSVDTransposeMode)
    int SVDGetTransposeMode(SlepcSVD,SlepcSVDTransposeMode*)
    int SVDSetDimensions(SlepcSVD,PetscInt,PetscInt,PetscInt)
    int SVDGetDimensions(SlepcSVD,PetscInt*,PetscInt*,PetscInt*)
    int SVDSetTolerances(SlepcSVD,PetscReal,PetscInt)
    int SVDGetTolerances(SlepcSVD,PetscReal*,PetscInt*)
    int SVDSetWhichSingularTriplets(SlepcSVD,SlepcSVDWhich)
    int SVDGetWhichSingularTriplets(SlepcSVD,SlepcSVDWhich*)

    int SVDSetUp(SlepcSVD)
    int SVDSolve(SlepcSVD)
    int SVDGetIterationNumber(SlepcSVD,PetscInt*)
    int SVDGetConvergedReason(SlepcSVD,SlepcSVDConvergedReason*)
    int SVDGetConverged(SlepcSVD,PetscInt*)
    int SVDGetSingularTriplet(SlepcSVD,PetscInt,PetscReal*OUTPUT,PetscVec,PetscVec)
    int SVDComputeResidualNorms(SlepcSVD,PetscInt,PetscReal*,PetscReal*)
    int SVDComputeRelativeError(SlepcSVD,PetscInt,PetscReal*)
    int SVDGetOperationCounters(SlepcSVD,PetscInt*,PetscInt*)

    int SVDCrossSetEPS(SlepcSVD,SlepcEPS)
    int SVDCrossGetEPS(SlepcSVD,SlepcEPS*)

    int SVDCyclicSetExplicitMatrix(SlepcSVD,PetscTruth)
    int SVDCyclicGetExplicitMatrix(SlepcSVD,PetscTruth*)
    int SVDCyclicSetEPS(SlepcSVD,SlepcEPS)
    int SVDCyclicGetEPS(SlepcSVD,SlepcEPS*)

    int SVDLanczosSetOneSide(SlepcSVD,PetscTruth)

    int SVDTRLanczosSetOneSide(SlepcSVD,PetscTruth)
