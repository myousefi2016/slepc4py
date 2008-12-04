cdef extern from "slepcsvd.h":

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
    int SVDSetDimensions(SlepcSVD,int,int)
    int SVDGetDimensions(SlepcSVD,int*,int*)
    int SVDSetTolerances(SlepcSVD,PetscReal,int)
    int SVDGetTolerances(SlepcSVD,PetscReal*,int*)
    int SVDSetWhichSingularTriplets(SlepcSVD,SlepcSVDWhich)
    int SVDGetWhichSingularTriplets(SlepcSVD,SlepcSVDWhich*)

    int SVDSetUp(SlepcSVD)
    int SVDSolve(SlepcSVD)
    int SVDGetIterationNumber(SlepcSVD, int*)
    int SVDGetConvergedReason(SlepcSVD,SlepcSVDConvergedReason*)
    int SVDGetConverged(SlepcSVD,int*)
    int SVDGetSingularTriplet(SlepcSVD,int,PetscReal*OUTPUT,PetscVec,PetscVec)
    int SVDComputeResidualNorms(SlepcSVD,int,PetscReal*,PetscReal*)
    int SVDComputeRelativeError(SlepcSVD,int,PetscReal*)
    int SVDGetOperationCounters(SlepcSVD,int*,int*)
