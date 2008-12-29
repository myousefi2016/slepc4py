cdef extern from "slepceps.h":

    ctypedef char* SlepcEPSType "const char*"
    SlepcEPSType EPSPOWER
    SlepcEPSType EPSSUBSPACE
    SlepcEPSType EPSARNOLDI
    SlepcEPSType EPSLANCZOS
    SlepcEPSType EPSKRYLOVSCHUR
    SlepcEPSType EPSLAPACK
    SlepcEPSType EPSARPACK
    SlepcEPSType EPSBLZPACK
    SlepcEPSType EPSTRLAN
    SlepcEPSType EPSBLOPEX
    SlepcEPSType EPSPRIMME

    ctypedef enum SlepcEPSProblemType "EPSProblemType":
        EPS_HEP
        EPS_GHEP
        EPS_NHEP
        EPS_GNHEP
        EPS_PGNHEP

    ctypedef enum SlepcEPSExtraction "EPSExtraction":
        EPS_RITZ
        EPS_HARMONIC
        EPS_REFINED
        EPS_REFINED_HARMONIC

    ctypedef enum SlepcEPSWhich "EPSWhich":
        EPS_LARGEST_MAGNITUDE, EPS_SMALLEST_MAGNITUDE,
        EPS_LARGEST_REAL,      EPS_SMALLEST_REAL,
        EPS_LARGEST_IMAGINARY, EPS_SMALLEST_IMAGINARY

    ctypedef enum SlepcEPSClass "EPSClass":
         EPS_ONE_SIDE
         EPS_TWO_SIDE

    ctypedef enum SlepcEPSConvergedReason "EPSConvergedReason":
        EPS_CONVERGED_ITERATING
        EPS_CONVERGED_TOL
        EPS_DIVERGED_ITS
        EPS_DIVERGED_BREAKDOWN
        EPS_DIVERGED_NONSYMMETRIC

    ctypedef enum SlepcEPSPowerShiftType "EPSPowerShiftType":
        EPSPOWER_SHIFT_CONSTANT
        EPSPOWER_SHIFT_RAYLEIGH
        EPSPOWER_SHIFT_WILKINSON

    int EPSView(SlepcEPS,PetscViewer)
    int EPSDestroy(SlepcEPS)
    int EPSCreate(MPI_Comm,SlepcEPS*)
    int EPSSetType(SlepcEPS,SlepcEPSType)
    int EPSGetType(SlepcEPS,SlepcEPSType*)
    int EPSSetOptionsPrefix(SlepcEPS,char[])
    int EPSAppendOptionsPrefix(SlepcEPS,char [])
    int EPSGetOptionsPrefix(SlepcEPS,char*[])
    int EPSSetFromOptions(SlepcEPS)

    int EPSSetProblemType(SlepcEPS,SlepcEPSProblemType)
    int EPSGetProblemType(SlepcEPS,SlepcEPSProblemType*)
    int EPSSetExtraction(SlepcEPS,SlepcEPSExtraction)
    int EPSGetExtraction(SlepcEPS,SlepcEPSExtraction*)
    int EPSSetClass(SlepcEPS,SlepcEPSClass)
    int EPSGetClass(SlepcEPS,SlepcEPSClass*)
    int EPSSetWhichEigenpairs(SlepcEPS,SlepcEPSWhich)
    int EPSGetWhichEigenpairs(SlepcEPS,SlepcEPSWhich*)

    int EPSSetTolerances(SlepcEPS,PetscReal,int)
    int EPSGetTolerances(SlepcEPS,PetscReal*,int*)
    int EPSSetDimensions(SlepcEPS,int,int)
    int EPSGetDimensions(SlepcEPS,int*,int*)

    int EPSSetIP(SlepcEPS,SlepcIP)
    int EPSGetIP(SlepcEPS,SlepcIP*)
    int EPSSetST(SlepcEPS,SlepcST)
    int EPSGetST(SlepcEPS,SlepcST*)

    int EPSSetOperators(SlepcEPS,PetscMat,PetscMat)
    int EPSGetOperators(SlepcEPS,PetscMat*,PetscMat*)

    int EPSSetInitialVector(SlepcEPS,PetscVec)
    int EPSGetInitialVector(SlepcEPS,PetscVec*)
    int EPSSetLeftInitialVector(SlepcEPS,PetscVec)
    int EPSGetLeftInitialVector(SlepcEPS,PetscVec*)

    int EPSSetUp(SlepcEPS)
    int EPSSolve(SlepcEPS)

    int EPSGetIterationNumber(SlepcEPS,int*)
    int EPSGetConvergedReason(SlepcEPS,SlepcEPSConvergedReason*)
    int EPSGetConverged(SlepcEPS,int*)
    int EPSGetValue(SlepcEPS,int,PetscScalar*,PetscScalar*)
    int EPSGetRightVector(SlepcEPS,int,PetscVec,PetscVec)
    int EPSGetLeftVector(SlepcEPS,int,PetscVec,PetscVec)
    int EPSGetEigenpair(SlepcEPS,int,PetscScalar*,PetscScalar*,PetscVec,PetscVec)

    int EPSGetErrorEstimate(SlepcEPS,int,PetscReal*)
    int EPSGetErrorEstimateLeft(SlepcEPS,int,PetscReal*)
    int EPSComputeRelativeError(SlepcEPS,int,PetscReal*)
    int EPSComputeRelativeErrorLeft(SlepcEPS,int,PetscReal*)
    int EPSComputeResidualNorm(SlepcEPS,int,PetscReal*)
    int EPSComputeResidualNormLeft(SlepcEPS,int,PetscReal*)
