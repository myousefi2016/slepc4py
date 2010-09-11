cdef extern from "slepceps.h" nogil:

    ctypedef char* SlepcEPSType "const char*"
    SlepcEPSType EPSPOWER
    SlepcEPSType EPSSUBSPACE
    SlepcEPSType EPSARNOLDI
    SlepcEPSType EPSLANCZOS
    SlepcEPSType EPSKRYLOVSCHUR
    SlepcEPSType EPSGD
    SlepcEPSType EPSJD
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
        EPS_HARMONIC_RELATIVE
        EPS_HARMONIC_RIGHT
        EPS_HARMONIC_LARGEST
        EPS_REFINED
        EPS_REFINED_HARMONIC

    ctypedef enum SlepcEPSWhich "EPSWhich":
        EPS_LARGEST_MAGNITUDE
        EPS_LARGEST_REAL
        EPS_LARGEST_IMAGINARY
        EPS_SMALLEST_MAGNITUDE
        EPS_SMALLEST_REAL
        EPS_SMALLEST_IMAGINARY
        EPS_TARGET_MAGNITUDE
        EPS_TARGET_REAL
        EPS_TARGET_IMAGINARY
        EPS_WHICH_USER

    ctypedef enum SlepcEPSBalance "EPSBalance":
        EPS_BALANCE_NONE
        EPS_BALANCE_ONESIDE
        EPS_BALANCE_TWOSIDE
        EPS_BALANCE_USER

    ctypedef enum SlepcEPSConvergedReason "EPSConvergedReason":
        EPS_CONVERGED_ITERATING
        EPS_CONVERGED_TOL
        EPS_DIVERGED_ITS
        EPS_DIVERGED_BREAKDOWN
        EPS_DIVERGED_NONSYMMETRIC

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
    int EPSIsGeneralized(SlepcEPS,PetscTruth*)
    int EPSIsHermitian(SlepcEPS,PetscTruth*)
    int EPSSetExtraction(SlepcEPS,SlepcEPSExtraction)
    int EPSGetExtraction(SlepcEPS,SlepcEPSExtraction*)
    int EPSSetBalance(SlepcEPS,SlepcEPSBalance,PetscInt,PetscReal)
    int EPSGetBalance(SlepcEPS,SlepcEPSBalance*,PetscInt*,PetscReal*)
    int EPSSetWhichEigenpairs(SlepcEPS,SlepcEPSWhich)
    int EPSGetWhichEigenpairs(SlepcEPS,SlepcEPSWhich*)
    int EPSSetLeftVectorsWanted(SlepcEPS,PetscTruth)
    int EPSGetLeftVectorsWanted(SlepcEPS,PetscTruth*)
    int EPSSetTarget(SlepcEPS,PetscScalar)
    int EPSGetTarget(SlepcEPS,PetscScalar*)

    int EPSSetTolerances(SlepcEPS,PetscReal,PetscInt)
    int EPSGetTolerances(SlepcEPS,PetscReal*,PetscInt*)
    int EPSSetDimensions(SlepcEPS,PetscInt,PetscInt,PetscInt)
    int EPSGetDimensions(SlepcEPS,PetscInt*,PetscInt*,PetscInt*)

    int EPSSetIP(SlepcEPS,SlepcIP)
    int EPSGetIP(SlepcEPS,SlepcIP*)
    int EPSSetST(SlepcEPS,SlepcST)
    int EPSGetST(SlepcEPS,SlepcST*)

    int EPSSetOperators(SlepcEPS,PetscMat,PetscMat)
    int EPSGetOperators(SlepcEPS,PetscMat*,PetscMat*)

    int EPSSetTrackAll(SlepcEPS,PetscTruth)
    int EPSGetTrackAll(SlepcEPS,PetscTruth*)

    int EPSSetDeflationSpace(SlepcEPS,PetscInt,PetscVec*)
    int EPSRemoveDeflationSpace(SlepcEPS)

    int EPSSetInitialSpace(SlepcEPS,PetscInt,PetscVec*)
    int EPSSetInitialSpaceLeft(SlepcEPS,PetscInt,PetscVec*)

    int EPSMonitorCancel(SlepcEPS)

    int EPSSetUp(SlepcEPS)
    int EPSSolve(SlepcEPS)

    int EPSGetIterationNumber(SlepcEPS,PetscInt*)
    int EPSGetConvergedReason(SlepcEPS,SlepcEPSConvergedReason*)
    int EPSGetConverged(SlepcEPS,PetscInt*)
    int EPSGetEigenvalue(SlepcEPS,PetscInt,PetscScalar*,PetscScalar*)
    int EPSGetEigenvector(SlepcEPS,PetscInt,PetscVec,PetscVec)
    int EPSGetEigenvectorLeft(SlepcEPS,PetscInt,PetscVec,PetscVec)
    int EPSGetEigenpair(SlepcEPS,PetscInt,PetscScalar*,PetscScalar*,PetscVec,PetscVec)
    int EPSGetInvariantSubspace(SlepcEPS,PetscVec*)
    int EPSGetInvariantSubspaceLeft(SlepcEPS,PetscVec*)

    int EPSGetErrorEstimate(SlepcEPS,PetscInt,PetscReal*)
    int EPSGetErrorEstimateLeft(SlepcEPS,PetscInt,PetscReal*)
    int EPSComputeRelativeError(SlepcEPS,PetscInt,PetscReal*)
    int EPSComputeRelativeErrorLeft(SlepcEPS,PetscInt,PetscReal*)
    int EPSComputeResidualNorm(SlepcEPS,PetscInt,PetscReal*)
    int EPSComputeResidualNormLeft(SlepcEPS,PetscInt,PetscReal*)
    int EPSGetOperationCounters(SlepcEPS,PetscInt*,PetscInt*,PetscInt*)


    ctypedef enum SlepcEPSPowerShiftType "EPSPowerShiftType":
        EPS_POWER_SHIFT_CONSTANT
        EPS_POWER_SHIFT_RAYLEIGH
        EPS_POWER_SHIFT_WILKINSON
    int EPSPowerSetShiftType(SlepcEPS,SlepcEPSPowerShiftType)
    int EPSPowerGetShiftType(SlepcEPS,SlepcEPSPowerShiftType*)

    int EPSArnoldiSetDelayed(SlepcEPS,PetscTruth)
    int EPSArnoldiGetDelayed(SlepcEPS,PetscTruth*)

    ctypedef enum SlepcEPSLanczosReorthogType "EPSLanczosReorthogType":
        EPS_LANCZOS_REORTHOG_LOCAL
        EPS_LANCZOS_REORTHOG_FULL
        EPS_LANCZOS_REORTHOG_SELECTIVE
        EPS_LANCZOS_REORTHOG_PERIODIC
        EPS_LANCZOS_REORTHOG_PARTIAL
        EPS_LANCZOS_REORTHOG_DELAYED
    int EPSLanczosSetReorthog(SlepcEPS,SlepcEPSLanczosReorthogType)
    int EPSLanczosGetReorthog(SlepcEPS,SlepcEPSLanczosReorthogType*)


cdef extern from * nogil:
    int VecDuplicate(PetscVec,PetscVec*)
