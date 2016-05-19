cdef extern from * nogil:

    ctypedef char* SlepcEPSType "const char*"
    SlepcEPSType EPSPOWER
    SlepcEPSType EPSSUBSPACE
    SlepcEPSType EPSARNOLDI
    SlepcEPSType EPSLANCZOS
    SlepcEPSType EPSKRYLOVSCHUR
    SlepcEPSType EPSGD
    SlepcEPSType EPSJD
    SlepcEPSType EPSRQCG
    SlepcEPSType EPSLOBPCG
    SlepcEPSType EPSCISS
    SlepcEPSType EPSLAPACK
    SlepcEPSType EPSARPACK
    SlepcEPSType EPSBLZPACK
    SlepcEPSType EPSTRLAN
    SlepcEPSType EPSBLOPEX
    SlepcEPSType EPSPRIMME
    SlepcEPSType EPSFEAST

    ctypedef enum SlepcEPSProblemType "EPSProblemType":
        EPS_HEP
        EPS_GHEP
        EPS_NHEP
        EPS_GNHEP
        EPS_PGNHEP
        EPS_GHIEP

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
        EPS_ALL
        EPS_WHICH_USER

    ctypedef enum SlepcEPSBalance "EPSBalance":
        EPS_BALANCE_NONE
        EPS_BALANCE_ONESIDE
        EPS_BALANCE_TWOSIDE
        EPS_BALANCE_USER

    ctypedef enum SlepcEPSErrorType "EPSErrorType":
        EPS_ERROR_ABSOLUTE
        EPS_ERROR_RELATIVE
        EPS_ERROR_BACKWARD

    ctypedef enum SlepcEPSConv "EPSConv":
        EPS_CONV_ABS
        EPS_CONV_REL
        EPS_CONV_NORM
        EPS_CONV_USER

    ctypedef enum SlepcEPSConvergedReason "EPSConvergedReason":
        EPS_CONVERGED_TOL
        EPS_CONVERGED_USER
        EPS_DIVERGED_ITS
        EPS_DIVERGED_BREAKDOWN
        EPS_DIVERGED_SYMMETRY_LOST
        EPS_CONVERGED_ITERATING

    int EPSView(SlepcEPS,PetscViewer)
    int EPSDestroy(SlepcEPS*)
    int EPSReset(SlepcEPS)
    int EPSCreate(MPI_Comm,SlepcEPS*)
    int EPSSetType(SlepcEPS,SlepcEPSType)
    int EPSGetType(SlepcEPS,SlepcEPSType*)
    int EPSSetOptionsPrefix(SlepcEPS,char[])
    int EPSAppendOptionsPrefix(SlepcEPS,char [])
    int EPSGetOptionsPrefix(SlepcEPS,char*[])
    int EPSSetFromOptions(SlepcEPS)

    int EPSSetProblemType(SlepcEPS,SlepcEPSProblemType)
    int EPSGetProblemType(SlepcEPS,SlepcEPSProblemType*)
    int EPSIsGeneralized(SlepcEPS,PetscBool*)
    int EPSIsHermitian(SlepcEPS,PetscBool*)
    int EPSIsPositive(SlepcEPS,PetscBool*)
    int EPSSetExtraction(SlepcEPS,SlepcEPSExtraction)
    int EPSGetExtraction(SlepcEPS,SlepcEPSExtraction*)
    int EPSSetBalance(SlepcEPS,SlepcEPSBalance,PetscInt,PetscReal)
    int EPSGetBalance(SlepcEPS,SlepcEPSBalance*,PetscInt*,PetscReal*)
    int EPSSetWhichEigenpairs(SlepcEPS,SlepcEPSWhich)
    int EPSGetWhichEigenpairs(SlepcEPS,SlepcEPSWhich*)
    int EPSSetTarget(SlepcEPS,PetscScalar)
    int EPSGetTarget(SlepcEPS,PetscScalar*)
    int EPSSetInterval(SlepcEPS,PetscReal,PetscReal)
    int EPSGetInterval(SlepcEPS,PetscReal*,PetscReal*)
    int EPSSetTolerances(SlepcEPS,PetscReal,PetscInt)
    int EPSGetTolerances(SlepcEPS,PetscReal*,PetscInt*)
    int EPSSetDimensions(SlepcEPS,PetscInt,PetscInt,PetscInt)
    int EPSGetDimensions(SlepcEPS,PetscInt*,PetscInt*,PetscInt*)

    int EPSSetBV(SlepcEPS,SlepcBV)
    int EPSGetBV(SlepcEPS,SlepcBV*)
    int EPSSetDS(SlepcEPS,SlepcDS)
    int EPSGetDS(SlepcEPS,SlepcDS*)
    int EPSSetST(SlepcEPS,SlepcST)
    int EPSGetST(SlepcEPS,SlepcST*)
    int EPSSetRG(SlepcEPS,SlepcRG)
    int EPSGetRG(SlepcEPS,SlepcRG*)

    int EPSSetOperators(SlepcEPS,PetscMat,PetscMat)
    int EPSGetOperators(SlepcEPS,PetscMat*,PetscMat*)

    int EPSSetConvergenceTest(SlepcEPS,SlepcEPSConv)
    int EPSGetConvergenceTest(SlepcEPS,SlepcEPSConv*)

    int EPSSetTrueResidual(SlepcEPS,PetscBool)
    int EPSGetTrueResidual(SlepcEPS,PetscBool*)

    int EPSSetTrackAll(SlepcEPS,PetscBool)
    int EPSGetTrackAll(SlepcEPS,PetscBool*)

    int EPSSetDeflationSpace(SlepcEPS,PetscInt,PetscVec*)
    int EPSSetInitialSpace(SlepcEPS,PetscInt,PetscVec*)

    int EPSMonitorCancel(SlepcEPS)

    int EPSSetUp(SlepcEPS)
    int EPSSolve(SlepcEPS)

    int EPSGetIterationNumber(SlepcEPS,PetscInt*)
    int EPSGetConvergedReason(SlepcEPS,SlepcEPSConvergedReason*)
    int EPSGetConverged(SlepcEPS,PetscInt*)
    int EPSGetEigenvalue(SlepcEPS,PetscInt,PetscScalar*,PetscScalar*)
    int EPSGetEigenvector(SlepcEPS,PetscInt,PetscVec,PetscVec)
    int EPSGetEigenpair(SlepcEPS,PetscInt,PetscScalar*,PetscScalar*,PetscVec,PetscVec)
    int EPSGetInvariantSubspace(SlepcEPS,PetscVec*)

    int EPSGetErrorEstimate(SlepcEPS,PetscInt,PetscReal*)
    int EPSComputeError(SlepcEPS,PetscInt,SlepcEPSErrorType,PetscReal*)
    int EPSErrorView(SlepcEPS,SlepcEPSErrorType,PetscViewer)

    ctypedef enum SlepcEPSPowerShiftType "EPSPowerShiftType":
        EPS_POWER_SHIFT_CONSTANT
        EPS_POWER_SHIFT_RAYLEIGH
        EPS_POWER_SHIFT_WILKINSON
    int EPSPowerSetShiftType(SlepcEPS,SlepcEPSPowerShiftType)
    int EPSPowerGetShiftType(SlepcEPS,SlepcEPSPowerShiftType*)

    int EPSArnoldiSetDelayed(SlepcEPS,PetscBool)
    int EPSArnoldiGetDelayed(SlepcEPS,PetscBool*)

    int EPSKrylovSchurSetRestart(SlepcEPS,PetscReal)
    int EPSKrylovSchurGetRestart(SlepcEPS,PetscReal*)
    int EPSKrylovSchurSetLocking(SlepcEPS,PetscBool);
    int EPSKrylovSchurGetLocking(SlepcEPS,PetscBool*);
    int EPSKrylovSchurSetPartitions(SlepcEPS,PetscInt);
    int EPSKrylovSchurGetPartitions(SlepcEPS,PetscInt*);
    int EPSKrylovSchurSetDetectZeros(SlepcEPS,PetscBool);
    int EPSKrylovSchurGetDetectZeros(SlepcEPS,PetscBool*);
    int EPSKrylovSchurSetDimensions(SlepcEPS,PetscInt,PetscInt,PetscInt);
    int EPSKrylovSchurGetDimensions(SlepcEPS,PetscInt*,PetscInt*,PetscInt*);
    int EPSKrylovSchurGetSubcommInfo(SlepcEPS,PetscInt*,PetscInt*,PetscVec*);
    int EPSKrylovSchurGetSubcommPairs(SlepcEPS,PetscInt,PetscScalar*,PetscVec);
    int EPSKrylovSchurGetSubcommMats(SlepcEPS,PetscMat*,PetscMat*);
    int EPSKrylovSchurUpdateSubcommMats(SlepcEPS,PetscScalar,PetscScalar,PetscMat,PetscScalar,PetscScalar,PetscMat,PetscMatStructure,PetscBool);
    int EPSKrylovSchurSetSubintervals(SlepcEPS,PetscReal*);

    ctypedef enum SlepcEPSLanczosReorthogType "EPSLanczosReorthogType":
        EPS_LANCZOS_REORTHOG_LOCAL
        EPS_LANCZOS_REORTHOG_FULL
        EPS_LANCZOS_REORTHOG_SELECTIVE
        EPS_LANCZOS_REORTHOG_PERIODIC
        EPS_LANCZOS_REORTHOG_PARTIAL
        EPS_LANCZOS_REORTHOG_DELAYED
    int EPSLanczosSetReorthog(SlepcEPS,SlepcEPSLanczosReorthogType)
    int EPSLanczosGetReorthog(SlepcEPS,SlepcEPSLanczosReorthogType*)

    int EPSGDSetKrylovStart(SlepcEPS,PetscBool)
    int EPSGDGetKrylovStart(SlepcEPS,PetscBool*)
    int EPSGDSetBlockSize(SlepcEPS,PetscInt)
    int EPSGDGetBlockSize(SlepcEPS,PetscInt*)
    int EPSGDSetRestart(SlepcEPS,PetscInt,PetscInt)
    int EPSGDGetRestart(SlepcEPS,PetscInt*,PetscInt*)
    int EPSGDSetInitialSize(SlepcEPS,PetscInt)
    int EPSGDGetInitialSize(SlepcEPS,PetscInt*)
    int EPSGDSetBOrth(SlepcEPS,PetscBool)
    int EPSGDGetBOrth(SlepcEPS,PetscBool*)
    int EPSGDSetWindowSizes(SlepcEPS,PetscInt,PetscInt)
    int EPSGDGetWindowSizes(SlepcEPS,PetscInt*,PetscInt*)
    int EPSGDSetDoubleExpansion(SlepcEPS,PetscBool)
    int EPSGDGetDoubleExpansion(SlepcEPS,PetscBool*)

    int EPSJDSetKrylovStart(SlepcEPS,PetscBool)
    int EPSJDGetKrylovStart(SlepcEPS,PetscBool*)
    int EPSJDSetBlockSize(SlepcEPS,PetscInt)
    int EPSJDGetBlockSize(SlepcEPS,PetscInt*)
    int EPSJDSetRestart(SlepcEPS,PetscInt,PetscInt)
    int EPSJDGetRestart(SlepcEPS,PetscInt*,PetscInt*)
    int EPSJDSetInitialSize(SlepcEPS,PetscInt)
    int EPSJDGetInitialSize(SlepcEPS,PetscInt*)
    int EPSJDSetFix(SlepcEPS,PetscReal)
    int EPSJDGetFix(SlepcEPS,PetscReal*)
    int EPSJDSetConstantCorrectionTolerance(SlepcEPS,PetscBool)
    int EPSJDGetConstantCorrectionTolerance(SlepcEPS,PetscBool*)
    int EPSJDSetBOrth(SlepcEPS,PetscBool)
    int EPSJDGetBOrth(SlepcEPS,PetscBool*)
    int EPSJDGetWindowSizes(SlepcEPS,PetscInt*,PetscInt*)
    int EPSJDSetWindowSizes(SlepcEPS,PetscInt,PetscInt)

    int EPSRQCGSetReset(SlepcEPS,PetscInt)
    int EPSRQCGGetReset(SlepcEPS,PetscInt*)

    int EPSCISSSetRegion(SlepcEPS,PetscScalar,PetscReal,PetscReal)
    int EPSCISSGetRegion(SlepcEPS,PetscScalar*,PetscReal*,PetscReal*)
    int EPSCISSSetSizes(SlepcEPS,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscBool)
    int EPSCISSGetSizes(SlepcEPS,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscBool*)
    int EPSCISSSetThreshold(SlepcEPS,PetscReal,PetscReal)
    int EPSCISSGetThreshold(SlepcEPS,PetscReal*,PetscReal*)
    int EPSCISSSetRefinement(SlepcEPS,PetscInt,PetscInt,PetscInt)
    int EPSCISSGetRefinement(SlepcEPS,PetscInt*,PetscInt*,PetscInt*)

cdef extern from * nogil:
    int VecDuplicate(PetscVec,PetscVec*)
    int MatCreateVecs(PetscMat,PetscVec*,PetscVec*)
