cdef extern from * nogil:

    ctypedef char* SlepcDSType "const char*"
    SlepcDSType DSHEP
    SlepcDSType DSNHEP
    SlepcDSType DSGHEP
    SlepcDSType DSGHIEP
    SlepcDSType DSGNHEP
    SlepcDSType DSSVD
    SlepcDSType DSPEP
    SlepcDSType DSNEP

    ctypedef enum SlepcDSStateType "DSStateType":
        DS_STATE_RAW
        DS_STATE_INTERMEDIATE
        DS_STATE_CONDENSED
        DS_STATE_TRUNCATED

    ctypedef enum SlepcDSMatType "DSMatType":
        DS_MAT_A
        DS_MAT_B
        DS_MAT_C
        DS_MAT_T
        DS_MAT_D
        DS_MAT_Q
        DS_MAT_Z
        DS_MAT_X
        DS_MAT_Y
        DS_MAT_U
        DS_MAT_VT
        DS_MAT_W
        DS_NUM_MAT

    int DSCreate(MPI_Comm,SlepcDS*)
    int DSView(SlepcDS,PetscViewer)
    int DSDestroy(SlepcDS*)
    int DSReset(SlepcDS)
    int DSSetType(SlepcDS,SlepcDSType)
    int DSGetType(SlepcDS,SlepcDSType*)

    int DSSetOptionsPrefix(SlepcDS,char[])
    int DSGetOptionsPrefix(SlepcDS,char*[])
    int DSAppendOptionsPrefix(SlepcDS,char[])
    int DSSetFromOptions(SlepcDS)

    int DSAllocate(SlepcDS,PetscInt)
    int DSGetLeadingDimension(SlepcDS,PetscInt*)
    int DSSetState(SlepcDS,SlepcDSStateType)
    int DSGetState(SlepcDS,SlepcDSStateType*)
    int DSSetDimensions(SlepcDS,PetscInt,PetscInt,PetscInt,PetscInt)
    int DSGetDimensions(SlepcDS,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*)
    int DSTruncate(SlepcDS,PetscInt)
    int DSSetMethod(SlepcDS,PetscInt)
    int DSGetMethod(SlepcDS,PetscInt*)
    int DSSetCompact(SlepcDS,PetscBool)
    int DSGetCompact(SlepcDS,PetscBool*)
    int DSSetExtraRow(SlepcDS,PetscBool)
    int DSGetExtraRow(SlepcDS,PetscBool*)
    int DSSetRefined(SlepcDS,PetscBool)
    int DSGetRefined(SlepcDS,PetscBool*)
    int DSGetArray(SlepcDS,SlepcDSMatType,PetscScalar *a[])
    int DSRestoreArray(SlepcDS,SlepcDSMatType,PetscScalar *a[])
    int DSGetArrayReal(SlepcDS,SlepcDSMatType,PetscReal *a[])
    int DSRestoreArrayReal(SlepcDS,SlepcDSMatType,PetscReal *a[])
    int DSVectors(SlepcDS,SlepcDSMatType,PetscInt*,PetscReal*)
    int DSSolve(SlepcDS,PetscScalar*,PetscScalar*)
    int DSSort(SlepcDS,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt*)
    int DSUpdateExtraRow(SlepcDS)
    int DSCond(SlepcDS,PetscReal*)
    int DSTranslateHarmonic(SlepcDS,PetscScalar,PetscReal,PetscBool,PetscScalar*,PetscReal*)
    int DSTranslateRKS(SlepcDS,PetscScalar)
    int DSNormalize(SlepcDS,SlepcDSMatType,PetscInt)
