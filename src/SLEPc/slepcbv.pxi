cdef extern from * nogil:

    ctypedef char* SlepcBVType "const char*"
    SlepcBVType BVMAT
    SlepcBVType BVSVEC
    SlepcBVType BVVECS
    SlepcBVType BVCONTIGUOUS

    ctypedef enum SlepcBVOrthogType "BVOrthogType":
        BV_ORTHOG_CGS
        BV_ORTHOG_MGS

    ctypedef enum SlepcBVOrthogRefineType "BVOrthogRefineType":
        BV_ORTHOG_REFINE_IFNEEDED
        BV_ORTHOG_REFINE_NEVER
        BV_ORTHOG_REFINE_ALWAYS

    ctypedef enum SlepcBVOrthogBlockType "BVOrthogBlockType":
        BV_ORTHOG_BLOCK_GS
        BV_ORTHOG_BLOCK_CHOL

    int BVCreate(MPI_Comm,SlepcBV*)
    int BVDuplicate(SlepcBV,SlepcBV*)
    int BVCopy(SlepcBV,SlepcBV)
    int BVView(SlepcBV,PetscViewer)
    int BVDestroy(SlepcBV*)
    int BVSetType(SlepcBV,SlepcBVType)
    int BVGetType(SlepcBV,SlepcBVType*)
    int BVSetSizes(SlepcBV,PetscInt,PetscInt,PetscInt)
    int BVSetSizesFromVec(SlepcBV,PetscVec,PetscInt)
    int BVGetSizes(SlepcBV,PetscInt*,PetscInt*,PetscInt*)

    int BVSetOptionsPrefix(SlepcBV,char[])
    int BVGetOptionsPrefix(SlepcBV,char*[])
    int BVAppendOptionsPrefix(SlepcBV,char[])
    int BVSetFromOptions(SlepcBV)

    int BVSetOrthogonalization(SlepcBV,SlepcBVOrthogType,SlepcBVOrthogRefineType,PetscReal,SlepcBVOrthogBlockType)
    int BVGetOrthogonalization(SlepcBV,SlepcBVOrthogType*,SlepcBVOrthogRefineType*,PetscReal*,SlepcBVOrthogBlockType*)

    int BVSetRandom(SlepcBV)

    int BVSetMatrix(SlepcBV,PetscMat,PetscBool)
    int BVGetMatrix(SlepcBV,PetscMat*,PetscBool*)
    int BVApplyMatrix(SlepcBV,PetscVec,PetscVec)

    int BVSetActiveColumns(SlepcBV,PetscInt,PetscInt)
    int BVGetActiveColumns(SlepcBV,PetscInt*,PetscInt*)

    int BVInsertVec(SlepcBV,PetscInt,PetscVec)
    int BVInsertVecs(SlepcBV,PetscInt,PetscInt*,PetscVec*,PetscBool)
    int BVGetColumn(SlepcBV,PetscInt,PetscVec*)
    int BVRestoreColumn(SlepcBV,PetscInt,PetscVec*)

    int BVDot(SlepcBV,SlepcBV,PetscMat)
    int BVDotVec(SlepcBV,PetscVec,PetscScalar*)

    int BVMatProject(SlepcBV,PetscMat,SlepcBV,PetscMat)
    int BVMatMult(SlepcBV,PetscMat,SlepcBV)
    int BVMatMultHermitianTranspose(SlepcBV,PetscMat,SlepcBV)
    int BVMultVec(SlepcBV,PetscScalar,PetscScalar,PetscVec,PetscScalar*)

    int BVScaleColumn(SlepcBV,PetscInt,PetscScalar)
    int BVScale(SlepcBV,PetscScalar)

    int BVNormColumn(SlepcBV,PetscInt,PetscNormType,PetscReal*)
    int BVNorm(SlepcBV,PetscNormType,PetscReal*)

    int BVOrthogonalizeVec(SlepcBV,PetscVec,PetscScalar*,PetscReal*,PetscBool*)
    int BVOrthogonalize(SlepcBV,PetscMat)


cdef inline int BV_Sizes(
    object size,
    PetscInt *_n,
    PetscInt *_N,
    ) except -1:
    # unpack and get local and global sizes
    cdef PetscInt n=PETSC_DECIDE, N=PETSC_DECIDE
    cdef object on, oN
    try:
        on, oN = size
    except (TypeError, ValueError):
        on = None; oN = size
    if on is not None: n = asInt(on)
    if oN is not None: N = asInt(oN)
    if n==PETSC_DECIDE and N==PETSC_DECIDE: raise ValueError(
        "local and global sizes cannot be both 'DECIDE'")
    # return result to the caller
    if _n != NULL: _n[0] = n
    if _N != NULL: _N[0] = N
    return 0
