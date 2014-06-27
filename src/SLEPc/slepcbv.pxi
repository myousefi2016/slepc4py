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

    int BVCreate(MPI_Comm,SlepcBV*)
    int BVView(SlepcBV,PetscViewer)
    int BVDestroy(SlepcBV*)
    int BVReset(SlepcBV)
    int BVSetType(SlepcBV,SlepcBVType)
    int BVGetType(SlepcBV,SlepcBVType*)
    int BVSetSizes(SlepcBV,PetscInt,PetscInt,PetscInt)
    int BVSetSizesFromVec(SlepcBV,PetscVec,PetscInt)
    int BVGetSizes(SlepcBV,PetscInt*,PetscInt*,PetscInt*)

    int BVSetOptionsPrefix(SlepcBV,char[])
    int BVGetOptionsPrefix(SlepcBV,char*[])
    int BVAppendOptionsPrefix(SlepcBV,char[])
    int BVSetFromOptions(SlepcBV)

    int BVSetOrthogonalization(SlepcBV,SlepcBVOrthogType,SlepcBVOrthogRefineType,PetscReal)
    int BVGetOrthogonalization(SlepcBV,SlepcBVOrthogType*,SlepcBVOrthogRefineType*,PetscReal*)

    int BVSetMatrix(SlepcBV,PetscMat,PetscBool)
    int BVGetMatrix(SlepcBV,PetscMat*,PetscBool*)
    int BVApplyMatrix(SlepcBV,PetscVec,PetscVec)

    int BVOrthogonalizeVec(SlepcBV,PetscVec,PetscScalar*,PetscReal*,PetscBool*)
