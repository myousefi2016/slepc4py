cdef extern from * nogil:

    ctypedef char* SlepcSTType "const char*"
    SlepcSTType STSHELL
    SlepcSTType STSHIFT
    SlepcSTType STSINVERT
    SlepcSTType STCAYLEY
    SlepcSTType STPRECOND

    ctypedef enum SlepcSTMatMode "STMatMode":
        ST_MATMODE_COPY
        ST_MATMODE_INPLACE
        ST_MATMODE_SHELL

    ctypedef enum  PetscMatStructure "MatStructure":
        MAT_SAME_NONZERO_PATTERN      "SAME_NONZERO_PATTERN"
        MAT_DIFFERENT_NONZERO_PATTERN "DIFFERENT_NONZERO_PATTERN"
        MAT_SUBSET_NONZERO_PATTERN    "SUBSET_NONZERO_PATTERN"
        MAT_SAME_PRECONDITIONER       "SAME_PRECONDITIONER"

    int STView(SlepcST,PetscViewer)
    int STDestroy(SlepcST*)
    int STReset(SlepcST)
    int STCreate(MPI_Comm,SlepcST*)
    int STGetType(SlepcST,SlepcSTType*)
    int STSetType(SlepcST,SlepcSTType)
    int STGetOptionsPrefix(SlepcST,char*[])
    int STSetOptionsPrefix(SlepcST,char[])
    int STAppendOptionsPrefix(SlepcST,char[])
    int STSetFromOptions(SlepcST)

    int STGetShift(SlepcST,PetscScalar*)
    int STSetShift(SlepcST,PetscScalar)

    int STGetKSP(SlepcST,PetscKSP*)
    int STSetKSP(SlepcST,PetscKSP)

    int STGetNumMatrices(SlepcST,PetscInt*)
    int STGetOperators(SlepcST,PetscInt,PetscMat*)
    int STSetOperators(SlepcST,PetscInt,PetscMat*)
    int STSetMatStructure(SlepcST,PetscMatStructure)

    int STGetOperationCounters(SlepcST,PetscInt*,PetscInt*)
    int STResetOperationCounters(SlepcST)

    int STGetMatMode(SlepcST,SlepcSTMatMode*)
    int STSetMatMode(SlepcST,SlepcSTMatMode)

    int STSetUp(SlepcST)
    int STApply(SlepcST,PetscVec,PetscVec)
    int STApplyTranspose(SlepcST,PetscVec,PetscVec)

    int STCayleySetAntishift(SlepcST,PetscScalar)
