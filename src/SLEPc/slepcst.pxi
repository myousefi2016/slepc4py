cdef extern from "slepcst.h":

    ctypedef char* SlepcSTType "const char*"
    SlepcSTType STSHELL
    SlepcSTType STSHIFT
    SlepcSTType STSINV
    SlepcSTType STCAYLEY
    SlepcSTType STFOLD

    ctypedef enum SlepcSTMatMode "STMatMode":
        STMATMODE_COPY
        STMATMODE_INPLACE
        STMATMODE_SHELL

    int STView(SlepcST,PetscViewer)
    int STDestroy(SlepcST)
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

    int STGetOperators(SlepcST,PetscMat*,PetscMat*)
    int STSetOperators(SlepcST,PetscMat,PetscMat)

    int STGetOperationCounters(SlepcST,PetscInt*,PetscInt*)
    int STResetOperationCounters(SlepcST)

    int STGetMatMode(SlepcST,SlepcSTMatMode*)
    int STSetMatMode(SlepcST,SlepcSTMatMode)

    int STSetUp(SlepcST)
    int STApply(SlepcST,PetscVec,PetscVec)
    int STApplyTranspose(SlepcST,PetscVec,PetscVec)
