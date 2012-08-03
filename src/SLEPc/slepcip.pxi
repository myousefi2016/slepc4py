cdef extern from * nogil:

    ctypedef char* SlepcIPType "const char*"
    SlepcIPType IPBILINEAR
    SlepcIPType IPSESQUILINEAR
    SlepcIPType IPINDEFINITE

    ctypedef enum SlepcIPOrthogType "IPOrthogType":
        IP_ORTHOG_MGS
        IP_ORTHOG_CGS

    ctypedef enum SlepcIPOrthogRefineType "IPOrthogRefineType":
        IP_ORTHOG_REFINE_NEVER
        IP_ORTHOG_REFINE_IFNEEDED
        IP_ORTHOG_REFINE_ALWAYS

    int IPCreate(MPI_Comm,SlepcIP*)
    int IPView(SlepcIP,PetscViewer)
    int IPDestroy(SlepcIP*)
    int IPReset(SlepcIP)
    int IPSetType(SlepcIP,SlepcIPType)
    int IPGetType(SlepcIP,SlepcIPType*)

    int IPSetOptionsPrefix(SlepcIP,char[])
    int IPGetOptionsPrefix(SlepcIP,char*[])
    int IPAppendOptionsPrefix(SlepcIP,char[])
    int IPSetFromOptions(SlepcIP)

    int IPSetOrthogonalization(SlepcIP,SlepcIPOrthogType,SlepcIPOrthogRefineType,PetscReal)
    int IPGetOrthogonalization(SlepcIP,SlepcIPOrthogType*,SlepcIPOrthogRefineType*,PetscReal*)

    int IPSetMatrix(SlepcIP,PetscMat)
    int IPGetMatrix(SlepcIP,PetscMat*)
    int IPApplyMatrix(SlepcIP,PetscVec,PetscVec)

    int IPNorm(SlepcIP,PetscVec,PetscReal*)
    int IPInnerProduct(SlepcIP,PetscVec,PetscVec,PetscReal*)

    int IPOrthogonalize(SlepcIP,PetscInt,PetscVec*,PetscInt,PetscBool*,PetscVec*,PetscVec,PetscScalar*,PetscReal*,PetscBool*)
