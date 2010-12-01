cdef extern from * nogil:

    ctypedef enum SlepcIPOrthogonalizationType "IPOrthogonalizationType":
        IP_ORTH_MGS
        IP_ORTH_CGS

    ctypedef enum SlepcIPOrthogonalizationRefinementType "IPOrthogonalizationRefinementType":
        IP_ORTH_REFINE_NEVER
        IP_ORTH_REFINE_IFNEEDED
        IP_ORTH_REFINE_ALWAYS

    ctypedef enum SlepcIPBilinearForm "IPBilinearForm":
        IP_INNER_HERMITIAN
        IP_INNER_SYMMETRIC

    int IPCreate(MPI_Comm,SlepcIP*)
    int IPView(SlepcIP,PetscViewer)
    int IPDestroy(SlepcIP)

    int IPSetOptionsPrefix(SlepcIP,char[])
    int IPGetOptionsPrefix(SlepcIP,char*[])
    int IPAppendOptionsPrefix(SlepcIP,char[])
    int IPSetFromOptions(SlepcIP)

    int IPSetOrthogonalization(SlepcIP,SlepcIPOrthogonalizationType,SlepcIPOrthogonalizationRefinementType,PetscReal)
    int IPGetOrthogonalization(SlepcIP,SlepcIPOrthogonalizationType*,SlepcIPOrthogonalizationRefinementType*,PetscReal*)

    int IPSetBilinearForm(SlepcIP,PetscMat,SlepcIPBilinearForm)
    int IPGetBilinearForm(SlepcIP,PetscMat*,SlepcIPBilinearForm*)
    int IPApplyMatrix(SlepcIP,PetscVec,PetscVec)

    int IPNorm(SlepcIP,PetscVec,PetscReal*)
    int IPInnerProduct(SlepcIP,PetscVec,PetscVec,PetscReal*)

    int IPOrthogonalize(SlepcIP,PetscInt,PetscVec*,PetscInt,PetscBool*,PetscVec*,PetscVec,PetscScalar*,PetscReal*,PetscBool*)
