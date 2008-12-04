cdef extern from "slepcip.h":

    ctypedef enum SlepcIPOrthogonalizationType "IPOrthogonalizationType":
        IP_MGS_ORTH
        IP_CGS_ORTH

    ctypedef enum SlepcIPOrthogonalizationRefinementType "IPOrthogonalizationRefinementType":
        IP_ORTH_REFINE_NEVER
        IP_ORTH_REFINE_IFNEEDED
        IP_ORTH_REFINE_ALWAYS

    ctypedef enum SlepcIPBilinearForm "IPBilinearForm":
        IPINNER_HERMITIAN
        IPINNER_SYMMETRIC

    int IPCreate(MPI_Comm,SlepcIP*)
    int IPView(SlepcIP,PetscViewer)
    int IPDestroy(SlepcIP)

    int IPSetOptionsPrefix(SlepcIP,char[])
    int IPGetOptionsPrefix(SlepcIP,char*[])
    int IPAppendOptionsPrefix(SlepcIP,char[])
    int IPSetFromOptions(SlepcIP)
