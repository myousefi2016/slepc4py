cdef extern from * nogil:

    ctypedef char* SlepcRGType "const char*"
    SlepcRGType RGINTERVAL
    SlepcRGType RGPOLYGON
    SlepcRGType RGELLIPSE
    SlepcRGType RGRING

    int RGCreate(MPI_Comm,SlepcRG*)
    int RGView(SlepcRG,PetscViewer)
    int RGDestroy(SlepcRG*)
    int RGSetType(SlepcRG,SlepcRGType)
    int RGGetType(SlepcRG,SlepcRGType*)

    int RGSetOptionsPrefix(SlepcRG,char[])
    int RGGetOptionsPrefix(SlepcRG,char*[])
    int RGAppendOptionsPrefix(SlepcRG,char[])
    int RGSetFromOptions(SlepcRG)

    int RGIsTrivial(SlepcRG,PetscBool*)
    int RGSetComplement(SlepcRG,PetscBool)
    int RGGetComplement(SlepcRG,PetscBool*)

    int RGEllipseSetParameters(SlepcRG,PetscScalar,PetscReal,PetscReal)
    int RGEllipseGetParameters(SlepcRG,PetscScalar*,PetscReal*,PetscReal*)
    int RGIntervalSetEndpoints(SlepcRG,PetscReal,PetscReal,PetscReal,PetscReal)
    int RGIntervalGetEndpoints(SlepcRG,PetscReal*,PetscReal*,PetscReal*,PetscReal*)

