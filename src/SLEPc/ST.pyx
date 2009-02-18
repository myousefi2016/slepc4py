# --------------------------------------------------------------------

class STType:
    """
    ST Types
    """
    SHELL  = STSHELL
    SHIFT  = STSHIFT
    SINV   = STSINV
    CAYLEY = STCAYLEY
    FOLD   = STFOLD

class STMatMode:
    """
    ST Matrix Mode
    """
    COPY    = STMATMODE_COPY
    INPLACE = STMATMODE_INPLACE
    SHELL   = STMATMODE_SHELL

# --------------------------------------------------------------------

cdef class ST(Object):

    """
    ST
    """

    Type    = STType
    MatMode = STMatMode

    def __cinit__(self):
        self.obj = <PetscObject*> &self.st
        self.st = NULL

    def view(self, Viewer viewer=None):
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( STView(self.st, vwr) )

    def destroy(self):
        CHKERR( STDestroy(self.st) )
        self.st = NULL
        return self

    def create(self, comm=None):
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcST newst = NULL
        CHKERR( STCreate(ccomm, &newst) )
        self.dec_ref(); self.st = newst
        return self

    def setType(self, st_type):
        CHKERR( STSetType(self.st, str2cp(st_type)) )

    def getType(self):
        cdef SlepcSTType st_type = NULL
        CHKERR( STGetType(self.st, &st_type) )
        return cp2str(st_type)

    def setOptionsPrefix(self, prefix):
        CHKERR( STSetOptionsPrefix(self.st, str2cp(prefix)) )

    def getOptionsPrefix(self):
        cdef const_char_p prefix = NULL
        CHKERR( STGetOptionsPrefix(self.st, &prefix) )
        return cp2str(prefix)

    def setFromOptions(self):
        CHKERR( STSetFromOptions(self.st) )

    #

    def setShift(self, shift):
        cdef PetscScalar sval = asScalar(shift)
        CHKERR( STSetShift(self.st, sval) )

    def getShift(self):
        cdef PetscScalar sval = 0
        CHKERR( STGetShift(self.st, &sval) )
        return toScalar(sval)

    def setMatMode(self, mode):
        cdef SlepcSTMatMode cmode = mode
        CHKERR( STSetMatMode(self.st, cmode) )

    def getMatMode(self):
        cdef SlepcSTMatMode cmode = STMATMODE_INPLACE
        CHKERR( STGetMatMode(self.st, &cmode) )
        return cmode

    def setOperators(self, Mat A not None, Mat B=None):
        cdef PetscMat Bmat = NULL
        if B is not None: Bmat = B.mat
        CHKERR( STSetOperators(self.st, A.mat, Bmat) )

    def getOperators(self):
        cdef Mat A = Mat()
        cdef Mat B = Mat()
        CHKERR( STGetOperators(self.st, &A.mat, &B.mat) )
        A.inc_ref(); B.inc_ref(); return (A, B)

    def setKSP(self, KSP ksp not None):
        CHKERR( STSetKSP(self.st, ksp.ksp) )

    def getKSP(self):
        cdef KSP ksp = KSP()
        CHKERR( STGetKSP(self.st, &ksp.ksp) )
        ksp.inc_ref(); return ksp

    #

    def setUp(self):
        CHKERR( STSetUp(self.st) )

    def apply(self, Vec x not None, Vec y not None):
        CHKERR( STApply(self.st, x.vec, y.vec) )

    def applyTranspose(self, Vec x not None, Vec y not None):
        CHKERR( STApplyTranspose(self.st, x.vec, y.vec) )


# --------------------------------------------------------------------
