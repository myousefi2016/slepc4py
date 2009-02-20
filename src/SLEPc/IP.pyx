# --------------------------------------------------------------------

class IPOrthoType(object):
    """
    IP orthogonalization types
    """
    CGS = IP_CGS_ORTH
    MGS = IP_MGS_ORTH

class IPRefineType(object):
    """
    IP orthogonalization refinement types
    """
    NEVER    = IP_ORTH_REFINE_NEVER
    IFNEEDED = IP_ORTH_REFINE_IFNEEDED
    ALWAYS   = IP_ORTH_REFINE_ALWAYS

class IPBilinearForm(object):
    """
    IP bilinear form types
    """
    HERMITIAN = IPINNER_HERMITIAN
    SYMMETRIC = IPINNER_SYMMETRIC

# --------------------------------------------------------------------

cdef class IP(Object):

    """
    IP
    """

    OrthoType    = IPOrthoType
    RefineType   = IPRefineType
    BilinearForm = IPBilinearForm

    def __cinit__(self):
        self.obj = <PetscObject*> &self.ip
        self.ip = NULL

    def view(self, Viewer viewer=None):
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( IPView(self.ip, vwr) )

    def destroy(self):
        CHKERR( IPDestroy(self.ip) )
        self.ip = NULL
        return self

    def create(self, comm=None):
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcIP newip = NULL
        CHKERR( IPCreate(ccomm, &newip) )
        self.dec_ref(); self.ip = newip
        return self

    def setOptionsPrefix(self, prefix):
        CHKERR( IPSetOptionsPrefix(self.ip, str2cp(prefix)) )

    def getOptionsPrefix(self):
        cdef const_char_p prefix = NULL
        CHKERR( IPGetOptionsPrefix(self.ip, &prefix) )
        return cp2str(prefix)

    def setFromOptions(self):
        CHKERR( IPSetFromOptions(self.ip) )

    #

    def getOrthogonalization(self):
        cdef SlepcIPOrthogonalizationType otype= IP_CGS_ORTH
        cdef SlepcIPOrthogonalizationRefinementType rtype = IP_ORTH_REFINE_IFNEEDED
        cdef PetscReal rval = PETSC_DEFAULT
        CHKERR( IPGetOrthogonalization(self.ip, &otype, &rtype, &rval) )
        return (otype, rtype, rval)

    def setOrthogonalization(self, type=None, refine=None, eta=None):
        cdef SlepcIPOrthogonalizationType otype= IP_CGS_ORTH
        cdef SlepcIPOrthogonalizationRefinementType rtype = IP_ORTH_REFINE_IFNEEDED
        cdef PetscReal rval = PETSC_DEFAULT
        if type   is not None: otype= type
        if refine is not None: rtype= refine
        if eta    is not None: rval = eta
        CHKERR( IPSetOrthogonalization(self.ip, otype, rtype, rval) )

    #

    def getBilinearForm(self):
        cdef Mat mat = Mat()
        cdef PetscMat m = NULL
        cdef SlepcIPBilinearForm f = IPINNER_HERMITIAN
        CHKERR( IPGetBilinearForm(self.ip, &m, &f) )
        mat.mat = m; mat.inc_ref()
        return (mat, f)

    def setBilinearForm(self, Mat mat=None, form=None):
        cdef PetscMat m = NULL
        cdef SlepcIPBilinearForm f = IPINNER_HERMITIAN
        if mat  is not None: m = mat.mat
        if form is not None: f = form
        CHKERR( IPSetBilinearForm(self.ip, m, f) )

    def applyMatrix(self, Vec x not None, Vec y not None):
        CHKERR( IPApplyMatrix(self.ip, x.vec, y.vec) )

    #

    def norm(self, Vec x not None):
        cdef PetscReal rval = 0
        CHKERR( IPNorm(self.ip, x.vec, &rval) )
        return rval

    def innerProduct(self, Vec x not None, Vec y not None):
        cdef PetscReal rval = 0
        CHKERR( IPInnerProduct(self.ip, x.vec, y.vec, &rval) )
        return rval

    def orthogonalize(self, VS, Vec v not None, Vec work=None):
        cdef PetscInt i = 0
        cdef PetscInt n = 0
        cdef PetscTruth* which = NULL
        cdef PetscVec* V = NULL
        cdef PetscScalar* H = NULL, h = 0
        cdef PetscReal norm = 0
        cdef PetscTruth lindep = PETSC_FALSE
        cdef PetscVec w = NULL
        cdef PetscScalar* sw = NULL
        cdef object tmp1 = None, tmp2 = None
        if isinstance(VS, Vec):
            n = 1
            V = &((<Vec>VS).vec)
            H = &h
        else:
            n = len(VS)
            tmp1 = allocate(n*sizeof(Vec),<void**>&V)
            tmp2 = allocate(n*sizeof(PetscScalar),<void**>&H)
            for i in range(n):
                V[i] = (<Vec?>VS[i]).vec
                H[i] = 0
        if work is not None: w = work.vec
        CHKERR( IPOrthogonalize(self.ip,
                                n, which, V, v.vec,
                                H, &norm, &lindep,
                                w, sw) )
        cdef object coefs = None
        if isinstance(VS, Vec):
            coefs = toScalar(H[0])
        else:
            coefs = [toScalar(H[i]) for i in range(n)]
        return (coefs, norm, <bint>lindep)


# --------------------------------------------------------------------
