# --------------------------------------------------------------------

class IPOrthoType:
    """
    IP Orthogonalization Types.
    """
    CGS = IP_CGS_ORTH
    MGS = IP_MGS_ORTH

class IPRefineType:
    """
    IP Orthogonalization Refinement Types.
    """
    NEVER    = IP_ORTH_REFINE_NEVER
    IFNEEDED = IP_ORTH_REFINE_IFNEEDED
    ALWAYS   = IP_ORTH_REFINE_ALWAYS

class IPBilinearForm:
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

# --------------------------------------------------------------------
