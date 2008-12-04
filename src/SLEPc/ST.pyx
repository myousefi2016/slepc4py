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
        cdef MPI_Comm ccomm = def_Comm(comm, PETSC_COMM_DEFAULT)
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

# --------------------------------------------------------------------
