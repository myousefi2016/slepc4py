# --------------------------------------------------------------------

class SVDType(object):
    """
    SVD Types
    """
    CROSS     = SVDCROSS
    CYCLIC    = SVDCYCLIC
    LAPACK    = SVDLAPACK
    LANCZOS   = SVDLANCZOS
    TRLANCZOS = SVDTRLANCZOS

class SVDWhich(object):
    """
    SVD desired piece of spectrum
    """
    LARGEST  = SVD_LARGEST
    SMALLEST = SVD_SMALLEST

class SVDTransposeMode(object):
    """
    SVD handling of the transpose of the matrix
    """
    EXPLICIT = SVD_TRANSPOSE_EXPLICIT
    IMPLICIT = SVD_TRANSPOSE_IMPLICIT

class SVDConvergedReason(object):
    """
    SVD convergence reasons
    """
    CONVERGED_TOL       = SVD_CONVERGED_TOL
    DIVERGED_ITS        = SVD_DIVERGED_ITS
    DIVERGED_BREAKDOWN  = SVD_DIVERGED_BREAKDOWN
    CONVERGED_ITERATING = SVD_CONVERGED_ITERATING

# --------------------------------------------------------------------

cdef class SVD(Object):

    """
    SVD
    """

    Type            = SVDType
    Which           = SVDWhich
    TransposeMode   = SVDTransposeMode
    ConvergedReason = SVDConvergedReason

    def __cinit__(self):
        self.obj = <PetscObject*> &self.svd
        self.svd = NULL

    def view(self, Viewer viewer=None):
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( SVDView(self.svd, vwr) )

    def destroy(self):
        CHKERR( SVDDestroy(self.svd) )
        self.svd = NULL
        return self

    def create(self, comm=None):
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcSVD newsvd = NULL
        CHKERR( SVDCreate(ccomm, &newsvd) )
        self.dec_ref(); self.svd = newsvd
        return self

    def setType(self, svd_type):
        CHKERR( SVDSetType(self.svd, str2cp(svd_type)) )

    def getType(self):
        cdef SlepcSVDType svd_type = NULL
        CHKERR( SVDGetType(self.svd, &svd_type) )
        return cp2str(svd_type)

    def setOptionsPrefix(self, prefix):
        CHKERR( SVDSetOptionsPrefix(self.svd, str2cp(prefix)) )

    def getOptionsPrefix(self):
        cdef const_char_p prefix = NULL
        CHKERR( SVDGetOptionsPrefix(self.svd, &prefix) )
        return cp2str(prefix)

    def setFromOptions(self):
        CHKERR( SVDSetFromOptions(self.svd) )

    #

    def getTransposeMode(self):
        cdef SlepcSVDTransposeMode val
        CHKERR( SVDGetTransposeMode(self.svd, &val) )
        return val

    def setTransposeMode(self, mode):
        cdef SlepcSVDTransposeMode val = mode
        CHKERR( SVDSetTransposeMode(self.svd, val) )

    def getWhichSingularTriplets(self):
        cdef SlepcSVDWhich val
        CHKERR( SVDGetWhichSingularTriplets(self.svd, &val) )
        return val

    def setWhichSingularTriplets(self, which):
        cdef SlepcSVDWhich val = which
        CHKERR( SVDSetWhichSingularTriplets(self.svd, val) )
    #

    def getTolerances(self):
        cdef PetscReal ctol = 0
        cdef PetscInt cmaxit = 0
        CHKERR( SVDGetTolerances(self.svd, &ctol, &cmaxit) )
        return (ctol, cmaxit)

    def setTolerances(self, tol=None, max_it=None):
        cdef PetscReal ctol = PETSC_IGNORE
        cdef PetscInt cmaxit = PETSC_IGNORE
        if tol    is not None: ctol   = tol
        if max_it is not None: cmaxit = max_it
        CHKERR( SVDSetTolerances(self.svd, ctol, cmaxit) )

    def getDimensions(self):
        cdef PetscInt nev = 0
        cdef PetscInt ncv = 0
        cdef PetscInt mpd = 0
        CHKERR( SVDGetDimensions(self.svd, &nev, &ncv, &mpd) )
        return (nev, ncv, mpd)

    def setDimensions(self, nev=None, ncv=None, mpd=None):
        cdef PetscInt cnev = PETSC_IGNORE
        cdef PetscInt cncv = PETSC_IGNORE
        cdef PetscInt cmpd = PETSC_IGNORE
        if nev is not None: cnev = nev
        if ncv is not None: cncv = ncv
        if mpd is not None: cmpd = mpd
        CHKERR( SVDSetDimensions(self.svd, cnev, cncv, cmpd) )

    def getIP(self):
        cdef IP ip = IP()
        CHKERR( SVDGetIP(self.svd, &ip.ip) )
        ip.inc_ref(); return ip

    def setIP(self, IP ip not None):
        CHKERR( SVDSetIP(self.svd, ip.ip) )

    def getOperator(self):
        cdef Mat A = Mat()
        CHKERR( SVDGetOperator(self.svd, &A.mat) )
        A.inc_ref(); return A

    def setOperator(self, Mat A not None):
        CHKERR( SVDSetOperator(self.svd, A.mat) )

    #

    def getInitialVector(self):
        cdef Vec V = Vec()
        CHKERR( SVDGetInitialVector(self.svd, &V.vec) )
        V.inc_ref(); return V

    def setInitialVector(self, Vec V not None):
        CHKERR( SVDSetInitialVector(self.svd, V.vec) )

    #

    def setUp(self):
        CHKERR( SVDSetUp(self.svd) )

    def solve(self):
        CHKERR( SVDSolve(self.svd) )

    def getIterationNumber(self):
        cdef PetscInt its = 0
        CHKERR( SVDGetIterationNumber(self.svd, &its) )
        return its

    def getConvergedReason(self):
        cdef SlepcSVDConvergedReason reason
        reason = SVD_CONVERGED_ITERATING
        CHKERR( SVDGetConvergedReason(self.svd, &reason) )
        return reason

    def getConverged(self):
        cdef PetscInt nconv = 0
        CHKERR( SVDGetConverged(self.svd, &nconv) )
        return nconv

    def getValue(self, int i):
        cdef PetscReal sigma = 0
        CHKERR( SVDGetSingularTriplet(self.svd, i, &sigma, NULL, NULL) )
        return sigma

    def getVectors(self, int i, Vec U not None, Vec V not None):
        cdef PetscReal sigma = 0
        CHKERR( SVDGetSingularTriplet(self.svd, i, &sigma, U.vec, V.vec) )

    def getSingularTriplet(self, i, Vec U=None, Vec V=None):
        cdef PetscReal sigma=0
        cdef PetscVec Uvec = NULL
        cdef PetscVec Vvec = NULL
        if U is not None: Uvec = U.vec
        if V is not None: Vvec = V.vec
        CHKERR( SVDGetSingularTriplet(self.svd, i, &sigma, Uvec, Vvec) )
        return sigma

    #

    def computeRelativeError(self, int i):
        cdef PetscReal val = 0
        CHKERR( SVDComputeRelativeError(self.svd, i, &val) )
        return val

    def computeResidualNorms(self, int i):
        cdef PetscReal val1=0,
        cdef PetscReal val2=0
        CHKERR( SVDComputeResidualNorms(self.svd, i, &val1, &val2) )
        return (val1, val2)

    #

    property transpose_mode:
        def __get__(self):
            return self.getTransposeMode()
        def __set__(self, value):
            self.setTransposeMode(value)

    property which:
        def __get__(self):
            return self.getWhichSingularTriplets()
        def __set__(self, value):
            self.setWhichSingularTriplets(value)

    property tol:
        def __get__(self):
            return self.getTolerances()[0]
        def __set__(self, value):
            self.setTolerances(tol=value)

    property max_it:
        def __get__(self):
            return self.getTolerances()[1]
        def __set__(self, value):
            self.setTolerances(max_it=value)

    property ip:
        def __get__(self):
            return self.getIP()
        def __set__(self, value):
            self.setIP(value)

# --------------------------------------------------------------------
