# --------------------------------------------------------------------

class SVDType(object):
    """
    SVD types
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
        cdef PetscReal rval = 0
        cdef PetscInt  ival = 0
        CHKERR( SVDGetTolerances(self.svd, &rval, &ival) )
        return (rval, ival)

    def setTolerances(self, tol=None, max_it=None):
        cdef PetscReal rval = PETSC_IGNORE
        cdef PetscInt  ival = PETSC_IGNORE
        if tol    is not None: rval = tol
        if max_it is not None: ival = max_it
        CHKERR( SVDSetTolerances(self.svd, rval, ival) )

    def getDimensions(self):
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        cdef PetscInt ival3 = 0
        CHKERR( SVDGetDimensions(self.svd, &ival1, &ival2, &ival3) )
        return (ival1, ival2, ival3)

    def setDimensions(self, nev=None, ncv=None, mpd=None):
        cdef PetscInt ival1 = PETSC_IGNORE
        cdef PetscInt ival2 = PETSC_IGNORE
        cdef PetscInt ival3 = PETSC_IGNORE
        if nev is not None: ival1 = nev
        if ncv is not None: ival2 = ncv
        if mpd is not None: ival3 = mpd
        CHKERR( SVDSetDimensions(self.svd, ival1, ival2, ival2) )

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
        cdef PetscInt ival = 0
        CHKERR( SVDGetIterationNumber(self.svd, &ival) )
        return ival

    def getConvergedReason(self):
        cdef SlepcSVDConvergedReason reason
        reason = SVD_CONVERGED_ITERATING
        CHKERR( SVDGetConvergedReason(self.svd, &reason) )
        return reason

    def getConverged(self):
        cdef PetscInt ival = 0
        CHKERR( SVDGetConverged(self.svd, &ival) )
        return ival

    def getValue(self, int i):
        cdef PetscReal rval = 0
        CHKERR( SVDGetSingularTriplet(self.svd, i, &rval, NULL, NULL) )
        return rval

    def getVectors(self, int i, Vec U not None, Vec V not None):
        cdef PetscReal rval = 0
        CHKERR( SVDGetSingularTriplet(self.svd, i, &rval, U.vec, V.vec) )

    def getSingularTriplet(self, int i, Vec U=None, Vec V=None):
        cdef PetscReal rval = 0
        cdef PetscVec Uvec = NULL
        cdef PetscVec Vvec = NULL
        if U is not None: Uvec = U.vec
        if V is not None: Vvec = V.vec
        CHKERR( SVDGetSingularTriplet(self.svd, i, &rval, Uvec, Vvec) )
        return rval

    #

    def computeRelativeError(self, int i):
        cdef PetscReal rval = 0
        CHKERR( SVDComputeRelativeError(self.svd, i, &rval) )
        return rval

    def computeResidualNorms(self, int i):
        cdef PetscReal rval1 = 0
        cdef PetscReal rval2 = 0
        CHKERR( SVDComputeResidualNorms(self.svd, i, &rval1, &rval2) )
        return (rval1, rval2)

    def getOperationCounters(self):
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        CHKERR( SVDGetOperationCounters(self.svd, &ival1, &ival2) )
        return (ival1, ival2)

    #

    def setCrossEPS(self, EPS eps not None):
        CHKERR( SVDCrossSetEPS(self.svd, eps.eps) )

    def getCrossEPS(self):
        cdef EPS eps = EPS()
        CHKERR( SVDCrossGetEPS(self.svd, &eps.eps) )
        eps.inc_ref(); return eps

    def setCyclicEPS(self, EPS eps not None):
        CHKERR( SVDCyclicSetEPS(self.svd, eps.eps) )

    def getCyclicEPS(self):
        cdef EPS eps = EPS()
        CHKERR( SVDCyclicGetEPS(self.svd, &eps.eps) )
        eps.inc_ref(); return eps

    def setCyclicExplicitMatrix(self, flag=True):
        cdef PetscTruth tval = PETSC_FALSE
        if flag: tval = PETSC_TRUE
        CHKERR( SVDCyclicSetExplicitMatrix(self.svd, tval) )

    def getCyclicExplicitMatrix(self):
        cdef PetscTruth tval = PETSC_FALSE
        CHKERR( SVDCyclicGetExplicitMatrix(self.svd, &tval) )
        return <bint>tval

    def setLanczosOneSide(self, flag=True):
        cdef PetscTruth tval = PETSC_FALSE
        if flag: tval = PETSC_TRUE
        CHKERR( SVDLanczosSetOneSide(self.svd, tval) )

    def setTRLanczosOneSide(self, flag=True):
        cdef PetscTruth tval = PETSC_FALSE
        if flag: tval = PETSC_TRUE
        CHKERR( SVDLanczosSetOneSide(self.svd, tval) )

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
