# --------------------------------------------------------------------

class EPSType:
    """
    EPS Types
    """
    # provided implementations
    KRYLOVSCHUR = EPSKRYLOVSCHUR
    LANCZOS     = EPSLANCZOS
    ARNOLDI     = EPSARNOLDI
    SUBSPACE    = EPSSUBSPACE
    POWER       = EPSPOWER
    LAPACK      = EPSLAPACK
    # with external libraries
    ARPACK      = EPSARPACK
    BLZPACK     = EPSBLZPACK
    TRLAN       = EPSTRLAN
    BLOPEX      = EPSBLOPEX
    PRIMME      = EPSPRIMME

class EPSProblemType:
    """
    EPS problem type
    """
    HEP    = EPS_HEP
    NHEP   = EPS_NHEP
    GHEP   = EPS_GHEP
    GNHEP  = EPS_GNHEP
    PGNHEP = EPS_PGNHEP

class EPSExtraction:
    """
    EPS extraction technique
    """
    RITZ             = EPS_RITZ
    HARMONIC         = EPS_HARMONIC
    REFINED          = EPS_REFINED
    REFINED_HARMONIC = EPS_REFINED_HARMONIC

class EPSClass:
    """
    EPS class of method
    """
    ONE_SIDE  = EPS_ONE_SIDE
    TWO_SIDE  = EPS_TWO_SIDE

class EPSWhich:
    """
    EPS desired piece of spectrum
    """
    LARGEST_MAGNITUDE  = EPS_LARGEST_MAGNITUDE
    SMALLEST_MAGNITUDE = EPS_SMALLEST_MAGNITUDE
    LARGEST_REAL       = EPS_LARGEST_REAL
    SMALLEST_REAL      = EPS_SMALLEST_REAL
    LARGEST_IMAGINARY  = EPS_LARGEST_IMAGINARY
    SMALLEST_IMAGINARY = EPS_SMALLEST_IMAGINARY

class EPSConvergedReason:
    """
    EPS convergence reasons
    """
    CONVERGED_TOL         = EPS_CONVERGED_TOL
    DIVERGED_ITS          = EPS_DIVERGED_ITS
    DIVERGED_BREAKDOWN    = EPS_DIVERGED_BREAKDOWN
    DIVERGED_NONSYMMETRIC = EPS_DIVERGED_NONSYMMETRIC
    CONVERGED_ITERATING   = EPS_CONVERGED_ITERATING

class EPSPowerShiftType:
    """
    EPS type of power shift
    """
    CONSTANT  = EPSPOWER_SHIFT_CONSTANT
    RAYLEIGH  = EPSPOWER_SHIFT_RAYLEIGH
    WILKINSON = EPSPOWER_SHIFT_WILKINSON

# --------------------------------------------------------------------

cdef class EPS(Object):

    """
    EPS
    """

    Type            = EPSType
    ProblemType     = EPSProblemType
    Extraction      = EPSExtraction
    Class           = EPSClass
    Which           = EPSWhich
    ConvergedReason = EPSConvergedReason

    PowerShiftType  = EPSPowerShiftType

    def __cinit__(self):
        self.obj = <PetscObject*> &self.eps
        self.eps = NULL

    def view(self, Viewer viewer=None):
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( EPSView(self.eps, vwr) )

    def destroy(self):
        CHKERR( EPSDestroy(self.eps) )
        self.eps = NULL
        return self

    def create(self, comm=None):
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcEPS neweps = NULL
        CHKERR( EPSCreate(ccomm, &neweps) )
        self.dec_ref(); self.eps = neweps
        return self

    def setType(self, eps_type):
        CHKERR( EPSSetType(self.eps, str2cp(eps_type)) )

    def getType(self):
        cdef SlepcEPSType eps_type = NULL
        CHKERR( EPSGetType(self.eps, &eps_type) )
        return cp2str(eps_type)

    def setOptionsPrefix(self, prefix):
        CHKERR( EPSSetOptionsPrefix(self.eps, str2cp(prefix)) )

    def getOptionsPrefix(self):
        cdef const_char_p prefix = NULL
        CHKERR( EPSGetOptionsPrefix(self.eps, &prefix) )
        return cp2str(prefix)

    def setFromOptions(self):
        CHKERR( EPSSetFromOptions(self.eps) )

    #

    def getProblemType(self):
        cdef SlepcEPSProblemType val
        CHKERR( EPSGetProblemType(self.eps, &val) )
        return val

    def setProblemType(self, problem_type):
        cdef SlepcEPSProblemType val = problem_type
        CHKERR( EPSSetProblemType(self.eps, val) )

    def getClass(self):
        cdef SlepcEPSClass val
        CHKERR( EPSGetClass(self.eps, &val) )
        return val

    def setClass(self, klass):
        cdef SlepcEPSClass val = klass
        CHKERR( EPSSetClass(self.eps, val) )

    def getExtraction(self):
        cdef SlepcEPSExtraction val = EPS_RITZ
        CHKERR( EPSGetExtraction(self.eps, &val) )
        return val

    def setExtraction(self, extraction):
        cdef SlepcEPSExtraction val = extraction
        CHKERR( EPSSetExtraction(self.eps, val) )

    def getWhichEigenpairs(self):
        cdef SlepcEPSWhich val
        CHKERR( EPSGetWhichEigenpairs(self.eps, &val) )
        return val

    def setWhichEigenpairs(self, which):
        cdef SlepcEPSWhich val = which
        CHKERR( EPSSetWhichEigenpairs(self.eps, val) )

    #

    def getTolerances(self):
        cdef PetscReal ctol = 0
        cdef PetscInt cmaxit = 0
        CHKERR( EPSGetTolerances(self.eps, &ctol, &cmaxit) )
        return (ctol, cmaxit)

    def setTolerances(self, tol=None, max_it=None):
        cdef PetscReal ctol = PETSC_IGNORE
        cdef PetscInt cmaxit = PETSC_IGNORE
        if tol    is not None: ctol   = tol
        if max_it is not None: cmaxit = max_it
        CHKERR( EPSSetTolerances(self.eps, ctol, cmaxit) )

    def getDimensions(self):
        cdef PetscInt nev = 0
        cdef PetscInt ncv = 0
        cdef PetscInt mpd = 0
        CHKERR( EPSGetDimensions(self.eps, &nev, &ncv, &mpd) )
        return (nev, ncv, mpd)

    def setDimensions(self, nev=None, ncv=None, mpd=None):
        cdef PetscInt cnev = PETSC_IGNORE
        cdef PetscInt cncv = PETSC_IGNORE
        cdef PetscInt cmpd = PETSC_IGNORE
        if nev is not None: cnev = nev
        if ncv is not None: cncv = ncv
        if mpd is not None: cmpd = mpd
        CHKERR( EPSSetDimensions(self.eps, cnev, cncv, cmpd) )

    def getST(self):
        cdef ST st = ST()
        CHKERR( EPSGetST(self.eps, &st.st) )
        st.inc_ref(); return st

    def setST(self, ST st not None):
        CHKERR( EPSSetST(self.eps, st.st) )

    def getIP(self):
        cdef IP ip = IP()
        CHKERR( EPSGetIP(self.eps, &ip.ip) )
        ip.inc_ref(); return ip

    def setIP(self, IP ip not None):
        CHKERR( EPSSetIP(self.eps, ip.ip) )

    def getOperators(self):
        cdef Mat A = Mat()
        cdef Mat B = Mat()
        CHKERR( EPSGetOperators(self.eps, &A.mat, &B.mat) )
        A.inc_ref(); B.inc_ref(); return (A, B)

    def setOperators(self, Mat A not None, Mat B=None):
        cdef PetscMat Bmat = NULL
        if B is not None: Bmat = B.mat
        CHKERR( EPSSetOperators(self.eps, A.mat, Bmat) )

    #

    def getInitialVector(self):
        cdef Vec V = Vec()
        CHKERR( EPSGetInitialVector(self.eps, &V.vec) )
        V.inc_ref(); return V

    def setInitialVector(self, Vec V not None):
        CHKERR( EPSSetInitialVector(self.eps, V.vec) )

    def getInitialVectorLeft(self):
        cdef Vec V = Vec()
        CHKERR( EPSGetLeftInitialVector(self.eps, &V.vec) )
        V.inc_ref(); return V

    def setInitialVectorLeft(self, Vec W not None):
        CHKERR( EPSSetLeftInitialVector(self.eps, W.vec) )

    #

    def setUp(self):
        CHKERR( EPSSetUp(self.eps) )

    def solve(self):
        CHKERR( EPSSolve(self.eps) )

    def getIterationNumber(self):
        cdef PetscInt val = 0
        CHKERR( EPSGetIterationNumber(self.eps, &val) )
        return val

    def getConvergedReason(self):
        cdef SlepcEPSConvergedReason val = EPS_CONVERGED_ITERATING
        CHKERR( EPSGetConvergedReason(self.eps, &val) )
        return val

    def getConverged(self):
        cdef PetscInt val = 0
        CHKERR( EPSGetConverged(self.eps, &val) )
        return val

    def getValue(self, int i):
        cdef PetscScalar sval1 = 0
        cdef PetscScalar sval2 = 0
        CHKERR( EPSGetValue(self.eps, i, &sval1, &sval2) )
        return complex(toScalar(sval1), toScalar(sval2))

    def getVector(self, int i, Vec Vr not None, Vec Vi=None):
        cdef PetscVec vecr = NULL
        cdef PetscVec veci = NULL
        if Vr is not None: vecr = Vr.vec
        if Vi is not None: veci = Vi.vec
        CHKERR( EPSGetRightVector(self.eps, i, vecr, veci) )

    def getVectorLeft(self, int i, Vec Wr not None, Vec Wi=None):
        cdef PetscVec vecr = NULL
        cdef PetscVec veci = NULL
        if Wr is not None: vecr = Wr.vec
        if Wi is not None: veci = Wi.vec
        CHKERR( EPSGetLeftVector(self.eps, i, vecr, veci) )

    def getEigenpair(self, int i, Vec Vr=None, Vec Vi=None):
        cdef PetscScalar sval1 = 0
        cdef PetscScalar sval2 = 0
        cdef PetscVec vecr = NULL
        cdef PetscVec veci = NULL
        if Vr is not None: vecr = Vr.vec
        if Vi is not None: veci = Vi.vec
        CHKERR( EPSGetEigenpair(self.eps, i, &sval1, &sval2, vecr, veci) )
        return complex(toScalar(sval1), toScalar(sval2))

    #

    def getErrorEstimate(self, int i):
        cdef PetscReal val = 0
        CHKERR( EPSGetErrorEstimate(self.eps, i, &val) )
        return val

    def getErrorEstimateLeft(self, int i):
        cdef PetscReal val = 0
        CHKERR( EPSGetErrorEstimateLeft(self.eps, i, &val) )
        return val

    def computeRelativeError(self, int i):
        cdef PetscReal val = 0
        CHKERR( EPSComputeRelativeError(self.eps, i, &val) )
        return val

    def computeRelativeErrorLeft(self, int i):
        cdef PetscReal val = 0
        CHKERR( EPSComputeRelativeErrorLeft(self.eps, i, &val) )
        return val

    def computeResidualNorm(self, int i):
        cdef PetscReal val = 0
        CHKERR( EPSComputeResidualNorm(self.eps, i, &val) )
        return val

    def computeResidualNormLeft(self, int i):
        cdef PetscReal val = 0
        CHKERR( EPSComputeResidualNormLeft(self.eps, i, &val) )
        return val
    #

    property problem_type:
        def __get__(self):
            return self.getProblemType()
        def __set__(self, value):
            self.setProblemType(value)

    property extraction:
        def __get__(self):
            return self.getExtraction()
        def __set__(self, value):
            self.setExtraction(value)

    property which:
        def __get__(self):
            return self.getWhichEigenpairs()
        def __set__(self, value):
            self.setWhichEigenpairs(value)

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

    property st:
        def __get__(self):
            return self.getST()
        def __set__(self, value):
            self.setST(value)

    property ip:
        def __get__(self):
            return self.getIP()
        def __set__(self, value):
            self.setIP(value)

# --------------------------------------------------------------------
