# --------------------------------------------------------------------

class EPSType(object):
    """
    EPS types
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

class EPSProblemType(object):
    """
    EPS problem type
    """
    HEP    = EPS_HEP
    NHEP   = EPS_NHEP
    GHEP   = EPS_GHEP
    GNHEP  = EPS_GNHEP
    PGNHEP = EPS_PGNHEP

class EPSExtraction(object):
    """
    EPS extraction technique
    """
    RITZ             = EPS_RITZ
    HARMONIC         = EPS_HARMONIC
    REFINED          = EPS_REFINED
    REFINED_HARMONIC = EPS_REFINED_HARMONIC

class EPSClass(object):
    """
    EPS class of method
    """
    ONE_SIDE  = EPS_ONE_SIDE
    TWO_SIDE  = EPS_TWO_SIDE

class EPSWhich(object):
    """
    EPS desired piece of spectrum
    """
    LARGEST_MAGNITUDE  = EPS_LARGEST_MAGNITUDE
    SMALLEST_MAGNITUDE = EPS_SMALLEST_MAGNITUDE
    LARGEST_REAL       = EPS_LARGEST_REAL
    SMALLEST_REAL      = EPS_SMALLEST_REAL
    LARGEST_IMAGINARY  = EPS_LARGEST_IMAGINARY
    SMALLEST_IMAGINARY = EPS_SMALLEST_IMAGINARY

class EPSConvergedReason(object):
    """
    EPS convergence reasons
    """
    CONVERGED_TOL         = EPS_CONVERGED_TOL
    DIVERGED_ITS          = EPS_DIVERGED_ITS
    DIVERGED_BREAKDOWN    = EPS_DIVERGED_BREAKDOWN
    DIVERGED_NONSYMMETRIC = EPS_DIVERGED_NONSYMMETRIC
    CONVERGED_ITERATING   = EPS_CONVERGED_ITERATING

class EPSPowerShiftType(object):
    """
    EPS Power shift type
    """
    CONSTANT  = EPSPOWER_SHIFT_CONSTANT
    RAYLEIGH  = EPSPOWER_SHIFT_RAYLEIGH
    WILKINSON = EPSPOWER_SHIFT_WILKINSON

class EPSLanczosReorthogType(object):
    """
    EPS Lanczos reorthogonalization type
    """
    LOCAL     =  EPSLANCZOS_REORTHOG_LOCAL
    FULL      =  EPSLANCZOS_REORTHOG_FULL
    SELECTIVE =  EPSLANCZOS_REORTHOG_SELECTIVE
    PERIODIC  =  EPSLANCZOS_REORTHOG_PERIODIC
    PARTIAL   =  EPSLANCZOS_REORTHOG_PARTIAL
    DELAYED   =  EPSLANCZOS_REORTHOG_DELAYED

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

    PowerShiftType      = EPSPowerShiftType
    LanczosReorthogType = EPSLanczosReorthogType

    def __cinit__(self):
        self.obj = <PetscObject*> &self.eps
        self.eps = NULL

    def view(self, Viewer viewer=None):
        """
        Prints the EPS data structure.

        Parameters
        ----------
        viewer: Viewer, optional
                Visualization context; if not provided, the standard output is used.
        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( EPSView(self.eps, vwr) )

    def destroy(self):
        """
        Destroys the EPS object.
        """
        CHKERR( EPSDestroy(self.eps) )
        self.eps = NULL
        return self

    def create(self, comm=None):
        """
        Creates the EPS object.

        Parameters
        ----------
        comm: MPI_Comm, optional
              MPI communicator; if not provided, it defaults to all processes.
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcEPS neweps = NULL
        CHKERR( EPSCreate(ccomm, &neweps) )
        self.dec_ref(); self.eps = neweps
        return self

    def setType(self, eps_type):
        """
        Selects the particular solver to be used in the EPS object.

        Parameters
        ----------
        eps_type: EPS.Type enumerate
                  The solver to be used.

        Options Database Keys
        ---------------------
        -eps_type <method>: Sets the method; use -help for a list 
                            of available methods.
    
        Notes
        -----
        See `EPSType` for available methods. The default is KRYLOVSCHUR.
        Normally, it is best to use `setFromOptions()` and then set the EPS
        type from the options database rather than by using this routine.
        Using the options database provides the user with maximum flexibility
        in evaluating the different available methods.
        """
        CHKERR( EPSSetType(self.eps, str2cp(eps_type)) )

    def getType(self):
        """
        Gets the EPS type of this object.

        Returns
        -------
        type: EPS.Type enumerate
              The solver currently being used.
        """
        cdef SlepcEPSType eps_type = NULL
        CHKERR( EPSGetType(self.eps, &eps_type) )
        return cp2str(eps_type)

    def setOptionsPrefix(self, prefix):
        """
        Sets the prefix used for searching for all EPS options in the database.

        Parameters
        ----------
        prefix: string
                The prefix string to prepend to all EPS option requests.

        Notes
        -----
        A hyphen (-) must NOT be given at the beginning of the prefix name.
        The first character of all runtime options is AUTOMATICALLY the
        hyphen.

        For example, to distinguish between the runtime options for two
        different EPS contexts, one could call::

            E1.setOptionsPrefix("eig1_")
            E2.setOptionsPrefix("eig2_")
        """
        CHKERR( EPSSetOptionsPrefix(self.eps, str2cp(prefix)) )

    def getOptionsPrefix(self):
        """
        Gets the prefix used for searching for all EPS options in the database.

        Returns
        -------
        prefix: string
                The prefix string set for this EPS object.
        """
        cdef const_char_p prefix = NULL
        CHKERR( EPSGetOptionsPrefix(self.eps, &prefix) )
        return cp2str(prefix)

    def setFromOptions(self):
        """
        Sets EPS options from the options database. This routine must be 
        called before `setUp()` if the user is to be allowed to set the 
        solver type.

        Notes
        -----
        To see all options, run your program with the -help option.
        """
        CHKERR( EPSSetFromOptions(self.eps) )

    #

    def getProblemType(self):
        """
        Gets the problem type from the EPS object.

        Returns
        -------
        problem_type: EPS.ProblemType enumerate
                      The problem type that was previously set.
        """
        cdef SlepcEPSProblemType val = EPS_NHEP
        CHKERR( EPSGetProblemType(self.eps, &val) )
        return val

    def setProblemType(self, problem_type):
        """
        Specifies the type of the eigenvalue problem.

        Parameters
        ----------
        problem_type: EPS.ProblemType enumerate
               The problem type to be set.

        Notes
        -----
        Allowed values are: Hermitian (HEP), non-Hermitian (NHEP), generalized
        Hermitian (GHEP), generalized non-Hermitian (GNHEP), and generalized
        non-Hermitian with positive semi-definite B (PGNHEP).

        This function must be used to instruct SLEPc to exploit symmetry. If no
        problem type is specified, by default a non-Hermitian problem is assumed
        (either standard or generalized). If the user knows that the problem is
        Hermitian (i.e. A=A^H) or generalized Hermitian (i.e. A=A^H, B=B^H, and 
        B positive definite) then it is recommended to set the problem type so
        that eigensolver can exploit these properties. 
        """
        cdef SlepcEPSProblemType val = problem_type
        CHKERR( EPSSetProblemType(self.eps, val) )

    def isGeneralized(self):
        """
        Tells whether the EPS object corresponds to a generalized eigenvalue problem.

        Returns
        -------
        flag: boolean
              True if two matrices were set with `setOperators()`.
        """
        cdef PetscTruth tval = PETSC_FALSE
        CHKERR( EPSIsGeneralized(self.eps, &tval) )
        return <bint> tval

    def isHermitian(self):
        """
        Tells whether the EPS object corresponds to a Hermitian eigenvalue problem.

        Returns
        -------
        flag: boolean
              True if the problem type set with `setProblemType()` was Hermitian.
        """
        cdef PetscTruth tval = PETSC_FALSE
        CHKERR( EPSIsHermitian(self.eps, &tval) )
        return <bint> tval

    def getClass(self):
        """
        Gets the eigensolver class from the EPS object.

        Returns
        -------
        klass: EPS.Class enumerate
               The class of eigensolver, either ONE_SIDE or TWO_SIDE.
        """
        cdef SlepcEPSClass val = EPS_ONE_SIDE
        CHKERR( EPSGetClass(self.eps, &val) )
        return val

    def setClass(self, klass):
        """
        Specifies the eigensolver class: either one-sided or two-sided.

        Parameters
        ----------
        klass: EPS.Class enumerate
               The class of eigensolver, either ONE_SIDE or TWO_SIDE.

        Notes
        -----
        Allowed solver classes are: one-sided (ONE_SIDE) and two-sided (TWO_SIDE).
        One-sided eigensolvers are the standard ones, which allow the computation of
        eigenvalues and (right) eigenvectors, whereas two-sided eigensolvers compute
        left eigenvectors as well.
        """
        cdef SlepcEPSClass val = klass
        CHKERR( EPSSetClass(self.eps, val) )

    def getExtraction(self):
        """
        Gets the extraction type used by the EPS object.

        Returns
        -------
        extraction: EPS.Extraction enumerate
                    The method of extraction.
        """
        cdef SlepcEPSExtraction val = EPS_RITZ
        CHKERR( EPSGetExtraction(self.eps, &val) )
        return val

    def setExtraction(self, extraction):
        """
        Sets the extraction type used by the EPS object.

        Parameters
        ----------
        extraction: EPS.Extraction enumerate
                    The extraction method to be used by the solver.

        Notes
        -----
        Not all eigensolvers support all types of extraction. See the SLEPc
        documentation for details.

        By default, a standard Rayleigh-Ritz extraction is used. Other extractions
        may be useful when computing interior eigenvalues.

        Harmonic-type extractions are used in combination with a 'target'. See
        `setTarget()`.
        """
        cdef SlepcEPSExtraction val = extraction
        CHKERR( EPSSetExtraction(self.eps, val) )

    def getWhichEigenpairs(self):
        """
        Returns which portion of the spectrum is to be sought.

        Returns
        -------
        which: EPS.Which enumerate
               The portion of the spectrum to be sought by the solver.
        """
        cdef SlepcEPSWhich val = EPS_LARGEST_MAGNITUDE
        CHKERR( EPSGetWhichEigenpairs(self.eps, &val) )
        return val

    def setWhichEigenpairs(self, which):
        """
        Specifies which portion of the spectrum is to be sought.

        Parameters
        ----------
        which: EPS.Which enumerate
               The portion of the spectrum to be sought by the solver.

        Notes
        -----
        Not all eigensolvers implemented in EPS account for all the possible
        values. Also, some values make sense only for certain types of 
        problems. If SLEPc is compiled for real numbers EPS_LARGEST_IMAGINARY
        and SMALLEST_IMAGINARY use the absolute value of the imaginary part 
        for eigenvalue selection.
        """
        cdef SlepcEPSWhich val = which
        CHKERR( EPSSetWhichEigenpairs(self.eps, val) )

    def getTarget(self):
        """
        Gets the value of the target.

        Returns
        -------
        target: float (real or complex)
                The value of the target.

        Notes
        -----
        If the target was not set by the user, then zero is returned.
        """
        cdef PetscScalar sval = 0
        CHKERR( EPSGetTarget(self.eps, &sval) )
        return toScalar(sval)

    def setTarget(self, target):
        """
        Sets the value of the target.

        Parameters
        ----------
        target: float (real or complex)
                The value of the target.

        Notes
        -----
        The target is a scalar value used to determine the portion of the
        spectrum of interest.

        If the target is not specified, then eigenvalues are computed according 
        to the which parameter, see `setWhichEigenpairs()`.
   
        If the target is specified, then the sought-after eigenvalues are those
        closest to the target.
        """
        cdef PetscScalar sval = asScalar(target)
        CHKERR( EPSSetTarget(self.eps, sval) )

    #

    def getTolerances(self):
        """
        Gets the tolerance and maximum iteration count used by the default 
        EPS convergence tests.

        Returns
        -------
        tol: float
             The convergence tolerance.
        max_it: int
             The maximum number of iterations
        """
        cdef PetscReal rval = 0
        cdef PetscInt  ival = 0
        CHKERR( EPSGetTolerances(self.eps, &rval, &ival) )
        return (rval, ival)

    def setTolerances(self, tol=None, max_it=None):
        """
        Sets the tolerance and maximum iteration count used by the default 
        EPS convergence tests.

        Parameters
        ----------
        tol: float, optional
             The convergence tolerance.
        max_it: int, optional
             The maximum number of iterations

        Notes
        -----
        Use PETSC_DECIDE for maxits to assign a reasonably good value, which
        is dependent on the solution method.
        """
        cdef PetscReal rval = PETSC_IGNORE
        cdef PetscInt  ival = PETSC_IGNORE
        if tol    is not None: rval = tol
        if max_it is not None: ival = max_it
        CHKERR( EPSSetTolerances(self.eps, rval, ival) )

    def getDimensions(self):
        """
        Gets the number of eigenvalues to compute and the dimension of the subspace.

        Returns
        -------
        nev: int
             Number of eigenvalues to compute.
        ncv: int
             Maximum dimension of the subspace to be used by the solver.
        mpd: int
             Maximum dimension allowed for the projected problem.
        """
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        cdef PetscInt ival3 = 0
        CHKERR( EPSGetDimensions(self.eps, &ival1, &ival2, &ival3) )
        return (ival1, ival2, ival3)

    def setDimensions(self, nev=None, ncv=None, mpd=None):
        """
        Sets the number of eigenvalues to compute and the dimension of the subspace.

        Parameters
        ----------
        nev: int, optional
             Number of eigenvalues to compute.
        ncv: int, optional
             Maximum dimension of the subspace to be used by the solver.
        mpd: int, optional
             Maximum dimension allowed for the projected problem.

        Notes
        -----
        Use PETSC_DECIDE for `ncv` and `mpd` to assign a reasonably good value, 
        which is dependent on the solution method.

        The parameters `ncv` and `mpd` are intimately related, so that the user 
        is advised to set one of them at most. Normal usage is the following:

         - In cases where `nev` is small, the user sets `ncv` (a reasonable default
           is `2*nev`).
         - In cases where `nev` is large, the user sets `mpd`.

        The value of `ncv` should always be between `nev` and `(nev+mpd)`, typically
        `ncv=nev+mpd`. If `nev` is not too large, `mpd=nev` is a reasonable choice, 
        otherwise a smaller value should be used.
        """
        cdef PetscInt ival1 = PETSC_IGNORE
        cdef PetscInt ival2 = PETSC_IGNORE
        cdef PetscInt ival3 = PETSC_IGNORE
        if nev is not None: ival1 = nev
        if ncv is not None: ival2 = ncv
        if mpd is not None: ival3 = mpd
        CHKERR( EPSSetDimensions(self.eps, ival1, ival2, ival3) )

    def getST(self):
        """
        Obtain the spectral transformation (`ST`) object associated to the 
        eigensolver object.

        Returns
        -------
        st: ST
            The spectral transformation.
        """
        cdef ST st = ST()
        CHKERR( EPSGetST(self.eps, &st.st) )
        st.inc_ref(); return st

    def setST(self, ST st not None):
        """
        Associates a spectral transformation object to the eigensolver.

        Parameters
        ----------
        st: ST
            The spectral transformation.
        """
        CHKERR( EPSSetST(self.eps, st.st) )

    def getIP(self):
        """
        Obtain the inner product object associated to the eigensolver object.

        Returns
        -------
        ip: IP
            The inner product context.
        """
        cdef IP ip = IP()
        CHKERR( EPSGetIP(self.eps, &ip.ip) )
        ip.inc_ref(); return ip

    def setIP(self, IP ip not None):
        """
        Associates an inner product object to the eigensolver.

        Parameters
        ----------
        ip: IP
            The inner product context.
        """
        CHKERR( EPSSetIP(self.eps, ip.ip) )

    def getOperators(self):
        """
        Gets the matrices associated with the eigenvalue problem.

        Returns
        -------
        A: PETSc.Mat
           The matrix associated with the eigensystem.
        B: PETSc.Mat
           The second matrix in the case of generalized eigenproblems.
        """
        cdef Mat A = Mat()
        cdef Mat B = Mat()
        CHKERR( EPSGetOperators(self.eps, &A.mat, &B.mat) )
        A.inc_ref(); B.inc_ref(); return (A, B)

    def setOperators(self, Mat A not None, Mat B=None):
        """
        Sets the matrices associated with the eigenvalue problem.

        Parameters
        ----------
        A: PETSc.Mat
           The matrix associated with the eigensystem.
        B: PETSc.Mat, optional
           The second matrix in the case of generalized eigenproblems;
           if not provided, a standard eigenproblem is assumed.
        """
        cdef PetscMat Bmat = NULL
        if B is not None: Bmat = B.mat
        CHKERR( EPSSetOperators(self.eps, A.mat, Bmat) )

    def attachDeflationSpace(self, space, ortho=False):
        """
        Add vectors to the basis of the deflation space.

        Parameters
        ----------
        space: a PETSc.Vec or an array of PETSc.Vec
               Set of basis vectors to be added to the deflation space.
        ortho: boolean, optional
               True if basis vectors can be assumed to be mutullay orthonormal.

        Notes
        -----
        When a deflation space is given, the eigensolver seeks the eigensolution
        in the restriction of the problem to the orthogonal complement of this
        space. This can be used for instance in the case that an invariant 
        subspace is known beforehand (such as the nullspace of the matrix).

        The basis vectors can be provided all at once or incrementally with
        several calls to `attachDeflationSpace()`.

        Set `ortho` to True if all the vectors passed in are known to be 
        mutually orthonormal.
        """
        cdef PetscInt i = 0, nds = 0
        cdef PetscVec* vds = NULL
        cdef PetscTruth tval = PETSC_FALSE
        cdef object tmp = None
        if isinstance(space, Vec): space = [space]
        nds = len(space)
        tmp = allocate(nds*sizeof(Vec),<void**>&vds)
        if ortho: tval = PETSC_TRUE
        for i in range(nds): vds[i] = (<Vec?>space[i]).vec
        CHKERR( EPSAttachDeflationSpace(self.eps, nds, vds, tval) )

    def removeDeflationSpace(self):
        """
        Removes the deflation space previously set with `attachDeflationSpace()`.
        """
        CHKERR( EPSRemoveDeflationSpace(self.eps) )

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
        cdef PetscInt ival = 0
        CHKERR( EPSGetIterationNumber(self.eps, &ival) )
        return ival

    def getConvergedReason(self):
        cdef SlepcEPSConvergedReason val = EPS_CONVERGED_ITERATING
        CHKERR( EPSGetConvergedReason(self.eps, &val) )
        return val

    def getConverged(self):
        cdef PetscInt ival = 0
        CHKERR( EPSGetConverged(self.eps, &ival) )
        return ival

    def getInvariantSubspace(self):
        cdef PetscInt i = 0, ncv = 0
        cdef PetscVec v = NULL, *isp = NULL
        CHKERR( EPSGetConverged(self.eps, &ncv) )
        CHKERR( EPSGetInitialVector(self.eps, &v) )
        cdef Vec V = None
        cdef list subspace = []
        cdef object tmp = allocate(ncv*sizeof(Vec),<void**>&isp)
        for i in range(ncv):
            V = Vec(); subspace.append(V)
            CHKERR( VecDuplicate(v, &isp[i]) )
            V.vec = isp[i]
        CHKERR( EPSGetInvariantSubspace(self.eps, isp) )
        return subspace

    def getInvariantSubspaceLeft(self):
        cdef PetscInt i = 0, ncv = 0
        cdef PetscVec w = NULL, *isp = NULL
        CHKERR( EPSGetConverged(self.eps, &ncv) )
        CHKERR( EPSGetLeftInitialVector(self.eps, &w) )
        cdef Vec W = None
        cdef list subspace = []
        cdef object tmp = allocate(ncv*sizeof(Vec),<void**>&isp)
        for i in range(ncv):
            W = Vec(); subspace.append(W)
            CHKERR( VecDuplicate(w, &isp[i]) )
            W.vec = isp[i]
        CHKERR( EPSGetLeftInvariantSubspace(self.eps, isp) )
        return subspace

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
        cdef PetscReal rval = 0
        CHKERR( EPSGetErrorEstimate(self.eps, i, &rval) )
        return rval

    def getErrorEstimateLeft(self, int i):
        cdef PetscReal rval = 0
        CHKERR( EPSGetErrorEstimateLeft(self.eps, i, &rval) )
        return rval

    def computeRelativeError(self, int i):
        cdef PetscReal rval = 0
        CHKERR( EPSComputeRelativeError(self.eps, i, &rval) )
        return rval

    def computeRelativeErrorLeft(self, int i):
        cdef PetscReal rval = 0
        CHKERR( EPSComputeRelativeErrorLeft(self.eps, i, &rval) )
        return rval

    def computeResidualNorm(self, int i):
        cdef PetscReal rval = 0
        CHKERR( EPSComputeResidualNorm(self.eps, i, &rval) )
        return rval

    def computeResidualNormLeft(self, int i):
        cdef PetscReal rval = 0
        CHKERR( EPSComputeResidualNormLeft(self.eps, i, &rval) )
        return rval

    def getOperationCounters(self):
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        cdef PetscInt ival3 = 0
        CHKERR( EPSGetOperationCounters(self.eps, &ival1, &ival2, &ival3) )
        return (ival1, ival2, ival3)

    #

    def setPowerShiftType(self, shift):
        cdef SlepcEPSPowerShiftType val = shift
        CHKERR( EPSPowerSetShiftType(self.eps, val) )

    def getPowerShiftType(self):
        cdef SlepcEPSPowerShiftType val = EPSPOWER_SHIFT_CONSTANT
        CHKERR( EPSPowerGetShiftType(self.eps, &val) )
        return val

    def setArnoldiDelayed(self, delayed):
        cdef PetscTruth val = PETSC_FALSE
        if delayed: val = PETSC_TRUE
        CHKERR( EPSArnoldiSetDelayed(self.eps, val) )

    def getArnoldiDelayed(self):
        cdef PetscTruth val = PETSC_FALSE
        CHKERR( EPSArnoldiGetDelayed(self.eps, &val) )
        return val

    def setLanczosReorthogType(self, reorthog):
        cdef SlepcEPSLanczosReorthogType val = reorthog
        CHKERR( EPSLanczosSetReorthog(self.eps, val) )

    def getLanczosReorthogType(self):
        cdef SlepcEPSLanczosReorthogType val = EPSLANCZOS_REORTHOG_LOCAL
        CHKERR( EPSLanczosGetReorthog(self.eps, &val) )
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

    property target:
        def __get__(self):
            return self.getTarget()
        def __set__(self, value):
            self.setTarget(value)

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

del EPSType
del EPSProblemType
del EPSExtraction
del EPSClass
del EPSWhich
del EPSConvergedReason
del EPSPowerShiftType
del EPSLanczosReorthogType

# --------------------------------------------------------------------
