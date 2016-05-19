# -----------------------------------------------------------------------------

class EPSType(object):
    """
    EPS type

    Native sparse eigensolvers.

    - `POWER`:        Power Iteration, Inverse Iteration, RQI.
    - `SUBSPACE`:     Subspace Iteration.
    - `ARNOLDI`:      Arnoldi.
    - `LANCZOS`:      Lanczos.
    - `KRYLOVSCHUR`:  Krylov-Schur (default).
    - `GD`:           Generalized Davidson.
    - `JD`:           Jacobi-Davidson.
    - `RQCG`:         Rayleigh Quotient Conjugate Gradient.
    - `LOBPCG`:       Locally Optimal Block Preconditioned Conjugate Gradient.
    - `CISS`:         Contour Integral Spectrum Slicing.
    - `LAPACK`:       Wrappers to dense eigensolvers in Lapack.

    Wrappers to sparse eigensolvers
    (should be enabled during installation of SLEPc)

    - `ARPACK`:
    - `BLZPACK`:
    - `TRLAN`:
    - `BLOPEX`:
    - `PRIMME`:
    - `FEAST`:
    """
    # provided implementations
    POWER        = S_(EPSPOWER)
    SUBSPACE     = S_(EPSSUBSPACE)
    ARNOLDI      = S_(EPSARNOLDI)
    LANCZOS      = S_(EPSLANCZOS)
    KRYLOVSCHUR  = S_(EPSKRYLOVSCHUR)
    GD           = S_(EPSGD)
    JD           = S_(EPSJD)
    RQCG         = S_(EPSRQCG)
    LOBPCG       = S_(EPSLOBPCG)
    CISS         = S_(EPSCISS)
    LAPACK       = S_(EPSLAPACK)
    # with external libraries
    ARPACK       = S_(EPSARPACK)
    BLZPACK      = S_(EPSBLZPACK)
    TRLAN        = S_(EPSTRLAN)
    BLOPEX       = S_(EPSBLOPEX)
    PRIMME       = S_(EPSPRIMME)
    FEAST        = S_(EPSFEAST)

class EPSProblemType(object):
    """
    EPS problem type

    - `HEP`:    Hermitian eigenproblem.
    - `NHEP`:   Non-Hermitian eigenproblem.
    - `GHEP`:   Generalized Hermitian eigenproblem.
    - `GNHEP`:  Generalized Non-Hermitian eigenproblem.
    - `PGNHEP`: Generalized Non-Hermitian eigenproblem
                with positive definite ``B``.
    - `GHIEP`:  Generalized Hermitian-indefinite eigenproblem.
    """
    HEP    = EPS_HEP
    NHEP   = EPS_NHEP
    GHEP   = EPS_GHEP
    GNHEP  = EPS_GNHEP
    PGNHEP = EPS_PGNHEP
    GHIEP  = EPS_GHIEP

class EPSExtraction(object):
    """
    EPS extraction technique

    - `RITZ`:              Standard Rayleigh-Ritz extraction.
    - `HARMONIC`:          Harmonic extraction.
    - `HARMONIC_RELATIVE`: Harmonic extraction relative to the eigenvalue.
    - `HARMONIC_RIGHT`:    Harmonic extraction for rightmost eigenvalues.
    - `HARMONIC_LARGEST`:  Harmonic extraction for largest magnitude (without target).
    - `REFINED`:           Refined extraction.
    - `REFINED_HARMONIC`:  Refined harmonic extraction.
    """
    RITZ              = EPS_RITZ
    HARMONIC          = EPS_HARMONIC
    HARMONIC_RELATIVE = EPS_HARMONIC_RELATIVE
    HARMONIC_RIGHT    = EPS_HARMONIC_RIGHT
    HARMONIC_LARGEST  = EPS_HARMONIC_LARGEST
    REFINED           = EPS_REFINED
    REFINED_HARMONIC  = EPS_REFINED_HARMONIC

class EPSBalance(object):
    """
    EPS type of balancing used for non-Hermitian problems

    - `NONE`:     None.
    - `ONESIDE`:  One-sided eigensolver, only right eigenvectors.
    - `TWOSIDE`:  Two-sided eigensolver, right and left eigenvectors.
    - `USER`:     User-defined.
    """
    NONE    = EPS_BALANCE_NONE
    ONESIDE = EPS_BALANCE_ONESIDE
    TWOSIDE = EPS_BALANCE_TWOSIDE
    USER    = EPS_BALANCE_USER

class EPSErrorType(object):
    """
    EPS error type to assess accuracy of computed solutions

    - `ABSOLUTE`:  Absolute error.
    - `RELATIVE`:  Relative error.
    - `BACKWARD`:  Backward error.
    """
    ABSOLUTE = EPS_ERROR_ABSOLUTE
    RELATIVE = EPS_ERROR_RELATIVE
    BACKWARD = EPS_ERROR_BACKWARD

class EPSWhich(object):
    """
    EPS desired piece of spectrum

    - `LARGEST_MAGNITUDE`:  Largest magnitude (default).
    - `LARGEST_REAL`:       Largest real parts.
    - `LARGEST_IMAGINARY`:  Largest imaginary parts in magnitude.
    - `SMALLEST_MAGNITUDE`: Smallest magnitude.
    - `SMALLEST_REAL`:      Smallest real parts.
    - `SMALLEST_IMAGINARY`: Smallest imaginary parts in magnitude.
    - `TARGET_MAGNITUDE`:   Closest to target (in magnitude).
    - `TARGET_REAL`:        Real part closest to target.
    - `TARGET_IMAGINARY`:   Imaginary part closest to target.
    - `ALL`:                All eigenvalues in an interval.
    - `USER`:               User defined ordering.
    """
    LARGEST_MAGNITUDE  = EPS_LARGEST_MAGNITUDE
    LARGEST_REAL       = EPS_LARGEST_REAL
    LARGEST_IMAGINARY  = EPS_LARGEST_IMAGINARY
    SMALLEST_MAGNITUDE = EPS_SMALLEST_MAGNITUDE
    SMALLEST_REAL      = EPS_SMALLEST_REAL
    SMALLEST_IMAGINARY = EPS_SMALLEST_IMAGINARY
    TARGET_MAGNITUDE   = EPS_TARGET_MAGNITUDE
    TARGET_REAL        = EPS_TARGET_REAL
    TARGET_IMAGINARY   = EPS_TARGET_IMAGINARY
    ALL                = EPS_ALL
    USER               = EPS_WHICH_USER

class EPSConv(object):
    """
    EPS convergence test

    - `ABS`:
    - `REL`:
    - `NORM`:
    - `USER`:
    """
    ABS  = EPS_CONV_ABS
    REL  = EPS_CONV_REL
    NORM = EPS_CONV_NORM
    USER = EPS_CONV_USER

class EPSConvergedReason(object):
    """
    EPS convergence reasons

    - `CONVERGED_TOL`:
    - `CONVERGED_USER`:
    - `DIVERGED_ITS`:
    - `DIVERGED_BREAKDOWN`:
    - `DIVERGED_SYMMETRY_LOST`:
    - `CONVERGED_ITERATING`:
    """
    CONVERGED_TOL          = EPS_CONVERGED_TOL
    CONVERGED_USER         = EPS_CONVERGED_USER
    DIVERGED_ITS           = EPS_DIVERGED_ITS
    DIVERGED_BREAKDOWN     = EPS_DIVERGED_BREAKDOWN
    DIVERGED_SYMMETRY_LOST = EPS_DIVERGED_SYMMETRY_LOST
    CONVERGED_ITERATING    = EPS_CONVERGED_ITERATING
    ITERATING              = EPS_CONVERGED_ITERATING

class EPSPowerShiftType(object):
    """
    EPS Power shift type.

    - `CONSTANT`:
    - `RAYLEIGH`:
    - `WILKINSON`:
    """
    CONSTANT  = EPS_POWER_SHIFT_CONSTANT
    RAYLEIGH  = EPS_POWER_SHIFT_RAYLEIGH
    WILKINSON = EPS_POWER_SHIFT_WILKINSON

class EPSLanczosReorthogType(object):
    """
    EPS Lanczos reorthogonalization type

    - `LOCAL`:
    - `FULL`:
    - `SELECTIVE`:
    - `PERIODIC`:
    - `PARTIAL`:
    - `DELAYED`:
    """
    LOCAL     =  EPS_LANCZOS_REORTHOG_LOCAL
    FULL      =  EPS_LANCZOS_REORTHOG_FULL
    SELECTIVE =  EPS_LANCZOS_REORTHOG_SELECTIVE
    PERIODIC  =  EPS_LANCZOS_REORTHOG_PERIODIC
    PARTIAL   =  EPS_LANCZOS_REORTHOG_PARTIAL
    DELAYED   =  EPS_LANCZOS_REORTHOG_DELAYED

# -----------------------------------------------------------------------------

cdef class EPS(Object):

    """
    EPS
    """

    Type            = EPSType
    ProblemType     = EPSProblemType
    Extraction      = EPSExtraction
    Balance         = EPSBalance
    ErrorType       = EPSErrorType
    Which           = EPSWhich
    Conv            = EPSConv
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
        viewer: Viewer, optional.
                Visualization context; if not provided, the standard
                output is used.
        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( EPSView(self.eps, vwr) )

    def destroy(self):
        """
        Destroys the EPS object.
        """
        CHKERR( EPSDestroy(&self.eps) )
        self.eps = NULL
        return self

    def reset(self):
        """
        Resets the EPS object.
        """
        CHKERR( EPSReset(self.eps) )

    def create(self, comm=None):
        """
        Creates the EPS object.

        Parameters
        ----------
        comm: MPI_Comm, optional
              MPI communicator; if not provided, it defaults to all
              processes.
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcEPS neweps = NULL
        CHKERR( EPSCreate(ccomm, &neweps) )
        SlepcCLEAR(self.obj); self.eps = neweps
        return self

    def setType(self, eps_type):
        """
        Selects the particular solver to be used in the EPS object.

        Parameters
        ----------
        eps_type: `EPS.Type` enumerate
                  The solver to be used.

        Notes
        -----
        See `EPS.Type` for available methods. The default is
        `EPS.Type.KRYLOVSCHUR`.  Normally, it is best to use
        `setFromOptions()` and then set the EPS type from the options
        database rather than by using this routine.  Using the options
        database provides the user with maximum flexibility in
        evaluating the different available methods.
        """
        cdef SlepcEPSType cval = NULL
        eps_type = str2bytes(eps_type, &cval)
        CHKERR( EPSSetType(self.eps, cval) )

    def getType(self):
        """
        Gets the EPS type of this object.

        Returns
        -------
        type: `EPS.Type` enumerate
              The solver currently being used.
        """
        cdef SlepcEPSType eps_type = NULL
        CHKERR( EPSGetType(self.eps, &eps_type) )
        return bytes2str(eps_type)

    def getOptionsPrefix(self):
        """
        Gets the prefix used for searching for all EPS options in the
        database.

        Returns
        -------
        prefix: string
                The prefix string set for this EPS object.
        """
        cdef const_char *prefix = NULL
        CHKERR( EPSGetOptionsPrefix(self.eps, &prefix) )
        return bytes2str(prefix)

    def setOptionsPrefix(self, prefix):
        """
        Sets the prefix used for searching for all EPS options in the
        database.

        Parameters
        ----------
        prefix: string
                The prefix string to prepend to all EPS option
                requests.

        Notes
        -----
        A hyphen (-) must NOT be given at the beginning of the prefix
        name.  The first character of all runtime options is
        AUTOMATICALLY the hyphen.

        For example, to distinguish between the runtime options for
        two different EPS contexts, one could call::

            E1.setOptionsPrefix("eig1_")
            E2.setOptionsPrefix("eig2_")
        """
        cdef const_char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( EPSSetOptionsPrefix(self.eps, cval) )

    def appendOptionsPrefix(self, prefix):
        """
        Appends to the prefix used for searching for all EPS options
        in the database.

        Parameters
        ----------
        prefix: string
                The prefix string to prepend to all EPS option requests.
        """
        cdef const_char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( EPSAppendOptionsPrefix(self.eps, cval) )

    def setFromOptions(self):
        """
        Sets EPS options from the options database. This routine must
        be called before `setUp()` if the user is to be allowed to set
        the solver type.

        Notes
        -----
        To see all options, run your program with the ``-help``
        option.
        """
        CHKERR( EPSSetFromOptions(self.eps) )

    #

    def getProblemType(self):
        """
        Gets the problem type from the EPS object.

        Returns
        -------
        problem_type: `EPS.ProblemType` enumerate
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
        problem_type: `EPS.ProblemType` enumerate
               The problem type to be set.

        Notes
        -----
        Allowed values are: Hermitian (HEP), non-Hermitian (NHEP),
        generalized Hermitian (GHEP), generalized non-Hermitian
        (GNHEP), and generalized non-Hermitian with positive
        semi-definite B (PGNHEP).

        This function must be used to instruct SLEPc to exploit
        symmetry. If no problem type is specified, by default a
        non-Hermitian problem is assumed (either standard or
        generalized). If the user knows that the problem is Hermitian
        (i.e. ``A=A^H``) or generalized Hermitian (i.e. ``A=A^H``,
        ``B=B^H``, and ``B`` positive definite) then it is recommended
        to set the problem type so that eigensolver can exploit these
        properties.
        """
        cdef SlepcEPSProblemType val = problem_type
        CHKERR( EPSSetProblemType(self.eps, val) )

    def isGeneralized(self):
        """
        Tells whether the EPS object corresponds to a generalized
        eigenvalue problem.

        Returns
        -------
        flag: boolean
              True if two matrices were set with `setOperators()`.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSIsGeneralized(self.eps, &tval) )
        return <bint> tval

    def isHermitian(self):
        """
        Tells whether the EPS object corresponds to a Hermitian
        eigenvalue problem.

        Returns
        -------
        flag: boolean
              True if the problem type set with `setProblemType()` was
              Hermitian.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSIsHermitian(self.eps, &tval) )
        return <bint> tval

    def isPositive(self):
        """
        Tells whether the EPS object corresponds to an eigenvalue problem
        type that requires a positive (semi-) definite matrix B.

        Returns
        -------
        flag: boolean
              True if the problem type set with `setProblemType()` was
              positive.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSIsPositive(self.eps, &tval) )
        return <bint> tval

    def getBalance(self):
        """
        Gets the balancing type used by the EPS object,
        and the associated parameters.

        Returns
        -------
        balance: `EPS.Balance` enumerate
                 The balancing method
        iterations: integer
                    Number of iterations of the balancing algorithm
        cutoff: real
                Cutoff value
        """
        cdef SlepcEPSBalance val = EPS_BALANCE_ONESIDE
        cdef PetscInt ival = 0
        cdef PetscReal rval = 0
        CHKERR( EPSGetBalance(self.eps, &val, &ival, &rval) )
        return (val, toInt(ival), toReal(rval))

    def setBalance(self, balance=None, iterations=None, cutoff=None):
        """
        Specifies the balancing technique to be employed by the
        eigensolver, and some parameters associated to it.

        Parameters
        ----------
        balance: `EPS.Balance` enumerate
                 The balancing method
        iterations: integer
                    Number of iterations of the balancing algorithm
        cutoff: real
                Cutoff value
        """
        cdef SlepcEPSBalance val = <SlepcEPSBalance>PETSC_DEFAULT
        cdef PetscInt  ival = PETSC_DEFAULT
        cdef PetscReal rval = PETSC_DEFAULT
        if balance    is not None: val  = balance
        if iterations is not None: ival = asInt(iterations)
        if cutoff     is not None: rval = asReal(cutoff)
        CHKERR( EPSSetBalance(self.eps, val, ival, rval) )

    def getExtraction(self):
        """
        Gets the extraction type used by the EPS object.

        Returns
        -------
        extraction: `EPS.Extraction` enumerate
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
        extraction: `EPS.Extraction` enumerate
                    The extraction method to be used by the solver.

        Notes
        -----
        Not all eigensolvers support all types of extraction. See the
        SLEPc documentation for details.

        By default, a standard Rayleigh-Ritz extraction is used. Other
        extractions may be useful when computing interior eigenvalues.

        Harmonic-type extractions are used in combination with a
        *target*. See `setTarget()`.
        """
        cdef SlepcEPSExtraction val = extraction
        CHKERR( EPSSetExtraction(self.eps, val) )

    def getWhichEigenpairs(self):
        """
        Returns which portion of the spectrum is to be sought.

        Returns
        -------
        which: `EPS.Which` enumerate
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
        which: `EPS.Which` enumerate
               The portion of the spectrum to be sought by the solver.

        Notes
        -----
        Not all eigensolvers implemented in EPS account for all the
        possible values. Also, some values make sense only for certain
        types of problems. If SLEPc is compiled for real numbers
        `EPS.Which.LARGEST_IMAGINARY` and
        `EPS.Which.SMALLEST_IMAGINARY` use the absolute value of the
        imaginary part for eigenvalue selection.
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
        The target is a scalar value used to determine the portion of
        the spectrum of interest. It is used in combination with
        `setWhichEigenpairs()`.
        """
        cdef PetscScalar sval = asScalar(target)
        CHKERR( EPSSetTarget(self.eps, sval) )

    def getInterval(self):
        """
        Gets the computational interval for spectrum slicing.

        Returns
        -------
        inta: float
                The left end of the interval.
        intb: float
                The right end of the interval.

        Notes
        -----
        If the interval was not set by the user, then zeros are returned.
        """
        cdef PetscReal inta = 0
        cdef PetscReal intb = 0
        CHKERR( EPSGetInterval(self.eps, &inta, &intb) )
        return (toReal(inta), toReal(intb))

    def setInterval(self, inta, intb):
        """
        Defines the computational interval for spectrum slicing.

        Parameters
        ----------
        inta: float
                The left end of the interval.
        intb: float
                The right end of the interval.

        Notes
        -----
        Spectrum slicing is a technique employed for computing all
        eigenvalues of symmetric eigenproblems in a given interval.
        This function provides the interval to be considered. It must
        be used in combination with `EPS.Which.ALL`, see
        `setWhichEigenpairs()`.
        """
        cdef PetscReal rval1 = asReal(inta)
        cdef PetscReal rval2 = asReal(intb)
        CHKERR( EPSSetInterval(self.eps, rval1, rval2) )

    #

    def getTolerances(self):
        """
        Gets the tolerance and maximum iteration count used by the
        default EPS convergence tests.

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
        return (toReal(rval), toInt(ival))

    def setTolerances(self, tol=None, max_it=None):
        """
        Sets the tolerance and maximum iteration count used by the
        default EPS convergence tests.

        Parameters
        ----------
        tol: float, optional
             The convergence tolerance.
        max_it: int, optional
             The maximum number of iterations

        Notes
        -----
        Use `DECIDE` for maxits to assign a reasonably good value,
        which is dependent on the solution method.
        """
        cdef PetscReal rval = PETSC_DEFAULT
        cdef PetscInt  ival = PETSC_DEFAULT
        if tol    is not None: rval = asReal(tol)
        if max_it is not None: ival = asInt(max_it)
        CHKERR( EPSSetTolerances(self.eps, rval, ival) )

    def getConvergenceTest(self):
        """
        Return the method used to compute the error estimate 
        used in the convergence test. 

        Returns
        -------
        conv: EPS.Conv
            The method used to compute the error estimate 
            used in the convergence test. 
        """
        cdef SlepcEPSConv conv = EPS_CONV_REL
        CHKERR( EPSGetConvergenceTest(self.eps, &conv) )
        return conv

    def setConvergenceTest(self, conv):
        """
        Specifies how to compute the error estimate 
        used in the convergence test. 

        Parameters
        ----------
        conv: EPS.Conv
            The method used to compute the error estimate 
            used in the convergence test.
        """
        cdef SlepcEPSConv tconv = conv
        CHKERR( EPSSetConvergenceTest(self.eps, tconv) )

    def getTrueResidual(self):
        """
        Returns the flag indicating whether true residual must be
        computed explicitly or not.

        Returns
        -------
        trueres: bool
            Whether the solver compute all residuals or not.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSGetTrueResidual(self.eps, &tval) )
        return <bint> tval

    def setTrueResidual(self, trueres):
        """
        Specifies if the solver must compute the true residual 
        explicitly or not.

        Parameters
        ----------
        trueres: bool
            Whether compute the true residual or not.
        """
        cdef PetscBool tval = trueres
        CHKERR( EPSSetTrueResidual(self.eps, tval) )

    def getTrackAll(self):
        """
        Returns the flag indicating whether all residual norms must be
        computed or not.

        Returns
        -------
        trackall: bool
            Whether the solver compute all residuals or not.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSGetTrackAll(self.eps, &tval) )
        return <bint> tval

    def setTrackAll(self, trackall):
        """
        Specifies if the solver must compute the residual of all
        approximate eigenpairs or not.

        Parameters
        ----------
        trackall: bool
            Whether compute all residuals or not.
        """
        cdef PetscBool tval = trackall
        CHKERR( EPSSetTrackAll(self.eps, tval) )

    def getDimensions(self):
        """
        Gets the number of eigenvalues to compute and the dimension of
        the subspace.

        Returns
        -------
        nev: int
             Number of eigenvalues to compute.
        ncv: int
             Maximum dimension of the subspace to be used by the
             solver.
        mpd: int
             Maximum dimension allowed for the projected problem.
        """
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        cdef PetscInt ival3 = 0
        CHKERR( EPSGetDimensions(self.eps, &ival1, &ival2, &ival3) )
        return (toInt(ival1), toInt(ival2), toInt(ival3))

    def setDimensions(self, nev=None, ncv=None, mpd=None):
        """
        Sets the number of eigenvalues to compute and the dimension of
        the subspace.

        Parameters
        ----------
        nev: int, optional
             Number of eigenvalues to compute.
        ncv: int, optional
             Maximum dimension of the subspace to be used by the
             solver.
        mpd: int, optional
             Maximum dimension allowed for the projected problem.

        Notes
        -----
        Use `DECIDE` for `ncv` and `mpd` to assign a reasonably good
        value, which is dependent on the solution method.

        The parameters `ncv` and `mpd` are intimately related, so that
        the user is advised to set one of them at most. Normal usage
        is the following:

        + In cases where `nev` is small, the user sets `ncv`
          (a reasonable default is 2 * `nev`).

        + In cases where `nev` is large, the user sets `mpd`.

        The value of `ncv` should always be between `nev` and (`nev` +
        `mpd`), typically `ncv` = `nev` + `mpd`. If `nev` is not too
        large, `mpd` = `nev` is a reasonable choice, otherwise a
        smaller value should be used.
        """
        cdef PetscInt ival1 = PETSC_DEFAULT
        cdef PetscInt ival2 = PETSC_DEFAULT
        cdef PetscInt ival3 = PETSC_DEFAULT
        if nev is not None: ival1 = asInt(nev)
        if ncv is not None: ival2 = asInt(ncv)
        if mpd is not None: ival3 = asInt(mpd)
        CHKERR( EPSSetDimensions(self.eps, ival1, ival2, ival3) )

    def getST(self):
        """
        Obtain the spectral transformation (`ST`) object associated to
        the eigensolver object.

        Returns
        -------
        st: ST
            The spectral transformation.
        """
        cdef ST st = ST()
        CHKERR( EPSGetST(self.eps, &st.st) )
        PetscINCREF(st.obj)
        return st

    def setST(self, ST st not None):
        """
        Associates a spectral transformation object to the
        eigensolver.

        Parameters
        ----------
        st: ST
            The spectral transformation.
        """
        CHKERR( EPSSetST(self.eps, st.st) )

    def getBV(self):
        """
        Obtain the basis vector objects associated to the eigensolver.

        Returns
        -------
        bv: BV
            The basis vectors context.
        """
        cdef BV bv = BV()
        CHKERR( EPSGetBV(self.eps, &bv.bv) )
        PetscINCREF(bv.obj)
        return bv

    def setBV(self, BV bv not None):
        """
        Associates a basis vectors object to the eigensolver.

        Parameters
        ----------
        bv: BV
            The basis vectors context.
        """
        CHKERR( EPSSetBV(self.eps, bv.bv) )

    def getDS(self):
        """
        Obtain the direct solver associated to the eigensolver.

        Returns
        -------
        ds: DS
            The direct solver context.
        """
        cdef DS ds = DS()
        CHKERR( EPSGetDS(self.eps, &ds.ds) )
        PetscINCREF(ds.obj)
        return ds

    def setDS(self, DS ds not None):
        """
        Associates a direct solver object to the eigensolver.

        Parameters
        ----------
        ds: DS
            The direct solver context.
        """
        CHKERR( EPSSetDS(self.eps, ds.ds) )

    def getRG(self):
        """
        Obtain the region object associated to the eigensolver.

        Returns
        -------
        rg: RG
            The region context.
        """
        cdef RG rg = RG()
        CHKERR( EPSGetRG(self.eps, &rg.rg) )
        PetscINCREF(rg.obj)
        return rg

    def setRG(self, RG rg not None):
        """
        Associates a region object to the eigensolver.

        Parameters
        ----------
        rg: RG
            The region context.
        """
        CHKERR( EPSSetRG(self.eps, rg.rg) )

    def getOperators(self):
        """
        Gets the matrices associated with the eigenvalue problem.

        Returns
        -------
        A: Mat
           The matrix associated with the eigensystem.
        B: Mat
           The second matrix in the case of generalized eigenproblems.
        """
        cdef Mat A = Mat()
        cdef Mat B = Mat()
        CHKERR( EPSGetOperators(self.eps, &A.mat, &B.mat) )
        PetscINCREF(A.obj)
        PetscINCREF(B.obj)
        return (A, B)

    def setOperators(self, Mat A not None, Mat B=None):
        """
        Sets the matrices associated with the eigenvalue problem.

        Parameters
        ----------
        A: Mat
           The matrix associated with the eigensystem.
        B: Mat, optional
           The second matrix in the case of generalized eigenproblems;
           if not provided, a standard eigenproblem is assumed.
        """
        cdef PetscMat Bmat = NULL
        if B is not None: Bmat = B.mat
        CHKERR( EPSSetOperators(self.eps, A.mat, Bmat) )

    def setDeflationSpace(self, space):
        """
        Add vectors to the basis of the deflation space.

        Parameters
        ----------
        space: a Vec or an array of Vec
               Set of basis vectors to be added to the deflation
               space.

        Notes
        -----
        When a deflation space is given, the eigensolver seeks the
        eigensolution in the restriction of the problem to the
        orthogonal complement of this space. This can be used for
        instance in the case that an invariant subspace is known
        beforehand (such as the nullspace of the matrix).

        The vectors do not need to be mutually orthonormal, since they
        are explicitly orthonormalized internally.

        These vectors do not persist from one `solve()` call to the other,
        so the deflation space should be set every time.
        """
        if isinstance(space, Vec): space = [space]
        cdef PetscVec* vs = NULL
        cdef Py_ssize_t i = 0, ns = len(space)
        cdef tmp = allocate(<size_t>ns*sizeof(Vec),<void**>&vs)
        for i in range(ns): vs[i] = (<Vec?>space[i]).vec
        CHKERR( EPSSetDeflationSpace(self.eps, <PetscInt>ns, vs) )

    #

    def setInitialSpace(self, space):
        """
        Sets the initial space from which the eigensolver starts to
        iterate.

        Parameters
        ----------
        space: Vec or sequence of Vec
           The initial space

        Notes
        -----
        Some solvers start to iterate on a single vector (initial vector).
        In that case, the other vectors are ignored.

        In contrast to `setDeflationSpace()`, these vectors do not persist
        from one `solve()` call to the other, so the initial space should be
        set every time.

        The vectors do not need to be mutually orthonormal, since they are
        explicitly orthonormalized internally.

        Common usage of this function is when the user can provide a rough
        approximation of the wanted eigenspace. Then, convergence may be faster.
        """
        if isinstance(space, Vec): space = [space]
        cdef PetscVec *vs = NULL
        cdef Py_ssize_t i = 0, ns = len(space)
        cdef tmp = allocate(<size_t>ns*sizeof(Vec),<void**>&vs)
        for i in range(ns): vs[i] = (<Vec?>space[i]).vec
        CHKERR( EPSSetInitialSpace(self.eps, <PetscInt>ns, vs) )

    #

    def cancelMonitor(self):
        """
        Clears all monitors for an EPS object.
        """
        CHKERR( EPSMonitorCancel(self.eps) )

    #

    def setUp(self):
        """
        Sets up all the internal data structures necessary for the
        execution of the eigensolver.

        Notes
        -----
        This function need not be called explicitly in most cases,
        since `solve()` calls it. It can be useful when one wants to
        measure the set-up time separately from the solve time.
        """
        CHKERR( EPSSetUp(self.eps) )

    def solve(self):
        """
        Solves the eigensystem.
        """
        CHKERR( EPSSolve(self.eps) )

    def getIterationNumber(self):
        """
        Gets the current iteration number. If the call to `solve()` is
        complete, then it returns the number of iterations carried out
        by the solution method.

        Returns
        -------
        its: int
             Iteration number.
        """
        cdef PetscInt ival = 0
        CHKERR( EPSGetIterationNumber(self.eps, &ival) )
        return toInt(ival)

    def getConvergedReason(self):
        """
        Gets the reason why the `solve()` iteration was stopped.

        Returns
        -------
        reason: `EPS.ConvergedReason` enumerate
                Negative value indicates diverged, positive value
                converged.
        """
        cdef SlepcEPSConvergedReason val = EPS_CONVERGED_ITERATING
        CHKERR( EPSGetConvergedReason(self.eps, &val) )
        return val

    def getConverged(self):
        """
        Gets the number of converged eigenpairs.

        Returns
        -------
        nconv: int
               Number of converged eigenpairs.

        Notes
        -----
        This function should be called after `solve()` has finished.
        """
        cdef PetscInt ival = 0
        CHKERR( EPSGetConverged(self.eps, &ival) )
        return toInt(ival)

    def getEigenvalue(self, int i):
        """
        Gets the i-th eigenvalue as computed by `solve()`.

        Parameters
        ----------
        i: int
           Index of the solution to be obtained.

        Returns
        -------
        e: scalar (possibly complex)
           The computed eigenvalue.

        Notes
        -----
        The index ``i`` should be a value between ``0`` and
        ``nconv-1`` (see `getConverged()`). Eigenpairs are indexed
        according to the ordering criterion established with
        `setWhichEigenpairs()`.
        """
        cdef PetscScalar sval1 = 0
        cdef PetscScalar sval2 = 0
        CHKERR( EPSGetEigenvalue(self.eps, i, &sval1, &sval2) )
        return complex(toScalar(sval1), toScalar(sval2))

    def getEigenvector(self, int i, Vec Vr not None, Vec Vi=None):
        """
        Gets the i-th eigenvector as computed by `solve()`.

        Parameters
        ----------
        i: int
           Index of the solution to be obtained.
        Vr: Vec
            Placeholder for the returned eigenvector (real part).
        Vi: Vec, optional
            Placeholder for the returned eigenvector (imaginary part).

        Notes
        -----
        The index ``i`` should be a value between ``0`` and
        ``nconv-1`` (see `getConverged()`). Eigenpairs are indexed
        according to the ordering criterion established with
        `setWhichEigenpairs()`.
        """
        cdef PetscVec vecr = NULL
        cdef PetscVec veci = NULL
        if Vr is not None: vecr = Vr.vec
        if Vi is not None: veci = Vi.vec
        CHKERR( EPSGetEigenvector(self.eps, i, vecr, veci) )

    def getEigenpair(self, int i, Vec Vr=None, Vec Vi=None):
        """
        Gets the i-th solution of the eigenproblem as computed by
        `solve()`.  The solution consists of both the eigenvalue and
        the eigenvector.

        Parameters
        ----------
        i: int
           Index of the solution to be obtained.
        Vr: Vec
            Placeholder for the returned eigenvector (real part).
        Vi: Vec
            Placeholder for the returned eigenvector (imaginary part).

        Returns
        -------
        e: scalar (possibly complex)
           The computed eigenvalue.

        Notes
        -----
        The index ``i`` should be a value between ``0`` and
        ``nconv-1`` (see `getConverged()`). Eigenpairs are indexed
        according to the ordering criterion established with
        `setWhichEigenpairs()`.
        """
        cdef PetscScalar sval1 = 0
        cdef PetscScalar sval2 = 0
        cdef PetscVec vecr = NULL
        cdef PetscVec veci = NULL
        if Vr is not None: vecr = Vr.vec
        if Vi is not None: veci = Vi.vec
        CHKERR( EPSGetEigenpair(self.eps, i, &sval1, &sval2, vecr, veci) )
        return complex(toScalar(sval1), toScalar(sval2))

    def getInvariantSubspace(self):
        """
        Gets an orthonormal basis of the computed invariant subspace.

        Returns
        -------
        subspace: list of Vec
           Basis of the invariant subspace.

        Notes
        -----
        This function should be called after `solve()` has finished.

        The returned vectors span an invariant subspace associated
        with the computed eigenvalues. An invariant subspace ``X`` of
        ``A` satisfies ``A x`` in ``X`` for all ``x`` in ``X`` (a
        similar definition applies for generalized eigenproblems).
        """
        cdef PetscInt i = 0, ncv = 0
        cdef PetscVec v = NULL, *isp = NULL
        cdef list subspace = []
        CHKERR( EPSGetConverged(self.eps, &ncv) )
        if ncv == 0: return subspace
        cdef PetscMat A = NULL
        CHKERR( EPSGetOperators(self.eps, &A, NULL) )
        CHKERR( MatCreateVecs(A, &v, NULL) )
        cdef Vec V = None
        cdef object tmp = allocate(ncv*sizeof(Vec),<void**>&isp)
        for i in range(ncv):
            if i == 0: isp[0] = v
            if i >= 1: CHKERR( VecDuplicate(v, &isp[i]) )
            V = Vec(); V.vec = isp[i]; subspace.append(V)
        CHKERR( EPSGetInvariantSubspace(self.eps, isp) )
        return subspace

    #

    def getErrorEstimate(self, int i):
        """
        Returns the error estimate associated to the i-th computed
        eigenpair.

        Parameters
        ----------
        i: int
           Index of the solution to be considered.

        Returns
        -------
        e: real
           Error estimate.

        Notes
        -----
        This is the error estimate used internally by the
        eigensolver. The actual error bound can be computed with
        `computeError()`.
        """
        cdef PetscReal rval = 0
        CHKERR( EPSGetErrorEstimate(self.eps, i, &rval) )
        return toReal(rval)

    def computeError(self, int i, etype=None):
        """
        Computes the error (based on the residual norm) associated with the i-th
        computed eigenpair.

        Parameters
        ----------
        i: int
           Index of the solution to be considered.
        etype: `EPS.ErrorType` enumerate
           The error type to compute.

        Returns
        -------
        e: real
           The error bound, computed in various ways from the residual norm
           ``||Ax-kBx||_2`` where ``k`` is the eigenvalue and
           ``x`` is the eigenvector.

        Notes
        -----
        The index ``i`` should be a value between ``0`` and
        ``nconv-1`` (see `getConverged()`).
        """
        cdef SlepcEPSErrorType et = EPS_ERROR_RELATIVE
        cdef PetscReal rval = 0
        if etype is not None: et = etype
        CHKERR( EPSComputeError(self.eps, i, et, &rval) )
        return toReal(rval)

    def errorView(self, etype=None, Viewer viewer=None):
        """
        Displays the errors associated with the computed solution
        (as well as the eigenvalues).

        Parameters
        ----------
        etype: `EPS.ErrorType` enumerate, optional
           The error type to compute.
        viewer: Viewer, optional.
                Visualization context; if not provided, the standard
                output is used.

        Notes
        -----
        By default, this function checks the error of all eigenpairs and prints
        the eigenvalues if all of them are below the requested tolerance.
        If the viewer has format ``ASCII_INFO_DETAIL`` then a table with
        eigenvalues and corresponding errors is printed.

        """
        cdef SlepcEPSErrorType et = EPS_ERROR_RELATIVE
        if etype is not None: et = etype
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( EPSErrorView(self.eps, et, vwr) )

    #

    def setPowerShiftType(self, shift):
        """
        Sets the type of shifts used during the power iteration. This
        can be used to emulate the Rayleigh Quotient Iteration (RQI)
        method.

        Parameters
        ----------
        shift: `EPS.PowerShiftType` enumerate
               The type of shift.

        Notes
        -----
        This call is only relevant if the type was set to
        `EPS.Type.POWER` with `setType()`.

        By default, shifts are constant
        (`EPS.PowerShiftType.CONSTANT`) and the iteration is the
        simple power method (or inverse iteration if a
        shift-and-invert transformation is being used).

        A variable shift can be specified
        (`EPS.PowerShiftType.RAYLEIGH` or
        `EPS.PowerShiftType.WILKINSON`). In this case, the iteration
        behaves rather like a cubic converging method as RQI.
        """
        cdef SlepcEPSPowerShiftType val = shift
        CHKERR( EPSPowerSetShiftType(self.eps, val) )

    def getPowerShiftType(self):
        """
        Gets the type of shifts used during the power iteration.

        Returns
        -------
        shift: `EPS.PowerShiftType` enumerate
               The type of shift.
        """
        cdef SlepcEPSPowerShiftType val = EPS_POWER_SHIFT_CONSTANT
        CHKERR( EPSPowerGetShiftType(self.eps, &val) )
        return val

    def setArnoldiDelayed(self, delayed):
        """
        Activates or deactivates delayed reorthogonalization in the
        Arnoldi iteration.

        Parameters
        ----------
        delayed: boolean
                 True if delayed reorthogonalization is to be used.

        Notes
        -----
        This call is only relevant if the type was set to
        `EPS.Type.ARNOLDI` with `setType()`.

        Delayed reorthogonalization is an aggressive optimization for
        the Arnoldi eigensolver than may provide better scalability,
        but sometimes makes the solver converge less than the default
        algorithm.
        """
        cdef PetscBool val = PETSC_FALSE
        if delayed: val = PETSC_TRUE
        CHKERR( EPSArnoldiSetDelayed(self.eps, val) )

    def getArnoldiDelayed(self):
        """
        Gets the type of reorthogonalization used during the Arnoldi
        iteration.

        Returns
        -------
        delayed: boolean
                 True if delayed reorthogonalization is to be used.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSArnoldiGetDelayed(self.eps, &tval) )
        return <bint> tval

    def setLanczosReorthogType(self, reorthog):
        """
        Sets the type of reorthogonalization used during the Lanczos
        iteration.

        Parameters
        ----------
        reorthog: `EPS.LanczosReorthogType` enumerate
                  The type of reorthogonalization.

        Notes
        -----
        This call is only relevant if the type was set to
        `EPS.Type.LANCZOS` with `setType()`.
        """
        cdef SlepcEPSLanczosReorthogType val = reorthog
        CHKERR( EPSLanczosSetReorthog(self.eps, val) )

    def getLanczosReorthogType(self):
        """
        Gets the type of reorthogonalization used during the Lanczos
        iteration.

        Returns
        -------
        reorthog: `EPS.LanczosReorthogType` enumerate
                  The type of reorthogonalization.
        """
        cdef SlepcEPSLanczosReorthogType val = \
            EPS_LANCZOS_REORTHOG_LOCAL
        CHKERR( EPSLanczosGetReorthog(self.eps, &val) )
        return val

    #

    def setKrylovSchurRestart(self, keep):
        """
        Sets the restart parameter for the Krylov-Schur method, in
        particular the proportion of basis vectors that must be kept
        after restart.

        Parameters
        ----------
        keep: float
              The number of vectors to be kept at restart.

        Notes
        -----
        Allowed values are in the range [0.1,0.9]. The default is 0.5.
        """
        cdef PetscReal val = keep
        CHKERR( EPSKrylovSchurSetRestart(self.eps, val) )

    def getKrylovSchurRestart(self):
        """
        Gets the restart parameter used in the Krylov-Schur method.

        Returns
        -------
        keep: float
              The number of vectors to be kept at restart.
        """
        cdef PetscReal val = 0
        CHKERR( EPSKrylovSchurGetRestart(self.eps, &val) )
        return val

    def setKrylovSchurLocking(self, lock):
        """
        Choose between locking and non-locking variants of the
        Krylov-Schur method.

        Parameters
        ----------
        lock: bool
              True if the locking variant must be selected.

        Notes
        -----
        The default is to lock converged eigenpairs when the method restarts.
        This behaviour can be changed so that all directions are kept in the
        working subspace even if already converged to working accuracy (the
        non-locking variant).
        """
        cdef PetscBool val = lock
        CHKERR( EPSKrylovSchurSetLocking(self.eps, val) )

    def getKrylovSchurLocking(self):
        """
        Gets the locking flag used in the Krylov-Schur method.

        Returns
        -------
        lock: bool
              The locking flag.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSKrylovSchurGetLocking(self.eps, &tval) )
        return <bint> tval

    def setKrylovSchurPartitions(self, npart):
        """
        Sets the number of partitions for the case of doing spectrum
        slicing for a computational interval with the communicator split
        in several sub-communicators.

        Parameters
        ----------
        npart: int
              The number of partitions.

        Notes
        -----
        By default, npart=1 so all processes in the communicator participate in
        the processing of the whole interval. If npart>1 then the interval is
        divided into npart subintervals, each of them being processed by a
        subset of processes.
        """
        cdef PetscInt val = npart
        CHKERR( EPSKrylovSchurSetPartitions(self.eps, val) )

    def getKrylovSchurPartitions(self):
        """
        Gets the number of partitions of the communicator in case of
        spectrum slicing.

        Returns
        -------
        npart: int
              The number of partitions.
        """
        cdef PetscInt val = 0
        CHKERR( EPSKrylovSchurGetPartitions(self.eps, &val) )
        return val

    def setKrylovSchurDetectZeros(self, detect):
        """
        Sets a flag to enforce detection of zeros during the factorizations
        throughout the spectrum slicing computation.

        Parameters
        ----------
        detect: bool
              True if zeros must checked for.

        Notes
        -----
        A zero in the factorization indicates that a shift coincides with
        an eigenvalue.

        This flag is turned off by default, and may be necessary in some cases,
        especially when several partitions are being used. This feature currently
        requires an external package for factorizations with support for zero
        detection, e.g. MUMPS.
        """
        cdef PetscBool val = detect
        CHKERR( EPSKrylovSchurSetDetectZeros(self.eps, val) )

    def getKrylovSchurDetectZeros(self):
        """
        Gets the flag that enforces zero detection in spectrum slicing.

        Returns
        -------
        detect: bool
              The zero detection flag.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( EPSKrylovSchurGetDetectZeros(self.eps, &tval) )
        return <bint> tval

    def setKrylovSchurDimensions(self, nev=None, ncv=None, mpd=None):
        """
        Sets the dimensions used for each subsolve step in case of doing
        spectrum slicing for a computational interval. The meaning of the
        parameters is the same as in `setDimensions()`.

        Parameters
        ----------
        nev: int, optional
             Number of eigenvalues to compute.
        ncv: int, optional
             Maximum dimension of the subspace to be used by the solver.
        mpd: int, optional
             Maximum dimension allowed for the projected problem.
        """
        cdef PetscInt ival1 = PETSC_DEFAULT
        cdef PetscInt ival2 = PETSC_DEFAULT
        cdef PetscInt ival3 = PETSC_DEFAULT
        if nev is not None: ival1 = asInt(nev)
        if ncv is not None: ival2 = asInt(ncv)
        if mpd is not None: ival3 = asInt(mpd)
        CHKERR( EPSKrylovSchurSetDimensions(self.eps, ival1, ival2, ival3) )

    def getKrylovSchurDimensions(self):
        """
        Gets the dimensions used for each subsolve step in case of doing
        spectrum slicing for a computational interval.

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
        CHKERR( EPSKrylovSchurGetDimensions(self.eps, &ival1, &ival2, &ival3) )
        return (toInt(ival1), toInt(ival2), toInt(ival3))

    def getKrylovSchurSubcommInfo(self):
        """
        Gets information related to the case of doing spectrum slicing
        for a computational interval with multiple communicators.

        Returns
        -------
        k: int
             Number of the subinterval for the calling process.
        n: int
             Number of eigenvalues found in the k-th subinterval.
        v: Vec
             A vector owned by processes in the subcommunicator with dimensions
             compatible for locally computed eigenvectors.

        Notes
        -----
        This function is only available for spectrum slicing runs.

        The returned Vec should be destroyed by the user.
        """
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        cdef Vec vec
        vec = Vec()
        CHKERR( EPSKrylovSchurGetSubcommInfo(self.eps, &ival1, &ival2, &vec.vec) )
        return (toInt(ival1), toInt(ival2), vec)

    def getKrylovSchurSubcommPairs(self, int i, Vec V):
        """
        Gets the i-th eigenpair stored internally in the multi-communicator
        to which the calling process belongs.

        Parameters
        ----------
        i: int
           Index of the solution to be obtained.
        V: Vec
           Placeholder for the returned eigenvector.

        Returns
        -------
        e: scalar
           The computed eigenvalue.

        Notes
        -----
        The index ``i`` should be a value between ``0`` and ``n-1``,
        where ``n`` is the number of vectors in the local subinterval,
        see `getKrylovSchurSubcommInfo()`.
        """
        cdef PetscScalar sval = 0
        cdef PetscVec vec = NULL
        if V is not None: vec = V.vec
        CHKERR( EPSKrylovSchurGetSubcommPairs(self.eps, i, &sval, vec) )
        return toScalar(sval)

    def getKrylovSchurSubcommMats(self):
        """
        Gets the eigenproblem matrices stored internally in the subcommunicator
        to which the calling process belongs.

        Returns
        -------
        A: Mat
           The matrix associated with the eigensystem.
        B: Mat
           The second matrix in the case of generalized eigenproblems.

        Notes
        -----
        This is the analog of `getOperators()`, but returns the matrices distributed
        differently (in the subcommunicator rather than in the parent communicator).

        These matrices should not be modified by the user.
        """
        cdef Mat A = Mat()
        cdef Mat B = Mat()
        CHKERR( EPSKrylovSchurGetSubcommMats(self.eps, &A.mat, &B.mat) )
        PetscINCREF(A.obj)
        PetscINCREF(B.obj)
        return (A, B)

    def updateKrylovSchurSubcommMats(self, s=1.0, a=1.0, Mat Au=None, t=1.0, b=1.0, Mat Bu=None, structure=None, globalup=False):
        """
        Update the eigenproblem matrices stored internally in the subcommunicator
        to which the calling process belongs.

        Parameters
        ----------
        s: float (real or complex)
           Scalar that multiplies the existing A matrix.
        a: float (real or complex)
           Scalar used in the axpy operation on A.
        Au: Mat, optional
           The matrix used in the axpy operation on A.
        t: float (real or complex)
           Scalar that multiplies the existing B matrix.
        b: float (real or complex)
           Scalar used in the axpy operation on B.
        Bu: Mat, optional
           The matrix used in the axpy operation on B.
        structure: `PETSc.Mat.Structure` enumerate
           Either same, different, or a subset of the non-zero sparsity pattern.
        globalup: bool
           Whether global matrices must be updated or not.

        Notes
        -----
        This function modifies the eigenproblem matrices at subcommunicator level,
        and optionally updates the global matrices in the parent communicator.
        The updates are expressed as ``A <-- s*A + a*Au``,  ``B <-- t*B + b*Bu``.

        It is possible to update one of the matrices, or both.

        The matrices `Au` and `Bu` must be equal in all subcommunicators.

        The `str` flag is passed to the `Mat.axpy()` operations to perform the updates.

        If `globalup` is True, communication is carried out to reconstruct the updated
        matrices in the parent communicator.
        """
        cdef PetscMat Amat = NULL
        if Au is not None: Amat = Au.mat
        cdef PetscMat Bmat = NULL
        if Bu is not None: Bmat = Bu.mat
        cdef PetscMatStructure vstr = matstructure(structure)
        cdef PetscBool tval = globalup
        CHKERR( EPSKrylovSchurUpdateSubcommMats(self.eps, s, a, Amat, t, b, Bmat, vstr, tval) )

    def setKrylovSchurSubintervals(self, subint):
        """
        Sets the subinterval boundaries for spectrum slicing with a computational interval.
        
        Parameters
        ----------
        subint: list of real values specifying subintervals

        Notes
        -----
        Logically Collective on EPS
        This function must be called after setKrylovSchurPartitions(). 
        For npart partitions, the argument subint must contain npart+1 
        real values sorted in ascending order: 
        subint_0, subint_1, ..., subint_npart, 
        where the first and last values must coincide with the interval 
        endpoints set with EPSSetInterval().
        The subintervals are then defined by two consecutive points: 
        [subint_0,subint_1], [subint_1,subint_2], and so on.
        """
        cdef PetscBool match = PETSC_FALSE
        CHKERR( PetscObjectTypeCompare(<PetscObject>self.eps, EPSKRYLOVSCHUR, &match) )
        if match == PETSC_FALSE: return
        cdef PetscReal *subintarray = NULL
        cdef Py_ssize_t i = 0, n = len(subint)
        cdef PetscInt nparts = 0
        CHKERR( EPSKrylovSchurGetPartitions(self.eps, &nparts) )
        assert n >= nparts
        cdef tmp = allocate(n*sizeof(PetscReal),<void**>&subintarray)
        for i in range(n): subintarray[i] = asReal(subint[i])
        CHKERR(EPSKrylovSchurSetSubintervals(self.eps, subintarray))
        
    def setRQCGReset(self, nrest):
        """
        Sets the reset parameter of the RQCG iteration. Every nrest iterations,
        the solver performs a Rayleigh-Ritz projection step.

        Parameters
        ----------
        nrest: integer
               The number of iterations between resets.
        """
        cdef PetscInt val = nrest
        CHKERR( EPSRQCGSetReset(self.eps, val) )

    def getRQCGReset(self):
        """
        Gets the reset parameter used in the RQCG method.

        Returns
        -------
        nrest: integer
               The number of iterations between resets.
        """
        cdef PetscInt val = 0
        CHKERR( EPSRQCGGetReset(self.eps, &val) )
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

    property bv:
        def __get__(self):
            return self.getBV()
        def __set__(self, value):
            self.setBV(value)

# -----------------------------------------------------------------------------

del EPSType
del EPSProblemType
del EPSExtraction
del EPSBalance
del EPSErrorType
del EPSWhich
del EPSConv
del EPSConvergedReason
del EPSPowerShiftType
del EPSLanczosReorthogType

# -----------------------------------------------------------------------------
