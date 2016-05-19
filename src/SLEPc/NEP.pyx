# -----------------------------------------------------------------------------

class NEPType(object):
    """
    NEP type

    Nonlinear eigensolvers.

    - `RII`:      Residual inverse iteration.
    - `SLP`:      Successive linear problems.
    - `NARNOLDI`: Nonlinear Arnoldi.
    - `CISS`:     Contour integral spectrum slice.
    - `INTERPOL`: Polynomial interpolation.
    - `NLEIGS`:   Fully rational Krylov method for nonlinear eigenproblems.
    """
    RII      = S_(NEPRII)
    SLP      = S_(NEPSLP)
    NARNOLDI = S_(NEPNARNOLDI)
    CISS     = S_(NEPCISS)
    INTERPOL = S_(NEPINTERPOL)
    NLEIGS   = S_(NEPNLEIGS)

class NEPErrorType(object):
    """
    NEP error type to assess accuracy of computed solutions

    - `ABSOLUTE`:  Absolute error.
    - `RELATIVE`:  Relative error.
    - `BACKWARD`:  Backward error.
    """
    ABSOLUTE = NEP_ERROR_ABSOLUTE
    RELATIVE = NEP_ERROR_RELATIVE
    BACKWARD = NEP_ERROR_BACKWARD

class NEPWhich(object):
    LARGEST_MAGNITUDE  = NEP_LARGEST_MAGNITUDE
    SMALLEST_MAGNITUDE = NEP_SMALLEST_MAGNITUDE
    LARGEST_REAL       = NEP_LARGEST_REAL
    SMALLEST_REAL      = NEP_SMALLEST_REAL
    LARGEST_IMAGINARY  = NEP_LARGEST_IMAGINARY
    SMALLEST_IMAGINARY = NEP_SMALLEST_IMAGINARY
    TARGET_MAGNITUDE   = NEP_TARGET_MAGNITUDE
    TARGET_REAL        = NEP_TARGET_REAL
    TARGET_IMAGINARY   = NEP_TARGET_IMAGINARY
    ALL                = NEP_ALL
    USER               = NEP_WHICH_USER

class NEPConvergedReason(object):
    CONVERGED_TOL          = NEP_CONVERGED_TOL
    CONVERGED_USER         = NEP_CONVERGED_USER
    DIVERGED_ITS           = NEP_DIVERGED_ITS
    DIVERGED_BREAKDOWN     = NEP_DIVERGED_BREAKDOWN
    DIVERGED_LINEAR_SOLVE  = NEP_DIVERGED_LINEAR_SOLVE
    CONVERGED_ITERATING    = NEP_CONVERGED_ITERATING
    ITERATING              = NEP_CONVERGED_ITERATING

class NEPRefine(object):
    """
    NEP refinement strategy

    - `NONE`:     No refinement.
    - `SIMPLE`:   Refine eigenpairs one by one.
    - `MULTIPLE`: Refine all eigenpairs simultaneously (invariant pair).
    """
    NONE     = NEP_REFINE_NONE
    SIMPLE   = NEP_REFINE_SIMPLE
    MULTIPLE = NEP_REFINE_MULTIPLE

class NEPRefineScheme(object):
    """
    Scheme for solving linear systems during iterative refinement

    - `SCHUR`:    Schur complement.
    - `MBE`:      Mixed block elimination.
    - `EXPLICIT`: Build the explicit matrix.
    """
    SCHUR    = NEP_REFINE_SCHEME_SCHUR
    MBE      = NEP_REFINE_SCHEME_MBE
    EXPLICIT = NEP_REFINE_SCHEME_EXPLICIT

# -----------------------------------------------------------------------------

cdef class NEP(Object):

    """
    NEP
    """

    Type            = NEPType
    ErrorType       = NEPErrorType
    Which           = NEPWhich
    ConvergedReason = NEPConvergedReason
    Refine          = NEPRefine
    RefineScheme    = NEPRefineScheme

    def __cinit__(self):
        self.obj = <PetscObject*> &self.nep
        self.nep = NULL

    def view(self, Viewer viewer=None):
        """
        Prints the NEP data structure.

        Parameters
        ----------
        viewer: Viewer, optional.
            Visualization context; if not provided, the standard
            output is used.
        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( NEPView(self.nep, vwr) )

    def destroy(self):
        """
        Destroys the NEP object.
        """
        CHKERR( NEPDestroy(&self.nep) )
        self.nep = NULL
        return self

    def reset(self):
        """
        Resets the NEP object.
        """
        CHKERR( NEPReset(self.nep) )

    def create(self, comm=None):
        """
        Creates the NEP object.

        Parameters
        ----------
        comm: Comm, optional.
            MPI communicator. If not provided, it defaults to all
            processes.
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcNEP newnep = NULL
        CHKERR( NEPCreate(ccomm, &newnep) )
        SlepcCLEAR(self.obj); self.nep = newnep
        return self

    def setType(self, nep_type):
        """
        Selects the particular solver to be used in the NEP object.

        Parameters
        ----------
        nep_type: `NEP.Type` enumerate
            The solver to be used.
        """
        cdef SlepcNEPType cval = NULL
        nep_type = str2bytes(nep_type, &cval)
        CHKERR( NEPSetType(self.nep, cval) )

    def getType(self):
        """
        Gets the NEP type of this object.

        Returns
        -------
        type: `NEP.Type` enumerate
            The solver currently being used.
        """
        cdef SlepcNEPType nep_type = NULL
        CHKERR( NEPGetType(self.nep, &nep_type) )
        return bytes2str(nep_type)

    def getOptionsPrefix(self):
        """
        Gets the prefix used for searching for all NEP options in the
        database.

        Returns
        -------
        prefix: string
            The prefix string set for this NEP object.
        """
        cdef const_char *prefix = NULL
        CHKERR( NEPGetOptionsPrefix(self.nep, &prefix) )
        return bytes2str(prefix)

    def setOptionsPrefix(self, prefix):
        """
        Sets the prefix used for searching for all NEP options in the
        database.

        Parameters
        ----------
        prefix: string
            The prefix string to prepend to all NEP option requests.
        """
        cdef const_char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( NEPSetOptionsPrefix(self.nep, cval) )

    def appendOptionsPrefix(self, prefix):
        """
        Appends to the prefix used for searching for all NEP options
        in the database.

        Parameters
        ----------
        prefix: string
            The prefix string to prepend to all NEP option requests.
        """
        cdef const_char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( NEPAppendOptionsPrefix(self.nep, cval) )

    def setFromOptions(self):
        """
        Sets NEP options from the options database. This routine must
        be called before `setUp()` if the user is to be allowed to set
        the solver type.
        """
        CHKERR( NEPSetFromOptions(self.nep) )

    def getWhichEigenpairs(self):
        """
        Returns which portion of the spectrum is to be sought.

        Returns
        -------
        which: `NEP.Which` enumerate
            The portion of the spectrum to be sought by the solver.
        """
        cdef SlepcNEPWhich val = NEP_LARGEST_MAGNITUDE
        CHKERR( NEPGetWhichEigenpairs(self.nep, &val) )
        return val

    def setWhichEigenpairs(self, which):
        """
        Specifies which portion of the spectrum is to be sought.

        Parameters
        ----------
        which: `NEP.Which` enumerate
            The portion of the spectrum to be sought by the solver.
        """
        cdef SlepcNEPWhich val = which
        CHKERR( NEPSetWhichEigenpairs(self.nep, val) )

    def getTolerances(self):
        """
        Gets the tolerance and maximum iteration count used by the
        default NEP convergence tests.

        Returns
        -------
        tol: float
            The convergence tolerance.
        maxit: int
            The maximum number of iterations.
        """
        cdef PetscReal rval = 0
        cdef PetscInt  ival = 0
        CHKERR( NEPGetTolerances(self.nep, &rval, &ival) )
        return (toReal(rval), toInt(ival))

    def setTolerances(self, tol=None, maxit=None):
        """
        Sets the tolerance and maximum iteration count used in convergence tests.

        Parameters
        ----------
        tol: float, optional
            The convergence tolerance.
        maxit: int, optional
            The maximum number of iterations.
        """
        cdef PetscReal rval = PETSC_DEFAULT
        cdef PetscInt  ival = PETSC_DEFAULT
        if tol   is not None: rval = asReal(tol)
        if maxit is not None: ival = asInt(maxit)
        CHKERR( NEPSetTolerances(self.nep, rval, ival) )

    def getRIILagPreconditioner(self):
        """
        Indicates how often the preconditioner is rebuilt.

        Returns
        -------
        lag: int
            The lag parameter.
        """
        cdef PetscInt ival = 0
        CHKERR( NEPRIIGetLagPreconditioner(self.nep, &ival) )
        return ival

    def setRIILagPreconditioner(self, lag):
        """
        Determines when the preconditioner is rebuilt in the
        nonlinear solve.

        Parameters
        ----------
        lag: int
            0 indicates NEVER rebuild, 1 means rebuild every time the Jacobian is
            computed within the nonlinear iteration, 2 means every second time
            the Jacobian is built, etc.
        """
        cdef PetscInt ival = lag
        CHKERR( NEPRIISetLagPreconditioner(self.nep, ival) )

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
        CHKERR( NEPGetTrackAll(self.nep, &tval) )
        return <bint>tval

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
        CHKERR( NEPSetTrackAll(self.nep, tval) )

    def getDimensions(self):
        """
        Gets the number of eigenvalues to compute and the dimension of
        the subspace.

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
        CHKERR( NEPGetDimensions(self.nep, &ival1, &ival2, &ival3) )
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
        """
        cdef PetscInt ival1 = PETSC_DEFAULT
        cdef PetscInt ival2 = PETSC_DEFAULT
        cdef PetscInt ival3 = PETSC_DEFAULT
        if nev is not None: ival1 = asInt(nev)
        if ncv is not None: ival2 = asInt(ncv)
        if mpd is not None: ival3 = asInt(mpd)
        CHKERR( NEPSetDimensions(self.nep, ival1, ival2, ival3) )

    def getBV(self):
        """
        Obtain the basis vectors object associated to the eigensolver.

        Returns
        -------
        bv: BV
            The basis vectors context.
        """
        cdef BV bv = BV()
        CHKERR( NEPGetBV(self.nep, &bv.bv) )
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
        CHKERR( NEPSetBV(self.nep, bv.bv) )

    def getRG(self):
        """
        Obtain the region object associated to the eigensolver.

        Returns
        -------
        rg: RG
            The region context.
        """
        cdef RG rg = RG()
        CHKERR( NEPGetRG(self.nep, &rg.rg) )
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
        CHKERR( NEPSetRG(self.nep, rg.rg) )

    #

    def setInitialSpace(self, space):
        """
        Sets the initial space from which the eigensolver starts to
        iterate.

        Parameters
        ----------
        space: Vec or sequence of Vec
           The initial space
        """
        if isinstance(space, Vec): space = [space]
        cdef PetscVec *vs = NULL
        cdef Py_ssize_t i = 0, ns = len(space)
        cdef tmp = allocate(<size_t>ns*sizeof(Vec),<void**>&vs)
        for i in range(ns): vs[i] = (<Vec?>space[i]).vec
        CHKERR( NEPSetInitialSpace(self.nep, <PetscInt>ns, vs) )

    #

    def cancelMonitor(self):
        """
        Clears all monitors for a NEP object.
        """
        CHKERR( NEPMonitorCancel(self.nep) )

    #

    def setUp(self):
        """
        Sets up all the internal data structures necessary for the
        execution of the eigensolver.
        """
        CHKERR( NEPSetUp(self.nep) )

    def solve(self):
        """
        Solves the eigensystem.
        """
        CHKERR( NEPSolve(self.nep) )

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
        CHKERR( NEPGetIterationNumber(self.nep, &ival) )
        return toInt(ival)

    def getConvergedReason(self):
        """
        Gets the reason why the `solve()` iteration was stopped.

        Returns
        -------
        reason: `NEP.ConvergedReason` enumerate
            Negative value indicates diverged, positive value
            converged.
        """
        cdef SlepcNEPConvergedReason val = NEP_CONVERGED_ITERATING
        CHKERR( NEPGetConvergedReason(self.nep, &val) )
        return val


    def getConverged(self):
        """
        Gets the number of converged eigenpairs.

        Returns
        -------
        nconv: int
            Number of converged eigenpairs.
        """
        cdef PetscInt ival = 0
        CHKERR( NEPGetConverged(self.nep, &ival) )
        return toInt(ival)

    def getEigenpair(self, int i, Vec Vr=None, Vec Vi=None):
        """
        Gets the i-th solution of the eigenproblem as computed by
        `solve()`.  The solution consists of both the eigenvalue and
        the eigenvector.

        Parameters
        ----------
        i: int
            Index of the solution to be obtained.
        Vr: Vec, optional
            Placeholder for the returned eigenvector (real part).
        Vi: Vec, optional
            Placeholder for the returned eigenvector (imaginary part).

        Returns
        -------
        e: scalar (possibly complex)
            The computed eigenvalue.
        """
        cdef PetscScalar sval1 = 0
        cdef PetscScalar sval2 = 0
        cdef PetscVec vecr = NULL
        cdef PetscVec veci = NULL
        if Vr is not None: vecr = Vr.vec
        if Vi is not None: veci = Vi.vec
        CHKERR( NEPGetEigenpair(self.nep, i, &sval1, &sval2, vecr, veci) )
        return complex(toScalar(sval1), toScalar(sval2))

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
        error: real
            Error estimate.
        """
        cdef PetscReal rval = 0
        CHKERR( NEPGetErrorEstimate(self.nep, i, &rval) )
        return toReal(rval)

    def computeError(self, int i, etype=None):
        """
        Computes the error (based on the residual norm) associated with the i-th
        computed eigenpair.

        Parameters
        ----------
        i: int
            Index of the solution to be considered.
        etype: `NEP.ErrorType` enumerate
            The error type to compute.

        Returns
        -------
        error: real
            The error bound, computed in various ways from the residual norm
            ``||T(lambda)x||_2`` where ``lambda`` is the eigenvalue and
            ``x`` is the eigenvector.
        """
        cdef SlepcNEPErrorType et = NEP_ERROR_RELATIVE
        cdef PetscReal rval = 0
        if etype is not None: et = etype
        CHKERR( NEPComputeError(self.nep, i, et, &rval) )
        return toReal(rval)

    def errorView(self, etype=None, Viewer viewer=None):
        """
        Displays the errors associated with the computed solution
        (as well as the eigenvalues).

        Parameters
        ----------
        etype: `NEP.ErrorType` enumerate, optional
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
        cdef SlepcNEPErrorType et = NEP_ERROR_RELATIVE
        if etype is not None: et = etype
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( NEPErrorView(self.nep, et, vwr) )

    def setFunction(self, function, Mat F, Mat P=None, args=None, kargs=None):
        """
        Sets the function to compute the nonlinear Function T(lambda)
        as well as the location to store the matrix.

        Parameters
        ----------
        function:
            Function evaluation routine
        F: Mat
            Function matrix
        P: Mat
            preconditioner matrix (usually same as the Function)
        """
        cdef PetscMat Fmat=NULL
        if F is not None: Fmat = F.mat
        cdef PetscMat Pmat=Fmat
        if P is not None: Pmat = P.mat
        CHKERR( NEPSetFunction(self.nep, Fmat, Pmat, NEP_Function, NULL) )
        if args is None: args = ()
        if kargs is None: kargs = {}
        self.set_attr('__function__', (function, args, kargs))

    def setJacobian(self, jacobian, Mat J, args=None, kargs=None):
        """
        Sets the function to compute Jacobian T'(lambda) as well
        as the location to store the matrix.

        Parameters
        ----------
        jacobian:
            Jacobian evaluation routine
        J: Mat
            Jacobian matrix
        """
        cdef PetscMat Jmat=NULL
        if J is not None: Jmat = J.mat
        CHKERR( NEPSetJacobian(self.nep, Jmat, NEP_Jacobian, NULL) )
        if args is None: args = ()
        if kargs is None: kargs = {}
        self.set_attr('__jacobian__', (jacobian, args, kargs))

    def setSplitOperator(self, A, f, structure=None):
        """
        Sets the operator of the nonlinear eigenvalue problem
        in split form.

        Parameters
        ----------
        A: Mat or sequence of Mat
            Coefficient matrices of the split form.
        f: sequence of FN
            Scalar functions of the split form.
        structure: `PETSc.Mat.Structure` enumerate, optional
            Structure flag for matrices.
        """
        if isinstance(A, Mat): A = [A]
        if isinstance(f, FN):  f = [f]
        cdef PetscMat *As = NULL
        cdef SlepcFN  *Fs = NULL
        cdef Py_ssize_t i = 0, n = len(A)
        cdef PetscMatStructure mstr = matstructure(structure)
        assert n == len(f)
        cdef tmp1 = allocate(<size_t>n*sizeof(Mat),<void**>&As)
        cdef tmp2 = allocate(<size_t>n*sizeof(FN),<void**>&Fs)
        for i in range(n):
            As[i] = (<Mat?>A[i]).mat
            Fs[i] = (<FN?>f[i]).fn
        CHKERR( NEPSetSplitOperator(self.nep, <PetscInt>n, As, Fs, mstr) )

# -----------------------------------------------------------------------------

del NEPType
del NEPErrorType
del NEPWhich
del NEPConvergedReason
del NEPRefine
del NEPRefineScheme

# -----------------------------------------------------------------------------
