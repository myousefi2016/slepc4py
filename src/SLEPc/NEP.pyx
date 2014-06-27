# -----------------------------------------------------------------------------

class NEPType(object):
    RII      = S_(NEPRII)
    SLP      = S_(NEPSLP)
    NARNOLDI = S_(NEPNARNOLDI)

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

class NEPConvergedReason(object):
    CONVERGED_FNORM_ABS      = NEP_CONVERGED_FNORM_ABS
    CONVERGED_FNORM_RELATIVE = NEP_CONVERGED_FNORM_RELATIVE
    CONVERGED_SNORM_RELATIVE = NEP_CONVERGED_SNORM_RELATIVE
    DIVERGED_LINEAR_SOLVE    = NEP_DIVERGED_LINEAR_SOLVE
    DIVERGED_FUNCTION_COUNT  = NEP_DIVERGED_FUNCTION_COUNT
    DIVERGED_MAX_IT          = NEP_DIVERGED_MAX_IT
    DIVERGED_BREAKDOWN       = NEP_DIVERGED_BREAKDOWN
    DIVERGED_FNORM_NAN       = NEP_DIVERGED_FNORM_NAN
    CONVERGED_ITERATING      = NEP_CONVERGED_ITERATING
    ITERATING                = NEP_CONVERGED_ITERATING

# -----------------------------------------------------------------------------

cdef class NEP(Object):

    """
    NEP
    """

    Type            = NEPType
    Which           = NEPWhich
    ConvergedReason = NEPConvergedReason

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
        abstol: float
            The absolute convergence tolerance.
        rtol: float
            The relative convergence tolerance.
        stol: float
            Convergence tolerance in terms of the norm of the change in the
            solution between steps, || delta x || < stol*|| x ||.
        maxit: int
            The maximum number of iterations.
        maxf: int
            The maximum number of function evaluations.
        """
        cdef PetscReal rval1 = 0
        cdef PetscReal rval2 = 0
        cdef PetscReal rval3 = 0
        cdef PetscInt  ival1 = 0
        cdef PetscInt  ival2 = 0
        CHKERR( NEPGetTolerances(self.nep, &rval1, &rval2, &rval3, &ival1, &ival2) )
        return (toReal(rval1), toReal(rval2), toReal(rval3), toInt(ival1), toInt(ival2))

    def setTolerances(self, abstol=None, rtol=None, stol=None, maxit=None, maxf=None):
        """
        Sets various parameters used in convergence tests.

        Parameters
        ----------
        abstol: float, optional
            The absolute convergence tolerance.
        rtol: float, optional
            The relative convergence tolerance.
        stol: float, optional
            Convergence tolerance in terms of the norm of the change in the
            solution between steps, || delta x || < stol*|| x ||.
        maxit: int, optional
            The maximum number of iterations.
        maxf: int, optional
            The maximum number of function evaluations.
        """
        cdef PetscReal rval1 = PETSC_DECIDE
        cdef PetscReal rval2 = PETSC_DECIDE
        cdef PetscReal rval3 = PETSC_DECIDE
        cdef PetscInt  ival1 = PETSC_DECIDE
        cdef PetscInt  ival2 = PETSC_DECIDE
        CHKERR( NEPGetTolerances(self.nep, &rval1, &rval2, &rval3, &ival1, &ival2) )
        if abstol is not None: rval1 = asReal(abstol)
        if rtol   is not None: rval2 = asReal(rtol)
        if stol   is not None: rval3 = asReal(stol)
        if maxit  is not None: ival1 = asInt(maxit)
        if maxf   is not None: ival2 = asInt(maxf)
        CHKERR( NEPSetTolerances(self.nep, rval1, rval2, rval3, ival1, ival2) )

    def getLagPreconditioner(self):
        """
        Indicates how often the preconditioner is rebuilt.

        Returns
        -------
        lag: int
            The lag parameter.
        """
        cdef PetscInt ival = 0
        CHKERR( NEPGetLagPreconditioner(self.nep, &ival) )
        return ival

    def setLagPreconditioner(self, lag):
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
        CHKERR( NEPSetLagPreconditioner(self.nep, ival) )

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
        cdef PetscInt ival1 = PETSC_DECIDE
        cdef PetscInt ival2 = PETSC_DECIDE
        cdef PetscInt ival3 = PETSC_DECIDE
        CHKERR( NEPGetDimensions(self.nep, &ival1, &ival2, &ival3) )
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

    def getEigenpair(self, int i, Vec V=None):
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
        cdef PetscScalar sval = 0
        cdef PetscVec vec = NULL
        if V is not None: vec = V.vec
        CHKERR( NEPGetEigenpair(self.nep, i, &sval, vec) )
        return toScalar(sval)

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

    def computeRelativeError(self, int i):
        """
        Computes the relative error bound associated with the i-th
        computed eigenpair.

        Parameters
        ----------
        i: int
            Index of the solution to be considered.

        Returns
        -------
        error: real
            The relative error bound.
        """
        cdef PetscReal rval = 0
        CHKERR( NEPComputeRelativeError(self.nep, i, &rval) )
        return toReal(rval)

    def computeResidualNorm(self, int i):
        """
        Computes the norm of the residual vector associated with the
        i-th computed eigenpair.

        Parameters
        ----------
        i: int
            Index of the solution to be considered.

        Returns
        -------
        norm: real
            The residual norm.
        """
        cdef PetscReal rval = 0
        CHKERR( NEPComputeResidualNorm(self.nep, i, &rval) )
        return toReal(rval)

    def getOperationCounters(self):
        """
        Gets the total number of operator applications, inner product
        operations and linear iterations used by the `NEP` object
        during the last `solve()` call.

        Returns
        -------
        ops: int
            number of operator applications.
        dots: int
            number of inner product operations.
        lits: int
            number of linear iterations.
        """
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        cdef PetscInt ival3 = 0
        CHKERR( NEPGetOperationCounters(self.nep, &ival1, &ival2, &ival3) )
        return (toInt(ival1), toInt(ival2), toInt(ival3))

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
        structure: `Mat.Structure` enumerate, optional
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
del NEPWhich
del NEPConvergedReason

# -----------------------------------------------------------------------------
