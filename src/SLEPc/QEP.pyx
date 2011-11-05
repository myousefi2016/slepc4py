# -----------------------------------------------------------------------------

class QEPType(object):
    LINEAR   = S_(QEPLINEAR)
    QARNOLDI = S_(QEPQARNOLDI)

class QEPProblemType(object):
    GENERAL    = QEP_GENERAL
    HERMITIAN  = QEP_HERMITIAN
    GYROSCOPIC = QEP_GYROSCOPIC

class QEPWhich(object):
    LARGEST_MAGNITUDE  = QEP_LARGEST_MAGNITUDE
    SMALLEST_MAGNITUDE = QEP_SMALLEST_MAGNITUDE
    LARGEST_REAL       = QEP_LARGEST_REAL
    SMALLEST_REAL      = QEP_SMALLEST_REAL
    LARGEST_IMAGINARY  = QEP_LARGEST_IMAGINARY
    SMALLEST_IMAGINARY = QEP_SMALLEST_IMAGINARY

class QEPConvergedReason(object):
    CONVERGED_TOL       = QEP_CONVERGED_TOL
    DIVERGED_ITS        = QEP_DIVERGED_ITS
    DIVERGED_BREAKDOWN  = QEP_DIVERGED_BREAKDOWN
    CONVERGED_ITERATING = QEP_CONVERGED_ITERATING
    ITERATING           = QEP_CONVERGED_ITERATING

# -----------------------------------------------------------------------------

cdef class QEP(Object):

    """
    QEP
    """

    Type            = QEPType
    ProblemType     = QEPProblemType
    Which           = QEPWhich
    ConvergedReason = QEPConvergedReason

    def __cinit__(self):
        self.obj = <PetscObject*> &self.qep
        self.qep = NULL

    def view(self, Viewer viewer=None):
        """
        Prints the QEP data structure.

        Parameters
        ----------
        viewer: Viewer, optional.
            Visualization context; if not provided, the standard
            output is used.
        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( QEPView(self.qep, vwr) )

    def destroy(self):
        """
        Destroys the QEP object.
        """
        CHKERR( QEPDestroy(&self.qep) )
        self.qep = NULL
        return self

    def reset(self):
        """
        Resets the QEP object.
        """
        CHKERR( QEPReset(self.qep) )

    def create(self, comm=None):
        """
        Creates the QEP object.

        Parameters
        ----------
        comm: Comm, optional.
            MPI communicator. If not provided, it defaults to all
            processes.
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcQEP newqep = NULL
        CHKERR( QEPCreate(ccomm, &newqep) )
        SlepcCLEAR(self.obj); self.qep = newqep
        return self

    def setType(self, qep_type):
        """
        Selects the particular solver to be used in the QEP object.

        Parameters
        ----------
        qep_type: `QEP.Type` enumerate
            T5he solver to be used.
        """
        cdef SlepcQEPType cval = NULL
        qep_type = str2bytes(qep_type, &cval)
        CHKERR( QEPSetType(self.qep, cval) )

    def getType(self):
        """
        Gets the QEP type of this object.

        Returns
        -------
        type: `QEP.Type` enumerate
            The solver currently being used.
        """
        cdef SlepcQEPType qep_type = NULL
        CHKERR( QEPGetType(self.qep, &qep_type) )
        return bytes2str(qep_type)

    def getOptionsPrefix(self):
        """
        Gets the prefix used for searching for all QEP options in the
        database.

        Returns
        -------
        prefix: string
            The prefix string set for this QEP object.
        """
        cdef const_char *prefix = NULL
        CHKERR( QEPGetOptionsPrefix(self.qep, &prefix) )
        return bytes2str(prefix)

    def setOptionsPrefix(self, prefix):
        """
        Sets the prefix used for searching for all QEP options in the
        database.

        Parameters
        ----------
        prefix: string
            The prefix string to prepend to all QEP option requests.
        """
        cdef const_char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( QEPSetOptionsPrefix(self.qep, cval) )

    def appendOptionsPrefix(self, prefix):
        """
        Appends to the prefix used for searching for all QEP options
        in the database.

        Parameters
        ----------
        prefix: string
            The prefix string to prepend to all QEP option requests.
        """
        cdef const_char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( QEPAppendOptionsPrefix(self.qep, cval) )

    def setFromOptions(self):
        """
        Sets QEP options from the options database. This routine must
        be called before `setUp()` if the user is to be allowed to set
        the solver type.
        """
        CHKERR( QEPSetFromOptions(self.qep) )

    def getProblemType(self):
        """
        Gets the problem type from the QEP object.

        Returns
        -------
        problem_type: `QEP.ProblemType` enumerate
            The problem type that was previously set.
        """
        cdef SlepcQEPProblemType val = QEP_GENERAL
        CHKERR( QEPGetProblemType(self.qep, &val) )
        return val

    def setProblemType(self, problem_type):
        """
        Specifies the type of the eigenvalue problem.

        Parameters
        ----------
        problem_type: `QEP.ProblemType` enumerate
            The problem type to be set.
        """
        cdef SlepcQEPProblemType val = problem_type
        CHKERR( QEPSetProblemType(self.qep, val) )

    def getWhichEigenpairs(self):
        """
        Returns which portion of the spectrum is to be sought.

        Returns
        -------
        which: `QEP.Which` enumerate
            The portion of the spectrum to be sought by the solver.
        """
        cdef SlepcQEPWhich val = QEP_LARGEST_MAGNITUDE
        CHKERR( QEPGetWhichEigenpairs(self.qep, &val) )
        return val

    def setWhichEigenpairs(self, which):
        """
        Specifies which portion of the spectrum is to be sought.

        Parameters
        ----------
        which: `QEP.Which` enumerate
            The portion of the spectrum to be sought by the solver.
        """
        cdef SlepcQEPWhich val = which
        CHKERR( QEPSetWhichEigenpairs(self.qep, val) )

    def getLeftVectorsWanted(self):
        """
        Returns the flag indicating whether left eigenvectors are
        required or not.

        Returns
        -------
        wanted: boolean
                Whether left eigenvectors are required or not.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( QEPGetLeftVectorsWanted(self.qep, &tval) )
        return <bint>tval

    def setLeftVectorsWanted(self, wanted):
        """
        Specifies the flag indicating whether left eigenvectors are
        required or not.

        Parameters
        ----------
        wanted: boolean
                Whether left eigenvectors are required or not.
        """
        cdef PetscBool tval = wanted
        CHKERR( QEPSetLeftVectorsWanted(self.qep, tval) )

    def getTolerances(self):
        """
        Gets the tolerance and maximum iteration count used by the
        default QEP convergence tests.

        Returns
        -------
        tol: float
            The convergence tolerance.
        max_it: int
            The maximum number of iterations
        """
        cdef PetscReal rval = 0
        cdef PetscInt  ival = 0
        CHKERR( QEPGetTolerances(self.qep, &rval, &ival) )
        return (toReal(rval), toInt(ival))

    def setTolerances(self, tol=None, max_it=None):
        """
        Sets the tolerance and maximum iteration count used by the
        default QEP convergence tests.

        Parameters
        ----------
        tol: float, optional
            The convergence tolerance.
        max_it: int, optional
            The maximum number of iterations
        """
        cdef PetscReal rval = PETSC_IGNORE
        cdef PetscInt  ival = PETSC_IGNORE
        if tol    is not None: rval = asReal(tol)
        if max_it is not None: ival = asInt(max_it)
        CHKERR( QEPSetTolerances(self.qep, rval, ival) )

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
        CHKERR( QEPGetTrackAll(self.qep, &tval) )
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
        CHKERR( QEPSetTrackAll(self.qep, tval) )

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
        CHKERR( QEPGetDimensions(self.qep, &ival1, &ival2, &ival3) )
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
        cdef PetscInt ival1 = PETSC_IGNORE
        cdef PetscInt ival2 = PETSC_IGNORE
        cdef PetscInt ival3 = PETSC_IGNORE
        if nev is not None: ival1 = asInt(nev)
        if ncv is not None: ival2 = asInt(ncv)
        if mpd is not None: ival3 = asInt(mpd)
        CHKERR( QEPSetDimensions(self.qep, ival1, ival2, ival3) )

    def getScaleFactor(self):
        """
        Gets the factor used for scaling the quadratic eigenproblem.

        Returns
        -------
        alpha: real
            The scaling factor.
        """
        cdef PetscReal rval = 0
        CHKERR( QEPGetScaleFactor(self.qep, &rval) )
        return toReal(rval)

    def setScaleFactor(self, alpha):
        """
        Sets the scaling factor to be used for scaling the quadratic problem
        before attempting to solve.

        Parameters
        ----------
        alpha: real
            The scaling factor.
        """
        cdef PetscReal rval = asReal(alpha)
        CHKERR( QEPSetScaleFactor(self.qep, rval) )

    def getIP(self):
        """
        Obtain the inner product associated to the eigensolver.

        Returns
        -------
        ip: IP
            The inner product context.
        """
        cdef IP ip = IP()
        CHKERR( QEPGetIP(self.qep, &ip.ip) )
        PetscINCREF(ip.obj)
        return ip

    def setIP(self, IP ip not None):
        """
        Associates an inner product to the eigensolver.

        Parameters
        ----------
        ip: IP
            The inner product context.
        """
        CHKERR( QEPSetIP(self.qep, ip.ip) )

    def getOperators(self):
        """
        Gets the matrices associated with the eigenvalue problem.

        Returns
        -------
        M: Mat
            The fist coefficient matrix.
        C: Mat
            The second coefficient matrix.
        K: Mat
            The third coefficient matrix.
        """
        cdef Mat M = Mat()
        cdef Mat C = Mat()
        cdef Mat K = Mat()
        CHKERR( QEPGetOperators(self.qep, &M.mat, &C.mat, &K.mat) )
        PetscINCREF(M.obj)
        PetscINCREF(C.obj)
        PetscINCREF(K.obj)
        return (M, C, K)

    def setOperators(self, Mat M not None, Mat C not None, Mat K not None):
        """
        Sets the matrices associated with the eigenvalue problem.

        Parameters
        ----------
        M: Mat
            The fist coefficient matrix.
        C: Mat
            The second coefficient matrix.
        K: Mat
            The third coefficient matrix.
        """
        CHKERR( QEPSetOperators(self.qep, M.mat, C.mat, K.mat) )

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
        CHKERR( QEPSetInitialSpace(self.qep, <PetscInt>ns, vs) )

    def setInitialSpaceLeft(self, space):
        """
        Sets the initial left space from which the solver starts to
        iterate.

        Parameters
        ----------
        space: Vec or sequence of Vec
           The initial left space
        """
        if isinstance(space, Vec): space = [space]
        cdef PetscVec *vs = NULL
        cdef Py_ssize_t i = 0, ns = len(space)
        cdef tmp = allocate(<size_t>ns*sizeof(Vec),<void**>&vs)
        for i in range(ns): vs[i] = (<Vec?>space[i]).vec
        CHKERR( QEPSetInitialSpaceLeft(self.qep, <PetscInt>ns, vs) )

    #

    def cancelMonitor(self):
        """
        Clears all monitors for a QEP object.
        """
        CHKERR( QEPMonitorCancel(self.qep) )

    #

    def setUp(self):
        """
        Sets up all the internal data structures necessary for the
        execution of the eigensolver.
        """
        CHKERR( QEPSetUp(self.qep) )

    def solve(self):
        """
        Solves the eigensystem.
        """
        CHKERR( QEPSolve(self.qep) )

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
        CHKERR( QEPGetIterationNumber(self.qep, &ival) )
        return toInt(ival)

    def getConvergedReason(self):
        """
        Gets the reason why the `solve()` iteration was stopped.

        Returns
        -------
        reason: `QEP.ConvergedReason` enumerate
            Negative value indicates diverged, positive value
            converged.
        """
        cdef SlepcQEPConvergedReason val = QEP_CONVERGED_ITERATING
        CHKERR( QEPGetConvergedReason(self.qep, &val) )
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
        CHKERR( QEPGetConverged(self.qep, &ival) )
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
        CHKERR( QEPGetEigenpair(self.qep, i, &sval1, &sval2, vecr, veci) )
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
        CHKERR( QEPGetErrorEstimate(self.qep, i, &rval) )
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
        CHKERR( QEPComputeRelativeError(self.qep, i, &rval) )
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
        CHKERR( QEPComputeResidualNorm(self.qep, i, &rval) )
        return toReal(rval)

    def getOperationCounters(self):
        """
        Gets the total number of operator applications, inner product
        operations and linear iterations used by the `QEP` object
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
        CHKERR( QEPGetOperationCounters(self.qep, &ival1, &ival2, &ival3) )
        return (toInt(ival1), toInt(ival2), toInt(ival3))

# -----------------------------------------------------------------------------

del QEPType
del QEPProblemType
del QEPWhich
del QEPConvergedReason

# -----------------------------------------------------------------------------
