# -----------------------------------------------------------------------------

class MFNType(object):
    """
    MFN type

    Action of a matrix function on a vector.

    - `KRYLOV`:  Restarted Krylov solver.
    - `EXPOKIT`: Implementation of the method in Expokit.
    """
    KRYLOV   = S_(MFNKRYLOV)
    EXPOKIT  = S_(MFNEXPOKIT)

class MFNConvergedReason(object):
    CONVERGED_TOL       = MFN_CONVERGED_TOL
    CONVERGED_ITS       = MFN_CONVERGED_ITS
    DIVERGED_ITS        = MFN_DIVERGED_ITS
    DIVERGED_BREAKDOWN  = MFN_DIVERGED_BREAKDOWN
    CONVERGED_ITERATING = MFN_CONVERGED_ITERATING
    ITERATING           = MFN_CONVERGED_ITERATING

# -----------------------------------------------------------------------------

cdef class MFN(Object):

    """
    MFN
    """

    Type            = MFNType
    ConvergedReason = MFNConvergedReason

    def __cinit__(self):
        self.obj = <PetscObject*> &self.mfn
        self.mfn = NULL

    def view(self, Viewer viewer=None):
        """
        Prints the MFN data structure.

        Parameters
        ----------
        viewer: Viewer, optional.
            Visualization context; if not provided, the standard
            output is used.
        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( MFNView(self.mfn, vwr) )

    def destroy(self):
        """
        Destroys the MFN object.
        """
        CHKERR( MFNDestroy(&self.mfn) )
        self.mfn = NULL
        return self

    def reset(self):
        """
        Resets the MFN object.
        """
        CHKERR( MFNReset(self.mfn) )

    def create(self, comm=None):
        """
        Creates the MFN object.

        Parameters
        ----------
        comm: Comm, optional.
            MPI communicator. If not provided, it defaults to all
            processes.
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcMFN newmfn = NULL
        CHKERR( MFNCreate(ccomm, &newmfn) )
        SlepcCLEAR(self.obj); self.mfn = newmfn
        return self

    def setType(self, mfn_type):
        """
        Selects the particular solver to be used in the MFN object.

        Parameters
        ----------
        mfn_type: `MFN.Type` enumerate
            The solver to be used.
        """
        cdef SlepcMFNType cval = NULL
        mfn_type = str2bytes(mfn_type, &cval)
        CHKERR( MFNSetType(self.mfn, cval) )

    def getType(self):
        """
        Gets the MFN type of this object.

        Returns
        -------
        type: `MFN.Type` enumerate
            The solver currently being used.
        """
        cdef SlepcMFNType mfn_type = NULL
        CHKERR( MFNGetType(self.mfn, &mfn_type) )
        return bytes2str(mfn_type)

    def getOptionsPrefix(self):
        """
        Gets the prefix used for searching for all MFN options in the
        database.

        Returns
        -------
        prefix: string
            The prefix string set for this MFN object.
        """
        cdef const_char *prefix = NULL
        CHKERR( MFNGetOptionsPrefix(self.mfn, &prefix) )
        return bytes2str(prefix)

    def setOptionsPrefix(self, prefix):
        """
        Sets the prefix used for searching for all MFN options in the
        database.

        Parameters
        ----------
        prefix: string
            The prefix string to prepend to all MFN option requests.
        """
        cdef const_char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( MFNSetOptionsPrefix(self.mfn, cval) )

    def appendOptionsPrefix(self, prefix):
        """
        Appends to the prefix used for searching for all MFN options
        in the database.

        Parameters
        ----------
        prefix: string
            The prefix string to prepend to all MFN option requests.
        """
        cdef const_char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( MFNAppendOptionsPrefix(self.mfn, cval) )

    def setFromOptions(self):
        """
        Sets MFN options from the options database. This routine must
        be called before `setUp()` if the user is to be allowed to set
        the solver type.
        """
        CHKERR( MFNSetFromOptions(self.mfn) )

    def getTolerances(self):
        """
        Gets the tolerance and maximum iteration count used by the
        default MFN convergence tests.

        Returns
        -------
        tol: float
            The convergence tolerance.
        max_it: int
            The maximum number of iterations
        """
        cdef PetscReal rval = 0
        cdef PetscInt  ival = 0
        CHKERR( MFNGetTolerances(self.mfn, &rval, &ival) )
        return (toReal(rval), toInt(ival))

    def setTolerances(self, tol=None, max_it=None):
        """
        Sets the tolerance and maximum iteration count used by the
        default MFN convergence tests.

        Parameters
        ----------
        tol: float, optional
            The convergence tolerance.
        max_it: int, optional
            The maximum number of iterations
        """
        cdef PetscReal rval = PETSC_DEFAULT
        cdef PetscInt  ival = PETSC_DEFAULT
        if tol    is not None: rval = asReal(tol)
        if max_it is not None: ival = asInt(max_it)
        CHKERR( MFNSetTolerances(self.mfn, rval, ival) )

    def getDimensions(self):
        """
        Gets the dimension of the subspace used by the solver.

        Returns
        -------
        ncv: int
            Maximum dimension of the subspace to be used by the solver.
        """
        cdef PetscInt ival = 0
        CHKERR( MFNGetDimensions(self.mfn, &ival) )
        return toInt(ival)

    def setDimensions(self, ncv):
        """
        Sets the dimension of the subspace to be used by the solver.

        Parameters
        ----------
        ncv: int
            Maximum dimension of the subspace to be used by the
            solver.
        """
        cdef PetscInt ival = asInt(ncv)
        CHKERR( MFNSetDimensions(self.mfn, ival) )

    def getFN(self):
        """
        Obtain the math function object associated to the MFN object.

        Returns
        -------
        fn: FN
            The math function context.
        """
        cdef FN fn = FN()
        CHKERR( MFNGetFN(self.mfn, &fn.fn) )
        PetscINCREF(fn.obj)
        return fn

    def setFN(self, FN fn not None):
        """
        Associates a math function object to the MFN object.

        Parameters
        ----------
        fn: FN
            The math function context.
        """
        CHKERR( MFNSetFN(self.mfn, fn.fn) )

    def getBV(self):
        """
        Obtain the basis vector object associated to the MFN object.

        Returns
        -------
        bv: BV
            The basis vectors context.
        """
        cdef BV bv = BV()
        CHKERR( MFNGetBV(self.mfn, &bv.bv) )
        PetscINCREF(bv.obj)
        return bv

    def setBV(self, BV bv not None):
        """
        Associates a basis vector object to the MFN object.

        Parameters
        ----------
        bv: BV
            The basis vectors context.
        """
        CHKERR( MFNSetBV(self.mfn, bv.bv) )

    def getOperator(self):
        """
        Gets the matrix associated with the MFN object.

        Returns
        -------
        A: Mat
            The matrix for which the matrix function is to be computed.
        """
        cdef Mat A = Mat()
        CHKERR( MFNGetOperator(self.mfn, &A.mat) )
        PetscINCREF(A.obj)
        return A

    def setOperator(self, Mat A not None):
        """
        Sets the matrix associated with the MFN object.

        Parameters
        ----------
        A: Mat
            The problem matrix.
        """
        CHKERR( MFNSetOperator(self.mfn, A.mat) )

    #

    def cancelMonitor(self):
        """
        Clears all monitors for a MFN object.
        """
        CHKERR( MFNMonitorCancel(self.mfn) )

    #

    def setUp(self):
        """
        Sets up all the internal data structures necessary for the
        execution of the eigensolver.
        """
        CHKERR( MFNSetUp(self.mfn) )

    def solve(self, Vec b not None, Vec x not None):
        """
        Solves the matrix function problem. Given a vector b, the
        vector x = f(alpha*A)*b is returned.

        Parameters
        ----------
        b: Vec
            The right hand side vector.
        x: Vec
            The solution.
        """
        CHKERR( MFNSolve(self.mfn, b.vec, x.vec) )

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
        CHKERR( MFNGetIterationNumber(self.mfn, &ival) )
        return toInt(ival)

    def getConvergedReason(self):
        """
        Gets the reason why the `solve()` iteration was stopped.

        Returns
        -------
        reason: `MFN.ConvergedReason` enumerate
            Negative value indicates diverged, positive value
            converged.
        """
        cdef SlepcMFNConvergedReason val = MFN_CONVERGED_ITERATING
        CHKERR( MFNGetConvergedReason(self.mfn, &val) )
        return val


# -----------------------------------------------------------------------------

del MFNType
del MFNConvergedReason

# -----------------------------------------------------------------------------
