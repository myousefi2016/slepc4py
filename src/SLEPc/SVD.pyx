# -----------------------------------------------------------------------------

class SVDType(object):
    """
    SVD types

    - `CROSS`:     Eigenproblem with the cross-product matrix.
    - `CYCLIC`:    Eigenproblem with the cyclic matrix.
    - `LAPACK`:    Wrappers to dense SVD solvers in Lapack.
    - `LANCZOS`:   Lanczos.
    - `TRLANCZOS`: Thick-restart Lanczos.
    """
    CROSS     = S_(SVDCROSS)
    CYCLIC    = S_(SVDCYCLIC)
    LAPACK    = S_(SVDLAPACK)
    LANCZOS   = S_(SVDLANCZOS)
    TRLANCZOS = S_(SVDTRLANCZOS)

class SVDErrorType(object):
    """
    SVD error type to assess accuracy of computed solutions

    - `ABSOLUTE`:  Absolute error.
    - `RELATIVE`:  Relative error.
    """
    ABSOLUTE = SVD_ERROR_ABSOLUTE
    RELATIVE = SVD_ERROR_RELATIVE

class SVDWhich(object):
    """
    SVD desired piece of spectrum

    - `LARGEST`:  largest singular values.
    - `SMALLEST`: smallest singular values.
    """
    LARGEST  = SVD_LARGEST
    SMALLEST = SVD_SMALLEST

class SVDConvergedReason(object):
    """
    SVD convergence reasons

    - `CONVERGED_TOL`:
    - `CONVERGED_USER`:
    - `DIVERGED_ITS`:
    - `DIVERGED_BREAKDOWN`:
    - `CONVERGED_ITERATING`:
    """
    CONVERGED_TOL       = SVD_CONVERGED_TOL
    CONVERGED_USER      = SVD_CONVERGED_USER
    DIVERGED_ITS        = SVD_DIVERGED_ITS
    DIVERGED_BREAKDOWN  = SVD_DIVERGED_BREAKDOWN
    CONVERGED_ITERATING = SVD_CONVERGED_ITERATING
    ITERATING           = SVD_CONVERGED_ITERATING

# -----------------------------------------------------------------------------

cdef class SVD(Object):

    """
    SVD
    """

    Type            = SVDType
    ErrorType       = SVDErrorType
    Which           = SVDWhich
    ConvergedReason = SVDConvergedReason

    def __cinit__(self):
        self.obj = <PetscObject*> &self.svd
        self.svd = NULL

    def view(self, Viewer viewer=None):
        """
        Prints the SVD data structure.

        Parameters
        ----------
        viewer: Viewer, optional
                Visualization context; if not provided, the standard
                output is used.
        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( SVDView(self.svd, vwr) )

    def destroy(self):
        """
        Destroys the SVD object.
        """
        CHKERR( SVDDestroy(&self.svd) )
        self.svd = NULL
        return self

    def reset(self):
        """
        Resets the SVD object.
        """
        CHKERR( SVDReset(self.svd) )

    def create(self, comm=None):
        """
        Creates the SVD object.

        Parameters
        ----------
        comm: Comm, optional
              MPI communicator; if not provided, it defaults to all
              processes.
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcSVD newsvd = NULL
        CHKERR( SVDCreate(ccomm, &newsvd) )
        SlepcCLEAR(self.obj); self.svd = newsvd
        return self

    def setType(self, svd_type):
        """
        Selects the particular solver to be used in the SVD object.

        Parameters
        ----------
        svd_type: `SVD.Type` enumerate
                  The solver to be used.

        Notes
        -----
        See `SVD.Type` for available methods. The default is CROSS.
        Normally, it is best to use `setFromOptions()` and then set
        the SVD type from the options database rather than by using
        this routine.  Using the options database provides the user
        with maximum flexibility in evaluating the different available
        methods.
        """
        cdef SlepcSVDType cval = NULL
        svd_type = str2bytes(svd_type, &cval)
        CHKERR( SVDSetType(self.svd, cval) )

    def getType(self):
        """
        Gets the SVD type of this object.

        Returns
        -------
        type: `SVD.Type` enumerate
              The solver currently being used.
        """
        cdef SlepcSVDType svd_type = NULL
        CHKERR( SVDGetType(self.svd, &svd_type) )
        return bytes2str(svd_type)

    def getOptionsPrefix(self):
        """
        Gets the prefix used for searching for all SVD options in the
        database.

        Returns
        -------
        prefix: string
                The prefix string set for this SVD object.
        """
        cdef const_char *prefix = NULL
        CHKERR( SVDGetOptionsPrefix(self.svd, &prefix) )
        return bytes2str(prefix)

    def setOptionsPrefix(self, prefix):
        """
        Sets the prefix used for searching for all SVD options in the
        database.

        Parameters
        ----------
        prefix: string
                The prefix string to prepend to all SVD option
                requests.

        Notes
        -----
        A hyphen (-) must NOT be given at the beginning of the prefix
        name.  The first character of all runtime options is
        AUTOMATICALLY the hyphen.

        For example, to distinguish between the runtime options for
        two different SVD contexts, one could call::

            S1.setOptionsPrefix("svd1_")
            S2.setOptionsPrefix("svd2_")
        """
        cdef const_char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( SVDSetOptionsPrefix(self.svd, cval) )

    def appendOptionsPrefix(self, prefix):
        """
        Appends to the prefix used for searching for all SVD options
        in the database.

        Parameters
        ----------
        prefix: string
                The prefix string to prepend to all SVD option requests.
        """
        cdef const_char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( SVDAppendOptionsPrefix(self.svd, cval) )

    def setFromOptions(self):
        """
        Sets SVD options from the options database. This routine must
        be called before `setUp()` if the user is to be allowed to set
        the solver type.

        Notes
        -----
        To see all options, run your program with the ``-help``
        option.
        """
        CHKERR( SVDSetFromOptions(self.svd) )

    #

    def getImplicitTranspose(self):
        """
        Gets the mode used to handle the transpose of the matrix
        associated with the singular value problem.

        Returns
        -------
        impl: boolean
              How to handle the transpose (implicitly or not).
        """
        cdef PetscBool val = PETSC_FALSE
        CHKERR( SVDGetImplicitTranspose(self.svd, &val) )
        return val

    def setImplicitTranspose(self, mode):
        """
        Indicates how to handle the transpose of the matrix
        associated with the singular value problem.

        Parameters
        ----------
        impl: boolean
              How to handle the transpose (implicitly or not).

        Notes
        -----
        By default, the transpose of the matrix is explicitly built
        (if the matrix has defined the MatTranspose operation).

        If this flag is set to true, the solver does not build the
        transpose, but handles it implicitly via MatMultTranspose().
        """
        cdef PetscBool val = mode
        CHKERR( SVDSetImplicitTranspose(self.svd, val) )

    def getWhichSingularTriplets(self):
        """
        Returns which singular triplets are to be sought.

        Returns
        -------
        which: `SVD.Which` enumerate
               The singular values to be sought (either largest or
               smallest).
        """
        cdef SlepcSVDWhich val = SVD_LARGEST
        CHKERR( SVDGetWhichSingularTriplets(self.svd, &val) )
        return val

    def setWhichSingularTriplets(self, which):
        """
        Specifies which singular triplets are to be sought.

        Parameters
        ----------
        which: `SVD.Which` enumerate
               The singular values to be sought (either largest or
               smallest).
        """
        cdef SlepcSVDWhich val = which
        CHKERR( SVDSetWhichSingularTriplets(self.svd, val) )
    #

    def getTolerances(self):
        """
        Gets the tolerance and maximum iteration count used by the
        default SVD convergence tests.

        Returns
        -------
        tol: float
             The convergence tolerance.
        max_it: int
             The maximum number of iterations
        """
        cdef PetscReal rval = 0
        cdef PetscInt  ival = 0
        CHKERR( SVDGetTolerances(self.svd, &rval, &ival) )
        return (toReal(rval), toInt(ival))

    def setTolerances(self, tol=None, max_it=None):
        """
        Sets the tolerance and maximum iteration count used by the
        default SVD convergence tests.

        Parameters
        ----------
        tol: float, optional
             The convergence tolerance.
        max_it: int, optional
             The maximum number of iterations

        Notes
        -----
        Use `DECIDE` for `max_it` to assign a reasonably good value,
        which is dependent on the solution method.
        """
        cdef PetscReal rval = PETSC_DEFAULT
        cdef PetscInt  ival = PETSC_DEFAULT
        if tol    is not None: rval = asReal(tol)
        if max_it is not None: ival = asInt(max_it)
        CHKERR( SVDSetTolerances(self.svd, rval, ival) )

    def getDimensions(self):
        """
        Gets the number of singular values to compute and the
        dimension of the subspace.

        Returns
        -------
        nsv: int
             Number of singular values to compute.
        ncv: int
             Maximum dimension of the subspace to be used by the
             solver.
        mpd: int
             Maximum dimension allowed for the projected problem.
        """
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        cdef PetscInt ival3 = 0
        CHKERR( SVDGetDimensions(self.svd, &ival1, &ival2, &ival3) )
        return (toInt(ival1), toInt(ival2), toInt(ival3))

    def setDimensions(self, nsv=None, ncv=None, mpd=None):
        """
        Sets the number of singular values to compute and the
        dimension of the subspace.

        Parameters
        ----------
        nsv: int, optional
             Number of singular values to compute.
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

         - In cases where `nsv` is small, the user sets `ncv`
           (a reasonable default is 2 * `nsv`).
         - In cases where `nsv` is large, the user sets `mpd`.

        The value of `ncv` should always be between `nsv` and (`nsv` +
        `mpd`), typically `ncv` = `nsv` + `mpd`. If `nsv` is not too
        large, `mpd` = `nsv` is a reasonable choice, otherwise a
        smaller value should be used.
        """
        cdef PetscInt ival1 = PETSC_DEFAULT
        cdef PetscInt ival2 = PETSC_DEFAULT
        cdef PetscInt ival3 = PETSC_DEFAULT
        if nsv is not None: ival1 = asInt(nsv)
        if ncv is not None: ival2 = asInt(ncv)
        if mpd is not None: ival3 = asInt(mpd)
        CHKERR( SVDSetDimensions(self.svd, ival1, ival2, ival3) )

    def getBV(self):
        """
        Obtain the basis vectors objects associated to the SVD object.

        Returns
        -------
        V: BV
            The basis vectors context for right singular vectors.
        U: BV
            The basis vectors context for left singular vectors.
        """
        cdef BV V = BV()
        cdef BV U = BV()
        CHKERR( SVDGetBV(self.svd, &V.bv, &U.bv) )
        PetscINCREF(V.obj)
        PetscINCREF(U.obj)
        return (V,U)

    def setBV(self, BV V not None,BV U=None):
        """
        Associates basis vectors objects to the SVD solver.

        Parameters
        ----------
        V: BV
            The basis vectors context for right singular vectors.
        U: BV
            The basis vectors context for left singular vectors.
        """
        cdef SlepcBV VBV = V.bv
        cdef SlepcBV UBV = NULL
        if U is not None: UBV = U.bv
        CHKERR( SVDSetBV(self.svd, VBV, UBV) )

    def getOperator(self):
        """
        Gets the matrix associated with the singular value problem.

        Returns
        -------
        A: Mat
           The matrix associated with the singular value problem.
        """
        cdef Mat A = Mat()
        CHKERR( SVDGetOperator(self.svd, &A.mat) )
        PetscINCREF(A.obj)
        return A

    def setOperator(self, Mat A not None):
        """
        Sets the matrix associated with the singular value problem.

        Parameters
        ----------
        A: Mat
           The matrix associated with the singular value problem.
        """
        CHKERR( SVDSetOperator(self.svd, A.mat) )

    #

    def setInitialSpace(self, space):
        """
        Sets the initial space from which the SVD solver starts to
        iterate.

        Parameters
        ----------
        space: an sequence of Vec
           The initial space.
        """
        if isinstance(space, Vec): space = [space]
        cdef PetscVec *vs = NULL
        cdef Py_ssize_t i = 0, ns = len(space)
        cdef tmp = allocate(<size_t>ns*sizeof(Vec),<void**>&vs)
        for i in range(ns): vs[i] = (<Vec?>space[i]).vec
        CHKERR( SVDSetInitialSpace(self.svd, <PetscInt>ns, vs) )

    #

    def cancelMonitor(self):
        """
        Clears all monitors for an SVD object.
        """
        CHKERR( SVDMonitorCancel(self.svd) )

    #

    def setUp(self):
        """
        Sets up all the internal data structures necessary for the
        execution of the singular value solver.

        Notes
        -----
        This function need not be called explicitly in most cases,
        since `solve()` calls it. It can be useful when one wants to
        measure the set-up time separately from the solve time.
        """
        CHKERR( SVDSetUp(self.svd) )

    def solve(self):
        """
        Solves the singular value problem.
        """
        CHKERR( SVDSolve(self.svd) )

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
        CHKERR( SVDGetIterationNumber(self.svd, &ival) )
        return toInt(ival)

    def getConvergedReason(self):
        """
        Gets the reason why the `solve()` iteration was stopped.

        Returns
        -------
        reason: `SVD.ConvergedReason` enumerate
                Negative value indicates diverged, positive value
                converged.
        """
        cdef SlepcSVDConvergedReason val = SVD_CONVERGED_ITERATING
        CHKERR( SVDGetConvergedReason(self.svd, &val) )
        return val

    def getConverged(self):
        """
        Gets the number of converged singular triplets.

        Returns
        -------
        nconv: int
               Number of converged singular triplets.

        Notes
        -----
        This function should be called after `solve()` has finished.
        """
        cdef PetscInt ival = 0
        CHKERR( SVDGetConverged(self.svd, &ival) )
        return toInt(ival)

    def getValue(self, int i):
        """
        Gets the i-th singular value as computed by `solve()`.

        Parameters
        ----------
        i: int
           Index of the solution to be obtained.

        Returns
        -------
        s: float
           The computed singular value.

        Notes
        -----
        The index ``i`` should be a value between ``0`` and
        ``nconv-1`` (see `getConverged()`. Singular triplets are
        indexed according to the ordering criterion established with
        `setWhichSingularTriplets()`.
        """
        cdef PetscReal rval = 0
        CHKERR( SVDGetSingularTriplet(self.svd, i, &rval, NULL, NULL) )
        return toReal(rval)

    def getVectors(self, int i, Vec U not None, Vec V not None):
        """
        Gets the i-th left and right singular vectors as computed by
        `solve()`.

        Parameters
        ----------
        i: int
           Index of the solution to be obtained.
        U: Vec
           Placeholder for the returned left singular vector.
        V: Vec
           Placeholder for the returned right singular vector.

        Notes
        -----
        The index ``i`` should be a value between ``0`` and
        ``nconv-1`` (see `getConverged()`. Singular triplets are
        indexed according to the ordering criterion established with
        `setWhichSingularTriplets()`.
        """
        cdef PetscReal dummy = 0
        CHKERR( SVDGetSingularTriplet(self.svd, i, &dummy, U.vec, V.vec) )

    def getSingularTriplet(self, int i, Vec U=None, Vec V=None):
        """
        Gets the i-th triplet of the singular value decomposition as
        computed by `solve()`. The solution consists of the singular
        value and its left and right singular vectors.

        Parameters
        ----------
        i: int
           Index of the solution to be obtained.
        U: Vec
           Placeholder for the returned left singular vector.
        V: Vec
           Placeholder for the returned right singular vector.

        Returns
        -------
        s: float
           The computed singular value.

        Notes
        -----
        The index ``i`` should be a value between ``0`` and
        ``nconv-1`` (see `getConverged()`. Singular triplets are
        indexed according to the ordering criterion established with
        `setWhichSingularTriplets()`.
        """
        cdef PetscReal rval = 0
        cdef PetscVec Uvec = NULL
        cdef PetscVec Vvec = NULL
        if U is not None: Uvec = U.vec
        if V is not None: Vvec = V.vec
        CHKERR( SVDGetSingularTriplet(self.svd, i, &rval, Uvec, Vvec) )
        return toReal(rval)

    #

    def computeError(self, int i, etype=None):
        """
        Computes the error (based on the residual norm) associated with the i-th
        singular triplet.

        Parameters
        ----------
        i: int
           Index of the solution to be considered.
        etype: `SVD.ErrorType` enumerate
           The error type to compute.

        Returns
        -------
        e: real
           The relative error bound, computed in various ways from the residual norm
           ``sqrt(n1^2+n2^2)`` where ``n1 = ||A*v-sigma*u||_2``,
           ``n2 = ||A^T*u-sigma*v||_2``, ``sigma`` is the singular
           value, ``u`` and ``v`` are the left and right singular
           vectors.

        Notes
        -----
        The index ``i`` should be a value between ``0`` and
        ``nconv-1`` (see `getConverged()`).
        """
        cdef SlepcSVDErrorType et = SVD_ERROR_RELATIVE
        cdef PetscReal rval = 0
        if etype is not None: et = etype
        CHKERR( SVDComputeError(self.svd, i, et, &rval) )
        return toReal(rval)

    def errorView(self, etype=None, Viewer viewer=None):
        """
        Displays the errors associated with the computed solution
        (as well as the eigenvalues).

        Parameters
        ----------
        etype: `SVD.ErrorType` enumerate, optional
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
        cdef SlepcSVDErrorType et = SVD_ERROR_RELATIVE
        if etype is not None: et = etype
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( SVDErrorView(self.svd, et, vwr) )

    #

    def setCrossEPS(self, EPS eps not None):
        """
        Associate an eigensolver object (`EPS`) to the singular value
        solver.

        Parameters
        ----------
        eps: EPS
             The eigensolver object.
        """
        CHKERR( SVDCrossSetEPS(self.svd, eps.eps) )

    def getCrossEPS(self):
        """
        Retrieve the eigensolver object (`EPS`) associated to the
        singular value solver.

        Returns
        -------
        eps: EPS
             The eigensolver object.
        """
        cdef EPS eps = EPS()
        CHKERR( SVDCrossGetEPS(self.svd, &eps.eps) )
        PetscINCREF(eps.obj)
        return eps

    def setCyclicEPS(self, EPS eps not None):
        """
        Associate an eigensolver object (`EPS`) to the singular value
        solver.

        Parameters
        ----------
        eps: EPS
             The eigensolver object.
        """
        CHKERR( SVDCyclicSetEPS(self.svd, eps.eps) )

    def getCyclicEPS(self):
        """
        Retrieve the eigensolver object (`EPS`) associated to the
        singular value solver.

        Returns
        -------
        eps: EPS
             The eigensolver object.
        """
        cdef EPS eps = EPS()
        CHKERR( SVDCyclicGetEPS(self.svd, &eps.eps) )
        PetscINCREF(eps.obj)
        return eps

    def setCyclicExplicitMatrix(self, flag=True):
        """
        Indicate if the eigensolver operator ``H(A) = [ 0 A ; A^T 0
        ]`` must be computed explicitly.

        Parameters
        ----------
        flag: boolean
              True if ``H(A)`` is built explicitly.
        """
        cdef PetscBool tval = PETSC_FALSE
        if flag: tval = PETSC_TRUE
        CHKERR( SVDCyclicSetExplicitMatrix(self.svd, tval) )

    def getCyclicExplicitMatrix(self):
        """
        Returns the flag indicating if ``H(A) = [ 0 A ; A^T 0 ]`` is
        built explicitly.

        Returns
        -------
        flag: boolean
              True if ``H(A)`` is built explicitly.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( SVDCyclicGetExplicitMatrix(self.svd, &tval) )
        return <bint>tval

    def setLanczosOneSide(self, flag=True):
        """
        Indicate if the variant of the Lanczos method to be used is
        one-sided or two-sided.

        Parameters
        ----------
        flag: boolean
              True if the method is one-sided.

        Notes
        -----
        By default, a two-sided variant is selected, which is
        sometimes slightly more robust. However, the one-sided variant
        is faster because it avoids the orthogonalization associated
        to left singular vectors. It also saves the memory required
        for storing such vectors.
        """
        cdef PetscBool tval = PETSC_FALSE
        if flag: tval = PETSC_TRUE
        CHKERR( SVDLanczosSetOneSide(self.svd, tval) )

    def setTRLanczosOneSide(self, flag=True):
        """
        Indicate if the variant of the thick-restart Lanczos method to
        be used is one-sided or two-sided.

        Parameters
        ----------
        flag: boolean
              True if the method is one-sided.

        Notes
        -----
        By default, a two-sided variant is selected, which is
        sometimes slightly more robust. However, the one-sided variant
        is faster because it avoids the orthogonalization associated
        to left singular vectors.
        """
        cdef PetscBool tval = PETSC_FALSE
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

    property bv:
        def __get__(self):
            return self.getBV()
        def __set__(self, value):
            self.setBV(value)

# -----------------------------------------------------------------------------

del SVDType
del SVDErrorType
del SVDWhich
del SVDConvergedReason

# -----------------------------------------------------------------------------
