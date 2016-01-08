# -----------------------------------------------------------------------------

class STType(object):
    """
    ST types

    - `SHELL`:  User-defined.
    - `SHIFT`:  Shift from origin.
    - `SINVERT`: Shift-and-invert.
    - `CAYLEY`: Cayley transform.
    - `PRECOND`: Preconditioner.
    """
    SHELL   = S_(STSHELL)
    SHIFT   = S_(STSHIFT)
    SINVERT = S_(STSINVERT)
    CAYLEY  = S_(STCAYLEY)
    PRECOND = S_(STPRECOND)

class STMatMode(object):
    """
    ST matrix mode

    - `COPY`:    A working copy of the matrix is created.
    - `INPLACE`: The operation is computed in-place.
    - `SHELL`:   The matrix ``A-sigma*B`` is handled as an
      implicit matrix.
    """
    COPY    = ST_MATMODE_COPY
    INPLACE = ST_MATMODE_INPLACE
    SHELL   = ST_MATMODE_SHELL

# -----------------------------------------------------------------------------

cdef class ST(Object):

    """
    ST
    """

    Type         = STType
    MatMode      = STMatMode

    def __cinit__(self):
        self.obj = <PetscObject*> &self.st
        self.st = NULL

    def view(self, Viewer viewer=None):
        """
        Prints the ST data structure.

        Parameters
        ----------
        viewer: Viewer, optional
                Visualization context; if not provided, the standard
                output is used.
        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( STView(self.st, vwr) )

    def destroy(self):
        """
        Destroys the ST object.
        """
        CHKERR( STDestroy(&self.st) )
        self.st = NULL
        return self

    def reset(self):
        """
        Resets the ST object.
        """
        CHKERR( STReset(self.st) )

    def create(self, comm=None):
        """
        Creates the ST object.

        Parameters
        ----------
        comm: Comm, optional
              MPI communicator; if not provided, it defaults to all
              processes.
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcST newst = NULL
        CHKERR( STCreate(ccomm, &newst) )
        SlepcCLEAR(self.obj); self.st = newst
        return self

    def setType(self, st_type):
        """
        Builds ST for a particular spectral transformation.

        Parameters
        ----------
        st_type: `ST.Type` enumerate
                 The spectral transformation to be used.

        Notes
        -----
        See `ST.Type` for available methods. The default is
        `ST.Type.SHIFT` with a zero shift.  Normally, it is best to
        use `setFromOptions()` and then set the ST type from the
        options database rather than by using this routine.  Using the
        options database provides the user with maximum flexibility in
        evaluating the different available methods.
        """
        cdef SlepcSTType cval = NULL
        st_type = str2bytes(st_type, &cval)
        CHKERR( STSetType(self.st, cval) )

    def getType(self):
        """
        Gets the ST type of this object.

        Returns
        -------
        type: `ST.Type` enumerate
              The spectral transformation currently being used.
        """
        cdef SlepcSTType st_type = NULL
        CHKERR( STGetType(self.st, &st_type) )
        return bytes2str(st_type)

    def setOptionsPrefix(self, prefix):
        """
        Sets the prefix used for searching for all ST options in the
        database.

        Parameters
        ----------
        prefix: string
                The prefix string to prepend to all ST option
                requests.

        Notes
        -----
        A hyphen (``-``) must NOT be given at the beginning of the
        prefix name.  The first character of all runtime options is
        AUTOMATICALLY the hyphen.
        """
        cdef const_char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( STSetOptionsPrefix(self.st, cval) )

    def getOptionsPrefix(self):
        """
        Gets the prefix used for searching for all ST options in the
        database.

        Returns
        -------
        prefix: string
                The prefix string set for this ST object.
        """
        cdef const_char *prefix = NULL
        CHKERR( STGetOptionsPrefix(self.st, &prefix) )
        return bytes2str(prefix)

    def setFromOptions(self):
        """
        Sets ST options from the options database. This routine must
        be called before `setUp()` if the user is to be allowed to set
        the solver type.

        Notes
        -----
        To see all options, run your program with the -help option.
        """
        CHKERR( STSetFromOptions(self.st) )

    #

    def setShift(self, shift):
        """
        Sets the shift associated with the spectral transformation.

        Parameters
        ----------
        shift: scalar (possibly complex)
               The value of the shift.

        Notes
        -----
        In some spectral transformations, changing the shift may have
        associated a lot of work, for example recomputing a
        factorization.
        """
        cdef PetscScalar sval = asScalar(shift)
        CHKERR( STSetShift(self.st, sval) )

    def getShift(self):
        """
        Gets the shift associated with the spectral transformation.

        Returns
        -------
        shift: scalar (possibly complex)
               The value of the shift.
        """
        cdef PetscScalar sval = 0
        CHKERR( STGetShift(self.st, &sval) )
        return toScalar(sval)

    def setTransform(self, flag):
        """
        Sets a flag to indicate whether the transformed matrices
        are computed or not.

        Parameters
        ----------
        flag: boolean
              This flag is intended for the case of polynomial
              eigenproblems solved via linearization.
              If this flag is False (default) the spectral transformation
              is applied to the linearization (handled by the eigensolver),
              otherwise it is applied to the original problem.
        """
        cdef PetscBool sval = flag
        CHKERR( STSetTransform(self.st, sval) )

    def getTransform(self):
        """
        Gets the flag indicating whether the transformed matrices
        are computed or not.

        Returns
        -------
        flag: boolean
               This flag is intended for the case of polynomial
               eigenproblems solved via linearization.
               If this flag is False (default) the spectral transformation
               is applied to the linearization (handled by the eigensolver),
               otherwise it is applied to the original problem.
        """
        cdef PetscBool sval = PETSC_FALSE
        CHKERR( STGetTransform(self.st, &sval) )
        return sval

    def setMatMode(self, mode):
        """
        Sets a flag to indicate how the matrix is being shifted in the
        shift-and-invert and Cayley spectral transformations.

        Parameters
        ----------
        mode: `ST.MatMode` enumerate
              The mode flag.

        Notes
        -----
        By default (`ST.MatMode.COPY`), a copy of matrix ``A`` is made
        and then this copy is shifted explicitly, e.g. ``A <- (A - s
        B)``.

        With `ST.MatMode.INPLACE`, the original matrix ``A`` is
        shifted at `setUp()` and unshifted at the end of the
        computations. With respect to the previous one, this mode
        avoids a copy of matrix ``A``. However, a backdraw is that the
        recovered matrix might be slightly different from the original
        one (due to roundoff).

        With `ST.MatMode.SHELL`, the solver works with an implicit
        shell matrix that represents the shifted matrix. This mode is
        the most efficient in creating the shifted matrix but it
        places serious limitations to the linear solves performed in
        each iteration of the eigensolver (typically, only interative
        solvers with Jacobi preconditioning can be used).

        In the case of generalized problems, in the two first modes
        the matrix ``A - s B`` has to be computed explicitly. The
        efficiency of this computation can be controlled with
        `setMatStructure()`.
        """
        cdef SlepcSTMatMode val = mode
        CHKERR( STSetMatMode(self.st, val) )

    def getMatMode(self):
        """
        Gets a flag that indicates how the matrix is being shifted in
        the shift-and-invert and Cayley spectral transformations.

        Returns
        -------
        mode: `ST.MatMode` enumerate
              The mode flag.
        """
        cdef SlepcSTMatMode val = ST_MATMODE_INPLACE
        CHKERR( STGetMatMode(self.st, &val) )
        return val

    def setOperators(self, operators):
        """
        Sets the matrices associated with the eigenvalue problem.

        Parameters
        ----------
        operators: sequence of Mat
           The matrices associated with the eigensystem.
        """
        operators = tuple(operators)
        cdef PetscMat *mats = NULL
        cdef Py_ssize_t k=0, n = len(operators)
        cdef tmp = allocate(<size_t>n*sizeof(PetscMat),<void**>&mats)
        for k from 0 <= k < n: mats[k] = (<Mat?>operators[k]).mat
        CHKERR( STSetOperators(self.st, <PetscInt>n, mats) )

    def getOperators(self):
        """
        Gets the matrices associated with the eigenvalue problem.

        Returns
        -------
        operators: tuple of Mat
           The matrices associated with the eigensystem.
        """
        cdef Mat A
        cdef PetscMat mat = NULL
        cdef PetscInt k=0, n=0
        CHKERR( STGetNumMatrices(self.st, &n) )
        cdef object operators = []
        for k from 0 <= k < n:
            CHKERR( STGetOperators(self.st, k, &mat) )
            A = Mat(); A.mat = mat; PetscINCREF(A.obj)
            operators.append(A)
        return tuple(operators)

    def setMatStructure(self, structure):
        """
        Sets an internal Mat.Structure attribute to indicate which is
        the relation of the sparsity pattern of the two matrices ``A``
        and ``B`` constituting the generalized eigenvalue
        problem. This function has no effect in the case of standard
        eigenproblems.

        Parameters
        ----------
        structure: `PETSc.Mat.Structure` enumerate
                   Either same, different, or a subset of the non-zero
                   sparsity pattern.

        Notes
        -----
        By default, the sparsity patterns are assumed to be
        different. If the patterns are equal or a subset then it is
        recommended to set this attribute for efficiency reasons (in
        particular, for internal *AXPY()* matrix operations).
        """
        cdef PetscMatStructure val = matstructure(structure)
        CHKERR( STSetMatStructure(self.st, val) )

    def setKSP(self, KSP ksp not None):
        """
        Sets the KSP object associated with the spectral
        transformation.

        Parameters
        ----------
        ksp: KSP
             The linear solver object.
        """
        CHKERR( STSetKSP(self.st, ksp.ksp) )

    def getKSP(self):
        """
        Gets the KSP object associated with the spectral
        transformation.

        Returns
        -------
        ksp: KSP
             The linear solver object.

        Notes
        -----
        On output, the internal value of KSP can be ``NULL`` if the
        combination of eigenproblem type and selected transformation
        does not require to solve a linear system of equations.
        """
        cdef KSP ksp = KSP()
        CHKERR( STGetKSP(self.st, &ksp.ksp) )
        PetscINCREF(ksp.obj)
        return ksp

    #

    def setUp(self):
        """
        Prepares for the use of a spectral transformation.
        """
        CHKERR( STSetUp(self.st) )

    def apply(self, Vec x not None, Vec y not None):
        """
        Applies the spectral transformation operator to a vector, for
        instance ``(A - sB)^-1 B`` in the case of the shift-and-invert
        tranformation and generalized eigenproblem.

        Parameters
        ----------
        x: Vec
           The input vector.
        y: Vec
           The result vector.
        """
        CHKERR( STApply(self.st, x.vec, y.vec) )

    def applyTranspose(self, Vec x not None, Vec y not None):
        """
        Applies the transpose of the operator to a vector, for
        instance ``B^T(A - sB)^-T`` in the case of the
        shift-and-invert tranformation and generalized eigenproblem.

        Parameters
        ----------
        x: Vec
           The input vector.
        y: Vec
           The result vector.
        """
        CHKERR( STApplyTranspose(self.st, x.vec, y.vec) )

    #

    def setCayleyAntishift(self, tau):
        """
        Sets the value of the anti-shift for the Cayley spectral
        transformation.

        Parameters
        ----------
        tau: scalar (possibly complex)
             The anti-shift.

        Notes
        -----

        In the generalized Cayley transform, the operator can be
        expressed as ``OP = inv(A - sigma B)*(A + tau B)``. This
        function sets the value of `tau`.  Use `setShift()` for
        setting ``sigma``.
        """
        cdef PetscScalar sval = asScalar(tau)
        CHKERR( STCayleySetAntishift(self.st, sval) )

    #

    property shift:
        def __get__(self):
            return self.getShift()
        def __set__(self, value):
            self.setShift(value)

    property mat_mode:
        def __get__(self):
            return self.getMatMode()
        def __set__(self, value):
            self.setMatMode(value)

    property ksp:
        def __get__(self):
            return self.getKSP()
        def __set__(self, value):
            self.setKSP(value)

# -----------------------------------------------------------------------------

del STType
del STMatMode

# -----------------------------------------------------------------------------
