# -----------------------------------------------------------------------------

class DSType(object):
    """
    DS type
    """
    HEP     = S_(DSHEP)
    NHEP    = S_(DSNHEP)
    GHEP    = S_(DSGHEP)
    GHIEP   = S_(DSGHIEP)
    GNHEP   = S_(DSGNHEP)
    SVD     = S_(DSSVD)
    PEP     = S_(DSPEP)
    NEP     = S_(DSNEP)

class DSStateType(object):
    """
    DS state types

    - `RAW`: Not processed yet.
    - `INTERMEDIATE`: Reduced to Hessenberg or tridiagonal form (or equivalent).
    - `CONDENSED`: Reduced to Schur or diagonal form (or equivalent).
    - `TRUNCATED`: Condensed form truncated to a smaller size.
    """
    RAW          = DS_STATE_RAW
    INTERMEDIATE = DS_STATE_INTERMEDIATE
    CONDENSED    = DS_STATE_CONDENSED
    TRUNCATED    = DS_STATE_TRUNCATED

class DSMatType(object):
    """
    To refer to one of the matrices stored internally in DS

    - `A`:  first matrix of eigenproblem/singular value problem.
    - `B`:  second matrix of a generalized eigenproblem.
    - `C`:  third matrix of a quadratic eigenproblem.
    - `T`:  tridiagonal matrix.
    - `D`:  diagonal matrix.
    - `Q`:  orthogonal matrix of (right) Schur vectors.
    - `Z`:  orthogonal matrix of left Schur vectors.
    - `X`:  right eigenvectors.
    - `Y`:  left eigenvectors.
    - `U`:  left singular vectors.
    - `VT`: right singular vectors.
    - `W`:  workspace matrix.
    """
    A  = DS_MAT_A
    B  = DS_MAT_B
    C  = DS_MAT_C
    T  = DS_MAT_T
    D  = DS_MAT_D
    Q  = DS_MAT_Q
    Z  = DS_MAT_Z
    X  = DS_MAT_X
    Y  = DS_MAT_Y
    U  = DS_MAT_U
    VT = DS_MAT_VT
    W  = DS_MAT_W

# -----------------------------------------------------------------------------

cdef class DS(Object):

    """
    DS
    """

    Type      = DSType
    StateType = DSStateType
    MatType   = DSMatType

    def __cinit__(self):
        self.obj = <PetscObject*> &self.ds
        self.ds = NULL

    def view(self, Viewer viewer=None):
        """
        Prints the DS data structure.

        Parameters
        ----------
        viewer: Viewer, optional
                Visualization context; if not provided, the standard
                output is used.
        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( DSView(self.ds, vwr) )

    def destroy(self):
        """
        Destroys the DS object.
        """
        CHKERR( DSDestroy(&self.ds) )
        self.ds = NULL
        return self

    def reset(self):
        """
        Resets the DS object.
        """
        CHKERR( DSReset(self.ds) )

    def create(self, comm=None):
        """
        Creates the DS object.

        Parameters
        ----------
        comm: Comm, optional
              MPI communicator; if not provided, it defaults to all
              processes.
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcDS newds = NULL
        CHKERR( DSCreate(ccomm, &newds) )
        SlepcCLEAR(self.obj); self.ds = newds
        return self

    def setType(self, ds_type):
        """
        Selects the type for the DS object.

        Parameters
        ----------
        ds_type: `DS.Type` enumerate
                  The direct solver type to be used.
        """
        cdef SlepcDSType cval = NULL
        ds_type = str2bytes(ds_type, &cval)
        CHKERR( DSSetType(self.ds, cval) )

    def getType(self):
        """
        Gets the DS type of this object.

        Returns
        -------
        type: `DS.Type` enumerate
              The direct solver type currently being used.
        """
        cdef SlepcDSType ds_type = NULL
        CHKERR( DSGetType(self.ds, &ds_type) )
        return bytes2str(ds_type)

    def setOptionsPrefix(self, prefix):
        """
        Sets the prefix used for searching for all DS options in the
        database.

        Parameters
        ----------
        prefix: string
                The prefix string to prepend to all DS option
                requests.

        Notes
        -----
        A hyphen (``-``) must NOT be given at the beginning of the
        prefix name.  The first character of all runtime options is
        AUTOMATICALLY the hyphen.
        """
        cdef const_char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( DSSetOptionsPrefix(self.ds, cval) )

    def getOptionsPrefix(self):
        """
        Gets the prefix used for searching for all DS options in the
        database.

        Returns
        -------
        prefix: string
                The prefix string set for this DS object.
        """
        cdef const_char *prefix = NULL
        CHKERR( DSGetOptionsPrefix(self.ds, &prefix) )
        return bytes2str(prefix)

    def setFromOptions(self):
        """
        Sets DS options from the options database.

        Notes
        -----
        To see all options, run your program with the ``-help``
        option.
        """
        CHKERR( DSSetFromOptions(self.ds) )

    #

    def allocate(self, ld):
        """
        Allocates memory for internal storage or matrices in DS.

        Parameters
        ----------
        ld: integer
            Leading dimension (maximum allowed dimension for the
            matrices, including the extra row if present).
        """
        cdef PetscInt val = ld
        CHKERR( DSAllocate(self.ds, val) )

    def getLeadingDimension(self):
        """
        Returns the leading dimension of the allocated matrices.

        Returns
        -------
        ld: integer
            Leading dimension (maximum allowed dimension for the matrices).
        """
        cdef PetscInt val = 0
        CHKERR( DSGetLeadingDimension(self.ds, &val) )
        return val

    def setState(self, state):
        """
        Change the state of the DS object.

        Parameters
        ----------
        state: `DS.StateType` enumerate
               The new state.

        Notes
        -----
        The state indicates that the dense system is in an initial
        state (raw), in an intermediate state (such as tridiagonal,
        Hessenberg or Hessenberg-triangular), in a condensed state
        (such as diagonal, Schur or generalized Schur), or in a
        truncated state.

        This function is normally used to return to the raw state when
        the condensed structure is destroyed.
        """
        cdef SlepcDSStateType val = state
        CHKERR( DSSetState(self.ds, val) )

    def getState(self):
        """
        Returns the current state.

        Returns
        -------
        state: `DS.StateType` enumerate
               The current state.
        """
        cdef SlepcDSStateType val = DS_STATE_RAW
        CHKERR( DSGetState(self.ds, &val) )
        return val

    def setDimensions(self, n=None, m=None, l=None, k=None):
        """
        Resize the matrices in the DS object.

        Parameters
        ----------
        n: int, optional
           The new size.
        m: int, optional
           The new column size (only for SVD).
        l: int, optional
           Number of locked (inactive) leading columns.
        k: int, optional
           Intermediate dimension (e.g., position of arrow).

        Notes
        -----
        The internal arrays are not reallocated.

        The value `m` is not used except in the case of DS.SVD.
        """
        cdef PetscInt ival1 = PETSC_DEFAULT
        cdef PetscInt ival2 = PETSC_DEFAULT
        cdef PetscInt ival3 = 0
        cdef PetscInt ival4 = 0
        if n is not None: ival1 = asInt(n)
        if m is not None: ival2 = asInt(m)
        if l is not None: ival3 = asInt(l)
        if k is not None: ival4 = asInt(k)
        CHKERR( DSSetDimensions(self.ds, ival1, ival2, ival3, ival4) )

    def getDimensions(self):
        """
        Returns the current dimensions.

        Returns
        -------
        n: int
           The new size.
        m: int
           The new column size (only for SVD).
        l: int
           Number of locked (inactive) leading columns.
        k: int
           Intermediate dimension (e.g., position of arrow).
        t: int
           Truncated length.
        """
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        cdef PetscInt ival3 = 0
        cdef PetscInt ival4 = 0
        cdef PetscInt ival5 = 0
        CHKERR( DSGetDimensions(self.ds, &ival1, &ival2, &ival3, &ival4, &ival5) )
        return (toInt(ival1), toInt(ival2), toInt(ival3), toInt(ival4), toInt(ival5))

    def setMethod(self, meth):
        """
        Selects the method to be used to solve the problem.

        Parameters
        ----------
        meth: int
              An index indentifying the method.
        """
        cdef PetscInt val = meth
        CHKERR( DSSetMethod(self.ds, val) )

    def getMethod(self):
        """
        Gets the method currently used in the DS.

        Returns
        -------
        meth: int
              Identifier of the method.
        """
        cdef PetscInt val = 0
        CHKERR( DSGetMethod(self.ds, &val) )
        return val

    def setCompact(self, comp):
        """
        Switch to compact storage of matrices.

        Parameters
        ----------
        comp: boolean
              A boolean flag.

        Notes
        -----
        Compact storage is used in some `DS` types such as
        `DS.Type.HEP` when the matrix is tridiagonal. This flag
        can be used to indicate whether the user provides the
        matrix entries via the compact form (the tridiagonal
        `DS.MatType.T`) or the non-compact one (`DS.MatType.A`).

        The default is ``False``.
        """
        cdef PetscBool val = PETSC_FALSE
        if comp: val = PETSC_TRUE
        CHKERR( DSSetCompact(self.ds, val) )

    def getCompact(self):
        """
        Gets the compact storage flag.

        Returns
        -------
        comp: boolean
              The flag.
        """
        cdef PetscBool val = PETSC_FALSE
        CHKERR( DSGetCompact(self.ds, &val) )
        return val

    def setExtraRow(self, ext):
        """
        Sets a flag to indicate that the matrix has one extra row.

        Parameters
        ----------
        ext: boolean
             A boolean flag.

        Notes
        -----
        In Krylov methods it is useful that the matrix representing
        the direct solver has one extra row, i.e., has dimension (n+1)
        x n. If this flag is activated, all transformations applied to
        the right of the matrix also affect this additional row. In
        that case, (n+1) must be less or equal than the leading
        dimension.

        The default is ``False``.
        """
        cdef PetscBool val = PETSC_FALSE
        if ext: val = PETSC_TRUE
        CHKERR( DSSetExtraRow(self.ds, val) )

    def getExtraRow(self):
        """
        Gets the extra row flag.

        Returns
        -------
        comp: boolean
              The flag.
        """
        cdef PetscBool val = PETSC_FALSE
        CHKERR( DSGetExtraRow(self.ds, &val) )
        return val

    def setRefined(self, ref):
        """
        Sets a flag to indicate that refined vectors must be computed.

        Parameters
        ----------
        ref: boolean
             A boolean flag.

        Notes
        -----
        Normally the vectors returned in `DS.MatType.X` are eigenvectors
        of the projected matrix. With this flag activated, `vectors()`
        will return the right singular vector of the smallest singular
        value of matrix At-theta*I, where At is the extended (n+1)xn
        matrix and theta is the Ritz value. This is used in the
        refined Ritz approximation.

        The default is ``False``.
        """
        cdef PetscBool val = PETSC_FALSE
        if ref: val = PETSC_TRUE
        CHKERR( DSSetRefined(self.ds, val) )

    def getRefined(self):
        """
        Gets the refined vectors flag.

        Returns
        -------
        comp: boolean
              The flag.
        """
        cdef PetscBool val = PETSC_FALSE
        CHKERR( DSGetRefined(self.ds, &val) )
        return val

    def truncate(self, n):
        """
        Truncates the system represented in the DS object.

        Parameters
        ----------
        n: integer
           The new size.
        """
        cdef PetscInt val = n
        CHKERR( DSTruncate(self.ds, val) )

    def updateExtraRow(self):
        """
        Performs all necessary operations so that the extra
        row gets up-to-date after a call to `solve()`.
        """
        CHKERR( DSUpdateExtraRow(self.ds) )

# -----------------------------------------------------------------------------

del DSType
del DSStateType
del DSMatType

# -----------------------------------------------------------------------------
