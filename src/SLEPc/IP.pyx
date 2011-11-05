# -----------------------------------------------------------------------------

class IPOrthogType(object):
    """
    IP orthogonalization types

    - `CGS`: Classical Gram-Schmidt.
    - `MGS`: Modified Gram-Schmidt.
    """
    CGS = IP_ORTHOG_CGS
    MGS = IP_ORTHOG_MGS

class IPOrthogRefineType(object):
    """
    IP orthogonalization refinement types

    - `NEVER`:    Never reorthogonalize.
    - `IFNEEDED`: Reorthogonalize if a criterion is satisfied.
    - `ALWAYS`:   Always reorthogonalize.
    """
    NEVER    = IP_ORTHOG_REFINE_NEVER
    IFNEEDED = IP_ORTHOG_REFINE_IFNEEDED
    ALWAYS   = IP_ORTHOG_REFINE_ALWAYS

# -----------------------------------------------------------------------------

cdef class IP(Object):

    """
    IP
    """

    OrthogType       = IPOrthogType
    OrthogRefineType = IPOrthogRefineType
    RefineType       = IPOrthogRefineType

    def __cinit__(self):
        self.obj = <PetscObject*> &self.ip
        self.ip = NULL

    def view(self, Viewer viewer=None):
        """
        Prints the IP data structure.

        Parameters
        ----------
        viewer: Viewer, optional
                Visualization context; if not provided, the standard
                output is used.
        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( IPView(self.ip, vwr) )

    def destroy(self):
        """
        Destroys the IP object.
        """
        CHKERR( IPDestroy(&self.ip) )
        self.ip = NULL
        return self

    def reset(self):
        """
        Resets the IP object.
        """
        CHKERR( IPReset(self.ip) )

    def create(self, comm=None):
        """
        Creates the IP object.

        Parameters
        ----------
        comm: Comm, optional
              MPI communicator; if not provided, it defaults to all
              processes.
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcIP newip = NULL
        CHKERR( IPCreate(ccomm, &newip) )
        SlepcCLEAR(self.obj); self.ip = newip
        return self

    def setOptionsPrefix(self, prefix):
        """
        Sets the prefix used for searching for all IP options in the
        database.

        Parameters
        ----------
        prefix: string
                The prefix string to prepend to all IP option
                requests.

        Notes
        -----
        A hyphen (``-``) must NOT be given at the beginning of the
        prefix name.  The first character of all runtime options is
        AUTOMATICALLY the hyphen.
        """
        cdef const_char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( IPSetOptionsPrefix(self.ip, cval) )

    def getOptionsPrefix(self):
        """
        Gets the prefix used for searching for all IP options in the
        database.

        Returns
        -------
        prefix: string
                The prefix string set for this IP object.
        """
        cdef const_char *prefix = NULL
        CHKERR( IPGetOptionsPrefix(self.ip, &prefix) )
        return bytes2str(prefix)

    def setFromOptions(self):
        """
        Sets IP options from the options database.

        Notes
        -----
        To see all options, run your program with the ``-help``
        option.
        """
        CHKERR( IPSetFromOptions(self.ip) )

    #

    def getOrthogonalization(self):
        """
        Gets the orthogonalization settings from the IP object.

        Returns
        -------
        type: `IP.OrthogType` enumerate
              The type of orthogonalization technique.
        refine: `IP.OrthogRefineType` enumerate
              The type of refinement.
        eta:  float
              Parameter for selective refinement (used when the the
              refinement type `IP.OrthogRefineType.IFNEEDED`).
        """
        cdef SlepcIPOrthogType val1 = IP_ORTHOG_CGS
        cdef SlepcIPOrthogRefineType val2 = IP_ORTHOG_REFINE_IFNEEDED
        cdef PetscReal rval = PETSC_DEFAULT
        CHKERR( IPGetOrthogonalization(self.ip, &val1, &val2, &rval) )
        return (val1, val2, toReal(rval))

    def setOrthogonalization(self, type=None, refine=None, eta=None):
        """
        Specifies the type of orthogonalization technique to be used
        (classical or modified Gram-Schmidt with or without
        refinement).


        Parameters
        ----------
        type: `IP.OrthogType` enumerate, optional
              The type of orthogonalization technique.
        refine: `IP.OrthogRefineType` enumerate, optional
              The type of refinement.
        eta:  float, optional
              Parameter for selective refinement.

        Notes
        -----
        The default settings work well for most problems.

        The parameter `eta` should be a real value between ``0`` and
        ``1`` (or `DEFAULT`).  The value of `eta` is used only when
        the refinement type is `IP.OrthogRefineType.IFNEEDED`.

        When using several processors, `IP.OrthogType.MGS` is likely to
        result in bad scalability.
        """
        cdef SlepcIPOrthogType val1 = IP_ORTHOG_CGS
        cdef SlepcIPOrthogRefineType val2 = IP_ORTHOG_REFINE_IFNEEDED
        cdef PetscReal rval = PETSC_DEFAULT
        if type   is not None: val1= type
        if refine is not None: val2= refine
        if eta    is not None: rval = asReal(eta)
        CHKERR( IPSetOrthogonalization(self.ip, val1, val2, rval) )

    #

    def getMatrix(self):
        """
        Retrieves the matrix representation of the inner produc

        Returns
        -------
        mat: the matrix of the inner product
        """
        cdef Mat mat = Mat()
        CHKERR( IPGetMatrix(self.ip, &mat.mat) )
        PetscINCREF(mat.obj)
        return mat

    def setMatrix(self, Mat mat):
        """
        Sets the bilinear form to be used for inner products.

        Parameters
        ----------
        mat:  Mat, optional
              The matrix of the inner product.
        """
        cdef PetscMat m = NULL
        if mat is not None: m = mat.mat
        CHKERR( IPSetMatrix(self.ip, m) )

    def applyMatrix(self, Vec x not None, Vec y not None):
        """
        Multiplies a vector with the matrix associated to the bilinear
        form.

        Parameters
        ----------
        x: Vec
           The input vector.
        y: Vec
           The result vector.

        Notes
        -----
        If the bilinear form has no associated matrix this function
        copies the vector.
        """
        CHKERR( IPApplyMatrix(self.ip, x.vec, y.vec) )

    #

    def norm(self, Vec x not None):
        """
        Computes the norm of a vector as the square root of the inner
        product ``(x,x)`` as defined by `innerProduct()`.

        Parameters
        ----------
        x: Vec
           The input vector.

        Returns
        -------
        norm: float
              The computed norm.

        Notes
        -----
        This function will usually compute the 2-norm of a vector,
        ``||x||_2``. But this behaviour may be different if using a
        non-standard inner product changed via
        `setMatrix()`. For example, if using the B-inner product
        for positive definite ``B`, ``(x,y)_B=y^H Bx``, then the
        computed norm is ``||x||_B = sqrt( x^H Bx )``.
        """
        cdef PetscReal rval = 0
        CHKERR( IPNorm(self.ip, x.vec, &rval) )
        return toReal(rval)

    def innerProduct(self, Vec x not None, Vec y not None):
        """
        Computes the inner product of two vectors.

        Parameters
        ----------
        x: Vec
           The first input vector.
        y: Vec
           The second input vector.

        Returns
        -------
        p: float
           The result of the inner product.

        Notes
        -----
        This function will usually compute the standard dot product,
        ``(x,y)=y^H x``.  However this behaviour may be different if
        changed via `setMatrix()`.  This allows use of other
        inner products such as the indefinite product ``y^T x`` for
        complex symmetric problems or the B-inner product for positive
        definite ``B``, ``(x,y)_B=y^H Bx``.
        """
        cdef PetscScalar sval = 0
        CHKERR( IPInnerProduct(self.ip, x.vec, y.vec, &sval) )
        return toScalar(sval)

    def orthogonalize(self, VS, Vec v not None):
        """
        Orthogonalize a vector with respect to a set of vectors.

        Parameters
        ----------
        VS: list of Vec
            Set of orthonormal vectors.
        v:  Vec
            Vector to be orthogonalized, modified on return.

        Returns
        -------
        H:  list of float
            Coefficients computed during orthogonalization.
        norm: float
            The norm of the resulting vector.
        lindep: boolean
            Flag indicating that refinement did not improve the
            quality of orthogonalization.

        Notes
        -----
        This function applies an orthogonal projector to project
        vector ``v`` onto the orthogonal complement of the span of the
        columns of ``VS``.

        On exit, ``v0 = [V v]*H``, where ``v0`` is the original vector
        ``v``.

        This routine does not normalize the resulting vector.
        """
        cdef PetscBool* which = NULL
        cdef PetscVec* V = NULL
        cdef PetscScalar* H = NULL, h = 0
        cdef PetscReal rval = 0
        cdef PetscBool tval = PETSC_FALSE
        cdef object tmp1 = None, tmp2 = None
        cdef Py_ssize_t i = 0, n = 0
        if isinstance(VS, Vec):
            n = 1
            V = &((<Vec>VS).vec)
            H = &h
        else:
            n = len(VS)
            tmp1 = allocate(<size_t>n*sizeof(PetscVec),<void**>&V)
            tmp2 = allocate(<size_t>n*sizeof(PetscScalar),<void**>&H)
            for i in range(n):
                V[i] = (<Vec?>VS[i]).vec
                H[i] = 0
        CHKERR( IPOrthogonalize(self.ip, 0, NULL,
                                <PetscInt>n, which, V,
                                v.vec, H, &rval, &tval) )
        cdef object coefs = None
        if isinstance(VS, Vec):
            coefs = toScalar(H[0])
        else:
            coefs = [toScalar(H[i]) for i in range(n)]
        return (coefs, toReal(rval), <bint>tval)

# -----------------------------------------------------------------------------

del IPOrthogType
del IPOrthogRefineType

# -----------------------------------------------------------------------------
