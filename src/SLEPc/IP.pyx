# --------------------------------------------------------------------

class IPOrthoType(object):
    """
    IP orthogonalization types

    - `CGS`: Classical Gram-Schmidt.
    - `MGS`: Modified Gram-Schmidt.
    """
    CGS = IP_CGS_ORTH
    MGS = IP_MGS_ORTH

class IPRefineType(object):
    """
    IP orthogonalization refinement types

    - `NEVER`:    Never reorthogonalize.
    - `IFNEEDED`: Reorthogonalize if a criterion is satisfied.
    - `ALWAYS`:   Always reorthogonalize.
    """
    NEVER    = IP_ORTH_REFINE_NEVER
    IFNEEDED = IP_ORTH_REFINE_IFNEEDED
    ALWAYS   = IP_ORTH_REFINE_ALWAYS

class IPBilinearForm(object):
    """
    IP bilinear form types

    - `HERMITIAN`:
    - `SYMMETRIC`:
    """
    HERMITIAN = IPINNER_HERMITIAN
    SYMMETRIC = IPINNER_SYMMETRIC

# --------------------------------------------------------------------

cdef class IP(Object):

    """
    IP
    """

    OrthoType    = IPOrthoType
    RefineType   = IPRefineType
    BilinearForm = IPBilinearForm

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
        CHKERR( IPDestroy(self.ip) )
        self.ip = NULL
        return self

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
        self.dec_ref(); self.ip = newip
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
        CHKERR( IPSetOptionsPrefix(self.ip, str2cp(prefix)) )

    def getOptionsPrefix(self):
        """
        Gets the prefix used for searching for all IP options in the
        database.

        Returns
        -------
        prefix: string
                The prefix string set for this IP object.
        """
        cdef const_char_p prefix = NULL
        CHKERR( IPGetOptionsPrefix(self.ip, &prefix) )
        return cp2str(prefix)

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
        type: IP.OrthoType enumerate
              The type of orthogonalization technique.
        refine: IP.RefineType enumerate
              The type of refinement.
        eta:  float
              Parameter for selective refinement (used when the the
              refinement type `IP.RefineType.IFNEEDED`).
        """
        cdef SlepcIPOrthogonalizationType val1 = IP_CGS_ORTH
        cdef SlepcIPOrthogonalizationRefinementType val2 = IP_ORTH_REFINE_IFNEEDED
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
        type: IP.OrthoType enumerate, optional
              The type of orthogonalization technique.
        refine: IP.RefineType enumerate, optional
              The type of refinement.
        eta:  float, optional
              Parameter for selective refinement.

        Notes
        -----
        The default settings work well for most problems.

        The parameter `eta` should be a real value between ``0`` and
        ``1`` (or `DEFAULT`).  The value of `eta` is used only when
        the refinement type is `IP.RefineType.IFNEEDED`.

        When using several processors, `IP.OrthoType.MGS` is likely to
        result in bad scalability.
        """
        cdef SlepcIPOrthogonalizationType val1 = IP_CGS_ORTH
        cdef SlepcIPOrthogonalizationRefinementType val2 = IP_ORTH_REFINE_IFNEEDED
        cdef PetscReal rval = PETSC_DEFAULT
        if type   is not None: val1= type
        if refine is not None: val2= refine
        if eta    is not None: rval = asReal(eta)
        CHKERR( IPSetOrthogonalization(self.ip, val1, val2, rval) )

    #

    def getBilinearForm(self):
        """
        Gets the bilinear form to be used for inner products.

        Returns
        -------
        form: IP.BilinearForm enumerate
              The type of bilinear form.
        """
        cdef Mat mat = Mat()
        cdef PetscMat m = NULL
        cdef SlepcIPBilinearForm val = IPINNER_HERMITIAN
        CHKERR( IPGetBilinearForm(self.ip, &m, &val) )
        mat.mat = m; mat.inc_ref()
        return (mat, val)

    def setBilinearForm(self, Mat mat=None, form=None):
        """
        Sets the bilinear form to be used for inner products.

        Parameters
        ----------
        mat:  Mat, optional
              The matrix of the bilinear form.
        form: IP.BilinearForm enumerate, optional
              The type of bilinear form.
        """
        cdef PetscMat m = NULL
        cdef SlepcIPBilinearForm val = IPINNER_HERMITIAN
        if mat  is not None: m = mat.mat
        if form is not None: val = form
        CHKERR( IPSetBilinearForm(self.ip, m, val) )

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
        `setBilinearForm()`. For example, if using the B-inner product
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
        changed via `setBilinearForm()`.  This allows use of other
        inner products such as the indefinite product ``y^T x`` for
        complex symmetric problems or the B-inner product for positive
        definite ``B``, ``(x,y)_B=y^H Bx``.
        """
        cdef PetscScalar sval = 0
        CHKERR( IPInnerProduct(self.ip, x.vec, y.vec, &sval) )
        return toScalar(sval)

    def orthogonalize(self, VS, Vec v not None, Vec work=None):
        """
        Orthogonalize a vector with respect to a set of vectors.

        Parameters
        ----------
        VS: list of Vec
            Set of orthonormal vectors.
        v:  Vec
            Vector to be orthogonalized, modified on return.
        work: Vec, optional
            Workspace.

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
        cdef PetscInt i = 0, n = 0
        cdef PetscTruth* which = NULL
        cdef PetscVec* V = NULL
        cdef PetscScalar* H = NULL, h = 0
        cdef PetscReal rval = 0
        cdef PetscTruth tval = PETSC_FALSE
        cdef PetscVec w = NULL
        cdef PetscScalar* sw = NULL
        cdef object tmp1 = None, tmp2 = None
        if isinstance(VS, Vec):
            n = 1
            V = &((<Vec>VS).vec)
            H = &h
        else:
            n = len(VS)
            tmp1 = allocate(n*sizeof(Vec),<void**>&V)
            tmp2 = allocate(n*sizeof(PetscScalar),<void**>&H)
            for i in range(n):
                V[i] = (<Vec?>VS[i]).vec
                H[i] = 0
        if work is not None: w = work.vec
        CHKERR( IPOrthogonalize(self.ip,
                                n, which, V, v.vec,
                                H, &rval, &tval,
                                w, sw) )
        cdef object coefs = None
        if isinstance(VS, Vec):
            coefs = toScalar(H[0])
        else:
            coefs = [toScalar(H[i]) for i in range(n)]
        return (coefs, toReal(rval), <bint>tval)

# --------------------------------------------------------------------

del IPOrthoType
del IPRefineType
del IPBilinearForm

# --------------------------------------------------------------------
