# -----------------------------------------------------------------------------

class BVType(object):
    """
    BV type
    """
    MAT        = S_(BVMAT)
    SVEC       = S_(BVSVEC)
    VECS       = S_(BVVECS)
    CONTIGUOUS = S_(BVCONTIGUOUS)

class BVOrthogType(object):
    """
    BV orthogonalization types

    - `CGS`: Classical Gram-Schmidt.
    - `MGS`: Modified Gram-Schmidt.
    """
    CGS = BV_ORTHOG_CGS
    MGS = BV_ORTHOG_MGS

class BVOrthogRefineType(object):
    """
    BV orthogonalization refinement types

    - `IFNEEDED`: Reorthogonalize if a criterion is satisfied.
    - `NEVER`:    Never reorthogonalize.
    - `ALWAYS`:   Always reorthogonalize.
    """
    IFNEEDED = BV_ORTHOG_REFINE_IFNEEDED
    NEVER    = BV_ORTHOG_REFINE_NEVER
    ALWAYS   = BV_ORTHOG_REFINE_ALWAYS

class BVOrthogBlockType(object):
    """
    BV block-orthogonalization types

    - `GS`:   Gram-Schmidt.
    - `CHOL`: Cholesky.
    """
    GS   = BV_ORTHOG_BLOCK_GS
    CHOL = BV_ORTHOG_BLOCK_CHOL

# -----------------------------------------------------------------------------

cdef class BV(Object):

    """
    BV
    """

    Type             = BVType
    OrthogType       = BVOrthogType
    OrthogRefineType = BVOrthogRefineType
    RefineType       = BVOrthogRefineType
    OrthogBlockType  = BVOrthogBlockType
    BlockType        = BVOrthogBlockType

    def __cinit__(self):
        self.obj = <PetscObject*> &self.bv
        self.bv = NULL

    def view(self, Viewer viewer=None):
        """
        Prints the BV data structure.

        Parameters
        ----------
        viewer: Viewer, optional
                Visualization context; if not provided, the standard
                output is used.
        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( BVView(self.bv, vwr) )

    def destroy(self):
        """
        Destroys the BV object.
        """
        CHKERR( BVDestroy(&self.bv) )
        self.bv = NULL
        return self

    def create(self, comm=None):
        """
        Creates the BV object.

        Parameters
        ----------
        comm: Comm, optional
              MPI communicator; if not provided, it defaults to all
              processes.
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcBV newbv = NULL
        CHKERR( BVCreate(ccomm, &newbv) )
        SlepcCLEAR(self.obj); self.bv = newbv
        return self

    def duplicate(self):
        """
        Duplicate the BV object with the same type and dimensions.
        """
        cdef BV bv = type(self)()
        CHKERR( BVDuplicate(self.bv, &bv.bv) )
        return bv

    def copy(self, BV result=None):
        if result is None:
            result = type(self)()
        if result.bv == NULL:
            CHKERR( BVDuplicate(self.bv, &result.bv) )
        CHKERR( BVCopy(self.bv, result.bv) )
        return result

    def setType(self, bv_type):
        """
        Selects the type for the BV object.

        Parameters
        ----------
        bv_type: `BV.Type` enumerate
                  The inner product type to be used.
        """
        cdef SlepcBVType cval = NULL
        bv_type = str2bytes(bv_type, &cval)
        CHKERR( BVSetType(self.bv, cval) )

    def getType(self):
        """
        Gets the BV type of this object.

        Returns
        -------
        type: `BV.Type` enumerate
              The inner product type currently being used.
        """
        cdef SlepcBVType bv_type = NULL
        CHKERR( BVGetType(self.bv, &bv_type) )
        return bytes2str(bv_type)

    def setSizes(self, sizes, m):
        """
        Sets the local and global sizes, and the number of columns.

        Parameters
        ----------
        sizes: int or two-tuple of int
              The global size ``N`` or a two-tuple ``(n, N)``
              with the local and global sizes.
        m: int
              The number of columns.

        Notes
        -----
        Either ``n`` or ``N`` (but not both) can be ``PETSc.DECIDE``
        or ``None`` to have it automatically set.
        """
        cdef PetscInt n=0, N=0
        cdef PetscInt ival = asInt(m)
        BV_Sizes(sizes, &n, &N)
        CHKERR( BVSetSizes(self.bv, n, N, ival) )

    def setSizesFromVec(self, Vec w not None, m):
        """
        Sets the local and global sizes, and the number of columns. Local and
        global sizes are specified indirectly by passing a template vector.

        Parameters
        ----------
        w: Vec
            The template vector.
        m: int
            The number of columns.
        """
        cdef PetscInt ival = asInt(m)
        CHKERR( BVSetSizesFromVec(self.bv, w.vec, ival) )

    def getSizes(self):
        """
        Returns the local and global sizes, and the number of columns.

        Returns
        -------
        sizes: two-tuple of int
                The local and global sizes ``(n, N)``.
        m: int
                The number of columns.
        """
        cdef PetscInt n=0, N=0, m=0
        CHKERR( BVGetSizes(self.bv, &n, &N, &m) )
        return ((toInt(n), toInt(N)), toInt(m))


    def setOptionsPrefix(self, prefix):
        """
        Sets the prefix used for searching for all BV options in the
        database.

        Parameters
        ----------
        prefix: string
                The prefix string to prepend to all BV option
                requests.

        Notes
        -----
        A hyphen (``-``) must NOT be given at the beginning of the
        prefix name.  The first character of all runtime options is
        AUTOMATICALLY the hyphen.
        """
        cdef const_char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( BVSetOptionsPrefix(self.bv, cval) )

    def getOptionsPrefix(self):
        """
        Gets the prefix used for searching for all BV options in the
        database.

        Returns
        -------
        prefix: string
                The prefix string set for this BV object.
        """
        cdef const_char *prefix = NULL
        CHKERR( BVGetOptionsPrefix(self.bv, &prefix) )
        return bytes2str(prefix)

    def setFromOptions(self):
        """
        Sets BV options from the options database.

        Notes
        -----
        To see all options, run your program with the ``-help``
        option.
        """
        CHKERR( BVSetFromOptions(self.bv) )

    #

    def getOrthogonalization(self):
        """
        Gets the orthogonalization settings from the BV object.

        Returns
        -------
        type: `BV.OrthogType` enumerate
              The type of orthogonalization technique.
        refine: `BV.OrthogRefineType` enumerate
              The type of refinement.
        eta:  float
              Parameter for selective refinement (used when the the
              refinement type `BV.OrthogRefineType.IFNEEDED`).
        block: `BV.OrthogBlockType` enumerate
              The type of block orthogonalization .
        """
        cdef SlepcBVOrthogType val1 = BV_ORTHOG_CGS
        cdef SlepcBVOrthogRefineType val2 = BV_ORTHOG_REFINE_IFNEEDED
        cdef SlepcBVOrthogBlockType val3 = BV_ORTHOG_BLOCK_GS
        cdef PetscReal rval = PETSC_DEFAULT
        CHKERR( BVGetOrthogonalization(self.bv, &val1, &val2, &rval, &val3) )
        return (val1, val2, toReal(rval), val3)

    def setOrthogonalization(self, type=None, refine=None, eta=None, block=None):
        """
        Specifies the method used for the orthogonalization of vectors
        (classical or modified Gram-Schmidt with or without refinement),
        and for the block-orthogonalization (simultaneous orthogonalization
        of a set of vectors).

        Parameters
        ----------
        type: `BV.OrthogType` enumerate, optional
              The type of orthogonalization technique.
        refine: `BV.OrthogRefineType` enumerate, optional
              The type of refinement.
        eta:  float, optional
              Parameter for selective refinement.
        block: `BV.OrthogBlockType` enumerate, optional
              The type of block orthogonalization.

        Notes
        -----
        The default settings work well for most problems.

        The parameter `eta` should be a real value between ``0`` and
        ``1`` (or `DEFAULT`).  The value of `eta` is used only when
        the refinement type is `BV.OrthogRefineType.IFNEEDED`.

        When using several processors, `BV.OrthogType.MGS` is likely to
        result in bad scalability.

        If the method set for block orthogonalization is GS, then the
        computation is done column by column with the vector orthogonalization.
        """
        cdef SlepcBVOrthogType val1 = BV_ORTHOG_CGS
        cdef SlepcBVOrthogRefineType val2 = BV_ORTHOG_REFINE_IFNEEDED
        cdef SlepcBVOrthogBlockType val3 = BV_ORTHOG_BLOCK_GS
        cdef PetscReal rval = PETSC_DEFAULT
        CHKERR( BVGetOrthogonalization(self.bv, &val1, &val2, &rval, &val3) )
        if type   is not None: val1 = type
        if refine is not None: val2 = refine
        if block  is not None: val3 = block
        if eta    is not None: rval = asReal(eta)
        CHKERR( BVSetOrthogonalization(self.bv, val1, val2, rval, val3) )

    #

    def getMatrix(self):
        """
        Retrieves the matrix representation of the inner product.

        Returns
        -------
        mat: the matrix of the inner product
        """
        cdef Mat mat = Mat()
        cdef PetscBool indef = PETSC_FALSE
        CHKERR( BVGetMatrix(self.bv, &mat.mat, &indef) )
        PetscINCREF(mat.obj)
        return mat, <bint>indef

    def setMatrix(self, Mat mat, bint indef):
        """
        Sets the bilinear form to be used for inner products.

        Parameters
        ----------
        mat:  Mat, optional
              The matrix of the inner product.
        indef: bool, optional
               Whether the matrix is indefinite
        """
        cdef PetscMat m = NULL if mat is None else mat.mat
        cdef PetscBool tval = PETSC_TRUE if indef else PETSC_FALSE
        CHKERR( BVSetMatrix(self.bv, m, tval) )

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
        CHKERR( BVApplyMatrix(self.bv, x.vec, y.vec) )

    def setActiveColumns(self, int l, int k):
        """
        Specify the columns that will be involved in operations.

        Parameters
        ----------
        l: int
            The leading number of columns.
        k: int
            The active number of columns.
        """
        CHKERR( BVSetActiveColumns(self.bv, l, k) )

    def getActiveColumns(self):
        """
        Returns the current active dimensions.

        Returns
        -------
        l: int
            The leading number of columns.
        k: int
            The active number of columns.
        """
        cdef PetscInt l=0, k=0
        CHKERR( BVGetActiveColumns(self.bv, &l, &k) )
        return (toInt(l), toInt(k))

    def scaleColumn(self, int j, alpha):
        """
        Scale column j by alpha

        Parameters
        ----------
        j: int
            column number to be scaled.
        alpha: float
            scaling factor.
        """
        cdef PetscScalar sval = asScalar(alpha)
        CHKERR( BVScaleColumn(self.bv, j, sval) )

    def scale(self, alpha):
        """
        Multiply the entries by a scalar value.

        Parameters
        ----------
        alpha: float
            scaling factor.

        Notes
        -----
        All active columns (except the leading ones) are scaled.
        """
        cdef PetscScalar sval = asScalar(alpha)
        CHKERR( BVScale(self.bv, sval) )

    def insertVec(self, int j, Vec w not None):
        """
        Insert a vector into the specified column.

        Parameters
        ----------
        j: int
            The column to be overwritten.
        w: Vec
            The vector to be copied.
        """
        CHKERR( BVInsertVec(self.bv, j, w.vec) )

    def insertVecs(self, int s, W not None, bint orth):
        """
        Insert a set of vectors into specified columns.

        Parameters
        ----------
        s: int
            The first column to be overwritten.
        W: Vec or sequence of Vec.
            Set of vectors to be copied.
        orth:
            Flag indicating if the vectors must be orthogonalized.

        Returns
        -------
        m: int
            Number of linearly independent vectors.

        Notes
        -----
        Copies the contents of vectors W into self(:,s:s+n), where n is the
        length of W. If orthogonalization flag is set then the vectors are
        copied one by one then orthogonalized against the previous one.  If any
        are linearly dependent then it is discared and the value of m is
        decreased.
        """
        if isinstance(W, Vec): W = [W]
        cdef PetscVec *ws = NULL
        cdef Py_ssize_t i = 0, ns = len(W)
        cdef tmp = allocate(<size_t>ns*sizeof(Vec),<void**>&ws)
        for i in range(ns): ws[i] = (<Vec?>W[i]).vec
        cdef PetscInt m = <PetscInt>ns
        cdef PetscBool tval = PETSC_TRUE if orth else PETSC_FALSE
        CHKERR( BVInsertVecs(self.bv, <PetscInt>s, &m, ws, tval) )
        return toInt(m)

    def dotVec(self, Vec v not None):
        """
        Computes multiple dot products of a vector against all the column
        vectors of a BV.

        Parameters
        ----------
        v: Vec
            A vector.

        Returns
        -------
        m: Vec
            A vector with the results.

        This is analogue to VecMDot(), but using BV to represent a collection
        of vectors. The result is m = X^H*y, so m_i is equal to x_j^H y. Note
        that here X is transposed as opposed to BVDot().

        If a non-standard inner product has been specified with BVSetMatrix(),
        then the result is m = X^H*B*y.
        """
        l, k = self.getActiveColumns()
        cdef PetscScalar* mval = NULL
        cdef tmp = allocate(<size_t>(k - l)*sizeof(PetscScalar),<void**>&mval)

        CHKERR( BVDotVec(self.bv, v.vec, mval) )

        v = Vec().create(COMM_SELF)
        v.setType('seq')
        v.setSizes((DECIDE,k-l))
        v.setArray([mval[i] for i in range(0, k - l)])
        v.ghostUpdate()

        return v

    def getColumn(self, int j):
        """
        Returns a Vec object that contains the entries of the requested column
        of the basis vectors object.

        Parameters
        ----------
        j: int
            The index of the requested column.

        Returns
        -------
        v: Vec
            The vector containing the jth column.

        Notes
        -----
        Modifying the returned Vec will change the BV entries as well.
        """
        cdef Vec v = Vec()
        CHKERR( BVGetColumn(self.bv, j, &v.vec) )
        return v

    def restoreColumn(self, int j, Vec v not None):
        """
        Restore a column obtained with BVGetColumn().

        Parameters
        ----------
        j: int
            The index of the requested column.

        v: Vec
            The vector obtained with BVGetColumn().

        Notes
        -----
        The arguments must match the corresponding call to BVGetColumn().
        """
        CHKERR( BVRestoreColumn(self.bv, j, &v.vec) )

    def dot(self, BV Y not None):
        """
        Computes the 'block-dot' product of two basis vectors objects.
            M = Y^H*X (m_ij = y_i^H x_j) or M = Y^H*B*X

        Parameters
        ----------
        Y: BV
            Left basis vectors, can be the same as self, giving M = X^H X.

        Returns
        -------
        M: Mat
            The resulting matrix.

        Notes
        -----
        This is the generalization of VecDot() for a collection of vectors, M =
        Y^H*X. The result is a matrix M whose entry m_ij is equal to y_i^H x_j
        (where y_i^H denotes the conjugate transpose of y_i).

        X and Y can be the same object.

        If a non-standard inner product has been specified with setMatrix(),
        then the result is M = Y^H*B*X. In this case, both X and Y must have
        the same associated matrix.

        Only rows (resp. columns) of M starting from ly (resp. lx) are
        computed, where ly (resp. lx) is the number of leading columns of Y
        (resp. X).
        """
        cdef BV X = self
        cdef PetscInt ky=0, kx=0
        CHKERR( BVGetActiveColumns(Y.bv, NULL, &ky) )
        CHKERR( BVGetActiveColumns(X.bv, NULL, &kx) )
        cdef Mat M = Mat().createDense((ky, kx), comm=COMM_SELF).setUp()
        CHKERR( BVDot(X.bv, Y.bv, M.mat) )
        return M

    def matProject(self, Mat A, BV Y not None):
        """
        Computes the projection of a matrix onto a subspace.

        M = Y^H A X

        Parameters
        ----------
        A: Mat or None
            Matrix to be projected.

        Y: BV
            Left basis vectors, can be the same as self, giving M = X^H A X.

        Returns
        -------
        M: Mat
            Projection of the matrix A onto the subspace.
        """
        cdef BV X = self
        cdef PetscInt ky=0, kx=0
        CHKERR( BVGetActiveColumns(Y.bv, NULL, &ky) )
        CHKERR( BVGetActiveColumns(X.bv, NULL, &kx) )
        cdef Mat M = Mat().createDense((ky, kx), comm=COMM_SELF).setUp()
        cdef PetscMat Amat = NULL if A is None else A.mat
        CHKERR( BVMatProject(X.bv, Amat, Y.bv, M.mat) )
        return M

    def matMult(self, Mat A not None, BV Y=None):
        """
        Computes the matrix-vector product for each column, Y = A*V.

        Parameters
        ----------
        A: Mat
            The matrix.

        Returns
        -------
        Y: BV
            The result.

        Notes
        -----
        Only active columns (excluding the leading ones) are processed.

        It is possible to choose whether the computation is done column by column
        or using dense matrices using the options database keys:

            -bv_matmult_vecs
            -bv_matmult_mat

        The default is bv_matmult_mat.
        """
        cdef MPI_Comm comm = PetscObjectComm(<PetscObject>self.bv)
        cdef SlepcBVType bv_type = NULL
        cdef PetscInt n=0, N=0, m=0
        cdef SlepcBVOrthogType val1 = BV_ORTHOG_CGS
        cdef SlepcBVOrthogRefineType val2 = BV_ORTHOG_REFINE_IFNEEDED
        cdef SlepcBVOrthogBlockType val3 = BV_ORTHOG_BLOCK_GS
        cdef PetscReal rval = PETSC_DEFAULT
        if Y is None: Y = BV()
        if Y.bv == NULL:
            CHKERR( BVGetType(self.bv, &bv_type) )
            CHKERR( MatGetLocalSize(A.mat, &n, NULL) )
            CHKERR( MatGetSize(A.mat, &N, NULL) )
            CHKERR( BVGetSizes(self.bv, NULL, NULL, &m) )
            CHKERR( BVGetOrthogonalization(self.bv, &val1, &val2, &rval, &val3) )
        if Y.bv == NULL:
            CHKERR( BVCreate(comm, &Y.bv) )
            CHKERR( BVSetType(Y.bv, bv_type) )
            CHKERR( BVSetSizes(Y.bv, n, N, m) )
            CHKERR( BVSetOrthogonalization(Y.bv, val1, val2, rval, val3) )
        CHKERR( BVMatMult(self.bv, A.mat, Y.bv) )
        return Y

    def matMultHermitianTranspose(self, Mat A not None, BV Y=None):
        """
        Computes the matrix-vector product with the conjugate transpose of a
        matrix for each column, Y=A^H*V.

        Parameters
        ----------
        A: Mat
            The matrix.

        Returns
        -------
        Y: BV
            The result.

        Notes
        -----
        Only active columns (excluding the leading ones) are processed.

        As opoosed to matMult(), this operation is always done by column by
        column, with a sequence of calls to MatMultHermitianTranspose().
        """
        cdef MPI_Comm comm = PetscObjectComm(<PetscObject>self.bv)
        cdef SlepcBVType bv_type = NULL
        cdef PetscInt n=0, N=0, m=0
        cdef SlepcBVOrthogType val1 = BV_ORTHOG_CGS
        cdef SlepcBVOrthogRefineType val2 = BV_ORTHOG_REFINE_IFNEEDED
        cdef SlepcBVOrthogBlockType val3 = BV_ORTHOG_BLOCK_GS
        cdef PetscReal rval = PETSC_DEFAULT
        if Y is None: Y = BV()
        if Y.bv == NULL:
            CHKERR( BVGetType(self.bv, &bv_type) )
            CHKERR( MatGetLocalSize(A.mat, &n, NULL) )
            CHKERR( MatGetSize(A.mat, &N, NULL) )
            CHKERR( BVGetSizes(self.bv, NULL, NULL, &m) )
            CHKERR( BVGetOrthogonalization(self.bv, &val1, &val2, &rval, &val3) )
        if Y.bv == NULL:
            CHKERR( BVCreate(comm, &Y.bv) )
            CHKERR( BVSetType(Y.bv, bv_type) )
            CHKERR( BVSetSizes(Y.bv, n, N, m) )
            CHKERR( BVSetOrthogonalization(Y.bv, val1, val2, rval, val3) )
        CHKERR( BVMatMultHermitianTranspose(self.bv, A.mat, Y.bv) )
        return Y

    def multVec(self, alpha, beta, Vec y not None, q):
        """
        Computes y = beta*y + alpha*X*q.

        Parameter
        ---------
        alpha: scalar
        beta: scalar
        q: scalar or sequence of scalars

        Return
        ------
        y: Vec
            The result.
        """
        cdef PetscScalar sval1 = asScalar(alpha)
        cdef PetscScalar sval2 = asScalar(beta)
        cdef PetscInt nq = 0
        cdef PetscScalar* qval = NULL
        cdef tmp = iarray_s(q, &nq, &qval)
        cdef PetscInt l=0, k=0
        CHKERR( BVGetActiveColumns(self.bv, &l, &k) )
        assert nq == k-l
        CHKERR( BVMultVec(self.bv, sval1, sval2, y.vec, qval) )

    def normColumn(self, int j, norm_type=None):
        """
        Computes the matrix norm of the BV.

        Parameters
        ----------
        j: int
            Index of column.
        norm_type: PETSc.NormType (int)
            The norm type.

        Returns
        -------
        norm: float

        Notes
        -----
        The norm of V[j] is computed (NORM_1, NORM_2, or NORM_INFINITY).

        If a non-standard inner product has been specified with BVSetMatrix(),
        then the returned value is ``sqrt(V[j]'* B*V[j])``, where B is the inner
        product matrix (argument 'type' is ignored).
        """
        cdef PetscNormType ntype = PETSC_NORM_2
        if norm_type is not None: ntype = norm_type
        cdef PetscReal norm = 0
        CHKERR( BVNormColumn(self.bv, j, ntype, &norm) )
        return toReal(norm)

    def norm(self, norm_type=None):
        """
        Computes the matrix norm of the BV.

        Parameters
        ----------
        norm_type: PETSC.NormType enumerate
            The norm type.

        Returns
        -------
        norm: float

        Notes
        -----
        All active columns (except the leading ones) are considered as a
        matrix. The allowed norms are NORM_1, NORM_FROBENIUS, and
        NORM_INFINITY.

        This operation fails if a non-standard inner product has been specified
        with BVSetMatrix().
        """
        cdef PetscNormType ntype = PETSC_NORM_FROBENIUS
        if norm_type is not None: ntype = norm_type
        cdef PetscReal norm = 0
        CHKERR( BVNorm(self.bv, ntype, &norm) )
        return toReal(norm)

    def setRandom(self):
        """
        Set the active columns of BV to random numbers.

        Notes
        -----
        All active columns (except the leading ones) are modified.
        """
        CHKERR( BVSetRandom(self.bv) )

    def orthogonalizeVec(self, Vec v not None):
        """
        Orthogonalize a vector with respect to a set of vectors.

        Parameters
        ----------
        v:  Vec
            Vector to be orthogonalized, modified on return.

        Returns
        -------
        norm: float
            The norm of the resulting vector.
        lindep: boolean
            Flag indicating that refinement did not improve the
            quality of orthogonalization.

        Notes
        -----
        This function applies an orthogonal projector to project
        vector ``v`` onto the orthogonal complement of the span of the
        columns of the BV.

        This routine does not normalize the resulting vector.
        """
        cdef PetscReal norm = 0
        cdef PetscBool ldep = PETSC_FALSE
        CHKERR( BVOrthogonalizeVec(self.bv, v.vec, NULL, &norm, &ldep) )
        return (toReal(norm), <bint>ldep)

    def orthogonalize(self, Mat R=None, **kargs):
        """
        Orthogonalize all columns (except leading ones),
        that is, compute the QR decomposition.

        Parameters
        ----------
        R: Mat or None
            A sequential dense matrix.

        Notes
        -----
        The output satisfies ``V0 = V*R`` (where V0 represent the input V) and ``V'*V = I``.
        """
        if kargs: self.setOrthogonalization(**kargs)
        cdef PetscMat Rmat = NULL if R is None else R.mat
        CHKERR( BVOrthogonalize(self.bv, Rmat) )

# -----------------------------------------------------------------------------

del BVType
del BVOrthogType
del BVOrthogRefineType
del BVOrthogBlockType

# -----------------------------------------------------------------------------
