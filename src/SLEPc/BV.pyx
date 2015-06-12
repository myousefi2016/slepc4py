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
        if type   is not None: val1= type
        if refine is not None: val2= refine
        if block  is not None: val3= block
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
        cdef PetscMat m = NULL
        if mat is not None: m = mat.mat
        cdef PetscBool cflag = PETSC_FALSE
        if indef: cflag = PETSC_TRUE
        CHKERR( BVSetMatrix(self.bv, m, cflag) )

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
        cdef PetscReal rval = 0
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( BVOrthogonalizeVec(self.bv, v.vec, NULL, &rval, &tval) )
        return (toReal(rval), <bint>tval)

# -----------------------------------------------------------------------------

del BVType
del BVOrthogType
del BVOrthogRefineType

# -----------------------------------------------------------------------------
