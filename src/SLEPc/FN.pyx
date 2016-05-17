# -----------------------------------------------------------------------------

class FNType(object):
    """
    FN type
    """
    COMBINE  = S_(FNCOMBINE)
    RATIONAL = S_(FNRATIONAL)
    EXP      = S_(FNEXP)
    LOG      = S_(FNLOG)
    PHI      = S_(FNPHI)
    SQRT     = S_(FNSQRT)
    INVSQRT  = S_(FNINVSQRT)

class FNCombineType(object):
    """
    FN type of combination of child functions

    - `ADD`:       Addition         f(x) = f1(x)+f2(x)
    - `MULTIPLY`:  Multiplication   f(x) = f1(x)*f2(x)
    - `DIVIDE`:    Division         f(x) = f1(x)/f2(x)
    - `COMPOSE`:   Composition      f(x) = f2(f1(x))
    """
    ADD      = FN_COMBINE_ADD
    MULTIPLY = FN_COMBINE_MULTIPLY
    DIVIDE   = FN_COMBINE_DIVIDE
    COMPOSE  = FN_COMBINE_COMPOSE

# -----------------------------------------------------------------------------

cdef class FN(Object):

    """
    FN
    """

    Type        = FNType
    CombineType = FNCombineType

    def __cinit__(self):
        self.obj = <PetscObject*> &self.fn
        self.fn = NULL

    def view(self, Viewer viewer=None):
        """
        Prints the FN data structure.

        Parameters
        ----------
        viewer: Viewer, optional
                Visualization context; if not provided, the standard
                output is used.
        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( FNView(self.fn, vwr) )

    def destroy(self):
        """
        Destroys the FN object.
        """
        CHKERR( FNDestroy(&self.fn) )
        self.fn = NULL
        return self

    def create(self, comm=None):
        """
        Creates the FN object.

        Parameters
        ----------
        comm: Comm, optional
              MPI communicator; if not provided, it defaults to all
              processes.
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcFN newfn = NULL
        CHKERR( FNCreate(ccomm, &newfn) )
        SlepcCLEAR(self.obj); self.fn = newfn
        return self

    def setType(self, fn_type):
        """
        Selects the type for the FN object.

        Parameters
        ----------
        fn_type: `FN.Type` enumerate
                  The inner product type to be used.
        """
        cdef SlepcFNType cval = NULL
        fn_type = str2bytes(fn_type, &cval)
        CHKERR( FNSetType(self.fn, cval) )

    def getType(self):
        """
        Gets the FN type of this object.

        Returns
        -------
        type: `FN.Type` enumerate
              The inner product type currently being used.
        """
        cdef SlepcFNType fn_type = NULL
        CHKERR( FNGetType(self.fn, &fn_type) )
        return bytes2str(fn_type)

    def setOptionsPrefix(self, prefix):
        """
        Sets the prefix used for searching for all FN options in the
        database.

        Parameters
        ----------
        prefix: string
                The prefix string to prepend to all FN option
                requests.

        Notes
        -----
        A hyphen (``-``) must NOT be given at the beginning of the
        prefix name.  The first character of all runtime options is
        AUTOMATICALLY the hyphen.
        """
        cdef const_char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( FNSetOptionsPrefix(self.fn, cval) )

    def getOptionsPrefix(self):
        """
        Gets the prefix used for searching for all FN options in the
        database.

        Returns
        -------
        prefix: string
                The prefix string set for this FN object.
        """
        cdef const_char *prefix = NULL
        CHKERR( FNGetOptionsPrefix(self.fn, &prefix) )
        return bytes2str(prefix)

    def setFromOptions(self):
        """
        Sets FN options from the options database.

        Notes
        -----
        To see all options, run your program with the ``-help``
        option.
        """
        CHKERR( FNSetFromOptions(self.fn) )

    #

    def evaluateFunction(self, x):
        """
        Computes the value of the function f(x) for a given x.

        Parameters
        ----------
        x: scalar
            Value where the function must be evaluated.

        Returns
        -------
        y: scalar
            The result of f(x).
        """
        cdef PetscScalar sval = 0
        CHKERR( FNEvaluateFunction(self.fn, x, &sval) )
        return toScalar(sval)

    def evaluateDerivative(self, x):
        """
        Computes the value of the derivative f'(x) for a given x.

        Parameters
        ----------
        x: scalar
            Value where the derivative must be evaluated.

        Returns
        -------
        y: scalar
            The result of f'(x).
        """
        cdef PetscScalar sval = 0
        CHKERR( FNEvaluateDerivative(self.fn, x, &sval) )
        return toScalar(sval)

    def setScale(self, alpha=None, beta=None):
        """
        Sets the scaling parameters that define the matematical function.

        Parameters
        ----------
        alpha: scalar (possibly complex)
               inner scaling (argument).
        beta: scalar (possibly complex)
               outer scaling (result).
        """
        cdef PetscScalar aval = 1.0
        cdef PetscScalar bval = 1.0
        if alpha is not None: aval = asScalar(alpha)
        if beta  is not None: bval = asScalar(beta)
        CHKERR( FNSetScale(self.fn, aval, bval) )

    def getScale(self):
        """
        Gets the scaling parameters that define the matematical function.

        Returns
        -------
        alpha: scalar (possibly complex)
               inner scaling (argument).
        beta: scalar (possibly complex)
               outer scaling (result).
        """
        cdef PetscScalar aval = 0, bval = 0
        CHKERR( FNGetScale(self.fn, &aval, &bval) )
        return (toScalar(aval), toScalar(bval))

    #

    def setRationalNumerator(self, alpha not None):
        """
        Sets the coefficients of the numerator of the rational function.

        Parameters
        ----------
        alpha: array of scalars
            Coefficients.
        """
        cdef PetscInt na = 0
        cdef PetscScalar *a = NULL
        cdef object tmp1 = iarray_s(alpha, &na, &a)
        CHKERR( FNRationalSetNumerator(self.fn, na, a) )

    def setRationalDenominator(self, alpha not None):
        """
        Sets the coefficients of the denominator of the rational function.

        Parameters
        ----------
        alpha: array of scalars
            Coefficients.
        """
        cdef PetscInt na = 0
        cdef PetscScalar *a = NULL
        cdef object tmp1 = iarray_s(alpha, &na, &a)
        CHKERR( FNRationalSetDenominator(self.fn, na, a) )

# -----------------------------------------------------------------------------

del FNType
del FNCombineType

# -----------------------------------------------------------------------------
