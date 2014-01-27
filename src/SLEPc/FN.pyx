# -----------------------------------------------------------------------------

class FNType(object):
    """
    FN type
    """
    RATIONAL = S_(FNRATIONAL)
    EXP      = S_(FNEXP)
    LOG      = S_(FNLOG)
    PHI      = S_(FNPHI)

# -----------------------------------------------------------------------------

cdef class FN(Object):

    """
    FN
    """

    Type             = FNType

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

    def setParameters(self, alpha=None, beta=None):
        """
        Sets the parameters that define the matematical function.

        Parameters
        ----------
        alpha: array of scalars
            First group of parameters.
        beta: array of scalars
            Second group of parameters.
        """
        cdef PetscInt na = 0, nb = 0
        cdef PetscScalar *a = NULL
        cdef PetscScalar *b = NULL
        cdef object tmp1, tmp2
        if alpha is not None:
            tmp1 = iarray_s(alpha, &na, &a)
        if beta is not None:
            tmp2 = iarray_s(beta, &nb, &b)
        CHKERR( FNSetParameters(self.fn, na, a, nb, b) )

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

# -----------------------------------------------------------------------------

del FNType

# -----------------------------------------------------------------------------
