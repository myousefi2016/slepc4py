# -----------------------------------------------------------------------------

class RGType(object):
    """
    RG type
    """
    INTERVAL   = S_(RGINTERVAL)
    POLYGON    = S_(RGPOLYGON)
    ELLIPSE    = S_(RGELLIPSE)
    RING       = S_(RGRING)

# -----------------------------------------------------------------------------

cdef class RG(Object):

    """
    RG
    """

    Type             = RGType

    def __cinit__(self):
        self.obj = <PetscObject*> &self.rg
        self.rg = NULL

    def view(self, Viewer viewer=None):
        """
        Prints the RG data structure.

        Parameters
        ----------
        viewer: Viewer, optional
                Visualization context; if not provided, the standard
                output is used.
        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( RGView(self.rg, vwr) )

    def destroy(self):
        """
        Destroys the RG object.
        """
        CHKERR( RGDestroy(&self.rg) )
        self.rg = NULL
        return self

    def create(self, comm=None):
        """
        Creates the RG object.

        Parameters
        ----------
        comm: Comm, optional
              MPI communicator; if not provided, it defaults to all
              processes.
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcRG newrg = NULL
        CHKERR( RGCreate(ccomm, &newrg) )
        SlepcCLEAR(self.obj); self.rg = newrg
        return self

    def setType(self, rg_type):
        """
        Selects the type for the RG object.

        Parameters
        ----------
        rg_type: `RG.Type` enumerate
                  The inner product type to be used.
        """
        cdef SlepcRGType cval = NULL
        rg_type = str2bytes(rg_type, &cval)
        CHKERR( RGSetType(self.rg, cval) )

    def getType(self):
        """
        Gets the RG type of this object.

        Returns
        -------
        type: `RG.Type` enumerate
              The inner product type currently being used.
        """
        cdef SlepcRGType rg_type = NULL
        CHKERR( RGGetType(self.rg, &rg_type) )
        return bytes2str(rg_type)

    def setOptionsPrefix(self, prefix):
        """
        Sets the prefix used for searching for all RG options in the
        database.

        Parameters
        ----------
        prefix: string
                The prefix string to prepend to all RG option
                requests.

        Notes
        -----
        A hyphen (``-``) must NOT be given at the beginning of the
        prefix name.  The first character of all runtime options is
        AUTOMATICALLY the hyphen.
        """
        cdef const_char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( RGSetOptionsPrefix(self.rg, cval) )

    def getOptionsPrefix(self):
        """
        Gets the prefix used for searching for all RG options in the
        database.

        Returns
        -------
        prefix: string
                The prefix string set for this RG object.
        """
        cdef const_char *prefix = NULL
        CHKERR( RGGetOptionsPrefix(self.rg, &prefix) )
        return bytes2str(prefix)

    def setFromOptions(self):
        """
        Sets RG options from the options database.

        Notes
        -----
        To see all options, run your program with the ``-help``
        option.
        """
        CHKERR( RGSetFromOptions(self.rg) )

    #

    def isTrivial(self):
        """
        Tells whether it is the trivial region (whole complex plane).

        Returns
        -------
        flag: boolean
             True if the region is equal to the whole complex plane, e.g.,
             an interval region with all four endpoints unbounded or an
             ellipse with infinite radius.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( RGIsTrivial(self.rg, &tval) )
        return <bint> tval

    def getComplement(self):
        """
        Returns the flag indicating whether the region is complemented or not.

        Returns
        -------
        flg: bool
            Whether the region is complemented or not.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( RGGetComplement(self.rg, &tval) )
        return <bint>tval

    def setComplement(self, comp):
        """
        Sets a flag to indicate that the region is the complement
        of the specified one.

        Parameters
        ----------
        comp: bool
            Activate/deactivate the complementation of the region.
        """
        cdef PetscBool tval = comp
        CHKERR( RGSetComplement(self.rg, tval) )

    #

    def setEllipseParameters(self, center, radius, vscale):
        """
        Sets the parameters defining the ellipse region.

        Parameters
        ----------
        center: float (real or complex)
              The center.
        radius: float
              The radius.
        vscale: float
              The vertical scale.
        """
        cdef PetscScalar sval = asScalar(center)
        cdef PetscReal val1 = radius
        cdef PetscReal val2 = vscale
        CHKERR( RGEllipseSetParameters(self.rg, sval, val1, val2) )

    def getEllipseParameters(self):
        """
        Gets the parameters that define the ellipse region.

        Returns
        -------
        center: float (real or complex)
              The center.
        radius: float
              The radius.
        vscale: float
              The vertical scale.
        """
        cdef PetscScalar sval = 0
        cdef PetscReal val1 = 0
        cdef PetscReal val2 = 0
        CHKERR( RGEllipseGetParameters(self.rg, &sval, &val1, &val2) )
        return (toScalar(sval), toReal(val1), toReal(val2))

    def setIntervalEndpoints(self, a, b, c, d):
        """
        Sets the parameters defining the interval region.

        Parameters
        ----------
        a: float
              The left endpoint in the real axis.
        b: float
              The right endpoint in the real axis.
        c: float
              The upper endpoint in the imaginary axis.
        d: float
              The lower endpoint in the imaginary axis.
        """
        cdef PetscReal va = a
        cdef PetscReal vb = b
        cdef PetscReal vc = c
        cdef PetscReal vd = d
        CHKERR( RGIntervalSetEndpoints(self.rg, va, vb, vc, vd) )

    def getIntervalEndpoints(self):
        """
        Gets the parameters that define the interval region.

        Returns
        -------
        a: float
              The left endpoint in the real axis.
        b: float
              The right endpoint in the real axis.
        c: float
              The upper endpoint in the imaginary axis.
        d: float
              The lower endpoint in the imaginary axis.
        """
        cdef PetscReal va = 0
        cdef PetscReal vb = 0
        cdef PetscReal vc = 0
        cdef PetscReal vd = 0
        CHKERR( RGIntervalGetEndpoints(self.rg, &va, &vb, &vc, &vd) )
        return (toReal(va), toReal(vb), toReal(vc), toReal(vd))

# -----------------------------------------------------------------------------

del RGType

# -----------------------------------------------------------------------------
