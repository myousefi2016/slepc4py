# -----------------------------------------------------------------------------

class PEPType(object):
    """
    PEP type

    Polynomial eigensolvers.

    - `LINEAR`:       Linearization via EPS.
    - `QARNOLDI`:     Q-Arnoldi for quadratic problems.
    - `TOAR`:         Two-level orthogonal Arnoldi.
    - `STOAR`:        Symmetric TOAR.
    - `JD`:           Polynomial Jacobi-Davidson.
    """
    LINEAR   = S_(PEPLINEAR)
    QARNOLDI = S_(PEPQARNOLDI)
    TOAR     = S_(PEPTOAR)
    STOAR    = S_(PEPSTOAR)
    JD       = S_(PEPJD)

class PEPProblemType(object):
    """
    PEP problem type

    - `GENERAL`:      No structure.
    - `HERMITIAN`:    Hermitian structure.
    - `GYROSCOPIC`:   Hamiltonian structure.
    """
    GENERAL    = PEP_GENERAL
    HERMITIAN  = PEP_HERMITIAN
    GYROSCOPIC = PEP_GYROSCOPIC

class PEPWhich(object):
    """
    PEP desired part of spectrum

    - `LARGEST_MAGNITUDE`:  Largest magnitude (default).
    - `LARGEST_REAL`:       Largest real parts.
    - `LARGEST_IMAGINARY`:  Largest imaginary parts in magnitude.
    - `SMALLEST_MAGNITUDE`: Smallest magnitude.
    - `SMALLEST_REAL`:      Smallest real parts.
    - `SMALLEST_IMAGINARY`: Smallest imaginary parts in magnitude.
    - `TARGET_MAGNITUDE`:   Closest to target (in magnitude).
    - `TARGET_REAL`:        Real part closest to target.
    - `TARGET_IMAGINARY`:   Imaginary part closest to target.
    - `USER`:               User-defined criterion.
    """
    LARGEST_MAGNITUDE  = PEP_LARGEST_MAGNITUDE
    SMALLEST_MAGNITUDE = PEP_SMALLEST_MAGNITUDE
    LARGEST_REAL       = PEP_LARGEST_REAL
    SMALLEST_REAL      = PEP_SMALLEST_REAL
    LARGEST_IMAGINARY  = PEP_LARGEST_IMAGINARY
    SMALLEST_IMAGINARY = PEP_SMALLEST_IMAGINARY
    TARGET_MAGNITUDE   = PEP_TARGET_MAGNITUDE
    TARGET_REAL        = PEP_TARGET_REAL
    TARGET_IMAGINARY   = PEP_TARGET_IMAGINARY
    USER               = PEP_WHICH_USER

class PEPBasis(object):
    MONOMIAL   = PEP_BASIS_MONOMIAL
    CHEBYSHEV1 = PEP_BASIS_CHEBYSHEV1
    CHEBYSHEV2 = PEP_BASIS_CHEBYSHEV2
    LEGENDRE   = PEP_BASIS_LEGENDRE
    LAGUERRE   = PEP_BASIS_LAGUERRE
    HERMITE    = PEP_BASIS_HERMITE

class PEPScale(object):
    """
    PEP scaling strategy

    - `NONE`:     No scaling.
    - `SCALAR`:   Parameter scaling.
    - `DIAGONAL`: Diagonal scaling.
    - `BOTH`:     Both parameter and diagonal scaling.
    """
    NONE     = PEP_SCALE_NONE
    SCALAR   = PEP_SCALE_SCALAR
    DIAGONAL = PEP_SCALE_DIAGONAL
    BOTH     = PEP_SCALE_BOTH

class PEPRefine(object):
    """
    PEP refinement strategy

    - `NONE`:     No refinement.
    - `SIMPLE`:   Refine eigenpairs one by one.
    - `MULTIPLE`: Refine all eigenpairs simultaneously (invariant pair).
    """
    NONE     = PEP_REFINE_NONE
    SIMPLE   = PEP_REFINE_SIMPLE
    MULTIPLE = PEP_REFINE_MULTIPLE

class PEPRefineScheme(object):
    """
    Scheme for solving linear systems during iterative refinement

    - `SCHUR`:    Schur complement.
    - `MBE`:      Mixed block elimination.
    - `EXPLICIT`: Build the explicit matrix.
    """
    SCHUR    = PEP_REFINE_SCHEME_SCHUR
    MBE      = PEP_REFINE_SCHEME_MBE
    EXPLICIT = PEP_REFINE_SCHEME_EXPLICIT

class PEPExtract(object):
    """
    Extraction strategy used to obtain eigenvectors of the PEP from the
    eigenvectors of the linearization

    - `NONE`:       Use the first block.
    - `NORM`:       Use the first or last block depending on norm of H.
    - `RESIDUAL`:   Use the block with smallest residual.
    - `STRUCTURED`: Combine all blocks in a certain way.
    """
    NONE       = PEP_EXTRACT_NONE
    NORM       = PEP_EXTRACT_NORM
    RESIDUAL   = PEP_EXTRACT_RESIDUAL
    STRUCTURED = PEP_EXTRACT_STRUCTURED

class PEPErrorType(object):
    """
    PEP error type to assess accuracy of computed solutions

    - `ABSOLUTE`:  Absolute error.
    - `RELATIVE`:  Relative error.
    - `BACKWARD`:  Backward error.
    """
    ABSOLUTE = PEP_ERROR_ABSOLUTE
    RELATIVE = PEP_ERROR_RELATIVE
    BACKWARD = PEP_ERROR_BACKWARD

class PEPConv(object):
    """
    PEP convergence test

    - `ABS`:
    - `REL`:
    - `NORM`:
    - `USER`:
    """
    ABS  = PEP_CONV_ABS
    REL  = PEP_CONV_REL
    NORM = PEP_CONV_NORM
    USER = PEP_CONV_USER

class PEPConvergedReason(object):
    """
    PEP convergence reasons

    - `CONVERGED_TOL`:
    - `CONVERGED_USER`:
    - `DIVERGED_ITS`:
    - `DIVERGED_BREAKDOWN`:
    - `DIVERGED_SYMMETRY_LOST`:
    - `CONVERGED_ITERATING`:
    """
    CONVERGED_TOL          = PEP_CONVERGED_TOL
    CONVERGED_USER         = PEP_CONVERGED_USER
    DIVERGED_ITS           = PEP_DIVERGED_ITS
    DIVERGED_BREAKDOWN     = PEP_DIVERGED_BREAKDOWN
    DIVERGED_SYMMETRY_LOST = PEP_DIVERGED_SYMMETRY_LOST
    CONVERGED_ITERATING    = PEP_CONVERGED_ITERATING
    ITERATING              = PEP_CONVERGED_ITERATING

# -----------------------------------------------------------------------------

cdef class PEP(Object):

    """
    PEP
    """

    Type            = PEPType
    ProblemType     = PEPProblemType
    Which           = PEPWhich
    Basis           = PEPBasis
    Scale           = PEPScale
    Refine          = PEPRefine
    RefineScheme    = PEPRefineScheme
    Extract         = PEPExtract
    ErrorType       = PEPErrorType
    Conv            = PEPConv
    ConvergedReason = PEPConvergedReason

    def __cinit__(self):
        self.obj = <PetscObject*> &self.pep
        self.pep = NULL

    def view(self, Viewer viewer=None):
        """
        Prints the PEP data structure.

        Parameters
        ----------
        viewer: Viewer, optional.
            Visualization context; if not provided, the standard
            output is used.
        """
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( PEPView(self.pep, vwr) )

    def destroy(self):
        """
        Destroys the PEP object.
        """
        CHKERR( PEPDestroy(&self.pep) )
        self.pep = NULL
        return self

    def reset(self):
        """
        Resets the PEP object.
        """
        CHKERR( PEPReset(self.pep) )

    def create(self, comm=None):
        """
        Creates the PEP object.

        Parameters
        ----------
        comm: Comm, optional.
            MPI communicator. If not provided, it defaults to all
            processes.
        """
        cdef MPI_Comm ccomm = def_Comm(comm, SLEPC_COMM_DEFAULT())
        cdef SlepcPEP newpep = NULL
        CHKERR( PEPCreate(ccomm, &newpep) )
        SlepcCLEAR(self.obj); self.pep = newpep
        return self

    def setType(self, pep_type):
        """
        Selects the particular solver to be used in the PEP object.

        Parameters
        ----------
        pep_type: `PEP.Type` enumerate
            The solver to be used.
        """
        cdef SlepcPEPType cval = NULL
        pep_type = str2bytes(pep_type, &cval)
        CHKERR( PEPSetType(self.pep, cval) )

    def getType(self):
        """
        Gets the PEP type of this object.

        Returns
        -------
        type: `PEP.Type` enumerate
            The solver currently being used.
        """
        cdef SlepcPEPType pep_type = NULL
        CHKERR( PEPGetType(self.pep, &pep_type) )
        return bytes2str(pep_type)

    def getOptionsPrefix(self):
        """
        Gets the prefix used for searching for all PEP options in the
        database.

        Returns
        -------
        prefix: string
            The prefix string set for this PEP object.
        """
        cdef const_char *prefix = NULL
        CHKERR( PEPGetOptionsPrefix(self.pep, &prefix) )
        return bytes2str(prefix)

    def setOptionsPrefix(self, prefix):
        """
        Sets the prefix used for searching for all PEP options in the
        database.

        Parameters
        ----------
        prefix: string
            The prefix string to prepend to all PEP option requests.
        """
        cdef const_char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( PEPSetOptionsPrefix(self.pep, cval) )

    def appendOptionsPrefix(self, prefix):
        """
        Appends to the prefix used for searching for all PEP options
        in the database.

        Parameters
        ----------
        prefix: string
            The prefix string to prepend to all PEP option requests.
        """
        cdef const_char *cval = NULL
        prefix = str2bytes(prefix, &cval)
        CHKERR( PEPAppendOptionsPrefix(self.pep, cval) )

    def setFromOptions(self):
        """
        Sets PEP options from the options database. This routine must
        be called before `setUp()` if the user is to be allowed to set
        the solver type.
        """
        CHKERR( PEPSetFromOptions(self.pep) )

    def getBasis(self):
        """
        Gets the type of polynomial basis used to 
        describe the polynomial eigenvalue problem.

        Returns
        -------
        basis: `PEP.Basis` enumerate
            the basis that was previously set.
        """
        cdef SlepcPEPBasis val = PEP_BASIS_MONOMIAL
        CHKERR( PEPGetBasis(self.pep, &val) )
        return val

    def setBasis(self, basis):
        """
        Specifies the type of polynomial basis used to 
        describe the polynomial eigenvalue problem.

        Parameters
        ----------
        basis: `PEP.Basis` enumerate
            the basis to be set.
        """
        cdef SlepcPEPBasis val = basis
        CHKERR( PEPSetBasis(self.pep, val) )

    def getProblemType(self):
        """
        Gets the problem type from the PEP object.

        Returns
        -------
        problem_type: `PEP.ProblemType` enumerate
            The problem type that was previously set.
        """
        cdef SlepcPEPProblemType val = PEP_GENERAL
        CHKERR( PEPGetProblemType(self.pep, &val) )
        return val

    def setProblemType(self, problem_type):
        """
        Specifies the type of the eigenvalue problem.

        Parameters
        ----------
        problem_type: `PEP.ProblemType` enumerate
            The problem type to be set.
        """
        cdef SlepcPEPProblemType val = problem_type
        CHKERR( PEPSetProblemType(self.pep, val) )

    def getWhichEigenpairs(self):
        """
        Returns which portion of the spectrum is to be sought.

        Returns
        -------
        which: `PEP.Which` enumerate
            The portion of the spectrum to be sought by the solver.
        """
        cdef SlepcPEPWhich val = PEP_LARGEST_MAGNITUDE
        CHKERR( PEPGetWhichEigenpairs(self.pep, &val) )
        return val

    def setWhichEigenpairs(self, which):
        """
        Specifies which portion of the spectrum is to be sought.

        Parameters
        ----------
        which: `PEP.Which` enumerate
            The portion of the spectrum to be sought by the solver.
        """
        cdef SlepcPEPWhich val = which
        CHKERR( PEPSetWhichEigenpairs(self.pep, val) )

    def getTolerances(self):
        """
        Gets the tolerance and maximum iteration count used by the
        default PEP convergence tests.

        Returns
        -------
        tol: float
            The convergence tolerance.
        max_it: int
            The maximum number of iterations
        """
        cdef PetscReal rval = 0
        cdef PetscInt  ival = 0
        CHKERR( PEPGetTolerances(self.pep, &rval, &ival) )
        return (toReal(rval), toInt(ival))

    def setTolerances(self, tol=None, max_it=None):
        """
        Sets the tolerance and maximum iteration count used by the
        default PEP convergence tests.

        Parameters
        ----------
        tol: float, optional
            The convergence tolerance.
        max_it: int, optional
            The maximum number of iterations
        """
        cdef PetscReal rval = PETSC_DEFAULT
        cdef PetscInt  ival = PETSC_DEFAULT
        if tol    is not None: rval = asReal(tol)
        if max_it is not None: ival = asInt(max_it)
        CHKERR( PEPSetTolerances(self.pep, rval, ival) )

    def getConvergenceTest(self):
        """
        Return the method used to compute the error estimate 
        used in the convergence test. 

        Returns
        -------
        conv: PEP.Conv
            The method used to compute the error estimate 
            used in the convergence test. 
        """
        cdef SlepcPEPConv conv = PEP_CONV_REL
        CHKERR( PEPGetConvergenceTest(self.pep, &conv) )
        return conv

    def setConvergenceTest(self, conv):
        """
        Specifies how to compute the error estimate 
        used in the convergence test. 

        Parameters
        ----------
        conv: PEP.Conv
            The method used to compute the error estimate 
            used in the convergence test.
        """
        cdef SlepcPEPConv tconv = conv
        CHKERR( PEPSetConvergenceTest(self.pep, tconv) )

    def getRefine(self):
        """
        Gets the refinement strategy used by the PEP object, 
        and the associated parameters. 

        Returns
        -------
        ref: PEP.Refine
            The refinement type.
        npart: int
            The number of partitions of the communicator.
        tol: real
            The convergence tolerance.
        its: int
            The maximum number of refinement iterations.
        scheme: PEP.RefineScheme
            Scheme for solving linear systems
        """
        cdef SlepcPEPRefine ref = PEP_REFINE_NONE
        cdef PetscInt npart = 1
        cdef PetscReal tol = PETSC_DEFAULT
        cdef PetscInt its = PETSC_DEFAULT
        cdef SlepcPEPRefineScheme scheme = PEP_REFINE_SCHEME_MBE
        CHKERR( PEPGetRefine(self.pep, &ref, &npart, &tol, &its, &scheme) )
        return (ref, toInt(npart), toReal(tol), toInt(its), scheme)

    def setRefine(self, ref, npart=None, tol=None, its=None, scheme=None):
        """
        Sets the refinement strategy used by the PEP object, 
        and the associated parameters. 

        Parameters
        ----------
        ref: PEP.Refine
            The refinement type.
        npart: int, optional
            The number of partitions of the communicator.
        tol: real, optional
            The convergence tolerance.
        its: int, optional
            The maximum number of refinement iterations.
        scheme: PEP.RefineScheme, optional
            Scheme for linear system solves
        """
        cdef SlepcPEPRefine tref = ref
        cdef PetscInt tnpart = 1
        cdef PetscReal ttol = PETSC_DEFAULT
        cdef PetscInt tits = PETSC_DEFAULT
        cdef SlepcPEPRefineScheme tscheme = PEP_REFINE_SCHEME_MBE
        if npart is not None: tnpart = asInt(npart)
        if tol is not None: ttol = asReal(tol)
        if its is not None: tits = asInt(its)
        if scheme is not None: tscheme = scheme
        CHKERR( PEPSetRefine(self.pep, tref, tnpart, ttol, tits, tscheme) )

    def getTrackAll(self):
        """
        Returns the flag indicating whether all residual norms must be
        computed or not.

        Returns
        -------
        trackall: bool
            Whether the solver compute all residuals or not.
        """
        cdef PetscBool tval = PETSC_FALSE
        CHKERR( PEPGetTrackAll(self.pep, &tval) )
        return <bint>tval

    def setTrackAll(self, trackall):
        """
        Specifies if the solver must compute the residual of all
        approximate eigenpairs or not.

        Parameters
        ----------
        trackall: bool
            Whether compute all residuals or not.
        """
        cdef PetscBool tval = trackall
        CHKERR( PEPSetTrackAll(self.pep, tval) )

    def getDimensions(self):
        """
        Gets the number of eigenvalues to compute and the dimension of
        the subspace.

        Returns
        -------
        nev: int
            Number of eigenvalues to compute.
        ncv: int
            Maximum dimension of the subspace to be used by the solver.
        mpd: int
            Maximum dimension allowed for the projected problem.
        """
        cdef PetscInt ival1 = 0
        cdef PetscInt ival2 = 0
        cdef PetscInt ival3 = 0
        CHKERR( PEPGetDimensions(self.pep, &ival1, &ival2, &ival3) )
        return (toInt(ival1), toInt(ival2), toInt(ival3))

    def setDimensions(self, nev=None, ncv=None, mpd=None):
        """
        Sets the number of eigenvalues to compute and the dimension of
        the subspace.

        Parameters
        ----------
        nev: int, optional
            Number of eigenvalues to compute.
        ncv: int, optional
            Maximum dimension of the subspace to be used by the
            solver.
        mpd: int, optional
            Maximum dimension allowed for the projected problem.
        """
        cdef PetscInt ival1 = PETSC_DEFAULT
        cdef PetscInt ival2 = PETSC_DEFAULT
        cdef PetscInt ival3 = PETSC_DEFAULT
        if nev is not None: ival1 = asInt(nev)
        if ncv is not None: ival2 = asInt(ncv)
        if mpd is not None: ival3 = asInt(mpd)
        CHKERR( PEPSetDimensions(self.pep, ival1, ival2, ival3) )

    def getST(self):
        """
        Obtain the spectral transformation (`ST`) object associated to
        the eigensolver object.

        Returns
        -------
        st: ST
            The spectral transformation.
        """
        cdef ST st = ST()
        CHKERR( PEPGetST(self.pep, &st.st) )
        PetscINCREF(st.obj)
        return st

    def setST(self, ST st not None):
        """
        Associates a spectral transformation object to the
        eigensolver.

        Parameters
        ----------
        st: ST
            The spectral transformation.
        """
        CHKERR( PEPSetST(self.pep, st.st) )

    def getScale(self, Vec Dl=None, Vec Dr=None):
        """
        Gets the strategy used for scaling the polynomial eigenproblem.

        Parameters
        ----------
        Dl: Vec, optional
            Placeholder for the returned left diagonal matrix.
        Dr: Vec, optional
            Placeholder for the returned right diagonal matrix.

        Returns
        -------
        scale: `PEP.Scale` enumerate
            The scaling strategy.
        alpha: real
            The scaling factor.
        its: integer
            The number of iteration of diagonal scaling.
        lbda: real
            Approximation of the wanted eigenvalues (modulus).
        """
        cdef SlepcPEPScale scale = PEP_SCALE_NONE
        cdef PetscReal alpha = 0
        cdef PetscInt its = 0
        cdef PetscReal lbda = 0
        cdef PetscVec vecl = NULL
        cdef PetscVec vecr = NULL
        CHKERR( PEPGetScale(self.pep, &scale, &alpha, &vecl, &vecr, &its, &lbda) )
        if Dl.vec != NULL:
            if vecl != NULL:
                CHKERR( VecCopy(vecl, Dl.vec) )
            else:
                CHKERR( VecSet(Dl.vec, 1.0) )
        if Dr.vec != NULL:
            if vecr != NULL:
                CHKERR( VecCopy(vecr, Dr.vec) )
            else:
                CHKERR( VecSet(Dr.vec, 1.0) )
        CHKERR( VecDestroy(&vecl) )
        CHKERR( VecDestroy(&vecr) )
        return (scale, toReal(alpha), toInt(its), toReal(lbda))

    def setScale(self, scale, alpha=None, Vec Dl=None, Vec Dr=None, its=None, lbda=None):
        """
        Sets the scaling strategy to be used for scaling the polynomial problem
        before attempting to solve.

        Parameters
        ----------
        scale: `PEP.Scale` enumerate
            The scaling strategy.
        alpha: real, optional
            The scaling factor.
        Dl: Vec, optional
            The left diagonal matrix.
        Dr: Vec, optional
            The right diagonal matrix.
        its: integer, optional
            The number of iteration of diagonal scaling.
        lbda: real, optional
            Approximation of the wanted eigenvalues (modulus).
        """
        cdef SlepcPEPScale senum = scale
        cdef PetscReal rval1 = PETSC_DEFAULT
        cdef PetscInt ival = PETSC_DEFAULT
        cdef PetscReal rval2 = PETSC_DEFAULT
        cdef PetscVec vecl = NULL
        cdef PetscVec vecr = NULL
        if alpha is not None: rval1 = asReal(alpha)
        if Dl is not None:    vecl = Dl.vec
        if Dr is not None:    vecr = Dr.vec
        if its is not None:   ival = asInt(its)
        if lbda is not None:  rval2 = asReal(lbda)
        CHKERR( PEPSetScale(self.pep, senum, rval1, vecl, vecr, ival, rval2) )

    def getBV(self):
        """
        Obtain the basis vectors object associated to the eigensolver.

        Returns
        -------
        bv: BV
            The basis vectors context.
        """
        cdef BV bv = BV()
        CHKERR( PEPGetBV(self.pep, &bv.bv) )
        PetscINCREF(bv.obj)
        return bv

    def setBV(self, BV bv not None):
        """
        Associates a basis vectors object to the eigensolver.

        Parameters
        ----------
        bv: BV
            The basis vectors context.
        """
        CHKERR( PEPSetBV(self.pep, bv.bv) )

    def getRG(self):
        """
        Obtain the region object associated to the eigensolver.

        Returns
        -------
        rg: RG
            The region context.
        """
        cdef RG rg = RG()
        CHKERR( PEPGetRG(self.pep, &rg.rg) )
        PetscINCREF(rg.obj)
        return rg

    def setRG(self, RG rg not None):
        """
        Associates a region object to the eigensolver.

        Parameters
        ----------
        rg: RG
            The region context.
        """
        CHKERR( PEPSetRG(self.pep, rg.rg) )

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
        CHKERR( PEPGetNumMatrices(self.pep, &n) )
        cdef object operators = []
        for k from 0 <= k < n:
            CHKERR( PEPGetOperators(self.pep, k, &mat) )
            A = Mat(); A.mat = mat; PetscINCREF(A.obj)
            operators.append(A)
        return tuple(operators)

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
        CHKERR( PEPSetOperators(self.pep, <PetscInt>n, mats) )

    #

    def setInitialSpace(self, space):
        """
        Sets the initial space from which the eigensolver starts to
        iterate.

        Parameters
        ----------
        space: Vec or sequence of Vec
           The initial space
        """
        if isinstance(space, Vec): space = [space]
        cdef PetscVec *vs = NULL
        cdef Py_ssize_t i = 0, ns = len(space)
        cdef tmp = allocate(<size_t>ns*sizeof(Vec),<void**>&vs)
        for i in range(ns): vs[i] = (<Vec?>space[i]).vec
        CHKERR( PEPSetInitialSpace(self.pep, <PetscInt>ns, vs) )

    #

    def cancelMonitor(self):
        """
        Clears all monitors for a PEP object.
        """
        CHKERR( PEPMonitorCancel(self.pep) )

    #

    def setUp(self):
        """
        Sets up all the internal data structures necessary for the
        execution of the eigensolver.
        """
        CHKERR( PEPSetUp(self.pep) )

    def solve(self):
        """
        Solves the eigensystem.
        """
        CHKERR( PEPSolve(self.pep) )

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
        CHKERR( PEPGetIterationNumber(self.pep, &ival) )
        return toInt(ival)

    def getConvergedReason(self):
        """
        Gets the reason why the `solve()` iteration was stopped.

        Returns
        -------
        reason: `PEP.ConvergedReason` enumerate
            Negative value indicates diverged, positive value
            converged.
        """
        cdef SlepcPEPConvergedReason val = PEP_CONVERGED_ITERATING
        CHKERR( PEPGetConvergedReason(self.pep, &val) )
        return val


    def getConverged(self):
        """
        Gets the number of converged eigenpairs.

        Returns
        -------
        nconv: int
            Number of converged eigenpairs.
        """
        cdef PetscInt ival = 0
        CHKERR( PEPGetConverged(self.pep, &ival) )
        return toInt(ival)

    def getEigenpair(self, int i, Vec Vr=None, Vec Vi=None):
        """
        Gets the i-th solution of the eigenproblem as computed by
        `solve()`.  The solution consists of both the eigenvalue and
        the eigenvector.

        Parameters
        ----------
        i: int
            Index of the solution to be obtained.
        Vr: Vec, optional
            Placeholder for the returned eigenvector (real part).
        Vi: Vec, optional
            Placeholder for the returned eigenvector (imaginary part).

        Returns
        -------
        e: scalar (possibly complex)
            The computed eigenvalue.
        """
        cdef PetscScalar sval1 = 0
        cdef PetscScalar sval2 = 0
        cdef PetscVec vecr = NULL
        cdef PetscVec veci = NULL
        if Vr is not None: vecr = Vr.vec
        if Vi is not None: veci = Vi.vec
        CHKERR( PEPGetEigenpair(self.pep, i, &sval1, &sval2, vecr, veci) )
        return complex(toScalar(sval1), toScalar(sval2))

    def getErrorEstimate(self, int i):
        """
        Returns the error estimate associated to the i-th computed
        eigenpair.

        Parameters
        ----------
        i: int
            Index of the solution to be considered.

        Returns
        -------
        error: real
            Error estimate.
        """
        cdef PetscReal rval = 0
        CHKERR( PEPGetErrorEstimate(self.pep, i, &rval) )
        return toReal(rval)

    def computeError(self, int i, etype=None):
        """
        Computes the error (based on the residual norm) associated with the i-th
        computed eigenpair.

        Parameters
        ----------
        i: int
           Index of the solution to be considered.
        etype: `PEP.ErrorType` enumerate
           The error type to compute.

        Returns
        -------
        error: real
           The error bound, computed in various ways from the
           residual norm ``||P(l)x||_2`` where ``l`` is the
           eigenvalue and ``x`` is the eigenvector.

        Notes
        -----
        The index ``i`` should be a value between ``0`` and
        ``nconv-1`` (see `getConverged()`).
        """
        cdef SlepcPEPErrorType et = PEP_ERROR_BACKWARD
        cdef PetscReal rval = 0
        if etype is not None: et = etype
        CHKERR( PEPComputeError(self.pep, i, et, &rval) )
        return toReal(rval)

    def errorView(self, etype=None, Viewer viewer=None):
        """
        Displays the errors associated with the computed solution
        (as well as the eigenvalues).

        Parameters
        ----------
        etype: `PEP.ErrorType` enumerate, optional
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
        cdef SlepcPEPErrorType et = PEP_ERROR_RELATIVE
        if etype is not None: et = etype
        cdef PetscViewer vwr = NULL
        if viewer is not None: vwr = viewer.vwr
        CHKERR( PEPErrorView(self.pep, et, vwr) )

    #

    def setLinearEPS(self, EPS eps not None):
        """
        Associate an eigensolver object (EPS) to the polynomial eigenvalue solver.

        Parameters
        ----------
        eps: EPS
            The linear eigensolver.
        """
        CHKERR( PEPLinearSetEPS(self.pep, eps.eps) )

    def getLinearEPS(self):
        """
        Retrieve the eigensolver object (EPS) associated to the polynomial
        eigenvalue solver.
 
        Returns
        -------
        eps: EPS
            The linear eigensolver.
        """
        cdef EPS eps = EPS()
        CHKERR( PEPLinearGetEPS(self.pep, &eps.eps) )
        PetscINCREF(eps.obj)
        return eps
        
    def setLinearCompanionForm(self, cform not None):
        """
        Choose between the two companion forms available for the linearization of
        a quadratic eigenproblem.

        Parameters
        ----------
        cform: integer
            1 or 2 (first or second companion form).
        """
        CHKERR( PEPLinearSetCompanionForm(self.pep, cform) )

    def getLinearCompanionForm(self):
        """
        Returns the number of the companion form that will be used for the
        linearization of a quadratic eigenproblem. 
 
        Returns
        -------
        cform: integer
            1 or 2 (first or second companion form).
        """
        cdef PetscInt cform = 0
        CHKERR( PEPLinearGetCompanionForm(self.pep, &cform) )
        return cform
        
    def setLinearExplicitMatrix(self, flag not None):
        """
        Indicate if the matrices A and B for the linearization of the problem
        must be built explicitly.

        Parameters
        ----------
        flag: boolean
            boolean flag indicating if the matrices are built explicitly .
        """
        cdef PetscBool sval = flag
        CHKERR( PEPLinearSetExplicitMatrix(self.pep, sval) )

    def getLinearExplicitMatrix(self):
        """
        Returns the flag indicating if the matrices A and B for the linearization
        are built explicitly.
 
        Returns
        -------
        flag: boolean
        """
        cdef PetscBool sval = PETSC_FALSE
        CHKERR( PEPLinearGetExplicitMatrix(self.pep, &sval) )
        return sval

# -----------------------------------------------------------------------------

del PEPType
del PEPProblemType
del PEPWhich
del PEPBasis
del PEPScale
del PEPRefine
del PEPRefineScheme
del PEPExtract
del PEPErrorType
del PEPConv
del PEPConvergedReason

# -----------------------------------------------------------------------------
