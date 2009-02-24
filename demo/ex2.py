import sys, petsc4py
petsc4py.init(sys.argv)

from petsc4py import PETSc
from slepc4py import SLEPc

from petsc4py.PETSc import Mat
from slepc4py.SLEPc import EPS


def solve_eigensystem(A, problem_type=EPS.ProblemType.HEP):
    # Create the results vectors
    xr, tmp = A.getVecs()
    xi, tmp = A.getVecs()

    # Setup the eigensolver
    E = EPS().create()
    E.setOperators(A,None)
    E.setDimensions(3,PETSc.DECIDE)
    E.setProblemType( problem_type )

    # Solve the eigensystem
    E.solve()

    print
    its = E.getIterationNumber()
    print "Number of iterations of the method: %i" % its
    sol_type = E.getType()
    print "Solution method: %s" % sol_type
    nev, ncv, mpd = E.getDimensions()
    print "Number of requested eigenvalues: %i" % nev
    tol, maxit = E.getTolerances()
    print "Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit)
    nconv = E.getConverged()
    print "Number of converged eigenpairs: %d" % nconv
    if nconv > 0:
        print
        print "        k          ||Ax-kx||/||kx|| "
        print "----------------- ------------------"
        for i in range(nconv):
            k = E.getEigenpair(i, xr, xi)
            error = E.computeRelativeError(i)
            if k.imag != 0.0: print " %9f%+9f j %12g" % (k.real, k.imag, error)
            else: print " %12f       %12g" % (k.real, error)
        print

def construct_operator():
    """
    Standard symmetric eigenproblem corresponding to the
    Laplacian operator in 2 dimensions.
    """
    # Problem size
    opts = PETSc.Options()
    m = n = opts.getInt('N', 40)
    # Create matrix for Laplacian operator
    A = Mat().create()
    A.setSizes([n*m, n*m])
    A.setFromOptions( )
    # Fill matrix
    Istart, Iend = A.getOwnershipRange()
    for I in range(Istart,Iend):
        v = -1.0; i = I/n; j = I-i*n;
        if i>0:
            J=I-n; A[I,J] = v
        if i<m-1:
            J=I+n; A[I,J] = v
        if j>0:
            J=I-1; A[I,J] = v
        if j<n-1:
            J=I+1; A[I,J] = v
        v=4.0; A[I,I] = v
    A.assemble()
    return A

if __name__ == '__main__':
    A = construct_operator()
    solve_eigensystem(A)

