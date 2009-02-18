import sys, petsc4py
petsc4py.init(sys.argv)

from petsc4py import PETSc
from slepc4py import SLEPc

from slepc4py.SLEPc import EPS
from petsc4py.PETSc import Mat

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


def laplace2d(x, f):
    # setup 5-points stencil
    u  = x[1:-1, 1:-1] # center
    uN = x[1:-1,  :-2] # north
    uS = x[1:-1, 2:  ] # south
    uW = x[ :-2, 1:-1] # west
    uE = x[2:,   1:-1] # east
    # apply Laplacian
    nx, ny = x.shape
    hx = 1.0/(nx-1) # x grid spacing
    hy = 1.0/(ny-1) # y grid spacing
    f[:,:] = x
    f[1:-1, 1:-1] = \
         (2*u - uE - uW) * (hy/hx) \
       + (2*u - uN - uS) * (hx/hy) \


class Laplacian2D(object):

    def __init__(self, M, N):
        self.M = M
        self.N = N

    def mult(self, A, x, y):
        m = self.N
        n = self.N
        xx = x[...].reshape(m,n)
        yy = y[...].reshape(m,n)
        laplace2d(xx, yy)


def construct_operator():
    """
    Standard symmetric eigenproblem corresponding to the
    Laplacian operator in 2 dimensions. Uses *shell* matrix.
    """
    # Problem size
    opts = PETSc.Options()
    m = n = opts.getInt('N', 40)
    # Create shell matrix
    context = Laplacian2D(m,n)
    A = Mat().createPython([m*n,m*n], context)
    return A

if __name__ == '__main__':
    A = construct_operator()
    solve_eigensystem(A)
