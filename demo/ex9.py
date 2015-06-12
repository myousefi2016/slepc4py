import sys, slepc4py
slepc4py.init(sys.argv)

from petsc4py import PETSc
from slepc4py import SLEPc

Print = PETSc.Sys.Print

def Laplacian2D(m, n):
    """
    Builds discretized Laplacian operator in 2 dimensions.
    """
    # Create matrix for 2D Laplacian operator
    A = PETSc.Mat().create()
    A.setSizes([m*n, m*n])
    A.setFromOptions( )
    A.setUp()
    # Fill matrix
    hx = 1.0/(m-1) # x grid spacing
    hy = 1.0/(n-1) # y grid spacing
    diagv = 2.0*hy/hx + 2.0*hx/hy
    offdx = -1.0*hy/hx
    offdy = -1.0*hx/hy
    Istart, Iend = A.getOwnershipRange()
    for I in xrange(Istart, Iend):
        A[I,I] = diagv
        i = I//n    # map row number to
        j = I - i*n # grid coordinates
        if i> 0  : J = I-n; A[I,J] = offdx
        if i< m-1: J = I+n; A[I,J] = offdx
        if j> 0  : J = I-1; A[I,J] = offdy
        if j< n-1: J = I+1; A[I,J] = offdy
    A.assemble()
    return A

def QuasiDiagonal(N):
    """
    Builds matrix diag(2)+[4 -1; -1 -1]
    """
    # Create matrix
    B = PETSc.Mat().create()
    B.setSizes([N, N])
    B.setFromOptions( )
    B.setUp()
    # Fill matrix
    Istart, Iend = B.getOwnershipRange()
    for I in xrange(Istart, Iend):
        B[I,I] = 2.0
    if Istart==0:
        B[0,0] = 6.0
        B[0,1] = -1.0
        B[1,0] = -1.0
        B[1,1] = 1.0
    B.assemble()
    return B

def solve_eigensystem(A, B, problem_type=SLEPc.EPS.ProblemType.GHEP):
    # Create the results vectors
    xr, tmp = A.getVecs()
    xi, tmp = A.getVecs()

    pc = PETSc.PC().create()
    # pc.setType(pc.Type.HYPRE)
    pc.setType(pc.Type.BJACOBI)
    
    ksp = PETSc.KSP().create()
    ksp.setType(ksp.Type.PREONLY)
    ksp.setPC( pc )
    
    
    F = SLEPc.ST().create()
    F.setType(F.Type.PRECOND)
    F.setKSP( ksp )
    F.setShift(0)

    # Setup the eigensolver
    E = SLEPc.EPS().create()
    E.setST(F)
    E.setOperators(A,B)
    E.setType(E.Type.LOBPCG)
    E.setDimensions(10,PETSc.DECIDE)
    E.setWhichEigenpairs(E.Which.SMALLEST_REAL)
    E.setProblemType( problem_type )
    E.setFromOptions()

    # Solve the eigensystem
    E.solve()

    Print("")
    its = E.getIterationNumber()
    Print("Number of iterations of the method: %i" % its)
    sol_type = E.getType()
    Print("Solution method: %s" % sol_type)
    nev, ncv, mpd = E.getDimensions()
    Print("Number of requested eigenvalues: %i" % nev)
    tol, maxit = E.getTolerances()
    Print("Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit))
    nconv = E.getConverged()
    Print("Number of converged eigenpairs: %d" % nconv)
    if nconv > 0:
        Print("")
        Print("        k          ||Ax-kx||/||kx|| ")
        Print("----------------- ------------------")
        for i in range(nconv):
            k = E.getEigenpair(i, xr, xi)
            error = E.computeError(i)
            if k.imag != 0.0:
              Print(" %9f%+9f j  %12g" % (k.real, k.imag, error))
            else:
              Print(" %12f       %12g" % (k.real, error))
        Print("")

def main():
    opts = PETSc.Options()
    N = opts.getInt('N', 10)
    m = opts.getInt('m', N)
    n = opts.getInt('n', m)
    Print("Symmetric-definite Eigenproblem, "
          "N=%d (%dx%d grid)" % (m*n, m, n))
    A = Laplacian2D(m,n)
    B = QuasiDiagonal(m*n)
    solve_eigensystem(A,B)

if __name__ == '__main__':
    main()
