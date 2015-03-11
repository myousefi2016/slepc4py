import sys, slepc4py
slepc4py.init(sys.argv)

from petsc4py import PETSc
from slepc4py import SLEPc

Print = PETSc.Sys.Print


def construct_operators(m,n):
    """
    Standard symmetric eigenproblem corresponding to the
    Laplacian operator in 2 dimensions.
    """
    Print("Quadratic Eigenproblem, N=%d (%dx%d grid)"% (m*n, m, n))
    # K is the 2-D Laplacian
    K = PETSc.Mat().create()
    K.setSizes([n*m, n*m])
    K.setFromOptions( )
    K.setUp()
    Istart, Iend = K.getOwnershipRange()
    for I in range(Istart,Iend):
        v = -1.0; i = I/n; j = I-i*n;
        if i>0:
            J=I-n; K[I,J] = v
        if i<m-1:
            J=I+n; K[I,J] = v
        if j>0:
            J=I-1; K[I,J] = v
        if j<n-1:
            J=I+1; K[I,J] = v
        v=4.0; K[I,I] = v
    K.assemble()
    # C is the zero matrix
    C = PETSc.Mat().create()
    C.setSizes([n*m, n*m])
    C.setFromOptions( )
    C.setUp()
    C.assemble()
    # M is the identity matrix
    M = PETSc.Mat().create()
    M.setSizes([n*m, n*m])
    M.setFromOptions( )
    M.setUp()
    M.assemble()
    M.shift(1)
    #
    return M, C, K

def solve_eigensystem(M, C, K):
    # Setup the eigensolver
    Q = SLEPc.PEP().create()
    A = [ ]
    A.append(K)
    A.append(C)
    A.append(M)
    Q.setOperators(A)
    Q.setDimensions(6)
    Q.setProblemType(SLEPc.PEP.ProblemType.GENERAL)
    Q.setFromOptions()
    # Solve the eigensystem
    Q.solve()
    # Create the results vectors
    xr, tmp = K.getVecs()
    xi, tmp = K.getVecs()

    its = Q.getIterationNumber()
    Print("Number of iterations of the method: %i" % its)
    sol_type = Q.getType()
    Print("Solution method: %s" % sol_type)
    nev, ncv, mpd = Q.getDimensions()
    Print("")
    Print("Number of requested eigenvalues: %i" % nev)
    tol, maxit = Q.getTolerances()
    Print("Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit))
    nconv = Q.getConverged()
    Print("Number of converged approximate eigenpairs: %d" % nconv)
    if nconv > 0:
        Print("")
        Print("          k           ||(k^2M+Ck+K)x||/||kx|| ")
        Print("-------------------- -------------------------")
        for i in range(nconv):
            k = Q.getEigenpair(i, xr, xi)
            error = Q.computeError(i)
            if k.imag != 0.0: 
                Print("%9f%+9f j    %12g" % (k.real, k.imag, error))
            else: 
                Print("%12f         %12g" % (k.real, error))
    Print("")

if __name__ == '__main__':
    opts = PETSc.Options()
    m = opts.getInt('m', 32)
    n = opts.getInt('n', m)
    M, C, K = construct_operators(m,n)
    solve_eigensystem(M, C, K)
    M = C = K = None

