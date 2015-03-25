import sys, slepc4py
slepc4py.init(sys.argv)

from petsc4py import PETSc
from slepc4py import SLEPc

Print = PETSc.Sys.Print


def build_matrix(m):
    """
    Markov model of a random walk on a triangular grid
    """
    N = m*(m+1)/2
    Print("Markov y=exp(t*A)*e_1, N=%d (m=%d)"% (N, m))
    A = PETSc.Mat().create()
    A.setSizes([N, N])
    A.setFromOptions( )
    A.setUp()
    Istart, Iend = A.getOwnershipRange()
    ix = 0
    cst = 0.5/(m-1)
    for i in range(1,m+1):
       jmax = m-i+1
       for j in range(1,jmax+1):
           ix = ix + 1
           if ix-1<Istart or ix>Iend:
               continue  # compute only owned rows
           if j!=jmax:
               pd = cst*(i+j-1)
               # north
               if i==1:
                   A[ix-1,ix] = 2*pd
               else:
                   A[ix-1,ix] = pd
               # east
               if j==1:
                   A[ix-1,ix+jmax-1] = 2*pd
               else:
                   A[ix-1,ix+jmax-1] = pd
           # south
           pu = 0.5 - cst*(i+j-3)
           if j>1:
               A[ix-1,ix-2] = pu
           # west
           if i>1:
               A[ix-1,ix-jmax-2] = pu

    A.assemble()
    return A

def solve_exp(t, A, b, x):
    # Setup the solver
    M = SLEPc.MFN().create()
    M.setOperator(A)
    f = M.getFN()
    f.setType(SLEPc.FN.Type.EXP)
    f.setScale(t)
    M.setTolerances(1e-7)
    M.setFromOptions()
    # Solve the problem
    M.solve(b,x)

    its = M.getIterationNumber()
    Print("Number of iterations of the method: %i" % its)
    sol_type = M.getType()
    Print("Solution method: %s" % sol_type)
    ncv = M.getDimensions()
    Print("")
    Print("Subspace dimension: %i" % ncv)
    tol, maxit = M.getTolerances()
    Print("Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit))
    Print("Computed vector at time t=%.4g has norm %g" % (t.real, x.norm()))
    Print("")

if __name__ == '__main__':
    opts = PETSc.Options()
    m = opts.getInt('m', 15)
    A = build_matrix(m)   # transition probability matrix
    x, b = A.getVecs()
    x.set(0)
    b.set(0)
    b[0] = 1
    b.assemble()
    t = 2
    solve_exp(t, A, b, x)    # compute x=exp(t*A)*b
    A = None
    b = x = None

