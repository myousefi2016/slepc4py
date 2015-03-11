import sys, slepc4py
slepc4py.init(sys.argv)

from petsc4py import PETSc
from slepc4py import SLEPc

opts = PETSc.Options()
n  = opts.getInt('n', 30)
mu = opts.getScalar('mu', 1e-6)

PETSc.Sys.Print( "Lauchli singular value decomposition, (%d x %d) mu=%g\n" % (n+1,n,mu) )

A = PETSc.Mat(); A.create()
A.setSizes([n+1, n])
A.setFromOptions( )
A.setUp()

rstart, rend = A.getOwnershipRange()

for i in xrange(rstart, rend):
  if i==0:
    for j in range(n):
      A[0,j] = 1.0
  else:
    A[i,i-1] = mu

A.assemble()

S = SLEPc.SVD(); S.create()

S.setOperator(A)
S.setType(S.Type.TRLANCZOS)
S.setFromOptions()

S.solve()

Print = PETSc.Sys.Print

Print( "******************************" )
Print( "*** SLEPc Solution Results ***" )
Print( "******************************\n" )

svd_type = S.getType()
Print( "Solution method: %s" % svd_type )

its = S.getIterationNumber()
Print( "Number of iterations of the method: %d" % its )

nsv, ncv, mpd = S.getDimensions()
Print( "Number of requested singular values: %d" % nsv )

tol, maxit = S.getTolerances()
Print( "Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit) )

nconv = S.getConverged()
Print( "Number of converged approximate singular triplets %d" % nconv )

if nconv > 0:
  v, u = A.getVecs()
  Print()
  Print("    sigma       residual norm ")
  Print("-------------  ---------------")
  for i in range(nconv):
    sigma = S.getSingularTriplet(i, u, v)
    error = S.computeError(i)
    Print( "   %6f     %12g" % (sigma, error) )
  Print()
