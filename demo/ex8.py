# ------------------------------------------------------------------------
#   Solve parabolic partial differential equation with time delay tau
#
#            u_t = u_xx + a*u(t) + b*u(t-tau)
#            u(0,t) = u(pi,t) = 0
#
#   with a = 20 and b(x) = -4.1+x*(1-exp(x-pi)).
#
#   Discretization leads to a DDE of dimension n
#
#            -u' = A*u(t) + B*u(t-tau)
#
#   which results in the nonlinear eigenproblem
#
#            (-lambda*I + A + exp(-tau*lambda)*B)*u = 0
# ------------------------------------------------------------------------

import sys, slepc4py
slepc4py.init(sys.argv)

from petsc4py import PETSc
from slepc4py import SLEPc
from numpy import exp
from math import pi

Print = PETSc.Sys.Print

opts = PETSc.Options()
n = opts.getInt('n', 128)
tau = opts.getReal('tau', 0.001)
a = 20
h = pi/(n+1)

# Setup the solver
nep = SLEPc.NEP().create()

# Create problem matrices
#   Identity matrix
Id = PETSc.Mat().create()
Id.setSizes([n, n])
Id.setFromOptions()
Id.setUp()
Id.assemble()
Id.shift(1.0)
Id.setOption(PETSc.Mat.Option.HERMITIAN, True)
#   A = 1/h^2*tridiag(1,-2,1) + a*I
A = PETSc.Mat().create()
A.setSizes([n, n])
A.setFromOptions()
A.setUp()
rstart, rend = A.getOwnershipRange()
vd = -2.0/(h*h)+a
vo = 1.0/(h*h)
if rstart == 0:
  A[0, :2] = [vd, vo]
  rstart += 1
if rend == n:
  A[n-1, -2:] = [vo, vd]
  rend -= 1
for i in range(rstart, rend):
  A[i, i-1:i+2] = [vo, vd, vo]
A.assemble()
#   B = diag(b(xi))
B = PETSc.Mat().create()
B.setSizes([n, n])
B.setFromOptions()
B.setUp()
rstart, rend = B.getOwnershipRange()
for i in range(rstart, rend):
  xi = (i+1)*h
  B[i, i] = -4.1+xi*(1.0-exp(xi-pi));
B.assemble()
B.setOption(PETSc.Mat.Option.HERMITIAN, True)

# Functions: f1=-lambda, f2=1.0, f3=exp(-tau*lambda)
f1 = SLEPc.FN().create()
f1.setType(SLEPc.FN.Type.RATIONAL)
f1.setRationalNumerator([-1, 0])
f2 = SLEPc.FN().create()
f2.setType(SLEPc.FN.Type.RATIONAL)
f2.setRationalNumerator([1])
f3 = SLEPc.FN().create()
f3.setType(SLEPc.FN.Type.EXP)
f3.setScale(-tau)

# Set the split operator. Note that A is passed first so that
#      SUBSET_NONZERO_PATTERN can be used
nep.setSplitOperator([A, Id, B], [f2, f1, f3], PETSc.Mat.Structure.SUBSET)

# Customize options
nep.setTolerances(tol=1e-9)
nep.setDimensions(1)
nep.setFromOptions()

# Solve the problem
nep.solve()

its = nep.getIterationNumber()
Print("Number of iterations of the method: %i" % its)
sol_type = nep.getType()
Print("Solution method: %s" % sol_type)
nev, ncv, mpd = nep.getDimensions()
Print("")
Print("Subspace dimension: %i" % ncv)
tol, maxit = nep.getTolerances()
Print("Stopping condition: tol=%.4g" % tol)
Print("")

nconv = nep.getConverged()
Print( "Number of converged eigenpairs %d" % nconv )

if nconv > 0:
  x, z = Id.getVecs()
  x.set(1.0)
  Print()
  Print("        k              ||T(k)x||")
  Print("----------------- ------------------")
  for i in range(nconv):
    k = nep.getEigenpair(i, x)
    res = nep.computeError(i)
    if k.imag != 0.0:
      Print( " %9f%+9f j %12g" % (k.real, k.imag, res) )
    else:
      Print( " %12f       %12g" % (k.real, res) )
  Print()


