import sys, slepc4py
slepc4py.init(sys.argv)

from petsc4py import PETSc
from slepc4py import SLEPc

Print = PETSc.Sys.Print

def construct_operator(m, n):
    """
    Standard symmetric eigenproblem corresponding to the
    Laplacian operator in 2 dimensions.
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
    for I in xrange(Istart, Iend) :
        A[I,I] = diagv
        i = I//n    # map row number to
        j = I - i*n # grid coordinates
        if i> 0  : J = I-n; A[I,J] = offdx
        if i< m-1: J = I+n; A[I,J] = offdx
        if j> 0  : J = I-1; A[I,J] = offdy
        if j< n-1: J = I+1; A[I,J] = offdy
    A.assemble()
    return A

def main():
    opts = PETSc.Options()
    n = opts.getInt('n', 32)
    m = opts.getInt('m', 32)
    Print("2-D Laplacian Eigenproblem solved with contour integral, "
          "N=%d (%dx%d grid)\n" % (m*n, m, n))
    A = construct_operator(m,n)

    # Solver object
    E = SLEPc.EPS().create()
    E.setOperators(A)
    E.setType(SLEPc.EPS.Type.CISS)

    # Define region of interest
    R = E.getRG()
    R.setType(SLEPc.RG.Type.ELLIPSE)
    R.setEllipseParameters(0.0,0.2,0.1)
    E.setFromOptions()

    # Compute solution
    E.solve()

    # Print solution
    vw = PETSc.Viewer.STDOUT()
    vw.pushFormat(PETSc.Viewer.Format.ASCII_INFO_DETAIL)
    E.errorView(viewer=vw)
    vw.popFormat()

if __name__ == '__main__':
    main()
