Overview
========

*SLEPc for Python* (slepc4py) is a Python package that provides
convenient access to the functionality of SLEPc.

SLEPc [1]_, [2]_ implements algorithms and tools for the numerical
solution of large, sparse eigenvalue problems on parallel
computers. It covers both standard and generalized eigenproblems
(either symmetric or non-symmetric) as well as the singular value
decomposition (SVD) and the quadratic eigenvalue problem (QEP).

SLEPc is intended for computing a subset of the spectrum of a matrix
(or matrix pair). One can for instance approximate the largest
magnitude eigenvalues, or the smallest ones, or even those eigenvalues
located near a given region of the complex plane. Interior eigenvalues
are harder to compute, so SLEPc provides different methodologies. One
such method is to use a spectral transformation. Cheaper alternatives
are also available.

.. [1] V. Hernandez, J. E. Roman, E. Romero, A. Tomas and
   V. Vidal. SLEPc Users Manual. DISC-II/24/02 - Revision 3.1
   D. Sistemas Informáticos y Computación, Universidad Politécnica de
   Valencia. 2010.

.. [2] Vicente Hernandez, Jose E. Roman and Vicente Vidal.
   SLEPc: A Scalable and Flexible Toolkit for the Solution of
   Eigenvalue Problems, ACM Trans. Math. Softw. 31(3), pp. 351-362,
   2005.

.. include:: links.txt


Features
--------

Currently, the following types of eigenproblems can be addressed:

* Standard eigenvalue problem, *Ax=kx*, either for Hermitian or
  non-Hermitian matrices.

* Generalized eigenvalue problem, *Ax=kBx*, either Hermitian
  positive-definite or not.

* Partial singular value decomposition of a rectangular matrix,
  *Au=sv*.

* Quadratic eigenvalue problem, *(k^2M+kC+K)x=0*.

For the eigenvalue problem, the following methods are available:

* Krylov eigensolvers, particularly Krylov-Schur, Arnoldi, and
  Lanczos.

* Davidson-type eigensolvers, including Generalized Davidson and
  Jacobi-Davidson.

* Subspace iteration and single vector iterations (inverse iteration,
  RQI).

For singular value computations, the following alternatives can be
used:

* Use an eigensolver via the cross-product matrix *A'A* or the cyclic
  matrix *[0 A; A' 0]*.

* Explicitly restarted Lanczos bidiagonalization.

* Implicitly restarted Lanczos bidiagonalization (thick-restart
  Lanczos).

For quadratic eigenvalue problems, the following methods are available:

* Use an eigensolver to solve the generalized eigenvalue problem 
  obtained after linearization.

* Q-Arnoldi, a memory efficient variant of Arnoldi for quadratic problems.

Computation of interior eigenvalues is supported by means of the
following methodologies:

* Spectral transformations, such as shift-and-invert. This technique
  implicitly uses the inverse of the shifted matrix *(A-tI)* in order
  to compute eigenvalues closest to a given target value, *t*.

* Harmonic extraction, a cheap alternative to shift-and-invert that
  also tries to approximate eigenvalues closest to a target, *t*, but
  without requiring a matrix inversion.

Other remarkable features include:

* High computational efficiency, by using NumPy and SLEPc under the
  hood.

* Data-structure neutral implementation, by using efficient sparse
  matrix storage provided by PETSc. Implicit matrix representation is
  also available by providing basic operations such as matrix-vector
  products as user-defined Python functions.

* Run-time flexibility, by specifying numerous setting at the command
  line.

* Ability to do the computation in parallel.


Components
----------

SLEPc provides the following components, which are mirrored by
slepc4py for its use from Python.

:EPS: The Eigenvalue Problem Solver is the component that provides all
      the functionality necessary to define and solve an
      eigenproblem. It provides mechanisms for completely specifying
      the problem: the problem type (e.g. standard symmetric), number
      of eigenvalues to compute, part of the spectrum of
      interest. Once the problem has been defined, a collection of
      solvers can be used to compute the required solutions.  The
      behaviour of the solvers can be tuned by means of a few
      parameters, such as the maximum dimension of the subspace to be
      used during the computation.

:SVD: This component is the analog of EPS for the case of Singular
      Value Decompositions. The user provides a rectangular matrix and
      specifies how many singular values and vectors are to be
      computed, whether the largest or smallest ones, as well as some
      other parameters for fine tuning the computation. Different
      solvers are available, as in the case of EPS.

:QEP: This component is the analog of EPS for the case of Quadratic
      Eigenvalue Problems. The user provides three square matrices that
      define the problem. Several parameters can be specified, as in
      the case of EPS. It is also possible to indicate whether the
      problem belongs to a special type, e.g., symmetric or gyroscopic.

:ST:  The Spectral Transformation is a component that provides
      convenient implementations of common spectral
      transformations. These are simple transformations that map
      eigenvalues to different positions, in such a way that
      convergence to wanted eigenvalues is enhanced. The most common
      spectral transformation is shift-and-invert, that allows for the
      computation of eigenvalues closest to a given target value.

:IP:  This component encapsulates the concept of an Inner Product in a
      vector space, which can be either the standard Hermitian inner
      product *x'y* or the positive definite product *x'By* for a
      given SPD matrix *B*. This component provides convenient access
      to common operations such as orthogonalization of vectors. The
      IP component is usually not required by end-users.
