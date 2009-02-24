================
SLEPc for Python
================

:Authors:      Lisandro Dalcín, Jose E. Roman
:Contact:      dalcinl@gmail.com, jroman@dsic.upv.es
:Organization: CIMEC_, GRyCAP_
:Address:      PTLC, (3000) Santa Fe, Argentina
:Web Site:     http://slepc4py.googlecode.com
:Date:         |today|
:Copyright:    This document has been placed in the public domain.

This document describes slepc4py_, a Python_ port to the SLEPc_
libraries.

SLEPc_ is a software package for the parallel solution of large-scale
eigenvalue problems. It can be used for computing eigenvalues and
eigenvectors of large, sparse matrices, or matrix pairs, and also for
computing singular values and vectors of a rectangular matrix.

SLEPc_ relies on PETSc_ for basic functionality such as the
representation of matrices and vectors, and the solution of linear
systems of equations. Thus, slepc4py_ must be used together with its
companion petsc4py_.


Contents
========

.. toctree::
   :maxdepth: 2

   overview
   tutorial
   install


.. include:: links.txt
