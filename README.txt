========================
README: SLEPc for Python
========================

:Author:       Lisandro Dalcin
:Organization: CIMEC_
:Address:      PTLC, (3000) Santa Fe, Argentina
:Contact:      dalcinl@gmail.com
:Web Site:     http://slepc4py.googlecode.com/
:Copyright:    This document has been placed in the public domain.

Thank you for downloading the *SLEPc for Python* project archive. As
this is a work in progress, please check the `project website`_ for
updates.  This project should be considered experimental, APIs are
subject to change at any time.

.. _CIMEC:            http://www.cimec.org.ar/
.. _project website:  http://slepc4py.googlecode.com/


- To build and install this package, you must meet the following
  requirements.

  + The latest PETSc_ distribution, built with *shared libraries*.

  + The latest SLEPc_ distribution, built with *shared libraries*.

  + Python_ 2.4/2.5/2.6.

  + NumPy_ 1.0.1 and above.

  + petsc4py_ 1.0.0 and above.

.. _PETSc:    http://www-unix.mcs.anl.gov/petsc/
.. _SLEPc:    http://www.grycap.upv.es/slepc/
.. _Python:   http://www.python.org/
.. _NumPy:    http://numpy.scipy.org/
.. _petsc4py: http://petsc4py.googlecode.com/


- This package uses standard `distutils`. For detailed instructions
  about requirements and the building/install process, read the file
  ``docs/install.txt``.


- The project documentation can be found in files ``docs/*.txt``.  It
  is written reStructuredText_ format. You can use Docutils_ to get
  HTML or LaTeX output. A basic ``Makefile`` is provided in ``docs/``
  directory. 
  
  + Try ``make html`` to obtain HTML output in ``docs/slepc4py.html``.

  + Try ``make pdf``  to obtain PDF output in ``docs/slepc4py.pdf``.

.. _Docutils:         http://docutils.sourceforge.net
.. _reStructuredText: http://docutils.sourceforge.net/rst.html
