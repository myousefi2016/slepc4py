# Author:    Lisandro Dalcin
# Contact:   dalcinl@gmail.com
# Copyright: This module has been placed in the public domain.

"""
SLEPc for Python
================

This package is an interface to SLEPc_ libraries.

SLEPc_ (the Scalable Library for Eigenvalue Problem Computations) is a
software library for the solution of large scale sparse eigenvalue
problems on parallel computers. It is an extension of PETSc_ and can
be used for either standard or generalized eigenproblems, with real or
complex arithmetic. It can also be used for computing a partial SVD of
a large, sparse, rectangular matrix.

.. _SLEPc: http://www.grycap.upv.es/slepc
.. _PETSc: http://www.mcs.anl.gov/petsc
"""

# --------------------------------------------------------------------

__author__    = 'Lisandro Dalcin'
__version__   = '1.0.0'
__credits__   = 'SLEPc Team <slepc-maint@grycap.upv.es>'

# --------------------------------------------------------------------

def init(args=None, arch=None):
    """
    Initializes SLEPc.

    :Parameters:
      - `args`: command-line arguments, usually the 'sys.argv' list.
      - `arch`: specific configuration to use.

    .. note:: This function should be called only once, typically at
       the very beginning of the bootstrap script of an application.
    """
    import petsc4py
    petsc4py.init(args, arch)
    import slepc4py.lib
    SLEPc = slepc4py.lib.ImportSLEPc(arch)
    SLEPc._initialize(args)


# --------------------------------------------------------------------

def get_include():
    """
    Return the directory in the package that contains header files.

    Extension modules that need to compile against slepc4py should use
    this function to locate the appropriate include directory. Using
    Python distutils (or perhaps NumPy distutils)::

      import petscc4py, slepc4py
      Extension('extension_name', ...
                include_dirs=[...,
                              petsc4py.get_include(),
                              slepc4py.get_include(),])
    """
    from os.path import dirname, join
    return join(dirname(__file__), 'include')

# --------------------------------------------------------------------

if __name__ == '__main__':
    import slepc4py.help
    slepc4py.help.help()

# --------------------------------------------------------------------
