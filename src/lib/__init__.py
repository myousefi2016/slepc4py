# Author:    Lisandro Dalcin
# Contact:   dalcinl@gmail.com
# Copyright: This module has been placed in the public domain.

# --------------------------------------------------------------------

"""
Extension modules for different SLEPc configurations.

SLEPc can be configured with different options (eg. debug/optimized,
single/double precisionm, C/C++ compilers, external packages). Each
configuration variant is associated to a name, frequently available as
an environmental variable named ``PETSC_ARCH``.

This package is a holds all the available variants of the SLEPc
extension module built agaist specific SLEPc configurations. It also
provides a convenience function using of the builtin ``imp`` module
for easily importing any of the available extension modules depending
on the value of a user-provided configuration name, the ``PETSC_ARCH``
environmental variable, or a configuration file.
"""

# --------------------------------------------------------------------

from petsc4py.lib import Import

def ImportSLEPc(arch=None):
    """
    Import the SLEPc extension module for a given configuration name.
    """
    return Import('slepc4py', 'SLEPc', __file__, arch,
                  'PETSC_ARCH', 'slepc.cfg')

# --------------------------------------------------------------------
