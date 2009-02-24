Installation
============


Requirements
------------

* PETSc_ 2.3.3/3.0.0 distribution, built with shared libraries.

* SLEPc_ 2.3.3/3.0.0 distribution, built with shared libraries.

* Python_ 2.4/2.5/2.6.

* NumPy_ package (latest version).

* petsc4py_ package.

.. include:: links.txt


Building
--------

Download and unpack the source distribution::

   $ wget -zxf http://slepc4py.googlecode.com/files/slepc4py-X.X.X.tar.gz
   $ tar -zxf slepc4py-X.X.X.tar.gz
   $ cd slepc4py-X.X.X

Some environmental configuration is needed to inform the location of
PETSc and SLEPc. You can set (using ``setenv``, ``export`` or what
applies to you shell or system) the environmental variables
``SLEPC_DIR``, ``PETSC_DIR``, and ``PETSC_ARCH`` indicating where you
have built/installed SLEPc and PETSc::

   $ export SLEPC_DIR=/usr/local/slepc/3.0.0
   $ export PETSC_DIR=/usr/local/petsc/3.0.0
   $ export PETSC_ARCH=linux-gnu

Alternatively, you can edit the file ``setup.cfg`` and provide the
required information below the ``[config]`` section::

   [config]
   slepc_dir  = /usr/local/slepc/3.0.0
   petsc_dir  = /usr/local/petsc/3.0.0
   petsc_arch = linux-gnu
   ...

Finally, you can build this distribution by typing::

   $ python setup.py build


Installing
----------

After building, this distribution is ready for installation. 

You can do a site-install type::

   $ python setup.py install

or, in case you need root privileges::

   $ su -c 'python setup.py install'

This will install the ``slepc4py`` package in the standard location
``<prefix>/lib/pythonX.X/site-packages``.

You can also do a user-install type::

   $ python setup.py install --home=$HOME

This will install the ``slepc4py`` package in the standard location
``$HOME/lib/python`` (or perhaps ``$HOME/lib64/python``). This
location should be listed in the ``PYTHONPATH`` environmental
variable.
