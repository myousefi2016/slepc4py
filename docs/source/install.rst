Installation
============

Requirements
------------

You need to have the following software properly installed in order to
build *SLEPc for Python*:

* Any MPI_ implementation [#]_ (e.g., MPICH_ or `Open MPI`_), 
  built with shared libraries.

* PETSc_ 2.3.3/3.0.0/3.1/3.2 release, built with shared libraries.

* SLEPc_ 2.3.3/3.0.0/3.1/3.2 release, built with shared libraries.

* Python_ 2.4 to 2.7 or 3.1 [#]_.

* NumPy_ package.

* petsc4py_ package.

.. [#] Unless you have appropriately configured and built SLEPc and 
       PETSc without MPI (configure option ``--with-mpi=0``).

.. [#] You may need to use a parallelized version of the Python
       interpreter with some MPI-1 implementations (e.g. MPICH1).

.. include:: links.txt


Using **pip** or **easy_install**
---------------------------------

If you already have a working SLEPc, set environment variables
:envvar:`SLEPC_DIR` and :envvar:`PETSC_DIR` (and perhaps
:envvar:`PETSC_ARCH`) to appropriate values::

    $ export SLEPC_DIR=/path/to/slepc
    $ export PETSC_DIR=/path/to/petsc
    $ export PETSC_ARCH=linux-gnu

.. note:: If you do not set these environment variables, the install
   process will attempt to download and install SLEPc for you.

Now you can use :program:`pip`::

    $ [sudo] pip install [--user] slepc4py

Alternatively, you can use *setuptools* :program:`easy_install`
(deprecated)::

    $ [sudo] easy_install slepc4py


Using **distutils**
-------------------

Downloading
^^^^^^^^^^^

The *SLEPc for Python* package is available for download at the
project website generously hosted by Google Code. You can use
:program:`wget` to get a release tarball::

   $ wget http://slepc4py.googlecode.com/files/slepc4py-X.X.X.tar.gz

Building
^^^^^^^^

After unpacking the release tarball::

   $ tar -zxf slepc4py-X.X.X.tar.gz
   $ cd slepc4py-X.X.X

the distribution is ready for building.

Some environmental configuration is needed to inform the location of
PETSc and SLEPc. You can set (using :command:`setenv`,
:command:`export` or what applies to you shell or system) the
environmental variables :envvar:`SLEPC_DIR``, :envvar:`PETSC_DIR`, and
:envvar:`PETSC_ARCH` indicating where you have built/installed SLEPc
and PETSc::

   $ export SLEPC_DIR=/usr/local/slepc/3.2
   $ export PETSC_DIR=/usr/local/petsc/3.2
   $ export PETSC_ARCH=linux-gnu

Alternatively, you can edit the file :file:`setup.cfg` and provide the
required information below the ``[config]`` section::

   [config]
   slepc_dir  = /usr/local/slepc/3.2
   petsc_dir  = /usr/local/petsc/3.2
   petsc_arch = linux-gnu
   ...

Finally, you can build the distribution by typing::

   $ python setup.py build

Installing
^^^^^^^^^^

After building, the distribution is ready for installation. 

You can do a site-install type::

   $ python setup.py install

or, in case you need root privileges::

   $ su -c 'python setup.py install'

This will install the :mod:`slepc4py` package in the standard location
:file:`{prefix}/lib/python{X}.{X}/site-packages`.

You can also do a user-install type. There are two options depending
on the target Python version.

* For Python 2.6 and up::

      $ python setup.py install --user

* For Python 2.5 and below (assuming your home directory is available
  through the :envvar:`HOME` environment variable)::

      $ python setup.py install --home=$HOME

  and then add :file:`$HOME/lib/python` or :file:`$HOME/lib64/python`
  to your :envvar:`PYTHONPATH` environment variable.
