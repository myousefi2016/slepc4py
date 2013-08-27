#!/usr/bin/env python

"""
SLEPc for Python
================

Python bindings for SLEPc libraries.

.. tip::

  You can also install `slepc4py-dev`_ with::

    $ pip install petsc==dev slepc==dev petsc4py==dev slepc4py==dev 

  .. _slepc4py-dev: hg+https://bitbucket.org/slepc/slepc4py.git#egg=slepc4py-dev
"""

## try:
##     import setuptools
## except ImportError:
##     pass

import sys, os

# --------------------------------------------------------------------
# Metadata
# --------------------------------------------------------------------

from conf.metadata import metadata

def name():
    return 'slepc4py'

def version():
    import re
    fh = open(os.path.join('src', '__init__.py'))
    try: data = fh.read()
    finally: fh.close()
    m = re.search(r"__version__\s*=\s*'(.*)'", data)
    return m.groups()[0]

name     = name()
version  = version()

url      = 'https://bitbucket.org/slepc/%(name)s/' % vars()
download = url + '%(name)s-%(version)s.tar.gz' % vars()

descr    = __doc__.strip().split('\n'); del descr[1:3]
devstat  = ['Development Status :: 5 - Production/Stable']
keywords = ['SLEPc', 'PETSc', 'MPI']

metadata['name'] = name
metadata['version'] = version
metadata['description'] = descr.pop(0)
metadata['long_description'] = '\n'.join(descr)
metadata['keywords'] += keywords
metadata['classifiers'] += devstat
metadata['url'] = url
metadata['download_url'] = download

metadata['provides'] = ['slepc4py']
metadata['requires'] = ['petsc4py']

# --------------------------------------------------------------------
# Extension modules
# --------------------------------------------------------------------

def get_ext_modules(Extension):
    from os   import walk, path
    from glob import glob
    depends = []
    for pth, dirs, files in walk('src'):
        depends += glob(path.join(pth, '*.h'))
        depends += glob(path.join(pth, '*.c'))
    try:
        import petsc4py
        petsc4py_includes = [petsc4py.get_include()]
    except ImportError:
        petsc4py_includes = []
    return [Extension('slepc4py.lib.SLEPc',
                      sources=['src/SLEPc.c',],
                      include_dirs=['src/include',
                                    ] + petsc4py_includes,
                      depends=depends)]

# --------------------------------------------------------------------
# Setup
# --------------------------------------------------------------------

from conf.slepcconf import setup, Extension
from conf.slepcconf import config, build, build_src, build_ext
from conf.slepcconf import test, sdist

CYTHON = '0.15'

def run_setup():
    if ('setuptools' in sys.modules):
        from os.path import exists, join
        metadata['zip_safe'] = False
        metadata['install_requires'] = ['petsc4py>=3.3.1']
        if not exists(join('src', 'slepc4py.SLEPc.c')):
            metadata['install_requires'] += ['Cython>='+CYTHON]
        SLEPC_DIR = os.environ.get('SLEPC_DIR')
        if not (SLEPC_DIR and os.path.isdir(SLEPC_DIR)):
            metadata['install_requires'].append('slepc>=3.3,<3.4')
    #
    setup(packages     = ['slepc4py',
                          'slepc4py.lib',],
          package_dir  = {'slepc4py'     : 'src',
                          'slepc4py.lib' : 'src/lib'},
          package_data = {'slepc4py'     : ['include/slepc4py/*.h',
                                            'include/slepc4py/*.i',
                                            'include/slepc4py/*.pxd',
                                            'include/slepc4py/*.pxi',
                                            'include/slepc4py/*.pyx',],
                          'slepc4py.lib' : ['slepc.cfg'],},
          ext_modules  = get_ext_modules(Extension),
          cmdclass     = {'config'     : config,
                          'build'      : build,
                          'build_src'  : build_src,
                          'build_ext'  : build_ext,
                          'test'       : test,
                          'sdist'      : sdist,
                          },
          **metadata)

def chk_cython(VERSION):
    import re
    from distutils import log
    from distutils.version import LooseVersion
    from distutils.version import StrictVersion
    warn = lambda msg='': sys.stderr.write(msg+'\n')
    #
    cython_zip = 'cython.zip'
    if os.path.isfile(cython_zip):
        path = os.path.abspath(cython_zip)
        if sys.path[0] != path:
            sys.path.insert(0, path)
            log.info("adding '%s' to sys.path", cython_zip)
    #
    try:
        import Cython
    except ImportError:
        warn("*"*80)
        warn()
        warn(" You need to generate C source files with Cython!!")
        warn(" Download and install Cython <http://www.cython.org>")
        warn()
        warn("*"*80)
        return False
    #
    try:
        CYTHON_VERSION = Cython.__version__
    except AttributeError:
        from Cython.Compiler.Version import version as CYTHON_VERSION
    REQUIRED = VERSION
    m = re.match(r"(\d+\.\d+(?:\.\d+)?).*", CYTHON_VERSION)
    if m:
        Version = StrictVersion
        AVAILABLE = m.groups()[0]
    else:
        Version = LooseVersion
        AVAILABLE = CYTHON_VERSION
    if (REQUIRED is not None and
        Version(AVAILABLE) < Version(REQUIRED)):
        warn("*"*80)
        warn()
        warn(" You need to install Cython %s (you have version %s)"
             % (REQUIRED, CYTHON_VERSION))
        warn(" Download and install Cython <http://www.cython.org>")
        warn()
        warn("*"*80)
        return False
    #
    return True

def run_cython(source, depends=(), includes=(),
               destdir_c=None, destdir_h=None, wdir=None,
               force=False, VERSION=None):
    from glob import glob
    from distutils import log
    from distutils import dep_util
    from distutils.errors import DistutilsError
    target = os.path.splitext(source)[0]+".c"
    cwd = os.getcwd()
    try:
        if wdir: os.chdir(wdir)
        alldeps = [source]
        for dep in depends:
            alldeps += glob(dep)
        if not (force or dep_util.newer_group(alldeps, target)):
            log.debug("skipping '%s' -> '%s' (up-to-date)",
                      source, target)
            return
    finally:
        os.chdir(cwd)
    if not chk_cython(VERSION):
        raise DistutilsError("requires Cython>=%s" % VERSION)
    log.info("cythonizing '%s' -> '%s'", source, target)
    from conf.cythonize import cythonize
    err = cythonize(source,
                    includes=includes,
                    destdir_c=destdir_c,
                    destdir_h=destdir_h,
                    wdir=wdir)
    if err:
        raise DistutilsError(
            "Cython failure: '%s' -> '%s'" % (source, target))

def build_sources(cmd):
    from os.path import exists, isdir, join
    if (exists(join('src', 'slepc4py.SLEPc.c')) and
        not (isdir('.hg') or isdir('.git')) and
        not cmd.force): return
    # slepc4py.SLEPc
    source = 'slepc4py.SLEPc.pyx'
    depends = ("include/*/*.pxd",
               "SLEPc/*.pyx",
               "SLEPc/*.pxi",)
    import petsc4py
    includes = ['include', petsc4py.get_include()]
    destdir_h = os.path.join('include', 'slepc4py')
    run_cython(source, depends, includes,
               destdir_c=None, destdir_h=destdir_h, wdir='src',
               force=cmd.force, VERSION=CYTHON)

build_src.run = build_sources

def run_testsuite(cmd):
    from distutils.errors import DistutilsError
    sys.path.insert(0, 'test')
    try:
        from runtests import main
    finally:
        del sys.path[0]
    err = main(cmd.args or [])
    if err:
        raise DistutilsError("test")

test.run = run_testsuite

# --------------------------------------------------------------------

def main():
    run_setup()

if __name__ == '__main__':
    main()

# --------------------------------------------------------------------
