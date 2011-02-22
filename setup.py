#!/usr/bin/env python

"""
SLEPc for Python
================

Python bindings for SLEPc libraries.
"""


## try:
##     import setuptools
## except ImportError:
##     pass

import sys
sys.dont_write_bytecode = True

# --------------------------------------------------------------------
# Metadata
# --------------------------------------------------------------------

from conf.metadata import metadata

def version():
    import os, re
    fh = open(os.path.join('src', '__init__.py'))
    try: data = fh.read()
    finally: fh.close()
    m = re.search(r"__version__\s*=\s*'(.*)'", data)
    return m.groups()[0]

name     = 'slepc4py'
version  = version()

url      = 'http://%(name)s.googlecode.com/' % vars()
download = url + 'files/%(name)s-%(version)s.tar.gz' % vars()

descr    = __doc__.strip().split('\n'); del descr[1:3]
devstat  = ['Development Status :: 3 - Alpha']
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
from conf.slepcconf import config, build, sdist
from conf.slepcconf import build_src, build_ext

def run_setup():
    import sys, os
    if (('distribute' in sys.modules) or
        ('setuptools' in sys.modules)):
        metadata['install_requires'] = ['petsc4py']
        if not os.environ.get('SLEPC_DIR'):
            metadata['install_requires'].append('slepc')
    if 'setuptools' in sys.modules:
        metadata['zip_safe'] = False
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
                          'sdist'      : sdist,
                          },
          **metadata)

def chk_cython(CYTHON_VERSION_REQUIRED):
    import sys, os
    from distutils.version import StrictVersion as Version
    warn = lambda msg='': sys.stderr.write(msg+'\n')
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
    CYTHON_VERSION = CYTHON_VERSION.split('+', 1)[0]
    for s in ('.alpha', 'alpha'):
        CYTHON_VERSION = CYTHON_VERSION.replace(s, 'a')
    for s in ('.beta',  'beta', '.rc', 'rc', '.c', 'c'):
        CYTHON_VERSION = CYTHON_VERSION.replace(s, 'b')
    if Version(CYTHON_VERSION) < Version(CYTHON_VERSION_REQUIRED):
        warn("*"*80)
        warn()
        warn(" You need to install Cython %s (you have version %s)"
             % (CYTHON_VERSION_REQUIRED, CYTHON_VERSION))
        warn(" Download and install Cython <http://www.cython.org>")
        warn()
        warn("*"*80)
        return False
    #
    return True

def run_cython(source, target, includes=(),
               depends=(), force=False,
               CYTHON_VERSION_REQUIRED=None):
    from distutils import log
    from distutils import dep_util
    from distutils.errors import DistutilsError
    depends = [source] + list(depends)
    if not (force or dep_util.newer_group(depends, target)):
        log.debug("skipping '%s' -> '%s' (up-to-date)",
                  source, target)
        return
    if (CYTHON_VERSION_REQUIRED and not
        chk_cython(CYTHON_VERSION_REQUIRED)):
        raise DistutilsError('requires Cython>=%s'
                             % CYTHON_VERSION_REQUIRED)
    log.info("cythonizing '%s' -> '%s'", source, target)
    from conf.cythonize import run as cythonize
    cythonize(source, includes=includes)

def build_sources(cmd):
    CYTHON_VERSION_REQUIRED = '0.13'
    import os, glob
    if not (os.path.isdir('.hg')  or
            os.path.isdir('.git') or
            cmd.force): return
    source = os.path.join('src', 'slepc4py.SLEPc.pyx')
    target = os.path.splitext(source)[0]+".c"
    depends = (glob.glob("src/include/*/*.pxd") +
               glob.glob("src/*/*.pyx") +
               glob.glob("src/*/*.pxi"))
    import petsc4py
    includes = [petsc4py.get_include()]
    run_cython(source, target, includes,
               depends, cmd.force,
               CYTHON_VERSION_REQUIRED)

build_src.run = build_sources

# --------------------------------------------------------------------

def main():
    run_setup()

if __name__ == '__main__':
    main()

# --------------------------------------------------------------------
