#!/usr/bin/env python

"""
SLEPc for Python
==============

Python bindings for SLEPc libraries.
"""


## try:
##     import setuptools
## except ImportError:
##     pass


# --------------------------------------------------------------------
# Metadata
# --------------------------------------------------------------------

from conf.metadata import metadata

name     = 'slepc4py'
version  = open('VERSION.txt').read().strip()
descr    = __doc__.strip().split('\n'); del descr[1:3]
keywords = ['SLEPc','PETSc','MPI']
devstat  = ['Development Status :: 3 - Alpha']
url      = 'http://%s.googlecode.com/' % name
download = 'http://%s.googlecode.com/files/%s-%s.tar.gz'
download = download % (name, name, version)

metadata['name'] = name
metadata['version'] = version
metadata['description'] = descr.pop(0)
metadata['long_description'] = '\n'.join(descr)
metadata['keywords'] += keywords
metadata['classifiers'] += devstat
metadata['url'] = url
metadata['download_url'] = download

# --------------------------------------------------------------------
# Extension modules
# --------------------------------------------------------------------


def get_ext_modules(Extension):
    from os   import walk, path
    from glob import glob
    depends = []
    for pth, dirs, files in walk('src'):
        depends += glob(path.join(pth, '*.h'))
    import petsc4py
    petsc4py_inc = petsc4py.get_include()
    return [Extension('slepc4py.lib.SLEPc',
                      sources=['src/SLEPc.c',],
                      include_dirs=['src/include',
                                    petsc4py_inc,],
                      depends=depends, language='c'),
            Extension('slepc4py.lib.SLEPc',
                      sources=['src/SLEPc.cpp',],
                      include_dirs=['src/include',
                                    petsc4py_inc,],
                      depends=depends, language='c++'),
            ]


# --------------------------------------------------------------------
# Setup
# --------------------------------------------------------------------

from conf.slepcconf import setup, Extension
from conf.slepcconf import config, build, build_py, build_ext

def main():
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
                          'build_py'   : build_py,
                          'build_ext'  : build_ext},
          **metadata)

# --------------------------------------------------------------------

if __name__ == '__main__':
    import sys, os
    C_SOURCE = os.path.join('src', 'slepc4py_SLEPc.c')
    def cython_help():
        if os.path.exists(C_SOURCE): return
        warn = lambda msg='': sys.stderr.write(msg+'\n')
        warn("*"*70)
        warn()
        warn("You need to generate C source files with Cython !!!")
        warn("Please execute in your shell:")
        warn()
        warn("$ python ./conf/cythonize.py")
        warn()
        warn("*"*70)
        warn()
    ## from distutils import log
    ## log.set_verbosity(log.DEBUG)
    cython_help()
    main()

# --------------------------------------------------------------------
