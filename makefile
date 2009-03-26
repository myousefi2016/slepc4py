.PHONY: default \
	src cython \
	config build test install \
	docs sphinx epydoc \
	sdist \
	clean distclean srcclean docsclean fullclean uninstall

PYTHON = python

default: build

src: src/SLEPc.c

config: 
	${PYTHON} setup.py config ${CONFIGOPT}

build: src
	${PYTHON} setup.py build ${BUILDOPT}

test:
	${MPIEXEC} ${VALGRIND} ${PYTHON} test/runtests.py < /dev/null

install: build
	${PYTHON} setup.py install ${INSTALLOPT} --home=${HOME}

docs: sphinx epydoc

clean:
	${PYTHON} setup.py clean --all

distclean: clean docsclean
	-${RM} -r build  _configtest.* *.py[co]
	-${RM} -r MANIFEST dist slepc4py.egg-info
	-${RM} `find . -name '*~'`
	-${RM} `find . -name '*.py[co]'`

srcclean:
	-${RM} src/slepc4py_SLEPc.c
	-${RM} src/include/slepc4py/slepc4py_SLEPc.h
	-${RM} src/include/slepc4py/slepc4py_SLEPc_api.h

docsclean:
	-${RM} -r docs/html

fullclean: distclean srcclean docsclean

uninstall:
	-${RM} -r ${HOME}/lib/python/slepc4py
	-${RM} -r ${HOME}/lib/python/slepc4py-*-py*.egg-info

CY_SRC_PXD = $(wildcard src/include/slepc4py/*.pxd)
CY_SRC_PXI = $(wildcard src/SLEPc/*.pxi)
CY_SRC_PYX = $(wildcard src/SLEPc/*.pyx)
src/SLEPc.c: src/slepc4py_SLEPc.c
src/slepc4py_SLEPc.c: ${CY_SRC_PXD} ${CY_SRC_PXI} ${CY_SRC_PYX}
	${PYTHON} ./conf/cythonize.py

cython:
	${PYTHON} ./conf/cythonize.py

SPHINXBUILD = sphinx-build
SPHINXOPTS  =
sphinx:
	mkdir -p build/doctrees docs/html/man
	${SPHINXBUILD} -d build/doctrees ${SPHINXOPTS} \
	docs/source docs/html/man

EPYDOCBUILD = ${PYTHON} ./conf/epydocify.py
EPYDOCOPTS  =
epydoc:
	${PYTHON} -c 'import petsc4py.PETSc'
	mkdir -p docs/html/api
	${EPYDOCBUILD} ${EPYDOCOPTS} --html -o docs/html/api 


sdist: src docs
	${PYTHON} setup.py sdist ${SDISTOPT}
