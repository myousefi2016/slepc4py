.PHONY: default src config build test install \
	uninstall clean distclean srcclean fullclean \
	sdist cython epydoc

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

uninstall:
	-${RM} -r ${HOME}/lib/python/slepc4py
	-${RM} -r ${HOME}/lib/python/slepc4py-*-py*.egg-info

clean:
	${PYTHON} setup.py clean --all
	-${RM} _configtest.* *.py[co]

distclean: clean 
	-${RM} -r build  *.py[co]
	-${RM} -r MANIFEST dist slepc4py.egg-info
	-${RM} `find . -name '*~'`
	-${RM} `find . -name '*.py[co]'`

srcclean:
	-${RM} src/slepc4py_SLEPc.c
	-${RM} src/include/slepc4py/slepc4py_SLEPc.h
	-${RM} src/include/slepc4py/slepc4py_SLEPc_api.h

fullclean: distclean srcclean


sdist:
	${PYTHON} setup.py sdist ${SDISTOPT}

CY_SRC_PXD = $(wildcard src/include/slepc4py/*.pxd)
CY_SRC_PXI = $(wildcard src/SLEPc/*.pxi)
CY_SRC_PYX = $(wildcard src/SLEPc/*.pyx)
src/SLEPc.c: src/slepc4py_SLEPc.c
src/slepc4py_SLEPc.c: ${CY_SRC_PXD} ${CY_SRC_PXI} ${CY_SRC_PYX}
	${PYTHON} ./conf/cythonize.py

cython:
	${PYTHON} ./conf/cythonize.py
epydoc:
	${PYTHON} ./conf/epydocify.py -o /tmp/slepc4py-api-doc
