.PHONY: default \
	src cython \
	config build test install \
	docs sphinx sphinx-html sphinx-pdf epydoc \
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
	-${RM} src/slepc4py.SLEPc.c
	-${RM} src/include/slepc4py/slepc4py.SLEPc.h
	-${RM} src/include/slepc4py/slepc4py.SLEPc_api.h

docsclean:
	-${RM} -r docs/html docs/*.pdf

fullclean: distclean srcclean docsclean

uninstall:
	-${RM} -r ${HOME}/lib/python/slepc4py
	-${RM} -r ${HOME}/lib/python/slepc4py-*-py*.egg-info

CY_SRC_PXD = $(wildcard src/include/slepc4py/*.pxd)
CY_SRC_PXI = $(wildcard src/SLEPc/*.pxi)
CY_SRC_PYX = $(wildcard src/SLEPc/*.pyx)
src/SLEPc.c: src/slepc4py.SLEPc.c
src/slepc4py.SLEPc.c: ${CY_SRC_PXD} ${CY_SRC_PXI} ${CY_SRC_PYX}
	${PYTHON} ./conf/cythonize.py

cython:
	${PYTHON} ./conf/cythonize.py

SPHINXBUILD = sphinx-build
SPHINXOPTS  =
sphinx: sphinx-html sphinx-pdf
sphinx-html:
	${PYTHON} -c 'import slepc4py.SLEPc'
	mkdir -p build/doctrees docs/html/man
	${SPHINXBUILD} -b html -d build/doctrees ${SPHINXOPTS} \
	docs/source docs/html/man
sphinx-pdf:
	${PYTHON} -c 'import slepc4py.SLEPc'
	mkdir -p build/doctrees build/latex
	${SPHINXBUILD} -b latex -d build/doctrees ${SPHINXOPTS} \
	docs/source build/latex
	${MAKE} -C build/latex all-pdf > /dev/null
	mv build/latex/*.pdf docs/

EPYDOCBUILD = ${PYTHON} ./conf/epydocify.py
EPYDOCOPTS  =
epydoc:
	${PYTHON} -c 'import slepc4py.SLEPc'
	mkdir -p docs/html/api
	${EPYDOCBUILD} ${EPYDOCOPTS} --html -o docs/html/api 


sdist: src docs
	${PYTHON} setup.py sdist ${SDISTOPT}
