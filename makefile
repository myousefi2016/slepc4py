.PHONY: default \
	config build test install sdist \
	docs rst2html sphinx sphinx-html sphinx-pdf epydoc epydoc-html \
	clean distclean srcclean docsclean fullclean uninstall

PYTHON = python

default: build

config: 
	${PYTHON} setup.py config ${CONFIGOPT}

build:
	${PYTHON} setup.py build ${BUILDOPT}

test:
	${MPIEXEC} ${VALGRIND} ${PYTHON} test/runtests.py < /dev/null

install: build
	${PYTHON} setup.py install ${INSTALLOPT} --home=${HOME}


sdist: docs
	${PYTHON} setup.py sdist ${SDISTOPT}

clean:
	${PYTHON} setup.py clean --all

distclean: clean
	-${RM} -r build  _configtest.* *.py[co]
	-${RM} -r MANIFEST dist slepc4py.egg-info
	-${RM} -r `find . -name '__pycache__'`
	-${RM} `find . -name '*.py[co]'`
	-${RM} `find . -name '*~'`

srcclean:
	-${RM} src/slepc4py.SLEPc.c
	-${RM} src/include/slepc4py/slepc4py.SLEPc.h
	-${RM} src/include/slepc4py/slepc4py.SLEPc_api.h

fullclean: distclean srcclean docsclean

uninstall:
	-${RM} -r ${HOME}/lib/python/slepc4py
	-${RM} -r ${HOME}/lib/python/slepc4py-*-py*.egg-info

# ----

docs: rst2html sphinx epydoc

docsclean:
	-${RM} docs/*.html docs/*.pdf
	-${RM} -r docs/usrman docs/apiref

RST2HTML = rst2html
RST2HTMLOPTS = --no-compact-lists --cloak-email-addresses
rst2html:
#	${RST2HTML} ${RST2HTMLOPTS} docs/index.rst > docs/index.html
	${RST2HTML} ${RST2HTMLOPTS} LICENSE.txt    > docs/LICENSE.html
	${RST2HTML} ${RST2HTMLOPTS} HISTORY.txt    > docs/HISTORY.html

SPHINXBUILD = sphinx-build
SPHINXOPTS  =
sphinx: sphinx-html sphinx-pdf
sphinx-html:
	${PYTHON} -c 'import slepc4py.SLEPc'
	mkdir -p build/doctrees docs/usrman
	${SPHINXBUILD} -b html -d build/doctrees ${SPHINXOPTS} \
	docs/source docs/usrman
sphinx-pdf:
	${PYTHON} -c 'import slepc4py.SLEPc'
	mkdir -p build/doctrees build/latex
	${SPHINXBUILD} -b latex -d build/doctrees ${SPHINXOPTS} \
	docs/source build/latex
	${MAKE} -C build/latex all-pdf > /dev/null
	mv build/latex/*.pdf docs/

EPYDOCBUILD = ${PYTHON} ./conf/epydocify.py
EPYDOCOPTS  =
epydoc: epydoc-html
epydoc-html:
	${PYTHON} -c 'import slepc4py.SLEPc'
	mkdir -p docs/apiref
	${EPYDOCBUILD} ${EPYDOCOPTS} -o docs/apiref

# ----
