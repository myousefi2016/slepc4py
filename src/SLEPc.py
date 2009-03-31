ARCH = None
from slepc4py.lib import ImportSLEPc
from slepc4py.lib import ImportPETSc
SLEPc = ImportSLEPc(ARCH)
PETSc = ImportPETSc(ARCH)
PETSc._initialize()
SLEPc._initialize()
del SLEPc, PETSc
del ImportSLEPc, ImportPETSc
del ARCH
