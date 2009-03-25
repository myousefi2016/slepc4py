ARCH = None
from slepc4py.lib import ImportSLEPc
SLEPc = ImportSLEPc(ARCH)
del ARCH, ImportSLEPc
SLEPc._initialize()
del SLEPc
