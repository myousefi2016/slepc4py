import sys, os, glob
import unittest

try:
    import slepc4py
except ImportError:
    from distutils.util import get_platform
    plat_specifier = ".%s-%s" % (get_platform(), sys.version[0:3])
    os.path.split(__file__)[0]
    path = os.path.join(os.path.split(__file__)[0], os.path.pardir,
                        'build', 'lib' + plat_specifier)
    sys.path.append(path)
    import slepc4py

args=['-malloc',
      '-malloc_debug',
      '-malloc_dump',
      #'-log_summary',
      ]
if '-slepc' in sys.argv:
    idx = sys.argv.index('-slepc')
    args.extend(sys.argv[idx+1:])
    del sys.argv[idx:]
    del idx

slepc4py.init(args)
from slepc4py import SLEPc

#version = SLEPc.Sys.getVersion()
exclude = {#'test_xxx'      : True,
           #'test_xxx'      : (2,3,3),
           }

def test_cases():
    from glob import glob
    directory = os.path.split(__file__)[0]
    pattern = os.path.join(directory, 'test_*.py')
    test_list = []
    for test_file in glob(pattern):
        filename = os.path.basename(test_file)
        modulename = os.path.splitext(filename)[0]
        if modulename in exclude:
            if exclude[modulename] is True or \
               exclude[modulename] == version:
                continue
        test = __import__(modulename)
        test_list.append(test)
    return test_list

for test in test_cases():
    try:
        unittest.main(test)
    except SystemExit:
        pass
