# --------------------------------------------------------------------

__all__ = ['setup',
           'Extension',
           'config',
           'build',
           'build_py',
           'build_ext',
           ]

# --------------------------------------------------------------------

import sys, os, platform

if not hasattr(sys, 'version_info') or \
       sys.version_info < (2, 4, 0,'final'):
    raise SystemExit("Python 2.4 or later is required "
                     "to build SLEPc for Python package.")

# --------------------------------------------------------------------

from core import PetscConfig
from core import setup, Extension, log
from core import config     as _config
from core import build      as _build
from core import build_py   as _build_py
from core import build_ext  as _build_ext

from distutils.errors import DistutilsError

# --------------------------------------------------------------------


class SlepcConfig(PetscConfig):

    def __init__(self, petsc_dir, petsc_arch, slepc_dir):
        PetscConfig.__init__(self, petsc_dir, petsc_arch)
        if not slepc_dir:
            raise DistutilsError("SLEPc not found")
        elif not os.path.isdir(slepc_dir):
            raise DistutilsError("invalid SLEPC_DIR")
        self.configdict['SLEPC_DIR'] = slepc_dir
        self.SLEPC_DIR = self['SLEPC_DIR']

    def configure_extension(self, extension):
        PetscConfig.configure_extension(self, extension)
        # define macros
        macros = [('SLEPC_DIR', self.SLEPC_DIR)]
        extension.define_macros.extend(macros)
        # includes and libraries
        if (os.path.exists(os.path.join(self.SLEPC_DIR, 'conf')) or
            os.path.exists(os.path.join(self.SLEPC_DIR, self.PETSC_ARCH, 'conf'))):
            SLEPC_INCLUDE = [
                os.path.join(self.SLEPC_DIR, self.PETSC_ARCH, 'include'),
                os.path.join(self.SLEPC_DIR, 'include'),
                ]
            SLEPC_LIB_DIR = [
                os.path.join(self.SLEPC_DIR, self.PETSC_ARCH, 'lib'),
                os.path.join(self.SLEPC_DIR, 'lib'),
                ]
        else:
            SLEPC_INCLUDE = [os.path.join(self.SLEPC_DIR, 'include'), self.SLEPC_DIR]
            SLEPC_LIB_DIR = [os.path.join(self.SLEPC_DIR, 'lib', self.PETSC_ARCH)]
        slepc_cfg = { }
        slepc_cfg['include_dirs'] = SLEPC_INCLUDE
        slepc_cfg['library_dirs'] = SLEPC_LIB_DIR
        slepc_cfg['libraries']    = ['slepc']
        slepc_cfg['runtime_library_dirs'] = slepc_cfg['library_dirs']
        self._configure_ext(extension, slepc_cfg)
        # extra configuration
        cflags = []
        extension.extra_compile_args.extend(cflags)
        lflags = []
        extension.extra_link_args.extend(lflags)

    def log_info(self):
        if not self.SLEPC_DIR: return
        log.info('SLEPC_DIR:   %s' % self.SLEPC_DIR)
        PetscConfig.log_info(self)


# --------------------------------------------------------------------

cmd_slepc_opts = [
    ('slepc-dir=', None,
     "define SLEPC_DIR, overriding environmental variable.")
    ]

class config(_config):

    user_options = _config.user_options + cmd_slepc_opts

    def initialize_options(self):
        _config.initialize_options(self)
        self.slepc_dir  = None

    def run(self):
        slepc_dir = config.get_slepc_dir(self.slepc_dir)
        if slepc_dir is None: return
        log.info('-' * 70)
        log.info('SLEPC_DIR:   %s' % slepc_dir)
        _config.run(self)

    @staticmethod
    def get_slepc_dir(slepc_dir):
        if not slepc_dir:
            log.warn("SLEPC_DIR not specified")
            return None
        slepc_dir = os.path.expandvars(slepc_dir)
        if '$SLEPC_DIR' in slepc_dir:
            log.warn("SLEPC_DIR not specified")
            return None
        slepc_dir = os.path.expanduser(slepc_dir)
        slepc_dir = os.path.abspath(slepc_dir)
        if not os.path.isdir(slepc_dir):
            log.warn('invalid SLEPC_DIR:  %s' % slepc_dir)
            return None
        return slepc_dir


class build(_build):

    user_options = _build.user_options + cmd_slepc_opts

    def initialize_options(self):
        _build.initialize_options(self)
        self.slepc_dir  = None

    def finalize_options(self):
        _build.finalize_options(self)
        self.set_undefined_options('config',
                                   ('slepc_dir', 'slepc_dir'),)
        self.slepc_dir = config.get_slepc_dir(self.slepc_dir)


class build_py(_build_py):

    config_file = 'slepc.cfg'

    def initialize_options(self):
        _build_py.initialize_options(self)
        self.slepc_dir = None

    def finalize_options(self):
        _build_py.finalize_options(self)
        self.set_undefined_options('build',
                                   ('slepc_dir',  'slepc_dir'),)

    def _config(self, py_file):
        SLEPC_DIR  = '$SLEPC_DIR'
        PETSC_DIR  = '$PETSC_DIR'
        PETSC_ARCH = '$PETSC_ARCH'
        config_py = open(py_file, 'r')
        config_data = config_py.read()
        config_py.close()
        #
        slepc_dir  = self.slepc_dir
        petsc_dir  = self.petsc_dir
        petsc_arch = self.petsc_arch
        pathsep    = os.path.pathsep
        #
        bmake_dir = os.path.join(petsc_dir, 'bmake')
        have_bmake = os.path.isdir(bmake_dir)
        #
        if '%(SLEPC_DIR)s' not in config_data:
            return # already configured
        if not slepc_dir:
            return # nothing known to put
        #
        if slepc_dir:
            SLEPC_DIR  = slepc_dir
        if petsc_dir:
            PETSC_DIR  = petsc_dir
        if petsc_arch:
            PETSC_ARCH = pathsep.join(petsc_arch)
        elif not have_bmake:
            PETSC_ARCH = 'default'
        log.info('writing %s' % py_file)
        config_py = open(py_file, 'w')
        config_py.write(config_data % vars())
        config_py.close()


class build_ext(_build_ext):

    def initialize_options(self):
        _build_ext.initialize_options(self)
        self.slepc_dir  = None

    def finalize_options(self):
        _build_ext.finalize_options(self)
        self.set_undefined_options('build',
                                   ('slepc_dir',  'slepc_dir'))

    def _get_config(self, petsc_dir, petsc_arch):
        return SlepcConfig(petsc_dir, petsc_arch, self.slepc_dir)


# --------------------------------------------------------------------
