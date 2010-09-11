# --------------------------------------------------------------------

__all__ = ['setup',
           'Extension',
           'config',
           'build',
           'build_ext',
           ]

# --------------------------------------------------------------------

import sys, os, platform

if not hasattr(sys, 'version_info') or \
       sys.version_info < (2, 4, 0, 'final'):
    raise SystemExit("Python 2.4 or later is required")

# --------------------------------------------------------------------

from conf.core import PetscConfig
from conf.core import setup, Extension, log
from conf.core import config     as _config
from conf.core import build      as _build
from conf.core import build_ext  as _build_ext

from distutils.errors import DistutilsError

# --------------------------------------------------------------------


class SlepcConfig(PetscConfig):

    def __init__(self,  slepc_dir, petsc_dir, petsc_arch):
        PetscConfig.__init__(self, petsc_dir, petsc_arch)
        if not slepc_dir:
            raise DistutilsError("SLEPc not found")
        elif not os.path.isdir(slepc_dir):
            raise DistutilsError("invalid SLEPC_DIR")
        self.configdict['SLEPC_DIR'] = slepc_dir
        self.SLEPC_DIR = self['SLEPC_DIR']

    def configure_extension(self, extension):
        PetscConfig.configure_extension(self, extension)
        SLEPC_DIR  = self.SLEPC_DIR
        PETSC_ARCH = self.PETSC_ARCH
        # define macros
        macros = [('SLEPC_DIR', SLEPC_DIR)]
        extension.define_macros.extend(macros)
        # includes and libraries
        if (os.path.exists(os.path.join(SLEPC_DIR, 'conf')) or
            os.path.exists(os.path.join(SLEPC_DIR, PETSC_ARCH, 'conf'))):
            SLEPC_INCLUDE = [
                os.path.join(SLEPC_DIR, PETSC_ARCH, 'include'),
                os.path.join(SLEPC_DIR, 'include'),
                ]
            SLEPC_LIB_DIR = [
                os.path.join(SLEPC_DIR, PETSC_ARCH, 'lib'),
                os.path.join(SLEPC_DIR, 'lib'),
                ]
        else:
            SLEPC_INCLUDE = [os.path.join(SLEPC_DIR, 'include'), SLEPC_DIR]
            SLEPC_LIB_DIR = [os.path.join(SLEPC_DIR, 'lib', PETSC_ARCH)]
        slepc_cfg = { }
        slepc_cfg['include_dirs'] = SLEPC_INCLUDE
        slepc_cfg['library_dirs'] = SLEPC_LIB_DIR
        slepc_cfg['libraries']    = ['slepc']
        slepc_cfg['runtime_library_dirs'] = slepc_cfg['library_dirs']
        self._configure_ext(extension, slepc_cfg, preppend=True)
        if self['BUILDSHAREDLIB'] == 'no':
            from petsc4py.lib import ImportPETSc
            PETSc = ImportPETSc(PETSC_ARCH)
            extension.extra_objects.append(PETSc.__file__)

        # extra configuration
        cflags = []
        extension.extra_compile_args.extend(cflags)
        lflags = []
        extension.extra_link_args.extend(lflags)

        slepcqep_h = os.path.join(SLEPC_DIR, 'include', 'slepcqep.h')
        if not os.path.exists(slepcqep_h):
            extension.include_dirs.append("src/include/compat")

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

    Configure = SlepcConfig

    user_options = _config.user_options + cmd_slepc_opts

    def initialize_options(self):
        _config.initialize_options(self)
        self.slepc_dir  = None

    def get_config_arch(self, arch):
        return config.Configure(self.slepc_dir, self.petsc_dir, arch)

    def run(self):
        self.slepc_dir = config.get_slepc_dir(self.slepc_dir)
        if self.slepc_dir is None: return
        log.info('-' * 70)
        log.info('SLEPC_DIR:   %s' % self.slepc_dir)
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


class build_ext(_build_ext):

    user_options = _build_ext.user_options + cmd_slepc_opts

    def initialize_options(self):
        _build_ext.initialize_options(self)
        self.slepc_dir  = None

    def finalize_options(self):
        _build_ext.finalize_options(self)
        self.set_undefined_options('build',
                                   ('slepc_dir',  'slepc_dir'))

    def get_config_arch(self, arch):
        return config.Configure(self.slepc_dir, self.petsc_dir, arch)

    def get_config_data(self, arch_list):
        template = """\
SLEPC_DIR  = %(SLEPC_DIR)s
PETSC_DIR  = %(PETSC_DIR)s
PETSC_ARCH = %(PETSC_ARCH)s
"""
        variables = {'SLEPC_DIR'  : self.slepc_dir,
                     'PETSC_DIR'  : self.petsc_dir,
                     'PETSC_ARCH' : os.path.pathsep.join(arch_list)}
        return template, variables

# --------------------------------------------------------------------
