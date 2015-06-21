# --------------------------------------------------------------------

__all__ = ['setup',
           'Extension',
           'config',
           'build',
           'build_src',
           'build_ext',
           'install',
           'clean',
           'test',
           'sdist',
           ]

# --------------------------------------------------------------------

import sys, os

from conf.baseconf import PetscConfig
from conf.baseconf import setup, Extension, log
from conf.baseconf import config     as _config
from conf.baseconf import build      as _build
from conf.baseconf import build_src  as _build_src
from conf.baseconf import build_ext  as _build_ext
from conf.baseconf import install    as _install
from conf.baseconf import clean      as _clean
from conf.baseconf import test       as _test
from conf.baseconf import sdist      as _sdist

from distutils.errors import DistutilsError

# --------------------------------------------------------------------


class SlepcConfig(PetscConfig):

    def __init__(self,  slepc_dir, petsc_dir, petsc_arch):
        PetscConfig.__init__(self, petsc_dir, petsc_arch)
        if not slepc_dir:
            raise DistutilsError("SLEPc not found")
        if not os.path.isdir(slepc_dir):
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
        SLEPC_INCLUDE = [
            os.path.join(SLEPC_DIR, PETSC_ARCH, 'include'),
            os.path.join(SLEPC_DIR, 'include'),
            ]
        SLEPC_LIB_DIR = [
            os.path.join(SLEPC_DIR, PETSC_ARCH, 'lib'),
            os.path.join(SLEPC_DIR, 'lib'),
            ]
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

    #@staticmethod
    def get_slepc_dir(slepc_dir):
        if not slepc_dir: return None
        slepc_dir = os.path.expandvars(slepc_dir)
        if not slepc_dir or '$SLEPC_DIR' in slepc_dir:
            try:
                import slepc
                slepc_dir = slepc.get_slepc_dir()
            except ImportError:
                log.warn("SLEPC_DIR not specified")
                return None
        slepc_dir = os.path.expanduser(slepc_dir)
        slepc_dir = os.path.abspath(slepc_dir)
        if not os.path.isdir(slepc_dir):
            log.warn('invalid SLEPC_DIR:  %s' % slepc_dir)
            return None
        return slepc_dir
    get_slepc_dir = staticmethod(get_slepc_dir)

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


class build_src(_build_src):
    pass


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


class install(_install):
    pass

class clean(_clean):
    pass

class test(_test):
    pass

class sdist(_sdist):
    pass

# --------------------------------------------------------------------
