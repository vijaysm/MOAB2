from distutils.core import setup, Extension, Command
from distutils.sysconfig import parse_makefile
from distutils.errors import DistutilsOptionError
from distutils.command.build_ext import build_ext
from unittest import TextTestRunner, TestLoader
import re
import os
import sys

def pair_fun(pre, post):
    def tmp(self):
        pre(self)
        post(self)
    return tmp

def add_imesh_defs(imesh_dir, self):
        defs = parse_makefile( os.path.join(imesh_dir, 'iMesh-Defs.inc') )

        lib_match = re.compile(r'(?:(?<=\s)|^)-([lL])\s*(\S*)')
        for match in lib_match.finditer( defs['IMESH_LIBS'] ):
            if match.group(1) == 'l':
                self.libraries.append( match.group(2) )
            elif match.group(1) == 'L':
                self.library_dirs.append( match.group(2) )

        inc_match = re.compile(r'(?:(?<=\s)|^)-(I)\s*(\S*)')
        for match in inc_match.finditer( defs['IMESH_INCLUDES'] ):
            self.include_dirs.append( match.group(2) )

def new_init(self):
    self.imesh_dir = None

def new_fin(self):
    if self.imesh_dir:
        add_imesh_defs(self.imesh_dir, self)

build_ext.user_options.append(('imesh-dir=', None,
                               'root directory for iMesh interface'))
build_ext.initialize_options = pair_fun(build_ext.initialize_options, new_init)
build_ext.finalize_options   = pair_fun(build_ext.finalize_options,   new_fin)


def rebase(path):
    return os.path.normpath( os.path.join( os.path.abspath(sys.argv[0]),
                                           '../'+path) )

iBase = Extension('itaps.iBase',
                  depends = map(rebase, ['common.h', 'iBase_Python.h']),
                  sources = map(rebase, ['iBase.c'])
                  )

iMesh = Extension('itaps.iMesh',
                  depends = map(rebase,
                                ['common.h', 'errors.h', 'iMesh_Python.h',
                                 'iBase_Python.h', 'iMesh_entSet.inl',
                                 'iMesh_iter.inl', 'iMesh_tag.inl',
                                 'numpy_extensions.h', 'numpy_extensions.inl']),
                  sources = map(rebase, ['iMesh.c'])
                  )

class TestCommand(Command):
    description = 'execute a variety of unit tests'
    user_options = [
        ('verbosity=', 'v', 'verbosity level'),
        ]

    def initialize_options(self):
        self.verbosity = 1

    def finalize_options(self):
        try:
            self.verbosity = int(self.verbosity)
        except ValueError:
            raise DistutilsOptionError('"verbosity" option must be an integer')

    def run(self):
        root = rebase('test')
        old = os.getcwd()
        tests = []
        regex = re.compile(r'^(.*).py$')

        os.chdir(root)
        for file in os.listdir('.'):
            match = regex.search(file)
            if match != None:
                tests.append( 'test.' + match.group(1) )

        TextTestRunner(verbosity=self.verbosity).run(
            TestLoader().loadTestsFromNames(tests) )
        os.chdir(old)


class DocCommand(Command):
    description = 'build documentation'
    user_options = [
        ('builder=', 'b', 'documentation builder'),
        ('target=',  't', 'target directory'),
        ]

    def initialize_options(self):
        self.builder = 'html'
        self.target = None

    def finalize_options(self):
        if self.target:
            self.target = os.path.abspath(self.target)
        else:
            self.target = '_build/' + self.builder

    def run(self):
        root = rebase('doc')
        old = os.getcwd()

        os.chdir(root)
        os.system('sphinx-build -b "%s" -d _build/doctrees . "%s"' %
                  (self.builder, self.target))
        os.chdir(old)        

class PerfCommand(Command):
    description = 'build/execute performance tests'
    user_options = [
        ('file=',  'F', 'file or directory containing test file(s)'),
        ('count=', 'c', 'number of times to test'),
        ]

    def initialize_options(self):
        self.file = None
        self.count = 20

    def finalize_options(self):
        if not self.file:
            raise DistutilsOptionError('"file" must be specified')

        try:
            self.count = int(self.count)
        except ValueError:
            raise DistutilsOptionError('"count" option must be an integer')

    def run(self):
        for cmd_name in self.get_sub_commands():
            self.run_command(cmd_name)

    sub_commands = [
        ('perf_build', None),
        ('perf_run',   None),
        ]

class PerfBuildCommand(Command):
    description = 'build performance tests'

    sep_by = " (separated by '%s')" % os.pathsep
    user_options = [
        ('include-dirs=', 'I',
         'list of directories to search for header files' + sep_by),
        ('libraries=', 'l',
         'external C libraries to link with'),
        ('library-dirs=', 'L',
         'directories to search for external C libraries' + sep_by),
        ('imesh-dir=', None, 'root directory for iMesh interface'),
        ]

    def initialize_options(self):
        self.include_dirs = []
        self.library_dirs = []
        self.libraries = []
        self.imesh_dir = None

    def finalize_options(self):
        if self.imesh_dir:
            add_imesh_defs(self.imesh_dir, self)

    def run(self):
        root = rebase('perf')
        old = os.getcwd()
        os.chdir(root)

        from distutils.ccompiler import new_compiler
        self.compiler = new_compiler()

        objs = self.compiler.compile(['perf.c'], include_dirs=self.include_dirs)
        self.compiler.link_executable(objs, 'perf',
                                      library_dirs=self.library_dirs,
                                      libraries=self.libraries)

        os.chdir(old)

class PerfRunCommand(Command):
    description = 'execute performance tests'
    user_options = [
        ('file=',  'F', 'file or directory containing test file(s)'),
        ('count=', 'c', 'number of times to test'),
        ]

    def initialize_options(self):
        self.file = None
        self.count = 20

    def finalize_options(self):
        self.set_undefined_options('perf',
                                   ('file', 'file'),
                                   ('count', 'count') )

        if not self.file:
            raise DistutilsOptionError('"file" must be specified')

        try:
            self.count = int(self.count)
        except ValueError:
            raise DistutilsOptionError('"count" option must be an integer')

    def run(self):
        root = rebase('perf')
        old = os.getcwd()
        os.chdir(root)
        os.system('python perf.py -c%d "%s"' % (self.count, self.file))
        os.chdir(old) 

setup(name = 'PyTAPS',
      version = '1.0',
      description = 'Python bindings for iBase and iMesh interfaces',
      author = 'Jim Porter',
      author_email = 'jvporter@wisc.edu',

      requires = ['numpy'],
      provides = ['itaps'],

      package_dir = {'itaps': rebase('pkg')},
      packages = ['itaps'],
      ext_modules = [iBase, iMesh],
      py_modules = ['itaps.helpers'],

      cmdclass = { 'test'       : TestCommand,
                   'doc'        : DocCommand,
                   'perf'       : PerfCommand,
                   'perf_build' : PerfBuildCommand,
                   'perf_run'   : PerfRunCommand
                   }
      )
