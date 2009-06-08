from distutils.core import setup, Extension, Command
from distutils.sysconfig import parse_makefile
from distutils.errors import DistutilsOptionError
from unittest import TextTestRunner, TestLoader
import re
import os
import sys

# this is an enormous hack because distutils won't let me specify libraries in
# 'extra_link_args' for unknown reasons
iMesh_libs = []
iMesh_libdirs = []
iMesh_incs = []
compileargs = ''

if 'IMESHPATH' in os.environ:
    defs = parse_makefile( os.path.join(os.environ['IMESHPATH'],
                                        'lib/iMesh-Defs.inc') )

    lib_match = re.compile(r'(?:(?<=\s)|^)-([lL])\s*(\S*)')
    for match in lib_match.finditer( defs['IMESH_LIBS'] ):
        if match.group(1) == 'l':
            iMesh_libs.append( match.group(2) )
        elif match.group(1) == 'L':
            iMesh_libdirs.append( match.group(2) )

    inc_match = re.compile(r'(?:(?<=\s|^)-(I)\s*(\S*)')
    for match in defs['IMESH_INCLUDES']:
        iMesh_incs.append( match.group(2) )

iBase = Extension('itaps.iBase',
                  include_dirs = iMesh_incs,
                  depends      = ['common.h', 'iBase_Python.h'],
                  sources      = ['iBase.c']
                  )

iMesh = Extension('itaps.iMesh',
                  include_dirs = iMesh_incs,
                  libraries    = iMesh_libs,
                  library_dirs = iMesh_libdirs,
                  depends      = ['common.h', 'iMesh_Python.h',
                                  'iBase_Python.h', 'errors.h'],
                  sources      = ['iMesh.c', 'iMesh_iter.c', 'iMesh_entSet.c',
                                  'iMesh_tag.c']
                  )

class TestCommand(Command):
    description = 'Execute a variety of unit tests'
    user_options = [
        ('verbosity=', 'v', 'verbosity level')
        ]

    def initialize_options(self):
        self.verbosity = 1

    def finalize_options(self):
        try:
            self.verbosity = int(self.verbosity)
        except ValueError:
            raise DistutilsOptionError('"verbosity" option must be an integer')

    def run(self):
        root = os.path.normpath(os.path.join(os.path.abspath(sys.argv[0]),
                                             '../test'))
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
    description = 'Build documentation'
    user_options = [
        ('builder=', 'b', 'documentation builder'),
        ('target=',  't', 'target directory')
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
        root = os.path.normpath(os.path.join(os.path.abspath(sys.argv[0]),
                                             '../doc'))
        old = os.getcwd()

        os.chdir(root)
        os.system('sphinx-build -b "%s" -d _build/doctrees . "%s"' %
                  (self.builder, self.target))
        os.chdir(old)        

class PerfCommand(Command):
    description = 'Execute performance tests'
    user_options = [
        ('file=',  'f', 'test file'),
        ('count=', 'c', 'number of times to test')
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
        root = os.path.normpath(os.path.join(os.path.abspath(sys.argv[0]),
                                             '../perf'))
        old = os.getcwd()

        os.chdir(root)
        os.system('make all')
        os.system('python perf.py -c%d "%s"' % (self.count, self.file))
        os.chdir(old)

setup(name = 'PyTAPS',
      version = '1.0b1',
      description = 'Python bindings for iBase and iMesh interfaces',
      author = 'Jim Porter',
      author_email = 'jvporter@wisc.edu',

      requires = ['numpy'],
      provides = ['itaps'],

      package_dir = {'itaps': 'pkg'},
      packages = ['itaps'],
      ext_modules = [iBase, iMesh],
      py_modules = ['itaps.helpers'],

      cmdclass = { 'test' : TestCommand,
                   'doc'  : DocCommand,
                   'perf' : PerfCommand }
      )
