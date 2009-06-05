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
compileargs = ''

if 'IMESHPATH' in os.environ:
    defs = parse_makefile( os.path.join(os.environ['IMESHPATH'],
                                        'lib/iMesh-Defs.inc') )

    compileargs = defs['IMESH_INCLUDES']

    lib_match = re.compile(r'(?:(?<=\s)|^)-([lL])\s*(\S*)')
    for match in lib_match.finditer( defs['IMESH_LIBS'] ):
        if match.group(1) == 'l':
            iMesh_libs.append( match.group(2) )
        elif match.group(1) == 'L':
            iMesh_libdirs.append( match.group(2) )

iBase = Extension('itaps.iBase',
                  extra_compile_args = [ compileargs ],
                  depends = ['common.h', 'iBase_Python.h'],
                  sources = ['iBase.c']
                  )

iMesh = Extension('itaps.iMesh',
                  extra_compile_args = [ compileargs ],
                  libraries = iMesh_libs,
                  library_dirs = iMesh_libdirs,
                  depends = ['common.h', 'iMesh_Python.h', 'iBase_Python.h',
                             'errors.h'],
                  sources = ['iMesh.c', 'iMesh_iter.c', 'iMesh_entSet.c',
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
        root = os.path.normpath(os.path.join(
                os.path.abspath(sys.argv[0]), '../test'))
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
        ('target=', 't', 'target directory')
        ]

    def initialize_options(self):
        self.builder = 'html'
        self.target = None
        self.root = os.path.normpath(os.path.join(
                os.path.abspath(sys.argv[0]), '../doc'))

    def finalize_options(self):
        if not self.target:
            self.target = '_build/' + self.builder
        self.target = os.path.abspath(os.path.join(self.root, self.target))

    def run(self):
        old = os.getcwd()

        os.chdir(self.root)
        os.system('sphinx-build -b "%s" -d _build/doctrees . "%s"' %
                  (self.builder, self.target))
        os.chdir(old)        


setup(name = 'PyTAPS',
      version = '1.0b1',
      description = 'Python bindings for iBase and iMesh interfaces',
      author = 'Jim Porter',
      author_email = 'jvporter@wisc.edu',

      package_dir = {'itaps': 'pkg'},
      packages = ['itaps'],
      ext_modules = [iBase, iMesh],
      py_modules = ['itaps.helpers'],

      cmdclass = { 'test' : TestCommand,
                   'doc'  : DocCommand }
      )
