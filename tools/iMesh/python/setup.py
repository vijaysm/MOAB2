from distutils.core import setup, Extension
from distutils.sysconfig import parse_makefile
import re
import os

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


setup(name = 'PyTAPS',
      version = '1.0b1',
      description = 'Python bindings for iBase and iMesh interfaces',
      author = 'Jim Porter',
      author_email = 'jvporter@wisc.edu',
      packages = ['itaps'],
      ext_modules = [iBase, iMesh],
      py_modules = ['itaps.helpers']
      )
