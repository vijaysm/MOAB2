from distutils.core import setup, Extension

iBase = Extension('itaps.iBase',
                  depends = ['common.h', 'iBase_Python.h'],
                  sources = ['iBase.c']
                  )

iMesh = Extension('itaps.iMesh',
                  libraries = ['iMesh', 'm', 'gcc_s', 'c', 'MOAB',
                               'netcdf_c++', 'netcdf', 'hdf5', 'stdc++'],
                  depends = ['common.h', 'iMesh_Python.h', 'iBase_Python.h',
                             'errors.h'],
                  sources = ['iMesh.c', 'iMesh_iter.c', 'iMesh_entSet.c',
                             'iMesh_tag.c']
                  )


setup(name = 'PyTAPS',
      version = '0.1',
      description = 'Python bindings for iBase and iMesh interfaces',
      author = 'Jim Porter',
      author_email = 'jvporter@wisc.edu',
      packages = ['itaps'],
      ext_modules = [iBase, iMesh],
      py_modules = ['itaps.helpers']
      )
