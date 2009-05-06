#pragma once

#ifndef _IBASE_MODULE

#if defined(PY_IBASE_UNIQUE_SYMBOL)
#define IBase_API PY_IBASE_UNIQUE_SYMBOL
#endif

#if defined(NO_IMPORT) || defined(NO_IMPORT_IBASE)
extern void **IBase_API;
#elif defined(PY_IBASE_UNIQUE_SYMBOL)
void **IBase_API;
#else
static void **IBase_API = NULL;
#endif

#define iBaseEntity_Type    (*(PyTypeObject*)IBase_API[0])
#define NPY_IBASEENT        (*(int*)         IBase_API[1])
#define iBaseEntitySet_Type (*(PyTypeObject*)IBase_API[2])
#define NPY_IBASEENTSET     (*(int*)         IBase_API[3])
#define iBaseTag_Type       (*(PyTypeObject*)IBase_API[4])
#define NPY_IBASETAG        (*(int*)         IBase_API[5])


#if !defined(NO_IMPORT_ARRAY) && !defined(NO_IMPORT)
static int import_iBase(void)
{
  PyObject *module = PyImport_ImportModule("itaps.iBase");
  PyObject *c_api = NULL;

  if(module == NULL)
    return -1;

  c_api = PyObject_GetAttrString(module,"_C_API");
  if(c_api == NULL)
  {
    Py_DECREF(module);
    return -1;
  }

  if(PyCObject_Check(c_api))
    IBase_API = (void **)PyCObject_AsVoidPtr(c_api);

  Py_DECREF(c_api);
  Py_DECREF(module);

  if(IBase_API == NULL)
    return -1;
  return 0;
}
#endif

#endif
