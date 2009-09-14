#ifndef PYTAPS_IBASE_PYTHON_H
#define PYTAPS_IBASE_PYTHON_H

#include <Python.h>
#include <iBase.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    PyObject_HEAD
    iBase_EntityHandle handle;
} iBaseEntity_Object;

#define iBaseEntity_New()                               \
    (iBaseEntity_Object*)PyObject_CallObject(           \
        (PyObject*)&iBaseEntity_Type,NULL)

#define iBaseEntity_Check(o)                            \
    PyObject_TypeCheck((o),&iBaseEntity_Type)

#define iBaseEntity_GetHandle(o)                        \
    ((iBaseEntity_Object*)(o))->handle

typedef struct
{
    PyObject_HEAD
    iBase_EntitySetHandle handle;
} iBaseEntitySet_Object;

#define iBaseEntitySet_New()                            \
    (iBaseEntitySet_Object*)PyObject_CallObject(        \
        (PyObject*)&iBaseEntitySet_Type,NULL)

#define iBaseEntitySet_Check(o)                         \
  PyObject_TypeCheck((o),&iBaseEntitySet_Type)

#define iBaseEntitySet_GetHandle(o)                     \
  ((iBaseEntitySet_Object*)(o))->handle

typedef struct
{
  PyObject_HEAD
  iBase_TagHandle handle;
} iBaseTag_Object;

#define iBaseTag_New()                                  \
    (iBaseTag_Object*)PyObject_CallObject(              \
        (PyObject*)&iBaseTag_Type,NULL)

#define iBaseTag_Check(o)                               \
    PyObject_TypeCheck((o),&iBaseTag_Type)

#define iBaseTag_GetHandle(o)                           \
    ((iBaseTag_Object*)(o))->handle




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


#if !defined(NO_IMPORT_IBASE) && !defined(NO_IMPORT)
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

#ifdef __cplusplus
} // extern "C"
#endif

#endif
