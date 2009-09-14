#ifndef PYTAPS_IMESH_PYTHON_H
#define PYTAPS_IMESH_PYTHON_H

#include <Python.h>
#include <iMesh.h>
#include "iBase_Python.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    PyObject_HEAD
    iMesh_Instance mesh;
} iMesh_Object;

typedef struct
{
    PyObject_HEAD
    iMesh_Object *mesh;
    int is_arr;
    union
    {
        iMesh_EntityIterator    iter;
        iMesh_EntityArrIterator arr_iter;
    };
} iMeshIter_Object;

typedef struct
{
    iBaseEntitySet_Object set;
    iMesh_Object *mesh;
} iMeshEntitySet_Object;

#define iMeshEntitySet_NewRaw()                         \
    (iMeshEntitySet_Object*)PyObject_CallObject(        \
        (PyObject*)&iMeshEntitySet_Type,NULL)

#define iMeshEntitySet_Check(o)                         \
    PyObject_TypeCheck((o),&iMeshEntitySet_Type)

#define iMeshEntitySet_GetMesh(o)                       \
    ( ((iMeshEntitySet_Object*)(o))->mesh )

typedef struct
{
    iBaseTag_Object tag;
    iMesh_Object *mesh;
} iMeshTag_Object;

#define iMeshTag_NewRaw()                               \
    (iMeshTag_Object*)PyObject_CallObject(              \
        (PyObject*)&iMeshTag_Type,NULL)

#define iMeshTag_Check(o)                               \
    PyObject_TypeCheck((o),&iMeshTag_Type)

#define iMeshTag_GetMesh(o)                             \
    ( ((iMeshTag_Object*)(o))->mesh )


#ifndef _IMESH_MODULE

#if defined(PY_IMESH_UNIQUE_SYMBOL)
#define IMesh_API PY_IMESH_UNIQUE_SYMBOL
#endif

#if defined(NO_IMPORT) || defined(NO_IMPORT_IMESH)
extern void **IMesh_API;
#elif defined(PY_IMESH_UNIQUE_SYMBOL)
void **IMesh_API;
#else
static void **IMesh_API = NULL;
#endif

#define iMesh_Type          (*(PyTypeObject*)IMesh_API[0])
#define iMeshIter_Type      (*(PyTypeObject*)IMesh_API[1])
#define iMeshEntitySet_Type (*(PyTypeObject*)IMesh_API[2])
#define NPY_IMESHENTSET     (*(int*)         IMesh_API[3])
#define iMeshTag_Type       (*(PyTypeObject*)IMesh_API[4])
#define NPY_IMESHTAG        (*(int*)         IMesh_API[5])



#if !defined(NO_IMPORT_IMESH) && !defined(NO_IMPORT)
static int import_iMesh(void)
{
    PyObject *module = PyImport_ImportModule("itaps.iMesh");
    PyObject *c_api = NULL;

    if(module == NULL)
        return -1;

    c_api = PyObject_GetAttrString(module,"_C_API");
    if(c_api == NULL)
    {
        Py_DECREF(module);
        return -2;
    }

    if(PyCObject_Check(c_api))
        IMesh_API = (void **)PyCObject_AsVoidPtr(c_api);

    Py_DECREF(c_api);
    Py_DECREF(module);

    if(IMesh_API == NULL)
        return -3;
    return 0;
}
#endif

#endif


#ifdef __cplusplus
} // extern "C"
#endif

#endif
