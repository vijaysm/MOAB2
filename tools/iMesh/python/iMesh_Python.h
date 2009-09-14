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


#ifdef __cplusplus
} // extern "C"
#endif

#endif
