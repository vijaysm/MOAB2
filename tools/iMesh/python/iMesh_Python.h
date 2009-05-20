#pragma once

#include "common.h"

#include <Python.h>
#include <iMesh.h>

#define PY_IBASE_UNIQUE_SYMBOL itaps_IBASE_API
#include "iBase_Python.h"

#define PY_ARRAY_UNIQUE_SYMBOL itaps_ARRAY_API
#include <numpy/arrayobject.h>

typedef struct
{
    PyObject_HEAD
    PyObject *base;
    void *memory;
} ArrDealloc_Object;


#define PyArray_NewFromMalloc(nd,dims,typenum,data) \
    PyArray_NewFromMallocBaseStrided(nd,dims,NULL,typenum,data,0)

#define PyArray_NewFromMallocBase(nd,dims,typenum,data,base) \
    PyArray_NewFromMallocBaseStrided(nd,dims,NULL,typenum,data,base)

#define PyArray_NewFromMallocStrided(nd,dims,strides,typenum,data) \
    PyArray_NewFromMallocBaseStrided(nd,dims,strides,typenum,data,0)

PyObject *
PyArray_NewFromMallocBaseStrided(int nd,npy_intp *dims,npy_intp *strides,
                                 int typenum,void *data,PyObject *base);

PyObject *
PyArray_TryFromObject(PyObject *obj,int typenum,int min_depth,int max_depth);


int checkError(iMesh_Instance mesh,int err);
enum iBase_TagValueType char_to_type(char c);
char type_to_char(enum iBase_TagValueType t);

typedef struct
{
    PyObject_HEAD
    iMesh_Instance mesh;
} iMeshObject;


typedef struct
{
    PyObject_HEAD
    iMesh_Instance mesh;
    int is_arr;
    union
    {
        iMesh_EntityIterator    iter;
        iMesh_EntityArrIterator arr_iter;
    };
} iMeshIter_Object;

extern PyTypeObject iMeshIter_Type;

typedef struct
{
    iBaseEntitySet_Object set;
    iMeshObject *mesh;
} iMeshEntitySet_Object;

extern PyTypeObject iMeshEntitySet_Type;
extern int NPY_IMESHENTSET;

#define iMeshEntitySet_New()                            \
    (iMeshEntitySet_Object*)PyObject_CallObject(        \
        (PyObject*)&iMeshEntitySet_Type,NULL)

#define iMeshEntitySet_Check(o)                         \
    PyObject_TypeCheck((o),&iMeshEntitySet_Type)

#define iMeshEntitySet_GetMesh(o)                       \
    ((iMeshEntitySet_Object*)(o))->mesh

typedef struct
{
    iBaseTag_Object tag;
    iMeshObject *mesh;
} iMeshTag_Object;

extern PyTypeObject iMeshTag_Type;
extern int NPY_IMESHTAG;

#define iMeshTag_New()                                  \
    (iMeshTag_Object*)PyObject_CallObject(              \
        (PyObject*)&iMeshTag_Type,NULL)

#define iMeshTag_Check(o)                               \
    PyObject_TypeCheck((o),&iMeshTag_Type)

#define iMeshTag_GetMesh(o)                             \
    ((iMeshTag_Object*)(o))->mesh
