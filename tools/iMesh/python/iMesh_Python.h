#pragma once

#include "common.h"

#include <Python.h>
#include <iMesh.h>

#define PY_IBASE_UNIQUE_SYMBOL itaps_IBASE_API
#include "iBase_Python.h"

#define PY_ARRAY_UNIQUE_SYMBOL itaps_ARRAY_API
#include <numpy/arrayobject.h>

typedef struct {
PyObject_HEAD
void *memory;
} _MyDeallocObject;

static void
_mydealloc_dealloc(_MyDeallocObject *self)
{
free(self->memory);
self->ob_type->tp_free((PyObject *)self);
}

static PyTypeObject _MyDeallocType = {
PyObject_HEAD_INIT(NULL)
0,                                          /*ob_size*/
"mydeallocator",                   /*tp_name*/
sizeof(_MyDeallocObject),    /*tp_basicsize*/
0,                                          /*tp_itemsize*/
_mydealloc_dealloc,             /*tp_dealloc*/
0,                         /*tp_print*/
0,                         /*tp_getattr*/
0,                         /*tp_setattr*/
0,                         /*tp_compare*/
0,                         /*tp_repr*/
0,                         /*tp_as_number*/
0,                         /*tp_as_sequence*/
0,                         /*tp_as_mapping*/
0,                         /*tp_hash */
0,                         /*tp_call*/
0,                         /*tp_str*/
0,                         /*tp_getattro*/
0,                         /*tp_setattro*/
0,                         /*tp_as_buffer*/
Py_TPFLAGS_DEFAULT,        /*tp_flags*/
"Internal deallocator object",           /* tp_doc */
};


/*#define PyArray_NewFromMallocData(nd, dims, typenum, data)    \
    PyArray_New(&PyArray_Type, nd, dims, typenum, NULL,         \
    data, 0, NPY_CARRAY|NPY_OWNDATA, NULL)*/

static PyObject *
PyArray_NewFromMallocData(int nd,npy_intp *dims,int typenum,void *data)
{
    PyObject *newobj;
    PyObject *arr =  PyArray_New(&PyArray_Type, nd, dims, typenum, NULL,
                                 data, 0, NPY_CARRAY, NULL);
    newobj = (PyObject*)PyObject_New(_MyDeallocObject, &_MyDeallocType);
    ((_MyDeallocObject *)newobj)->memory = data;
    PyArray_BASE(arr) = newobj;
    return arr;
}


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
