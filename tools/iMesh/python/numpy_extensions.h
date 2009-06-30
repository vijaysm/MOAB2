#ifndef PYTAPS_NUMPY_EXTENSIONS_H
#define PYTAPS_NUMPY_EXTENSIONS_H

#include <numpy/arrayobject.h>

typedef struct
{
    PyObject_HEAD
    PyObject *base;
    void *memory;
} ArrDealloc_Object;

static PyTypeObject ArrDealloc_Type;

#define PyArray_NewFromMalloc(nd,dims,typenum,data) \
    PyArray_NewFromMallocBaseStrided(nd,dims,NULL,typenum,data,0)

#define PyArray_NewFromMallocBase(nd,dims,typenum,data,base) \
    PyArray_NewFromMallocBaseStrided(nd,dims,NULL,typenum,data,base)

#define PyArray_NewFromMallocStrided(nd,dims,strides,typenum,data) \
    PyArray_NewFromMallocBaseStrided(nd,dims,strides,typenum,data,0)

static PyObject *
PyArray_NewFromMallocBaseStrided(int nd,npy_intp *dims,npy_intp *strides,
                                 int typenum,void *data,PyObject *base);

static PyObject *
PyArray_TryFromObject(PyObject *obj,int typenum,int min_depth,int max_depth);

#endif
