#ifndef PYTAPS_NUMPY_EXTENSIONS_H
#define PYTAPS_NUMPY_EXTENSIONS_H

#include <numpy/arrayobject.h>

typedef struct
{
    PyObject_HEAD
    PyObject *base;
    void *memory;
} ArrDealloc_Object;

static void
ArrDeallocObj_dealloc(ArrDealloc_Object *self)
{
    free(self->memory);
    Py_XDECREF(self->base);

    self->ob_type->tp_free((PyObject *)self);
}

static PyTypeObject ArrDealloc_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                                          /* ob_size */
    "arrdealloc",                               /* tp_name */
    sizeof(ArrDealloc_Object),                  /* tp_basicsize */
    0,                                          /* tp_itemsize */
    (destructor)ArrDeallocObj_dealloc,          /* tp_dealloc */
    0,                                          /* tp_print */
    0,                                          /* tp_getattr */
    0,                                          /* tp_setattr */
    0,                                          /* tp_compare */
    0,                                          /* tp_repr */
    0,                                          /* tp_as_number */
    0,                                          /* tp_as_sequence */
    0,                                          /* tp_as_mapping */
    0,                                          /* tp_hash */
    0,                                          /* tp_call */
    0,                                          /* tp_str */
    0,                                          /* tp_getattro */
    0,                                          /* tp_setattro */
    0,                                          /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                         /* tp_flags */
    "Internal deallocator object",              /* tp_doc */
};

#define PyArray_NewFromMalloc(nd,dims,typenum,data) \
    PyArray_NewFromMallocBaseStrided(nd,dims,NULL,typenum,data,0)

#define PyArray_NewFromMallocBase(nd,dims,typenum,data,base) \
    PyArray_NewFromMallocBaseStrided(nd,dims,NULL,typenum,data,base)

#define PyArray_NewFromMallocStrided(nd,dims,strides,typenum,data) \
    PyArray_NewFromMallocBaseStrided(nd,dims,strides,typenum,data,0)

static PyObject *
PyArray_NewFromMallocBaseStrided(int nd,npy_intp *dims,npy_intp *strides,
                                 int typenum,void *data,PyObject *base)
{
    ArrDealloc_Object *newobj;
    PyObject *arr = PyArray_New(&PyArray_Type,nd,dims,typenum,strides,data,0,
                                NPY_CARRAY,NULL);

    newobj = PyObject_New(ArrDealloc_Object,&ArrDealloc_Type);
    newobj->memory = data;
    Py_XINCREF(base);
    newobj->base = base;

    PyArray_BASE(arr) = (PyObject*)newobj;
    return arr;
}

static PyObject *
PyArray_TryFromObject(PyObject *obj,int typenum,int min_depth,int max_depth)
{
    PyObject *ret = PyArray_FromAny(obj,PyArray_DescrFromType(typenum),
                                    min_depth,max_depth,NPY_C_CONTIGUOUS,NULL);
    PyErr_Clear();
    return ret;
}

static PyObject *
PyArray_ToVectors(PyObject *obj,int typenum,int nd,npy_intp vec_dim,int index)
{
    assert(index < nd);
    PyObject *ret = PyArray_FromAny(obj,PyArray_DescrFromType(typenum),nd,nd,
                                    NPY_C_CONTIGUOUS,NULL);
    if(ret == NULL)
        return NULL;
    if(PyArray_DIM(ret,index) != vec_dim)
    {
        Py_DECREF(ret);
        PyErr_SetString(PyExc_ValueError,"Expected n-d vector"); /* TODO */
        return NULL;
    }
    return ret;
}

#endif
