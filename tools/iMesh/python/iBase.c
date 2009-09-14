#define _IBASE_MODULE
#include "iBase_Python.h"
#include "common.h"

#include <numpy/arrayobject.h>

SIMPLE_TYPE(iBaseEntity,"itaps.iBase.Entity","");
SIMPLE_TYPE(iBaseEntitySet,"itaps.iBase.EntitySet","");
SIMPLE_TYPE(iBaseTag,"itaps.iBase.Tag","");

static PyObject *
iBaseEntity_FromHandle(iBase_EntityHandle h)
{
    iBaseEntity_Object *o = iBaseEntity_New();
    o->handle = h;
    return (PyObject*)o;
}
static int NPY_IBASEENT;

static PyObject *
iBaseEntitySet_FromHandle(iBase_EntitySetHandle h)
{
    iBaseEntitySet_Object *o = iBaseEntitySet_New();
    o->handle = h;
    return (PyObject*)o;
}
static int NPY_IBASEENTSET;

static PyObject *
iBaseTag_FromHandle(iBase_TagHandle h)
{
    iBaseTag_Object *o = iBaseTag_New();
    o->handle = h;
    return (PyObject*)o;
}
static int NPY_IBASETAG;

ENUM_TYPE(iBaseType,           "iBase.Type",           "");
ENUM_TYPE(iBaseAdjCost,        "iBase.AdjCost",        "");
ENUM_TYPE(iBaseStorageOrder,   "iBase.StorageOrder",   "");
ENUM_TYPE(iBaseCreationStatus, "iBase.CreationStatus", "");

static PyMethodDef module_methods[] = {
    {0}
};


static PyObject *
iBaseEntArr_getitem(void *data,void *arr)
{
    return iBaseEntity_FromHandle(*(iBase_EntityHandle*)data);
}

static int
iBaseEntArr_setitem(PyObject *item,void *data,void *arr)
{
    if(!iBaseEntity_Check(item))
        return -1;
    *(iBase_EntityHandle*)data = iBaseEntity_GetHandle(item);
    return 0;
}

static void
iBaseEntArr_copyswapn(void *dest,npy_intp dstride,void *src,npy_intp sstride,
                      npy_intp n,int swap,void *arr)
{}

static void
iBaseEntArr_copyswap(void *dest,void *src,int swap,void *arr)
{}

static npy_bool
iBaseEntArr_nonzero(void *data,void *arr)
{
    return *(iBase_EntityHandle*)data != 0;
}

static PyObject *
iBaseEntObj_repr(iBaseEntity_Object *self)
{
    char out[64];
    snprintf(out,64,"<itaps.iBase.Entity %p>",self->handle);
    return Py_BuildValue("s",out);
}

static PyObject *
iBaseEntObj_richcompare(iBaseEntity_Object *lhs,iBaseEntity_Object *rhs,int op)
{
    if(!iBaseEntity_Check(lhs) || !iBaseEntity_Check(rhs))
        return Py_NotImplemented;

    switch(op)
    {
    case Py_EQ:
        return PyBool_FromLong(lhs->handle == rhs->handle);
    case Py_NE:
        return PyBool_FromLong(lhs->handle != rhs->handle);
    default:
        return Py_NotImplemented;
  }
}

static PyArray_ArrFuncs iBaseEntArr_funcs = {
    {0},
    iBaseEntArr_getitem,
    iBaseEntArr_setitem,
    iBaseEntArr_copyswapn,
    iBaseEntArr_copyswap,
    0,
    0,
    0,
    0,
    0,
    iBaseEntArr_nonzero
};


static PyObject *
iBaseEntSetArr_getitem(void *data,void *arr)
{
    return iBaseEntitySet_FromHandle(*(iBase_EntitySetHandle*)data);
}

static int
iBaseEntSetArr_setitem(PyObject *item,void *data,void *arr)
{
    if(!iBaseEntitySet_Check(item))
        return -1;
    *(iBase_EntitySetHandle*)data = iBaseEntitySet_GetHandle(item);
    return 0;
}

static void
iBaseEntSetArr_copyswapn(void *dest,npy_intp dstride,void *src,
                         npy_intp sstride,npy_intp n,int swap,void *arr)
{}

static void
iBaseEntSetArr_copyswap(void *dest,void *src,int swap,void *arr)
{}

static npy_bool
iBaseEntSetArr_nonzero(void *data,void *arr)
{
    return *(iBase_EntitySetHandle*)data != 0;
}

static PyObject *
iBaseEntSetObj_repr(iBaseEntitySet_Object *self)
{
    char out[64];
    snprintf(out,64,"<itaps.iBase.EntitySet %p>",self->handle);
    return Py_BuildValue("s",out);
}

static PyObject *
iBaseEntSetObj_richcompare(iBaseEntitySet_Object *lhs,
                           iBaseEntitySet_Object *rhs,int op)
{
    if(!iBaseEntitySet_Check(lhs) || !iBaseEntitySet_Check(rhs))
        return Py_NotImplemented;

    switch(op)
    {
    case Py_EQ:
        return PyBool_FromLong(lhs->handle == rhs->handle);
    case Py_NE:
        return PyBool_FromLong(lhs->handle != rhs->handle);
    default:
        return Py_NotImplemented;
    }
}

static PyArray_ArrFuncs iBaseEntSetArr_funcs = {
    {0},
    iBaseEntSetArr_getitem,
    iBaseEntSetArr_setitem,
    iBaseEntSetArr_copyswapn,
    iBaseEntSetArr_copyswap,
    0,
    0,
    0,
    0,
    0,
    iBaseEntSetArr_nonzero
};



static PyObject *
iBaseTagArr_getitem(void *data,void *arr)
{
    return iBaseTag_FromHandle(*(iBase_TagHandle*)data);
}

static int
iBaseTagArr_setitem(PyObject *item,void *data,void *arr)
{
    if(!iBaseTag_Check(item))
        return -1;
    *(iBase_TagHandle*)data = iBaseTag_GetHandle(item);
    return 0;
}

static void
iBaseTagArr_copyswapn(void *dest,npy_intp dstride,void *src,npy_intp sstride,
                      npy_intp n,int swap,void *arr)
{}

static void
iBaseTagArr_copyswap(void *dest,void *src,int swap,void *arr)
{}

static npy_bool
iBaseTagArr_nonzero(void *data,void *arr)
{
    return *(iBase_TagHandle*)data != 0;
}

static PyArray_ArrFuncs iBaseTagArr_funcs = {
    {0},
    iBaseTagArr_getitem,
    iBaseTagArr_setitem,
    iBaseTagArr_copyswapn,
    iBaseTagArr_copyswap,
    0,
    0,
    0,
    0,
    0,
    iBaseTagArr_nonzero
};



PyMODINIT_FUNC initiBase(void)
{
    PyArray_Descr *descr;
    PyObject *m = Py_InitModule("iBase",module_methods);
    import_array();

    /***** register C API *****/
    static void *IBase_API[6];
    PyObject *api_obj;

    /* Initialize the C API pointer array */
    IBase_API[0] = &iBaseEntity_Type;
    IBase_API[1] = &NPY_IBASEENT;
    IBase_API[2] = &iBaseEntitySet_Type;
    IBase_API[3] = &NPY_IBASEENTSET;
    IBase_API[4] = &iBaseTag_Type;
    IBase_API[5] = &NPY_IBASETAG;

    /* Create a CObject containing the API pointer array's address */
    api_obj = PyCObject_FromVoidPtr(IBase_API,NULL);

    if(api_obj != NULL)
        PyModule_AddObject(m, "_C_API", api_obj);

    /***** initialize type enum *****/
    REGISTER_SIMPLE(m,"Type",iBaseType);

    ADD_ENUM(iBaseType,"vertex", iBase_VERTEX);
    ADD_ENUM(iBaseType,"edge",   iBase_EDGE);
    ADD_ENUM(iBaseType,"face",   iBase_FACE);
    ADD_ENUM(iBaseType,"region", iBase_REGION);
    ADD_ENUM(iBaseType,"all",    iBase_ALL_TYPES);

    /***** initialize adjacency cost enum *****/
    REGISTER_SIMPLE(m,"AdjCost",iBaseAdjCost);

    ADD_ENUM(iBaseAdjCost,"unavailable",     iBase_UNAVAILABLE);
    ADD_ENUM(iBaseAdjCost,"all_order_1",     iBase_ALL_ORDER_1);
    ADD_ENUM(iBaseAdjCost,"all_order_logn",  iBase_ALL_ORDER_LOGN);
    ADD_ENUM(iBaseAdjCost,"all_order_n",     iBase_ALL_ORDER_N);
    ADD_ENUM(iBaseAdjCost,"some_order_1",    iBase_SOME_ORDER_1);
    ADD_ENUM(iBaseAdjCost,"some_order_logn", iBase_SOME_ORDER_LOGN);
    ADD_ENUM(iBaseAdjCost,"some_order_n",    iBase_SOME_ORDER_N);

    /***** initialize storage order enum *****/
    REGISTER_SIMPLE(m,"StorageOrder",iBaseStorageOrder);

    ADD_ENUM(iBaseStorageOrder,"blocked",     iBase_BLOCKED);
    ADD_ENUM(iBaseStorageOrder,"interleaved", iBase_INTERLEAVED);

    /***** initialize creation status enum *****/
    REGISTER_SIMPLE(m,"CreationStatus",iBaseCreationStatus);

    ADD_ENUM(iBaseCreationStatus,"new",        iBase_NEW);
    ADD_ENUM(iBaseCreationStatus,"exists",     iBase_ALREADY_EXISTED);
    ADD_ENUM(iBaseCreationStatus,"duplicated", iBase_CREATED_DUPLICATE);
    ADD_ENUM(iBaseCreationStatus,"failed",     iBase_CREATION_FAILED);

    /***** initialize iBaseEntity handle *****/
    iBaseEntity_Type.tp_repr = (reprfunc)iBaseEntObj_repr;
    iBaseEntity_Type.tp_richcompare = (richcmpfunc)iBaseEntObj_richcompare;

    REGISTER_SIMPLE(m,"Entity",iBaseEntity);

    /***** initialize iBaseEntitySet handle *****/
    iBaseEntitySet_Type.tp_repr = (reprfunc)iBaseEntSetObj_repr;
    iBaseEntitySet_Type.tp_richcompare = 
        (richcmpfunc)iBaseEntSetObj_richcompare;

    REGISTER_SIMPLE(m,"EntitySet",iBaseEntitySet);

    /***** initialize iBaseTag handle *****/
    REGISTER_SIMPLE(m,"Tag",iBaseTag);

    /***** initialize iBaseEntity array type *****/
    descr = PyArray_DescrNewFromType(NPY_INTP);
    descr->f = &iBaseEntArr_funcs;

    descr->typeobj = &iBaseEntity_Type;
    descr->kind = 'V';
    descr->type = 'j';
    descr->hasobject = NPY_USE_GETITEM|NPY_USE_SETITEM;
    descr->elsize = sizeof(iBase_EntityHandle);

    NPY_IBASEENT = PyArray_RegisterDataType(descr);

    /***** initialize iBaseEntitySet array type *****/
    descr = PyArray_DescrNewFromType(NPY_INTP);
    descr->f = &iBaseEntSetArr_funcs;

    descr->typeobj = &iBaseEntitySet_Type;
    descr->kind = 'V';
    descr->type = 'J';
    descr->hasobject = NPY_USE_GETITEM|NPY_USE_SETITEM;
    descr->elsize = sizeof(iBase_EntitySetHandle);

    NPY_IBASEENTSET = PyArray_RegisterDataType(descr);

    /***** initialize iBaseTag array type *****/
    descr = PyArray_DescrNewFromType(NPY_INTP);
    descr->f = &iBaseTagArr_funcs;

    descr->typeobj = &iBaseTag_Type;
    descr->kind = 'V';
    descr->type = 'T';
    descr->hasobject = NPY_USE_GETITEM|NPY_USE_SETITEM;
    descr->elsize = sizeof(iBase_TagHandle);

    NPY_IBASETAG = PyArray_RegisterDataType(descr);
}
