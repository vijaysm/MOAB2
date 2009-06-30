#include "iMesh_Python.h"
#include "iBase_Python.h"
#include "structmember.h"

static int
iMeshIterObj_init(iMeshIter_Object *self,PyObject *args,PyObject *kwds)
{
    static char *kwlist[] = {"set","type","topology","size",0};
    iMeshEntitySet_Object *set;
    int type,topo,array_size=1,err;

    if( !PyArg_ParseTupleAndKeywords(args,kwds,"O!ii|i",kwlist,
                                     &iMeshEntitySet_Type,&set,&type,&topo,
                                     &array_size) )
        return -1;

    self->mesh = set->mesh;
    Py_INCREF(self->mesh);

    if(array_size == 1)
    {
        self->is_arr = 0;
        iMesh_initEntIter(self->mesh->mesh,set->set.handle,type,topo,
                          &self->iter,&err);
    }
    else
    {
        self->is_arr = 1;
        iMesh_initEntArrIter(self->mesh->mesh,set->set.handle,type,topo,
                             array_size,&self->arr_iter,&err);
    }
    if(checkError(self->mesh->mesh,err))
        return -1;

    return 0;
}

static void
iMeshIterObj_dealloc(iMeshIter_Object *self)
{
    if(self->mesh && self->iter)
    {
        int err;
        if(self->is_arr)
            iMesh_endEntArrIter(self->mesh->mesh,self->arr_iter,&err);
        else
            iMesh_endEntIter(self->mesh->mesh,self->iter,&err);
    }

    Py_XDECREF(self->mesh);
    self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
iMeshIterObj_reset(iMeshIter_Object *self,PyObject *args)
{
    int err;
    if(self->is_arr)
        iMesh_resetEntArrIter(self->mesh->mesh,self->arr_iter,&err);
    else
        iMesh_resetEntIter(self->mesh->mesh,self->iter,&err);

    if(checkError(self->mesh->mesh,err))
        return NULL;
    Py_RETURN_NONE;
}

static PyObject *
iMeshIterObj_iternext(iMeshIter_Object *self)
{
    int has_data,err;

    if(self->is_arr)
    {
        iBase_EntityHandle *entities=0;
        int alloc=0,size;

        iMesh_getNextEntArrIter(self->mesh->mesh,self->arr_iter,&entities,
                                &alloc,&size,&has_data,&err);
        if(checkError(self->mesh->mesh,err))
            return NULL;
        if(!has_data)
            return NULL;

        npy_intp dims[] = {size};
        return PyArray_NewFromMalloc(1,dims,NPY_IBASEENT,entities);
    }
    else
    {
        iBase_EntityHandle entity;
        iMesh_getNextEntIter(self->mesh->mesh,self->iter,&entity,&has_data,
                             &err);
        if(checkError(self->mesh->mesh,err))
            return NULL;
        if(!has_data)
            return NULL;

        iBaseEntity_Object *o = iBaseEntity_New();
        o->handle = entity;
        return (PyObject*)o;
    }
}

static PyMethodDef iMeshIterObj_methods[] = {
    { "reset", (PyCFunction)iMeshIterObj_reset, METH_NOARGS,
      "Reset the iterator"
    },
    {0}
};

static PyMemberDef iMeshIterObj_members[] = {
    {"instance", T_OBJECT_EX, offsetof(iMeshIter_Object, mesh), READONLY,
     "base iMesh instance"},
    {0}
};

static PyTypeObject iMeshIter_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                                   /* ob_size */
    "itaps.iMesh.Iterator",              /* tp_name */
    sizeof(iMeshIter_Object),            /* tp_basicsize */
    0,                                   /* tp_itemsize */
    (destructor)iMeshIterObj_dealloc,    /* tp_dealloc */
    0,                                   /* tp_print */
    0,                                   /* tp_getattr */
    0,                                   /* tp_setattr */
    0,                                   /* tp_compare */
    0,                                   /* tp_repr */
    0,                                   /* tp_as_number */
    0,                                   /* tp_as_sequence */
    0,                                   /* tp_as_mapping */
    0,                                   /* tp_hash */
    0,                                   /* tp_call */
    0,                                   /* tp_str */
    0,                                   /* tp_getattro */
    0,                                   /* tp_setattro */
    0,                                   /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT |
    Py_TPFLAGS_HAVE_ITER |
    Py_TPFLAGS_BASETYPE,                 /* tp_flags */
    "iMesh iterator object",             /* tp_doc */
    0,                                   /* tp_traverse */
    0,                                   /* tp_clear */
    0,                                   /* tp_richcompare */
    0,                                   /* tp_weaklistoffset */
    PyObject_SelfIter,                   /* tp_iter */
    (iternextfunc)iMeshIterObj_iternext, /* tp_iternext */
    iMeshIterObj_methods,                /* tp_methods */
    iMeshIterObj_members,                /* tp_members */
    0,                                   /* tp_getset */
    0,                                   /* tp_base */
    0,                                   /* tp_dict */
    0,                                   /* tp_descr_get */
    0,                                   /* tp_descr_set */
    0,                                   /* tp_dictoffset */
    (initproc)iMeshIterObj_init,         /* tp_init */
    0,                                   /* tp_alloc */
    0,                                   /* tp_new */
};
