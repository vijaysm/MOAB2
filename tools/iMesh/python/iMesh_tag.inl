#define NO_IMPORT_ARRAY
#define NO_IMPORT_IBASE

#include "iMesh_Python.h"
#include "iBase_Python.h"
#include "errors.h"
#include "structmember.h"

static char typechars[] = {'i','d','E','b'};

static enum iBase_TagValueType
char_to_type(char c)
{
    int i;
    for(i=0; i<sizeof(typechars); i++)
    {
        if(typechars[i] == c)
            return i;
    }
    return -1;
}

static char
type_to_char(enum iBase_TagValueType t)
{
    return typechars[t];
}


static void
iMeshTagObj_dealloc(iMeshTag_Object *self)
{
    Py_XDECREF(self->mesh);
    ((PyObject*)self)->ob_type->tp_free((PyObject*)self);
}

static PyObject *
iMeshTagObj_getName(iMeshTag_Object *self,void *closure)
{
    char name[512];
    int err;
    iMesh_getTagName(self->mesh->mesh,self->tag.handle,name,&err,sizeof(name));
    if(checkError(self->mesh->mesh,err))
        return NULL;

    return Py_BuildValue("s",name);
}

static PyObject *
iMeshTagObj_getSizeValues(iMeshTag_Object *self,void *closure)
{
    int size,err;
    iMesh_getTagSizeValues(self->mesh->mesh,self->tag.handle,&size,&err);
    if(checkError(self->mesh->mesh,err))
        return NULL;

    return Py_BuildValue("i",size);
}

static PyObject *
iMeshTagObj_getSizeBytes(iMeshTag_Object *self,void *closure)
{
    int size,err;
    iMesh_getTagSizeBytes(self->mesh->mesh,self->tag.handle,&size,&err);
    if(checkError(self->mesh->mesh,err))
        return NULL;

    return Py_BuildValue("i",size);
}

static PyObject *
iMeshTagObj_getType(iMeshTag_Object *self,void *closure)
{
    int type,err;
    iMesh_getTagType(self->mesh->mesh,self->tag.handle,&type,&err);
    if(checkError(self->mesh->mesh,err))
        return NULL;

    return Py_BuildValue("c",type_to_char(type));
}

static PyObject *
iMeshTagObj_setData(iMeshTag_Object *self,PyObject *args)
{
    PyObject *obj;
    PyObject *data_obj;
    char typechar=0;
    int type;
    int err;

    if(!PyArg_ParseTuple(args,"OO|c",&obj,&data_obj,&typechar))
        return NULL;

    if(typechar == 0)
    {
        /* infer the type of the data */
        iMesh_getTagType(self->mesh->mesh,self->tag.handle,&type,&err);
        if(checkError(self->mesh->mesh,err))
            return NULL;
    }
    else
    {
        type = char_to_type(typechar);
        if(type == -1)
        {
            PyErr_SetString(PyExc_ValueError,ERR_TYPE_CODE);
            return NULL;
        }
    }
 
    PyObject *ents = PyArray_TryFromObject(obj,NPY_IBASEENT,1,1);
    if(ents)
    {
        int ent_size;
        iBase_EntityHandle *entities;
        int data_size;
        PyObject *data_arr=0;

        ent_size = PyArray_SIZE(ents);
        entities = PyArray_DATA(ents);

        if(type == iBase_INTEGER)
        {
            data_arr = PyArray_FROMANY(data_obj,NPY_INT,1,1,NPY_C_CONTIGUOUS);
            if(data_arr == NULL)
                return NULL;

            data_size = PyArray_SIZE(data_arr);
            int *data = PyArray_DATA(data_arr);
            iMesh_setIntArrData(self->mesh->mesh,entities,ent_size,
                                self->tag.handle,data,data_size,&err);
        }
        else if(type == iBase_DOUBLE)
        {
            data_arr = PyArray_FROMANY(data_obj,NPY_DOUBLE,1,1,
                                       NPY_C_CONTIGUOUS);
            if(data_arr == NULL)
                return NULL;

            data_size = PyArray_SIZE(data_arr);
            double *data = PyArray_DATA(data_arr);
            iMesh_setDblArrData(self->mesh->mesh,entities,ent_size,
                                self->tag.handle,data,data_size,&err);
        }
        else if(type == iBase_ENTITY_HANDLE)
        {
            data_arr = PyArray_FROMANY(data_obj,NPY_IBASEENT,1,1,
                                       NPY_C_CONTIGUOUS);
            if(data_arr == NULL)
                return NULL;

            data_size = PyArray_SIZE(data_arr);
            iBase_EntityHandle *data = PyArray_DATA(data_arr);
            iMesh_setEHArrData(self->mesh->mesh,entities,ent_size,
                               self->tag.handle,data,data_size,&err);
        }
        else /* iBase_BYTES */
        {
            data_arr = PyArray_FROMANY(data_obj,NPY_BYTE,1,1,NPY_C_CONTIGUOUS);
            if(data_arr == NULL)
                return NULL;

            data_size = PyArray_SIZE(data_arr);
            char *data = PyArray_DATA(data_arr);
            iMesh_setArrData(self->mesh->mesh,entities,ent_size,
                             self->tag.handle,data,data_size,&err);
        }

        Py_DECREF(ents);
        Py_XDECREF(data_arr);
    }
    else if(iBaseEntitySet_Check(obj))
    {
        iBase_EntitySetHandle set = iBaseEntitySet_GetHandle(obj);

        if(type == iBase_INTEGER)
        {
            PyObject *o = PyNumber_Int(data_obj);
            if(o == NULL)
                return NULL;
            iMesh_setEntSetIntData(self->mesh->mesh,set,self->tag.handle,
                                   PyInt_AsLong(o),&err);
            Py_DECREF(o);
        }
        else if(type == iBase_DOUBLE)
        {
            PyObject *o = PyNumber_Float(data_obj);
            if(o == NULL)
                return NULL;
            iMesh_setEntSetDblData(self->mesh->mesh,set,self->tag.handle,
                                   PyFloat_AsDouble(o),&err);
            Py_DECREF(o);
        }
        else if(type == iBase_ENTITY_HANDLE)
        {
            if(!iBaseEntity_Check(data_obj))
                return NULL;
            iMesh_setEntSetEHData(self->mesh->mesh,set,self->tag.handle,
                                  iBaseEntity_GetHandle(data_obj),&err);
        }
        else /* iBase_BYTES */
        {
            PyObject *data_arr = PyArray_FROMANY(data_obj,NPY_BYTE,1,1,
                                                 NPY_C_CONTIGUOUS);
            if(data_arr == NULL)
                return NULL;

            char *data = PyArray_DATA(data_arr);
            int data_size = PyArray_SIZE(data_arr);
            iMesh_setEntSetData(self->mesh->mesh,set,self->tag.handle,data,
                                data_size,&err);
            Py_DECREF(data_arr);
        }
    }
    else if(iBaseEntity_Check(obj))
    {
        iBase_EntityHandle entity = iBaseEntity_GetHandle(obj);

        if(type == iBase_INTEGER)
        {
            PyObject *o = PyNumber_Int(data_obj);
            if(o == NULL)
                return NULL;
            iMesh_setIntData(self->mesh->mesh,entity,self->tag.handle,
                             PyInt_AsLong(o),&err);
            Py_DECREF(o);
        }
        else if(type == iBase_DOUBLE)
        {
            PyObject *o = PyNumber_Float(data_obj);
            if(o == NULL)
                return NULL;
            iMesh_setDblData(self->mesh->mesh,entity,self->tag.handle,
                             PyFloat_AsDouble(o),&err);
            Py_DECREF(o);
        }
        else if(type == iBase_ENTITY_HANDLE)
        {
            if(!iBaseEntity_Check(data_obj))
                return NULL;
            iMesh_setEHData(self->mesh->mesh,entity,self->tag.handle,
                            iBaseEntity_GetHandle(data_obj),&err);
        }
        else /* iBase_BYTES */
        {
            PyObject *data_arr = PyArray_FROMANY(data_obj,NPY_BYTE,1,1,
                                                 NPY_C_CONTIGUOUS);
            if(data_arr == NULL)
                return NULL;

            char *data = PyArray_DATA(data_arr);
            int data_size = PyArray_SIZE(data_arr);
            iMesh_setData(self->mesh->mesh,entity,self->tag.handle,data,
                          data_size,&err);
            Py_DECREF(data_arr);
        }
    }
    else
    {
        PyErr_SetString(PyExc_ValueError,ERR_ANY_ENT);
        return NULL;
    }

    if(checkError(self->mesh->mesh,err))
        return NULL;
    Py_RETURN_NONE;
}

static PyObject *
iMeshTagObj_getData(iMeshTag_Object *self,PyObject *args)
{
    PyObject *obj;
    char typechar=0;
    int type;
    int err;

    if(!PyArg_ParseTuple(args,"O|c",&obj,&typechar))
        return NULL;

    if(typechar == 0)
    {
        /* infer the type of the data */
        iMesh_getTagType(self->mesh->mesh,self->tag.handle,&type,&err);
        if(checkError(self->mesh->mesh,err))
            return NULL;
    }
    else
    {
        type = char_to_type(typechar);
        if(type == -1)
        {
            PyErr_SetString(PyExc_ValueError,"invalid type code");
            return NULL;
        }
    }

    PyObject *ents = PyArray_TryFromObject(obj,NPY_IBASEENT,1,1);
    if(ents)
    {
        int ent_size = PyArray_SIZE(ents);
        iBase_EntityHandle *entities = PyArray_DATA(ents);
        PyObject *ret = 0;

        if(type == iBase_INTEGER)
        {
            int *data=0;
            int alloc=0,size;

            iMesh_getIntArrData(self->mesh->mesh,entities,ent_size,
                                self->tag.handle,&data,&alloc,&size,&err);
            if(!checkError(self->mesh->mesh,err))
            {
                npy_intp dims[] = {size};
                ret = PyArray_NewFromMalloc(1,dims,NPY_INT,data);
            }
        }
        else if(type == iBase_DOUBLE)
        {
            double *data=0;
            int alloc=0,size;

            iMesh_getDblArrData(self->mesh->mesh,entities,ent_size,
                                self->tag.handle,&data,&alloc,&size,&err);
            if(!checkError(self->mesh->mesh,err))
            {
                npy_intp dims[] = {size};
                ret = PyArray_NewFromMalloc(1,dims,NPY_DOUBLE,data);
            }
        }
        else if(type == iBase_ENTITY_HANDLE)
        {
            iBase_EntityHandle *data=0;
            int alloc=0,size;

            iMesh_getEHArrData(self->mesh->mesh,entities,ent_size,
                               self->tag.handle,&data,&alloc,&size,&err);
            if(!checkError(self->mesh->mesh,err))
            {
                npy_intp dims[] = {size};
                ret = PyArray_NewFromMalloc(1,dims,NPY_IBASEENT,data);
            }
        }
        else /* iBase_BYTES */
        {
            char *data=0;
            int alloc=0,size;

            iMesh_getArrData(self->mesh->mesh,entities,ent_size,
                             self->tag.handle,&data,&alloc,&size,&err);
            if(!checkError(self->mesh->mesh,err))
            {
                npy_intp dims[] = {size};
                ret = PyArray_NewFromMalloc(1,dims,NPY_BYTE,data);
            }
        }

        Py_DECREF(ents);
        return ret;
    }
    else if(iBaseEntitySet_Check(obj))
    {
        iBase_EntitySetHandle set = iBaseEntitySet_GetHandle(obj);

        if(type == iBase_INTEGER)
        {
            int data;
            iMesh_getEntSetIntData(self->mesh->mesh,set,self->tag.handle,&data,
                                   &err);
            if(checkError(self->mesh->mesh,err))
                return NULL;
            return Py_BuildValue("i",data);
        }
        else if(type == iBase_DOUBLE)
        {
            double data;
            iMesh_getEntSetDblData(self->mesh->mesh,set,self->tag.handle,&data,
                                   &err);
            if(checkError(self->mesh->mesh,err))
                return NULL;
            return Py_BuildValue("d",data);
        }
        else if(type == iBase_ENTITY_HANDLE)
        {
            iBaseEntity_Object *data = iBaseEntity_New();
            iMesh_getEntSetEHData(self->mesh->mesh,set,self->tag.handle,
                                  &data->handle,&err);
            if(checkError(self->mesh->mesh,err))
            {
                Py_DECREF(data);
                return NULL;
            }
            return (PyObject*)data;
        }
        else /* iBase_BYTES */
        {
            char *data=0;
            int alloc=0,size;
            iMesh_getEntSetData(self->mesh->mesh,set,self->tag.handle,&data,
                                &alloc,&size,&err);
            if(checkError(self->mesh->mesh,err))
                return NULL;
            npy_intp dims[] = {size};
            return PyArray_NewFromMalloc(1,dims,NPY_BYTE,data);
        }
    }
    else if(iBaseEntity_Check(obj))
    {
        iBase_EntityHandle entity = iBaseEntity_GetHandle(obj);

        if(type == iBase_INTEGER)
        {
            int data;
            iMesh_getIntData(self->mesh->mesh,entity,self->tag.handle,&data,
                             &err);
            if(checkError(self->mesh->mesh,err))
                return NULL;
            return Py_BuildValue("i",data);
        }
        else if(type == iBase_DOUBLE)
        {
            double data;
            iMesh_getDblData(self->mesh->mesh,entity,self->tag.handle,&data,
                             &err);
            if(checkError(self->mesh->mesh,err))
                return NULL;
            return Py_BuildValue("d",data);
        }
        else if(type == iBase_ENTITY_HANDLE)
        {
            iBaseEntity_Object *data = iBaseEntity_New();
            iMesh_getEHData(self->mesh->mesh,entity,self->tag.handle,
                            &data->handle,&err);
            if(checkError(self->mesh->mesh,err))
            {
                Py_DECREF(data);
                return NULL;
            }
            return (PyObject*)data;
        }
        else /* iBase_BYTES */
        {
            char *data=0;
            int alloc=0,size;
            iMesh_getData(self->mesh->mesh,entity,self->tag.handle,&data,
                          &alloc,&size,&err);
            if(checkError(self->mesh->mesh,err))
                return NULL;
            npy_intp dims[] = {size};
            return PyArray_NewFromMalloc(1,dims,NPY_BYTE,data);
        }
    }
    else
    {
        PyErr_SetString(PyExc_ValueError,ERR_ANY_ENT);
        return NULL;
    }
}

static PyObject *
iMeshTagObj_remove(iMeshTag_Object *self,PyObject *args)
{
    PyObject *obj;
    int err;

    if(!PyArg_ParseTuple(args,"O",&obj))
        return NULL;

    PyObject *ents = PyArray_TryFromObject(obj,NPY_IBASEENT,1,1);
    if(ents)
    {
        int ent_size = PyArray_SIZE(ents);
        iBase_EntityHandle *entities = PyArray_DATA(ents);
        iMesh_rmvArrTag(self->mesh->mesh,entities,ent_size,self->tag.handle,
                        &err);
        Py_DECREF(ents);
    }
    else if(iBaseEntitySet_Check(obj))
    {
        iBase_EntitySetHandle set = iBaseEntitySet_GetHandle(obj);
        iMesh_rmvEntSetTag(self->mesh->mesh,set,self->tag.handle,&err);
    }
    else if(iBaseEntity_Check(obj))
    {
        iBase_EntityHandle entity = iBaseEntity_GetHandle(obj);
        iMesh_rmvTag(self->mesh->mesh,entity,self->tag.handle,&err);
    }
    else
    {
        PyErr_SetString(PyExc_ValueError,ERR_ANY_ENT);
        return NULL;
    }

    if(checkError(self->mesh->mesh,err))
        return NULL;
    Py_RETURN_NONE;
}


static PyMethodDef iMeshTagObj_methods[] = {
    { "setData", (PyCFunction)iMeshTagObj_setData, METH_VARARGS,
      "Set tag values on an entity (or array/set of entities)"
    },
    { "getData", (PyCFunction)iMeshTagObj_getData, METH_VARARGS,
      "Get tag values on an entity (or array/set of entities)"
    },
    { "remove", (PyCFunction)iMeshTagObj_remove, METH_VARARGS,
      "Remove a tag value from an entity (or array/set of entities)"
    },
    {0}
};

static PyMemberDef iMeshTagObj_members[] = {
    {"instance", T_OBJECT_EX, offsetof(iMeshTag_Object, mesh), READONLY,
     "base iMesh instance"},
    {0}
};

static PyGetSetDef iMeshTagObj_getset[] = {
    { "name", (getter)iMeshTagObj_getName, 0,
      "Get the name for a given tag handle", 0
    },
    { "sizeValues", (getter)iMeshTagObj_getSizeValues, 0,
      "Get size of a tag in units of numbers of tag data type", 0
    },
    { "sizeBytes", (getter)iMeshTagObj_getSizeBytes, 0,
      "Get size of a tag in units of bytes", 0
    },
    { "type", (getter)iMeshTagObj_getType, 0,
      "Get the data type of the specified tag handle", 0
    },
    {0}
};


static PyTypeObject iMeshTag_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                                   /* ob_size */
    "itaps.iMesh.Tag",                   /* tp_name */
    sizeof(iMeshTag_Object),             /* tp_basicsize */
    0,                                   /* tp_itemsize */
    (destructor)iMeshTagObj_dealloc,     /* tp_dealloc */
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
    Py_TPFLAGS_BASETYPE,                 /* tp_flags */
    "iMesh tag object",                  /* tp_doc */
    0,                                   /* tp_traverse */
    0,                                   /* tp_clear */
    0,                                   /* tp_richcompare */
    0,                                   /* tp_weaklistoffset */
    0,                                   /* tp_iter */
    0,                                   /* tp_iternext */
    iMeshTagObj_methods,                 /* tp_methods */
    iMeshTagObj_members,                 /* tp_members */
    iMeshTagObj_getset,                  /* tp_getset */
    0,                                   /* tp_base */
    0,                                   /* tp_dict */
    0,                                   /* tp_descr_get */
    0,                                   /* tp_descr_set */
    0,                                   /* tp_dictoffset */
    0,                                   /* tp_init */
    0,                                   /* tp_alloc */
    0,                                   /* tp_new */
};

