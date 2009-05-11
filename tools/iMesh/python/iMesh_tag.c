#define NO_IMPORT_ARRAY
#define NO_IMPORT_IBASE

#include "iMesh_Python.h"
#include "iBase_Python.h"

static char typechars[] = {'i','d','E','b'};

enum iBase_TagValueType
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

char
type_to_char(enum iBase_TagValueType t)
{
    return typechars[t];
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


static PyMethodDef iMeshTagObj_methods[] = {
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


PyTypeObject iMeshTag_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                                   /* ob_size */
    "itaps.iMesh.tag",                   /* tp_name */
    sizeof(iMeshTag_Object),             /* tp_basicsize */
    0,                                   /* tp_itemsize */
    0,                                   /* tp_dealloc */
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
    0,                                   /* tp_members */
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

