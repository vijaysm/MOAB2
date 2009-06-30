#ifndef PYTAPS_IMESH_PYTHON_H
#define PYTAPS_IMESH_PYTHON_H

#include "common.h"

#include <Python.h>
#include <iMesh.h>

#define PY_IBASE_UNIQUE_SYMBOL itaps_IBASE_API
#include "iBase_Python.h"

#define PY_ARRAY_UNIQUE_SYMBOL itaps_ARRAY_API
#include "numpy_extensions.h"

static int checkError(iMesh_Instance mesh,int err);
static enum iBase_TagValueType char_to_type(char c);
static char type_to_char(enum iBase_TagValueType t);

static PyObject *
AdjacencyList_New(PyObject *adj,PyObject *offsets);
static PyObject *
IndexedAdjacencyList_New(PyObject *ents, PyObject *adj,PyObject *indices,
                         PyObject *offsets);

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

static PyTypeObject iMeshIter_Type;

typedef struct
{
    iBaseEntitySet_Object set;
    iMesh_Object *mesh;
} iMeshEntitySet_Object;

static PyTypeObject iMeshEntitySet_Type;
static int NPY_IMESHENTSET;

static iMeshEntitySet_Object * iMeshEntitySet_New(iMesh_Object *mesh);

#define iMeshEntitySet_NewRaw()                         \
    (iMeshEntitySet_Object*)PyObject_CallObject(        \
        (PyObject*)&iMeshEntitySet_Type,NULL)

#define iMeshEntitySet_Check(o)                         \
    PyObject_TypeCheck((o),&iMeshEntitySet_Type)

#define iMeshEntitySet_GetMesh(o)                       \
    ((iMeshEntitySet_Object*)(o))->mesh

typedef struct
{
    iBaseTag_Object tag;
    iMesh_Object *mesh;
} iMeshTag_Object;

static PyTypeObject iMeshTag_Type;
static int NPY_IMESHTAG;

static iMeshTag_Object * iMeshTag_New(iMesh_Object *mesh);

#define iMeshTag_NewRaw()                               \
    (iMeshTag_Object*)PyObject_CallObject(              \
        (PyObject*)&iMeshTag_Type,NULL)

#define iMeshTag_Check(o)                               \
    PyObject_TypeCheck((o),&iMeshTag_Type)

#define iMeshTag_GetMesh(o)                             \
    ((iMeshTag_Object*)(o))->mesh

#endif
