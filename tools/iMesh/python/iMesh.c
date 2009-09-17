#define _IMESH_MODULE
#include "iMesh_Python.h"
#include "errors.h"
#include "common.h"
#include "helpers.h"
#include "numpy_extensions.h"

static enum iBase_TagValueType char_to_type(char c);
static char type_to_char(enum iBase_TagValueType t);

static PyTypeObject iMesh_Type;
static PyTypeObject iMeshIter_Type;
static PyTypeObject iMeshEntitySet_Type;
static int NPY_IMESHENTSET;
static PyTypeObject iMeshTag_Type;
static int NPY_IMESHTAG;

static int
checkError(iMesh_Instance mesh,int err)
{
    if(err)
    {
        char descr[512];
        iMesh_getDescription(mesh,descr,&err,sizeof(descr));

        PyErr_SetString(PyExc_RuntimeError,descr);
        return 1;
    }
    else
        return 0;
}

static iMeshEntitySet_Object *
iMeshEntitySet_New(iMesh_Object *instance)
{
    iMeshEntitySet_Object *o = iMeshEntitySet_NewRaw();
    o->instance = instance;
    Py_INCREF(o->instance);
    return o;
}

static iMeshTag_Object *
iMeshTag_New(iMesh_Object *instance)
{
    iMeshTag_Object *o = iMeshTag_NewRaw();
    o->instance = instance;
    Py_INCREF(o->instance);
    return o;
}

static int
iMeshObj_init(iMesh_Object *self,PyObject *args,PyObject *kwds)
{
    static char *kwlist[] = {"options",0};
    const char *options = "";

    if( !PyArg_ParseTupleAndKeywords(args,kwds,"|s",kwlist,&options) )
        return -1;

    int err;
    iMesh_newMesh(options,&self->handle,&err,strlen(options));
    if(checkError(self->handle,err))
        return -1;
    return 0;
}

static void
iMeshObj_dealloc(iMesh_Object *self)
{
    if(self->handle)
    {
        int err;
        iMesh_dtor(self->handle,&err);
    }
    self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
iMeshObj_getRootSet(iMesh_Object *self,void *closure)
{
    iMeshEntitySet_Object *rootset = iMeshEntitySet_New(self);

    int err;
    iMesh_getRootSet(self->handle,&rootset->base.handle,&err);
    if(checkError(self->handle,err))
    {
        Py_DECREF((PyObject*)rootset);
        return NULL;
    }

    return (PyObject*)rootset;
}


static PyObject *
iMeshObj_getGeometricDimension(iMesh_Object *self,void *closure)
{
    int dim,err;
    iMesh_getGeometricDimension(self->handle,&dim,&err);
    if(checkError(self->handle,err))
        return NULL;

    return PyInt_FromLong(dim);
}

static int
iMeshObj_setGeometricDimension(iMesh_Object *self,PyObject *value,void *closure)
{
    if(value == NULL)
    {
        PyErr_SetString(PyExc_TypeError, "Cannot delete the "
                        "geometricDimension attribute");
        return -1;
    }
  
    int dim,err;
    if(!PyArg_Parse(value,"i",&dim))
        return -1;
    iMesh_setGeometricDimension(self->handle,dim,&err);
    if(checkError(self->handle,err))
        return -1;

    return 0;
}

static PyObject *
iMeshObj_getDfltStorage(iMesh_Object *self,void *closure)
{
    int order,err;
    iMesh_getDfltStorage(self->handle,&order,&err);
    if(checkError(self->handle,err))
        return NULL;

    return Py_BuildValue("i",order);
}

static PyObject *
iMeshObj_getAdjTable(iMesh_Object *self,void *closure)
{
    int *adjtable=0;
    int alloc=0,size,err;

    iMesh_getAdjTable(self->handle,&adjtable,&alloc,&size,&err);
    if(checkError(self->handle,err))
        return NULL;

    npy_intp dims[] = {4,4};
    return PyArray_NewFromMalloc(2,dims,NPY_INT,adjtable);
}

static PyObject *
iMeshObj_areEHValid(iMesh_Object *self,PyObject *args)
{
    int doReset,areInv,err;
    if(!PyArg_ParseTuple(args,"i",&doReset))
        return NULL;

    iMesh_areEHValid(self->handle,doReset,&areInv,&err);
    if(checkError(self->handle,err))
        return NULL;

    return Py_BuildValue("i",areInv);
}

static PyObject *
iMeshObj_getVtxCoords(iMesh_Object *self,PyObject *args)
{
    PyObject *obj;
    int storage_order=-1;
    int err;

    if(!PyArg_ParseTuple(args,"O|i",&obj,&storage_order))
        return NULL;

    PyObject *ents = PyArray_TryFromObject(obj,NPY_IBASEENT,1,1);
    if(ents)
    {
        if(storage_order == -1)
        {
            Py_DECREF(ents);
            PyErr_SetString(PyExc_ValueError,ERR_STORAGE_ORDER);
            return NULL;
        }

        int size;
        double *coords=0;
        int coords_alloc=0,coords_size;

        size = PyArray_SIZE(ents);
        iBase_EntityHandle *data = PyArray_DATA(ents);

        iMesh_getVtxArrCoords(self->handle,data,size,storage_order,&coords,
                              &coords_alloc,&coords_size,&err);
        Py_DECREF(ents);

        if(checkError(self->handle,err))
            return NULL;

        npy_intp outer = (storage_order == iBase_BLOCKED) ? 3:size;
        /* TODO: think about this */
        npy_intp dims[] = {outer, coords_size/outer};
        return PyArray_NewFromMalloc(2,dims,NPY_DOUBLE,coords);
    }
    else if(iBaseEntity_Check(obj))
    {
        double *v = malloc(3*sizeof(double));
        iMesh_getVtxCoord(self->handle,iBaseEntity_GetHandle(obj), v+0,v+1,v+2,
                          &err);
        if(checkError(self->handle,err))
        {
            free(v);
            return NULL;
        }

        npy_intp dims[] = {3};
        return PyArray_NewFromMalloc(1,dims,NPY_DOUBLE,v);
    }
    else
    {
        PyErr_SetString(PyExc_ValueError,ERR_ENT_OR_ENTARR);
        return NULL;
    }
}

static PyObject *
iMeshObj_getEntTopo(iMesh_Object *self,PyObject *args)
{
    PyObject *obj;
    int err;

    if(!PyArg_ParseTuple(args,"O",&obj))
        return NULL;

    PyObject *ents = PyArray_TryFromObject(obj,NPY_IBASEENT,1,1);
    if(ents)
    {
        int size;
        iBase_EntityHandle *data;
        int *topos=0;
        int topo_alloc=0,topo_size;

        size = PyArray_SIZE(ents);
        data = PyArray_DATA(ents);

        iMesh_getEntArrTopo(self->handle,data,size,&topos,&topo_alloc,
                            &topo_size,&err);
        Py_DECREF(ents);
        if(checkError(self->handle,err))
            return NULL;

        npy_intp dims[] = {topo_size};
        return PyArray_NewFromMalloc(1,dims,NPY_UINT,topos);
    }
    else if(iBaseEntity_Check(obj))
    {
        int topo;
        iBase_EntityHandle handle = ((iBaseEntity_Object*)obj)->handle;

        iMesh_getEntTopo(self->handle,handle,&topo,&err);
        if(checkError(self->handle,err))
            return NULL;

        return PyInt_FromLong(topo);
    }
    else
    {
        PyErr_SetString(PyExc_ValueError,ERR_ENT_OR_ENTARR);
        return NULL;
    }
}

static PyObject *
iMeshObj_getEntType(iMesh_Object *self,PyObject *args)
{
    PyObject *obj;
    int err;

    if(!PyArg_ParseTuple(args,"O",&obj))
        return NULL;

    PyObject *ents = PyArray_TryFromObject(obj,NPY_IBASEENT,1,1);
    if(ents)
    {
        int size;
        iBase_EntityHandle *data;
        int *types=0;
        int type_alloc=0,type_size;

        size = PyArray_SIZE(ents);
        data = PyArray_DATA(ents);
      
        iMesh_getEntArrType(self->handle,data,size,&types,&type_alloc,
                            &type_size,&err);
        Py_DECREF(ents);
        if(checkError(self->handle,err))
            return NULL;
    
        npy_intp dims[] = {type_size};
        return PyArray_NewFromMalloc(1,dims,NPY_UINT,types);
    }
    else if(iBaseEntity_Check(obj))
    {
        int type;
        iBase_EntityHandle handle = ((iBaseEntity_Object*)obj)->handle;
        iMesh_getEntType(self->handle,handle,&type,&err);
        if(checkError(self->handle,err))
            return NULL;
    
        return PyInt_FromLong(type);
    }
    else
    {
        PyErr_SetString(PyExc_ValueError,ERR_ENT_OR_ENTARR);
        return NULL;
    }
}

static PyObject *
iMeshObj_getEntAdj(iMesh_Object *self,PyObject *args)
{
    PyObject *obj;
    int type_req;
    int err;

    if(!PyArg_ParseTuple(args,"Oi",&obj,&type_req))
        return NULL;

    PyObject *ents = PyArray_TryFromObject(obj,NPY_IBASEENT,1,1);
    if(ents)
    {
        int size;
        iBase_EntityHandle *data;
        iBase_EntityHandle *adj=0;
        int adj_alloc=0,adj_size;
        int *offsets;
        int offsets_alloc=0,offsets_size;
    
        size = PyArray_SIZE(ents);
        data = PyArray_DATA(ents);

        iMesh_getEntArrAdj(self->handle,data,size,type_req,&adj,&adj_alloc,
                           &adj_size,&offsets,&offsets_alloc,&offsets_size,
                           &err);
        Py_DECREF(ents);
        if(checkError(self->handle,err))
            return NULL;

        npy_intp adj_dims[] = {adj_size};
        npy_intp off_dims[] = {offsets_size};

        return AdjacencyList_New(
            PyArray_NewFromMalloc(1,adj_dims,NPY_IBASEENT,adj),
            PyArray_NewFromMalloc(1,off_dims,NPY_INT,offsets) );
    }
    else if(iBaseEntity_Check(obj))
    {
        iBase_EntityHandle *adj=0;
        int adj_alloc=0,adj_size;
        iBase_EntityHandle handle = iBaseEntity_GetHandle(obj);

        iMesh_getEntAdj(self->handle,handle,type_req,&adj,&adj_alloc,&adj_size,
                        &err);
        if(checkError(self->handle,err))
            return NULL;

        npy_intp dims[] = {adj_size};
        return PyArray_NewFromMalloc(1,dims,NPY_IBASEENT,adj);
    }
    else
    {
        PyErr_SetString(PyExc_ValueError,ERR_ENT_OR_ENTARR);
        return NULL;
    }
}

static PyObject *
iMeshObj_getEnt2ndAdj(iMesh_Object *self,PyObject *args)
{
    PyObject *obj;
    int bridge_type,type_req;
    int err;

    if(!PyArg_ParseTuple(args,"Oii",&obj,&bridge_type,&type_req))
        return NULL;

    PyObject *ents = PyArray_TryFromObject(obj,NPY_IBASEENT,1,1);
    if(ents)
    {
        int size;
        iBase_EntityHandle *data;
        iBase_EntityHandle *adj=0;
        int adj_alloc=0,adj_size;
        int *offsets;
        int offsets_alloc=0,offsets_size;

        size = PyArray_SIZE(ents);
        data = PyArray_DATA(ents);

        iMesh_getEntArr2ndAdj(self->handle,data,size,bridge_type,type_req,&adj,
                              &adj_alloc,&adj_size,&offsets,&offsets_alloc,
                              &offsets_size,&err);
        Py_DECREF(ents);
        if(checkError(self->handle,err))
            return NULL;

        npy_intp adj_dims[] = {adj_size};
        npy_intp off_dims[] = {offsets_size};

        return AdjacencyList_New(
            PyArray_NewFromMalloc(1,adj_dims,NPY_IBASEENT,adj),
            PyArray_NewFromMalloc(1,off_dims,NPY_INT,offsets) );
    }
    else if(iBaseEntity_Check(obj))
    {
        iBase_EntityHandle *adj=0;
        int adj_alloc=0,adj_size;
        iBase_EntityHandle handle = iBaseEntity_GetHandle(obj);

        iMesh_getEnt2ndAdj(self->handle,handle,bridge_type,type_req,&adj,
                           &adj_alloc,&adj_size,&err);
        if(checkError(self->handle,err))
            return NULL;

        npy_intp dims[] = {adj_size};
        return PyArray_NewFromMalloc(1,dims,NPY_IBASEENT,adj);
    }
    else
    {
        PyErr_SetString(PyExc_ValueError,ERR_ENT_OR_ENTARR);
        return NULL;
    }
}

static PyObject *
iMeshObj_createEntSet(iMesh_Object *self,PyObject *args)
{
    int isList,err;
    PyObject *obj;
    iMeshEntitySet_Object *set;

    if(!PyArg_ParseTuple(args,"O!",&PyBool_Type,&obj))
        return NULL;

    set = iMeshEntitySet_New(self);

    isList = (obj == Py_True);
  
    iMesh_createEntSet(self->handle,isList,&set->base.handle,&err);
    if(checkError(self->handle,err))
    {
        Py_DECREF((PyObject*)set);
        return NULL;
    }

    return (PyObject*)set;  
}

static PyObject *
iMeshObj_destroyEntSet(iMesh_Object *self,PyObject *args)
{
    int err;
    iBaseEntitySet_Object *set;

    if(!PyArg_ParseTuple(args,"O!",&iBaseEntitySet_Type,&set))
        return NULL;

    iMesh_destroyEntSet(self->handle,set->handle,&err);
    if(checkError(self->handle,err))
        return NULL;

    Py_RETURN_NONE;
}

static PyObject *
iMeshObj_setVtxCoords(iMesh_Object *self,PyObject *args)
{
    PyObject *obj;
    int storage_order = -1;
    PyObject *data;
    PyObject *ents = 0;
    PyObject *verts = 0;
    int err;

    if(!PyArg_ParseTuple(args,"OO|i",&obj,&data,&storage_order))
        return NULL;

    ents = PyArray_TryFromObject(obj,NPY_IBASEENT,1,1);
    if(ents)
    {
        if(storage_order == -1)
        {
            Py_DECREF(ents);
            PyErr_SetString(PyExc_ValueError,ERR_STORAGE_ORDER);
            goto err;
        }

        verts = PyArray_ToVectors(data,NPY_DOUBLE,2,3,
                                  storage_order==iBase_INTERLEAVED);
        if(verts == NULL)
            goto err;

        int ent_size = PyArray_SIZE(ents);
        iBase_EntityHandle *entities = PyArray_DATA(ents);
        int coord_size = PyArray_SIZE(verts);
        double *coords = PyArray_DATA(verts);

        iMesh_setVtxArrCoords(self->handle,entities,ent_size,storage_order,
                              coords,coord_size,&err);
        Py_DECREF(ents);
        Py_DECREF(verts);
    }
    else if(iBaseEntity_Check(obj))
    {
        verts = PyArray_ToVectors(data,NPY_DOUBLE,1,3,0);
        if(verts == NULL)
            goto err;

        double *v = PyArray_DATA(verts);
        iBase_EntityHandle entity = iBaseEntity_GetHandle(obj);

        iMesh_setVtxCoord(self->handle,entity, v[0],v[1],v[2], &err);
        Py_DECREF(verts);
    }
    else
    {
        PyErr_SetString(PyExc_ValueError,ERR_ENT_OR_ENTARR);
        return NULL;
    }

    if(checkError(self->handle,err))
        return NULL;
    Py_RETURN_NONE;

err:
    Py_XDECREF(ents);
    Py_XDECREF(verts);
    return NULL;
}

static PyObject *
iMeshObj_createVtx(iMesh_Object *self,PyObject *args)
{
    int storage_order=-1;
    PyObject *data;
    PyObject *vertices;
    int err;

    if(!PyArg_ParseTuple(args,"O|i",&data,&storage_order))
        return NULL;

    if( (vertices = PyArray_TryFromObject(data,NPY_DOUBLE,2,2)) != NULL)
    {
        if(storage_order == -1)
        {
            Py_DECREF(vertices);
            PyErr_SetString(PyExc_ValueError,ERR_STORAGE_ORDER);
            return NULL;
        }

        int index = storage_order == iBase_BLOCKED;
        int count = PyArray_DIM(vertices,index); /* this is a bit odd! */
        int coord_size = PyArray_SIZE(vertices);
        double *coords = PyArray_DATA(vertices);
        iBase_EntityHandle *entities=0;
        int ent_alloc=0,ent_size;

        iMesh_createVtxArr(self->handle,count,storage_order,coords,coord_size,
                           &entities,&ent_alloc,&ent_size,&err);
        Py_DECREF(vertices);
        if(checkError(self->handle,err))
            return NULL;

        npy_intp dims[] = {ent_size};
        return PyArray_NewFromMalloc(1,dims,NPY_IBASEENT,entities);
    }
    else if( (vertices = PyArray_TryFromObject(data,NPY_DOUBLE,1,1)) != NULL)
    {
        if(PyArray_SIZE(vertices) != 3)
        {
            Py_DECREF(vertices);
            PyErr_SetString(PyExc_ValueError,ERR_ARR_SIZE);
            return NULL;
        }

        double *v = PyArray_DATA(vertices);

        iBaseEntity_Object *entity = iBaseEntity_New();
        iMesh_createVtx(self->handle, v[0],v[1],v[2], &entity->handle,&err);
        Py_DECREF(vertices);

        if(checkError(self->handle,err))
        {
            Py_DECREF((PyObject*)entity);
            return NULL;
        }
        return (PyObject*)entity;
    }
    else
    {
        PyErr_SetString(PyExc_ValueError,ERR_ARR_DIMS);
        return NULL;
    }
}

static PyObject *
iMeshObj_createEnt(iMesh_Object *self,PyObject *args)
{
    int topo,status,err;
    PyObject *obj;

    if(!PyArg_ParseTuple(args,"iO",&topo,&obj))
        return NULL;

    PyObject *ents = PyArray_FROMANY(obj,NPY_IBASEENT,1,1,NPY_C_CONTIGUOUS);
    if(ents == NULL)
        return NULL;

    int lower_size = PyArray_SIZE(ents);
    iBase_EntityHandle *lower = PyArray_DATA(ents);

    iBaseEntity_Object *entity = iBaseEntity_New();

    iMesh_createEnt(self->handle,topo,lower,lower_size,&entity->handle,&status,
                    &err);
    Py_DECREF(ents);
    if(checkError(self->handle,err))
    {
        Py_DECREF((PyObject*)entity);
        return NULL;
    }

    return Py_BuildValue("(Oi)",entity,status);
}

static PyObject *
iMeshObj_createEntArr(iMesh_Object *self,PyObject *args)
{
    int topo,err;
    PyObject *obj;

    if(!PyArg_ParseTuple(args,"iO",&topo,&obj))
        return NULL;

    PyObject *ents = PyArray_FROMANY(obj,NPY_IBASEENT,1,1,NPY_C_CONTIGUOUS);
    if(ents == NULL)
        return NULL;

    int lower_size = PyArray_SIZE(ents);
    iBase_EntityHandle *lower = PyArray_DATA(ents);

    iBase_EntityHandle *entities=0;
    int ent_alloc=0,ent_size;
    int *status;
    int stat_alloc=0,stat_size;

    iMesh_createEntArr(self->handle,topo,lower,lower_size,&entities,&ent_alloc,
                       &ent_size,&status,&stat_alloc,&stat_size,&err);
    Py_DECREF(ents);
    if(checkError(self->handle,err))
        return NULL;

    npy_intp ent_dims[] = {ent_size};
    npy_intp stat_dims[] = {stat_size};
    return Py_BuildValue("(OO)",
        PyArray_NewFromMalloc(1,ent_dims,NPY_IBASEENT,entities),
        PyArray_NewFromMalloc(1,stat_dims,NPY_INT,status)
        );
}


static PyObject *
iMeshObj_deleteEnt(iMesh_Object *self,PyObject *args)
{
    PyObject *obj;
    int err;

    if(!PyArg_ParseTuple(args,"O",&obj))
        return NULL;

    PyObject *ents = PyArray_TryFromObject(obj,NPY_IBASEENT,1,1);
    if(ents)
    {
        int size = PyArray_SIZE(ents);
        iBase_EntityHandle *entities = PyArray_DATA(ents);
        iMesh_deleteEntArr(self->handle,entities,size,&err);
        Py_DECREF(ents);
    }
    else if(iBaseEntity_Check(obj))
    {
        iBase_EntityHandle entity = iBaseEntity_GetHandle(obj);
        iMesh_deleteEnt(self->handle,entity,&err);
    }
    else
    {
        PyErr_SetString(PyExc_ValueError,ERR_ENT_OR_ENTARR);
        return NULL;
    }

    if(checkError(self->handle,err))
        return NULL;
    Py_RETURN_NONE;
}


static PyObject *
iMeshObj_createTag(iMesh_Object *self,PyObject *args)
{
    char *name;
    char typechar;
    int size,type,err;
    iMeshTag_Object *tag;

    if(!PyArg_ParseTuple(args,"sic",&name,&size,&typechar))
        return NULL;

    type = char_to_type(typechar);
    if(type == -1)
    {
        PyErr_SetString(PyExc_ValueError,ERR_TYPE_CODE);
        return NULL;
    }

    tag = iMeshTag_New(self);

    iMesh_createTag(self->handle,name,size,type,&tag->base.handle,&err,
                    strlen(name));
    if(checkError(self->handle,err))
    {
        Py_DECREF((PyObject*)tag);
        return NULL;
    }

    return (PyObject*)tag;
}

static PyObject *
iMeshObj_destroyTag(iMesh_Object *self,PyObject *args)
{
    int forced,err;
    iBaseTag_Object *tag;
    PyObject *obj;

    if(!PyArg_ParseTuple(args,"O!O!",&iBaseTag_Type,&tag,&PyBool_Type,&obj))
        return NULL;

    forced = (obj == Py_True);

    iMesh_destroyTag(self->handle,tag->handle,forced,&err);
    if(checkError(self->handle,err))
        return NULL;

    Py_RETURN_NONE;
}

static PyObject *
iMeshObj_getTagHandle(iMesh_Object *self,PyObject *args)
{
    char *name;
    iMeshTag_Object *tag;
    int err;

    if(!PyArg_ParseTuple(args,"s",&name))
        return NULL;

    tag = iMeshTag_New(self);

    iMesh_getTagHandle(self->handle,name,&tag->base.handle,&err,strlen(name));
    if(checkError(self->handle,err))
    {
        Py_DECREF((PyObject*)tag);
        return NULL;
    }

    return (PyObject*)tag;
}

static PyObject *
iMeshObj_getAllTags(iMesh_Object *self,PyObject *args)
{
    PyObject *ents;
    iBase_TagHandle *tags=0;
    int alloc=0,size;
    int err;

    if(!PyArg_ParseTuple(args,"O",&ents))
        return NULL;

    if(iBaseEntitySet_Check(ents))
    {
        iBase_EntitySetHandle set = iBaseEntitySet_GetHandle(ents);

        iMesh_getAllEntSetTags(self->handle,set,&tags,&alloc,&size,&err);
        if(checkError(self->handle,err))
            return NULL;
    }
    else if(iBaseEntity_Check(ents))
    {
        iBase_EntityHandle entity = iBaseEntity_GetHandle(ents);

        iMesh_getAllTags(self->handle,entity,&tags,&alloc,&size,&err);
        if(checkError(self->handle,err))
            return NULL;
    }
    else
    {
        PyErr_SetString(PyExc_ValueError,ERR_ENT_OR_ENTSET);
        return NULL;
    }

    npy_intp dims[] = {size};
    return PyArray_NewFromMallocBase(1,dims,NPY_IMESHTAG,tags,
                                     (PyObject*)self);
}


static PyMethodDef iMeshObj_methods[] = {
    { "areEHValid", (PyCFunction)iMeshObj_areEHValid, METH_VARARGS,
      "Return whether entity handles have changed since last reset or since "
      "instance construction"
    },
    { "getVtxCoords", (PyCFunction)iMeshObj_getVtxCoords, METH_VARARGS,
      "Get coordinates of specified vertex(ices)"
    },
    { "getEntTopo", (PyCFunction)iMeshObj_getEntTopo, METH_VARARGS,
      "Get the entity topology(ies) for the specified entity(ies)"
    },
    { "getEntType", (PyCFunction)iMeshObj_getEntType, METH_VARARGS,
      "Get the entity type(s) for the specified entity(ies)"
    },
    { "getEntAdj", (PyCFunction)iMeshObj_getEntAdj, METH_VARARGS,
      "Get entities of specified type adjacent to entity(ies)"
    },
    { "getEnt2ndAdj", (PyCFunction)iMeshObj_getEnt2ndAdj, METH_VARARGS,
      "Get \"2nd order\" adjacencies to entity(ies)"
    },
    { "createEntSet", (PyCFunction)iMeshObj_createEntSet, METH_VARARGS,
      "Create an entity set"
    },
    { "destroyEntSet", (PyCFunction)iMeshObj_destroyEntSet, METH_VARARGS,
      "Destroy an entity set"
    },
    { "setVtxCoords", (PyCFunction)iMeshObj_setVtxCoords, METH_VARARGS,
      "Set coordinates for a vertex or array of vertices"
    },
    { "createVtx", (PyCFunction)iMeshObj_createVtx, METH_VARARGS,
      "Create a new vertex or array of vertices at specified coordinates"
    },
    { "createEnt", (PyCFunction)iMeshObj_createEnt, METH_VARARGS,
      "Create a new entity with specified lower-order topology"
    },
    { "createEntArr", (PyCFunction)iMeshObj_createEntArr, METH_VARARGS,
      "Create an array of entities with specified lower-order topology"
    },
    { "deleteEnt", (PyCFunction)iMeshObj_deleteEnt, METH_VARARGS,
      "Delete specified entity(ies)"
    },
    { "createTag", (PyCFunction)iMeshObj_createTag, METH_VARARGS,
      "Create a tag with specified name, size, and type"
    },
    { "destroyTag", (PyCFunction)iMeshObj_destroyTag, METH_VARARGS,
      "Destroy a tag"
    },
    { "getTagHandle", (PyCFunction)iMeshObj_getTagHandle, METH_VARARGS,
      "Get the handle of an existing tag with the specified name"
    },
    { "getAllTags", (PyCFunction)iMeshObj_getAllTags, METH_VARARGS,
      "Get all the tags associated with a specified entity handle (or "
      "array/set of entities)"
    },
  {0}
};

static PyGetSetDef iMeshObj_getset[] = {
    { "rootSet", (getter)iMeshObj_getRootSet, 0,
      "root set", 0
    },
    { "geometricDimension", (getter)iMeshObj_getGeometricDimension,
      (setter)iMeshObj_setGeometricDimension,
      "geometric dimension", 0
    },
    { "dfltStorage",(getter)iMeshObj_getDfltStorage, 0,
      "default storage order", 0
    },
    { "adjTable",(getter)iMeshObj_getAdjTable, 0,
      "get adjacency table", 0
    },
    {0}
};

static PyObject * iMeshObj_getAttr(PyObject *self,PyObject *attr_name)
{
    PyObject *ret;

    ret = PyObject_GenericGetAttr(self,attr_name);
    if(ret)
        return ret;
    else
    {
        PyErr_Clear();
        PyObject *root = iMeshObj_getRootSet((iMesh_Object*)self,0);
        if(!root)
            return NULL;
        ret = PyObject_GetAttr(root,attr_name);
        Py_DECREF(root);
        return ret;
    }
}

static PyTypeObject iMesh_Type = {
    PyObject_HEAD_INIT(NULL)
    0,                            /* ob_size */
    "itaps.iMesh.Mesh",           /* tp_name */
    sizeof(iMesh_Object),         /* tp_basicsize */
    0,                            /* tp_itemsize */
    (destructor)iMeshObj_dealloc, /* tp_dealloc */
    0,                            /* tp_print */
    0,                            /* tp_getattr */
    0,                            /* tp_setattr */
    0,                            /* tp_compare */
    0,                            /* tp_repr */
    0,                            /* tp_as_number */
    0,                            /* tp_as_sequence */
    0,                            /* tp_as_mapping */
    0,                            /* tp_hash */
    0,                            /* tp_call */
    0,                            /* tp_str */
    iMeshObj_getAttr,             /* tp_getattro */
    0,                            /* tp_setattro */
    0,                            /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT |
    Py_TPFLAGS_BASETYPE,          /* tp_flags */
    "iMesh objects",              /* tp_doc */
    0,                            /* tp_traverse */
    0,                            /* tp_clear */
    0,                            /* tp_richcompare */
    0,                            /* tp_weaklistoffset */
    0,                            /* tp_iter */
    0,                            /* tp_iternext */
    iMeshObj_methods,             /* tp_methods */
    0,                            /* tp_members */
    iMeshObj_getset,              /* tp_getset */
    0,                            /* tp_base */
    0,                            /* tp_dict */
    0,                            /* tp_descr_get */
    0,                            /* tp_descr_set */
    0,                            /* tp_dictoffset */
    (initproc)iMeshObj_init,      /* tp_init */
    0,                            /* tp_alloc */
    0,                            /* tp_new */
};


static PyMethodDef module_methods[] = {
    {0}
};

static PyObject *
iMeshEntSetArr_getitem(void *data,void *arr)
{
    ArrDealloc_Object *b = (ArrDealloc_Object*)PyArray_BASE(arr);
    iMesh_Object *instance = (iMesh_Object*)b->base;
    iMeshEntitySet_Object *o = iMeshEntitySet_New(instance);

    o->base.handle = *(iBase_EntitySetHandle*)data;

    return (PyObject*)o;
}

static PyObject *
iMeshTagArr_getitem(void *data,void *arr)
{
    ArrDealloc_Object *b = (ArrDealloc_Object*)PyArray_BASE(arr);
    iMesh_Object *instance = (iMesh_Object*)b->base;
    iMeshTag_Object *o = iMeshTag_New(instance);

    o->base.handle = *(iBase_TagHandle*)data;

    return (PyObject*)o;
}

static PyArray_ArrFuncs iMeshEntSetArr_Funcs;
static int NPY_IMESHENTSET;

static PyArray_ArrFuncs iMeshTagArr_Funcs;
static int NPY_IMESHTAG;

ENUM_TYPE(iMeshTopology,"iMesh.Topology","");


PyMODINIT_FUNC initiMesh(void)
{
    PyObject *m;
    PyArray_Descr *descr;

    ArrDealloc_Type.tp_new = PyType_GenericNew;
    if (PyType_Ready(&ArrDealloc_Type) < 0)
        return;

    m = Py_InitModule("iMesh",module_methods);
    import_array();
    import_iBase();
    import_helpers();

    /***** register C API *****/
    static void *IMesh_API[6];
    PyObject *api_obj;

    /* Initialize the C API pointer array */
    IMesh_API[0] = &iMesh_Type;
    IMesh_API[1] = &iMeshIter_Type;
    IMesh_API[2] = &iMeshEntitySet_Type;
    IMesh_API[3] = &NPY_IMESHENTSET;
    IMesh_API[4] = &iMeshTag_Type;
    IMesh_API[5] = &NPY_IMESHTAG;

    /* Create a CObject containing the API pointer array's address */
    api_obj = PyCObject_FromVoidPtr(IMesh_API,NULL);

    if(api_obj != NULL)
        PyModule_AddObject(m, "_C_API", api_obj);

    REGISTER_SIMPLE(m,"Mesh",iMesh);

    /***** initialize topology enum *****/
    REGISTER_SIMPLE(m,"Topology",iMeshTopology);

    ADD_ENUM(iMeshTopology,"point",         iMesh_POINT);
    ADD_ENUM(iMeshTopology,"line_segment",  iMesh_LINE_SEGMENT);
    ADD_ENUM(iMeshTopology,"polygon",       iMesh_POLYGON);
    ADD_ENUM(iMeshTopology,"triangle",      iMesh_TRIANGLE);
    ADD_ENUM(iMeshTopology,"quadrilateral", iMesh_QUADRILATERAL);
    ADD_ENUM(iMeshTopology,"polyhedron",    iMesh_POLYHEDRON);
    ADD_ENUM(iMeshTopology,"tetrahedron",   iMesh_TETRAHEDRON);
    ADD_ENUM(iMeshTopology,"hexahedron",    iMesh_HEXAHEDRON);
    ADD_ENUM(iMeshTopology,"prism",         iMesh_PRISM);
    ADD_ENUM(iMeshTopology,"pyramid",       iMesh_PYRAMID);
    ADD_ENUM(iMeshTopology,"septahedron",   iMesh_SEPTAHEDRON);
    ADD_ENUM(iMeshTopology,"all",           iMesh_ALL_TOPOLOGIES);

    /***** initialize iterator type *****/
    iMeshIter_Type.tp_new = PyType_GenericNew;
    if(PyType_Ready(&iMeshIter_Type) < 0)
        return;
    PyModule_AddObject(m,"Iterator",(PyObject *)&iMeshIter_Type);

    /***** initialize entity set type *****/
    iMeshEntitySet_Type.tp_base = &iBaseEntitySet_Type;
    if(PyType_Ready(&iMeshEntitySet_Type) < 0)
        return;
    PyModule_AddObject(m,"EntitySet",(PyObject *)&iMeshEntitySet_Type);

    /***** initialize tag type *****/
    iMeshTag_Type.tp_base = &iBaseTag_Type;
    if(PyType_Ready(&iMeshTag_Type) < 0)
        return;
    PyModule_AddObject(m,"Tag",(PyObject *)&iMeshTag_Type);


    /***** initialize iMeshEntitySet array *****/
    descr = PyArray_DescrNewFromType(NPY_IBASEENTSET);
    memcpy(&iMeshEntSetArr_Funcs,descr->f,sizeof(PyArray_ArrFuncs));
    descr->f = &iMeshEntSetArr_Funcs;

    descr->typeobj = &iMeshEntitySet_Type;
    descr->type = 'M';
    descr->f->getitem = iMeshEntSetArr_getitem;

    NPY_IMESHENTSET = PyArray_RegisterDataType(descr);

    /***** initialize iMeshTag array *****/
    descr = PyArray_DescrNewFromType(NPY_IBASETAG);
    memcpy(&iMeshTagArr_Funcs,descr->f,sizeof(PyArray_ArrFuncs));
    descr->f = &iMeshTagArr_Funcs;

    descr->typeobj = &iMeshTag_Type;
    descr->type = 'M';
    descr->f->getitem = iMeshTagArr_getitem;

    NPY_IMESHTAG = PyArray_RegisterDataType(descr);
}

/* Include source files so that everything is in one translation unit */
#include "iMesh_entSet.inl"
#include "iMesh_iter.inl"
#include "iMesh_tag.inl"
