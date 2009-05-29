#include "errors.h"
#include "iMesh_Python.h"
#include "iBase_Python.h"

int checkError(iMesh_Instance mesh,int err)
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

PyObject *
PyArray_TryFromObject(PyObject *obj,int typenum,int min_depth,int max_depth)
{
    PyObject *ret = PyArray_FromAny(obj,PyArray_DescrFromType(typenum),
                                    min_depth,max_depth,NPY_C_CONTIGUOUS,NULL);
    PyErr_Clear();
    return ret;
}

/* TODO: these are never freed! */
static PyObject *g_helper_module;
static PyObject *g_adj_list;
static PyObject *g_ind_adj_list;

/* NOTE: steals references to adj and offsets */
PyObject *
AdjacencyList_New(PyObject *adj,PyObject *offsets)
{
    PyObject *res;

    if( (res = PyObject_CallFunction(g_adj_list,"OO",adj,offsets)) == NULL)
        PyErr_SetString(PyExc_RuntimeError,ERR_ADJ_LIST);

    Py_DECREF(adj);
    Py_DECREF(offsets);

    return res;
}

/* NOTE: steals references to adj and offsets */
PyObject *
IndexedAdjacencyList_New(PyObject *ents, PyObject *adj,PyObject *indices,
                         PyObject *offsets)
{
    PyObject *res;

    if( (res = PyObject_CallFunction(g_ind_adj_list,"OOOO",ents,adj,indices,
                                     offsets)) == NULL)
        PyErr_SetString(PyExc_RuntimeError,ERR_ADJ_LIST);

    Py_DECREF(ents);
    Py_DECREF(adj);
    Py_DECREF(indices);
    Py_DECREF(offsets);

    return res;
}

static int
iMeshObj_init(iMeshObject *self,PyObject *args,PyObject *kwds)
{
    static char *kwlist[] = {"options",0};
    const char *options = "";

    if( !PyArg_ParseTupleAndKeywords(args,kwds,"|s",kwlist,&options) )
        return -1;

    int err;
    iMesh_newMesh(options,&self->mesh,&err,strlen(options));
    if(checkError(self->mesh,err))
        return -1;

    return 0;
}

static void
iMeshObj_dealloc(iMeshObject *self)
{
    if(self->mesh)
    {
        int err;
        iMesh_dtor(self->mesh,&err);
    }

    self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
iMeshObj_load(iMeshObject *self,PyObject *args)
{
    const char *name = 0;
    const char *options = "";
    iBase_EntitySetHandle root;
    int err;

    if(!PyArg_ParseTuple(args,"s|s",&name,&options))
        return NULL;

    iMesh_getRootSet(self->mesh,&root,&err);
    if(checkError(self->mesh,err))
        return NULL;

    iMesh_load(self->mesh,root,name,options,&err,strlen(name),strlen(options));
    if(checkError(self->mesh,err))
        return NULL;

    Py_RETURN_NONE;
}

static PyObject *
iMeshObj_save(iMeshObject *self,PyObject *args)
{
    const char *name = 0;
    const char *options = "";
    iBase_EntitySetHandle root;
    int err;

    if(!PyArg_ParseTuple(args,"s|s",&name,&options))
        return NULL;

    iMesh_getRootSet(self->mesh,&root,&err);
    if(checkError(self->mesh,err))
        return NULL;

    iMesh_save(self->mesh,root,name,options,&err,strlen(name),strlen(options));
    if(checkError(self->mesh,err))
        return NULL;

    Py_RETURN_NONE;
}

static PyObject *
iMeshObj_getRootSet(iMeshObject *self,void *closure)
{
    iMeshEntitySet_Object *rootset = iMeshEntitySet_New();
    rootset->mesh = self; /* TODO: incref? */

    int err;
    iMesh_getRootSet(self->mesh,&rootset->set.handle,&err);
    if(checkError(self->mesh,err))
    {
        Py_DECREF((PyObject*)rootset);
        return NULL;
    }

    return (PyObject*)rootset;
}


static PyObject *
iMeshObj_getGeometricDimension(iMeshObject *self,void *closure)
{
    int dim,err;
    iMesh_getGeometricDimension(self->mesh,&dim,&err);
    if(checkError(self->mesh,err))
        return NULL;

    return Py_BuildValue("i",dim);
}

static int
iMeshObj_setGeometricDimension(iMeshObject *self,PyObject *value,void *closure)
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
    iMesh_setGeometricDimension(self->mesh,dim,&err);
    if(checkError(self->mesh,err))
        return -1;

    return 0;
}

static PyObject *
iMeshObj_getDfltStorage(iMeshObject *self,void *closure)
{
    int order,err;
    iMesh_getDfltStorage(self->mesh,&order,&err);
    if(checkError(self->mesh,err))
        return NULL;

    return Py_BuildValue("i",order);
}

static PyObject *
iMeshObj_getAdjTable(iMeshObject *self,void *closure)
{
    int *adjtable=0;
    int alloc=0,size,err;

    iMesh_getAdjTable(self->mesh,&adjtable,&alloc,&size,&err);
    if(checkError(self->mesh,err))
        return NULL;

    npy_intp dims[] = {4,4};
    return PyArray_NewFromMalloc(2,dims,NPY_INT,adjtable);
}

static PyObject *
iMeshObj_areEHValid(iMeshObject *self,PyObject *args)
{
    int doReset,areInv,err;
    if(!PyArg_ParseTuple(args,"i",&doReset))
        return NULL;

    iMesh_areEHValid(self->mesh,doReset,&areInv,&err);
    if(checkError(self->mesh,err))
        return NULL;

    return Py_BuildValue("i",areInv);
}

static PyObject *
iMeshObj_getVtxCoords(iMeshObject *self,PyObject *args)
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
            PyErr_SetString(PyExc_ValueError,ERR_STORAGE_ORDER);
            return NULL;
        }

        int size;
        double *coords=0;
        int coords_alloc=0,coords_size;

        size = PyArray_SIZE(ents);
        iBase_EntityHandle *data = PyArray_DATA(ents);

        iMesh_getVtxArrCoords(self->mesh,data,size,storage_order,&coords,
                              &coords_alloc,&coords_size,&err);
        Py_DECREF(ents);

        if(checkError(self->mesh,err))
            return NULL;

        npy_intp outer;
        if(storage_order == iBase_BLOCKED)
            outer = 3;
        else
            outer = size;
        /* TODO: think about this */
        npy_intp dims[] = {outer, coords_size/outer};
        return PyArray_NewFromMalloc(2,dims,NPY_DOUBLE,coords);
    }
    else if(iBaseEntity_Check(obj))
    {
        double *v = malloc(3*sizeof(double));
        iMesh_getVtxCoord(self->mesh,iBaseEntity_GetHandle(obj), v+0,v+1,v+2,
                          &err);
        if(checkError(self->mesh,err))
            return NULL;

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
iMeshObj_getEntTopo(iMeshObject *self,PyObject *args)
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

        iMesh_getEntArrTopo(self->mesh,data,size,&topos,&topo_alloc,&topo_size,
                            &err);
        Py_DECREF(ents);
        if(checkError(self->mesh,err))
            return NULL;

        npy_intp dims[] = {topo_size};
        return PyArray_NewFromMalloc(1,dims,NPY_UINT,topos);
    }
    else if(iBaseEntity_Check(obj))
    {
        int topo;
        iBase_EntityHandle handle = ((iBaseEntity_Object*)obj)->handle;

        iMesh_getEntTopo(self->mesh,handle,&topo,&err);
        if(checkError(self->mesh,err))
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
iMeshObj_getEntType(iMeshObject *self,PyObject *args)
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
      
        iMesh_getEntArrType(self->mesh,data,size,&types,&type_alloc,&type_size,
                            &err);
        Py_DECREF(ents);
        if(checkError(self->mesh,err))
            return NULL;
    
        npy_intp dims[] = {type_size};
        return PyArray_NewFromMalloc(1,dims,NPY_UINT,types);
    }
    else if(iBaseEntity_Check(obj))
    {
        int type;
        iBase_EntityHandle handle = ((iBaseEntity_Object*)obj)->handle;
        iMesh_getEntType(self->mesh,handle,&type,&err);
        if(checkError(self->mesh,err))
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
iMeshObj_getEntAdj(iMeshObject *self,PyObject *args)
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

        iMesh_getEntArrAdj(self->mesh,data,size,type_req,&adj,&adj_alloc,
                           &adj_size,&offsets,&offsets_alloc,&offsets_size,
                           &err);
        Py_DECREF(ents);
        if(checkError(self->mesh,err))
            return NULL;

        npy_intp adj_dims[] = {adj_size};
        npy_intp off_dims[] = {offsets_size};

        return AdjacencyList_New(
            PyArray_NewFromMalloc(1,adj_dims,NPY_IBASEENT,adj),
            PyArray_NewFromMalloc(1,off_dims,NPY_INT,offsets) );

        PyObject *pair = PyTuple_New(2);
        npy_intp dims[1];
 
        dims[0] = adj_size;
        PyTuple_SET_ITEM(pair, 0,
            PyArray_NewFromMalloc(1,dims,NPY_IBASEENT,adj));

        dims[0] = offsets_size;
        PyTuple_SET_ITEM(pair, 1,
            PyArray_NewFromMalloc(1,dims,NPY_INT,offsets));

        return pair;
    }
    else if(iBaseEntity_Check(obj))
    {
        iBase_EntityHandle *adj=0;
        int adj_alloc=0,adj_size;
        iBase_EntityHandle handle = iBaseEntity_GetHandle(obj);

        iMesh_getEntAdj(self->mesh,handle,type_req,&adj,&adj_alloc,&adj_size,
                        &err);
        if(checkError(self->mesh,err))
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
iMeshObj_getEnt2ndAdj(iMeshObject *self,PyObject *args)
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

        iMesh_getEntArr2ndAdj(self->mesh,data,size,bridge_type,type_req,&adj,
                              &adj_alloc,&adj_size,&offsets,&offsets_alloc,
                              &offsets_size,&err);
        Py_DECREF(ents);
        if(checkError(self->mesh,err))
            return NULL;

        npy_intp adj_dims[] = {adj_size};
        npy_intp off_dims[] = {offsets_size};

        return AdjacencyList_New(
            PyArray_NewFromMalloc(1,adj_dims,NPY_IBASEENT,adj),
            PyArray_NewFromMalloc(1,off_dims,NPY_INT,offsets) );

        PyObject *pair = PyTuple_New(2);
        npy_intp dims[1];
 
        dims[0] = adj_size;
        PyTuple_SET_ITEM(pair, 0,
            PyArray_NewFromMalloc(1,dims,NPY_IBASEENT,adj));

        dims[0] = offsets_size;
        PyTuple_SET_ITEM(pair, 1,
            PyArray_NewFromMalloc(1,dims,NPY_INT,offsets));

        return pair;
    }
    else if(iBaseEntity_Check(obj))
    {
        iBase_EntityHandle *adj=0;
        int adj_alloc=0,adj_size;
        iBase_EntityHandle handle = iBaseEntity_GetHandle(obj);

        iMesh_getEnt2ndAdj(self->mesh,handle,bridge_type,type_req,&adj,
                           &adj_alloc,&adj_size,&err);
        if(checkError(self->mesh,err))
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
iMeshObj_createEntSet(iMeshObject *self,PyObject *args)
{
    int isList,err;
    PyObject *obj;
    iMeshEntitySet_Object *set;

    if(!PyArg_ParseTuple(args,"O!",&PyBool_Type,&obj))
        return NULL;

    set = iMeshEntitySet_New();
    set->mesh = self;
    /*Py_INCREF(self); TODO?? */

    isList = (obj == Py_True);
  
    iMesh_createEntSet(self->mesh,isList,&set->set.handle,&err);
    if(checkError(self->mesh,err))
    {
        Py_DECREF((PyObject*)set);
        return NULL;
    }

    return (PyObject*)set;  
}

static PyObject *
iMeshObj_destroyEntSet(iMeshObject *self,PyObject *args)
{
    int err;
    iBaseEntitySet_Object *set;

    if(!PyArg_ParseTuple(args,"O!",&iBaseEntitySet_Type,&set))
        return NULL;

    iMesh_destroyEntSet(self->mesh,set->handle,&err);
    if(checkError(self->mesh,err))
        return NULL;

    Py_RETURN_NONE;
}

static PyObject *
iMeshObj_setVtxCoords(iMeshObject *self,PyObject *args)
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
            PyErr_SetString(PyExc_ValueError,ERR_STORAGE_ORDER);
            goto err;
        }

        verts = PyArray_FROMANY(data,NPY_DOUBLE,2,2,NPY_C_CONTIGUOUS);
        if(verts == NULL)
            goto err;

        int ent_size = PyArray_SIZE(ents);
        iBase_EntityHandle *entities = PyArray_DATA(ents);
        int coord_size = PyArray_SIZE(verts);
        double *coords = PyArray_DATA(verts);

        iMesh_setVtxArrCoords(self->mesh,entities,ent_size,storage_order,
                              coords,coord_size,&err);
        Py_DECREF(ents);
        Py_DECREF(verts);
    }
    else if(iBaseEntity_Check(obj))
    {
        verts = PyArray_FROMANY(data,NPY_DOUBLE,1,1,NPY_C_CONTIGUOUS);
        if(verts == NULL)
            goto err;

        if(PyArray_SIZE(verts) != 3)
        {
            PyErr_SetString(PyExc_ValueError,ERR_ARR_SIZE);
            goto err;
        }

        double *v = PyArray_DATA(verts);
        iBase_EntityHandle entity = iBaseEntity_GetHandle(obj);

        iMesh_setVtxCoord(self->mesh,entity, v[0],v[1],v[2], &err);
        Py_DECREF(verts);
    }
    else
    {
        PyErr_SetString(PyExc_ValueError,ERR_ENT_OR_ENTARR);
        return NULL;
    }

    if(checkError(self->mesh,err))
        return NULL;
    Py_RETURN_NONE;

err:
    Py_XDECREF(ents);
    Py_XDECREF(verts);
    return NULL;
}

static PyObject *
iMeshObj_createVtx(iMeshObject *self,PyObject *args)
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

        iMesh_createVtxArr(self->mesh,count,storage_order,coords,coord_size,
                           &entities,&ent_alloc,&ent_size,&err);
        Py_DECREF(vertices);
        if(checkError(self->mesh,err))
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
        iMesh_createVtx(self->mesh, v[0],v[1],v[2], &entity->handle,&err);
        Py_DECREF(vertices);

        if(checkError(self->mesh,err))
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
iMeshObj_createEnt(iMeshObject *self,PyObject *args)
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

    iMesh_createEnt(self->mesh,topo,lower,lower_size,&entity->handle,&status,
                    &err);
    Py_DECREF(ents);
    if(checkError(self->mesh,err))
    {
        Py_DECREF((PyObject*)entity);
        return NULL;
    }

    PyObject *pair = PyTuple_New(2);
    PyTuple_SET_ITEM(pair,0,(PyObject*)entity);
    PyTuple_SET_ITEM(pair,1,Py_BuildValue("i",status));
    return pair;
}

static PyObject *
iMeshObj_createEntArr(iMeshObject *self,PyObject *args)
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

    iMesh_createEntArr(self->mesh,topo,lower,lower_size,&entities,&ent_alloc,
                       &ent_size,&status,&stat_alloc,&stat_size,&err);
    Py_DECREF(ents);
    if(checkError(self->mesh,err))
        return NULL;

    PyObject *pair = PyTuple_New(2);
    npy_intp dims[1];

    dims[0] = ent_size;
    PyTuple_SET_ITEM(pair, 0,
        PyArray_NewFromMalloc(1,dims,NPY_IBASEENT,entities));

    dims[0] = stat_size;
    PyTuple_SET_ITEM(pair, 1,
        PyArray_NewFromMalloc(1,dims,NPY_INT,status));

    return pair;
}


static PyObject *
iMeshObj_deleteEnt(iMeshObject *self,PyObject *args)
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
        iMesh_deleteEntArr(self->mesh,entities,size,&err);
        Py_DECREF(ents);
    }
    else if(iBaseEntity_Check(obj))
    {
        iBase_EntityHandle entity = iBaseEntity_GetHandle(obj);
        iMesh_deleteEnt(self->mesh,entity,&err);
    }
    else
    {
        PyErr_SetString(PyExc_ValueError,ERR_ENT_OR_ENTARR);
        return NULL;
    }

    if(checkError(self->mesh,err))
        return NULL;
    Py_RETURN_NONE;
}


static PyObject *
iMeshObj_createTag(iMeshObject *self,PyObject *args)
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

    tag = iMeshTag_New();
    tag->mesh = self;
    /*Py_INCREF(self); TODO?? */

    iMesh_createTag(self->mesh,name,size,type,&tag->tag.handle,&err,
                    strlen(name));
    if(checkError(self->mesh,err))
    {
        Py_DECREF((PyObject*)tag);
        return NULL;
    }

    return (PyObject*)tag;
}

static PyObject *
iMeshObj_destroyTag(iMeshObject *self,PyObject *args)
{
    int forced,err;
    iBaseTag_Object *tag;
    PyObject *obj;

    if(!PyArg_ParseTuple(args,"O!O!",&iBaseTag_Type,&tag,&PyBool_Type,&obj))
        return NULL;

    forced = (obj == Py_True);

    iMesh_destroyTag(self->mesh,tag->handle,forced,&err);
    if(checkError(self->mesh,err))
        return NULL;

    Py_RETURN_NONE;
}

static PyObject *
iMeshObj_getTagHandle(iMeshObject *self,PyObject *args)
{
    char *name;
    iMeshTag_Object *tag;
    int err;

    if(!PyArg_ParseTuple(args,"s",&name))
        return NULL;

    tag = iMeshTag_New();
    tag->mesh = self;
    /*Py_INCREF(self); TODO?? */

    iMesh_getTagHandle(self->mesh,name,&tag->tag.handle,&err,strlen(name));
    if(checkError(self->mesh,err))
    {
        Py_DECREF((PyObject*)tag);
        return NULL;
    }

    return (PyObject*)tag;
}

static PyObject *
iMeshObj_getAllTags(iMeshObject *self,PyObject *args)
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

        iMesh_getAllEntSetTags(self->mesh,set,&tags,&alloc,&size,&err);
        if(checkError(self->mesh,err))
            return NULL;
    }
    else if(iBaseEntity_Check(ents))
    {
        iBase_EntityHandle entity = iBaseEntity_GetHandle(ents);

        iMesh_getAllTags(self->mesh,entity,&tags,&alloc,&size,&err);
        if(checkError(self->mesh,err))
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


static PyMethodDef iMesh_methods[] = {
    { "load", (PyCFunction)iMeshObj_load, METH_VARARGS,
      "Load a mesh from a file"
    },
    { "save", (PyCFunction)iMeshObj_save, METH_VARARGS,
      "Save the mesh to a file"
    },
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

static PyGetSetDef iMesh_getset[] = {
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

static PyTypeObject iMeshType = {
    PyObject_HEAD_INIT(NULL)
    0,                            /* ob_size */
    "itaps.iMesh.Mesh",           /* tp_name */
    sizeof(iMeshObject),          /* tp_basicsize */
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
    0,                            /* tp_getattro */
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
    iMesh_methods,                /* tp_methods */
    0,                            /* tp_members */
    iMesh_getset,                 /* tp_getset */
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
    iMeshObject *mesh = (iMeshObject*)b->base;
    iMeshEntitySet_Object *o = iMeshEntitySet_New();

    o->set.handle = *(iBase_EntitySetHandle*)data;
    o->mesh = mesh; /* TODO: incref? */

    return (PyObject*)o;
}

static PyObject *
iMeshTagArr_getitem(void *data,void *arr)
{
    ArrDealloc_Object *b = (ArrDealloc_Object*)PyArray_BASE(arr);
    iMeshObject *mesh = (iMeshObject*)b->base;
    iMeshTag_Object *o = iMeshTag_New();

    o->tag.handle = *(iBase_TagHandle*)data;
    o->mesh = mesh; /* TODO: incref? */

    return (PyObject*)o;
}

static PyArray_ArrFuncs iMeshEntSetArr_Funcs;
int NPY_IMESHENTSET;

static PyArray_ArrFuncs iMeshTagArr_Funcs;
int NPY_IMESHTAG;

ENUM_TYPE(Topology,"iMesh.Topology","");


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

PyObject *
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

    /***** import helper module *****/
    if( (g_helper_module = PyImport_ImportModule("itaps.helpers")) == NULL)
        return;
    if( (g_adj_list = PyObject_GetAttrString(g_helper_module,"AdjacencyList") )
        == NULL)
        return;
    if( (g_ind_adj_list = PyObject_GetAttrString(g_helper_module,
        "IndexedAdjacencyList")) == NULL)
        return;

    iMeshType.tp_new = PyType_GenericNew;
    if(PyType_Ready(&iMeshType) < 0)
        return;
    Py_INCREF(&iMeshType);
    PyModule_AddObject(m,"Mesh",(PyObject *)&iMeshType);

    /***** initialize topology enum *****/
    REGISTER_SIMPLE(m,Topology);

    ADD_ENUM(Topology,"point",         iMesh_POINT);
    ADD_ENUM(Topology,"line_segment",  iMesh_LINE_SEGMENT);
    ADD_ENUM(Topology,"polygon",       iMesh_POLYGON);
    ADD_ENUM(Topology,"triangle",      iMesh_TRIANGLE);
    ADD_ENUM(Topology,"quadrilateral", iMesh_QUADRILATERAL);
    ADD_ENUM(Topology,"polyhedron",    iMesh_POLYHEDRON);
    ADD_ENUM(Topology,"tetrahedron",   iMesh_TETRAHEDRON);
    ADD_ENUM(Topology,"hexahedron",    iMesh_HEXAHEDRON);
    ADD_ENUM(Topology,"prism",         iMesh_PRISM);
    ADD_ENUM(Topology,"pyramid",       iMesh_PYRAMID);
    ADD_ENUM(Topology,"septahedron",   iMesh_SEPTAHEDRON);
    ADD_ENUM(Topology,"all",           iMesh_ALL_TOPOLOGIES);

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
