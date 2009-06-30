static PyObject *
PyArray_TryFromObject(PyObject *obj,int typenum,int min_depth,int max_depth)
{
    PyObject *ret = PyArray_FromAny(obj,PyArray_DescrFromType(typenum),
                                    min_depth,max_depth,NPY_C_CONTIGUOUS,NULL);
    PyErr_Clear();
    return ret;
}

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
