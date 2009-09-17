#ifndef PYTAPS_HELPERS_H
#define PYTAPS_HELPERS_H

/* TODO: these are never freed! */
static PyObject *g_helper_module;
static PyObject *g_adj_list;
static PyObject *g_ind_adj_list;

/* NOTE: steals references to adj and offsets */
static PyObject *
AdjacencyList_New(PyObject *adj,PyObject *offsets)
{
    PyObject *res;

    if( (res = PyObject_CallFunction(g_adj_list,"OO",adj,offsets)) == NULL)
        PyErr_SetString(PyExc_RuntimeError,ERR_ADJ_LIST);

    Py_DECREF(adj);
    Py_DECREF(offsets);

    return res;
}

/* NOTE: steals references to ents, adj, indices, and offsets */
static PyObject *
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

static int import_helpers(void)
{
    if( (g_helper_module = PyImport_ImportModule("itaps.helpers")) == NULL)
        return -1;
    if( (g_adj_list = PyObject_GetAttrString(g_helper_module,"AdjacencyList") )
        == NULL)
        return -1;
    if( (g_ind_adj_list = PyObject_GetAttrString(g_helper_module,
        "IndexedAdjacencyList")) == NULL)
        return -1;

    return 0;

    /* Suppress warnings if above functions aren't used */
    (void)AdjacencyList_New;
    (void)IndexedAdjacencyList_New;

}

#endif
