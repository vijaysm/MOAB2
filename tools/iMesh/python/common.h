#ifndef PYTAPS_COMMON_H
#define PYTAPS_COMMON_H

#include <Python.h>

#define SIMPLE_TYPE(name,name_str,doc)                          \
    static PyTypeObject name ## _Type = {                       \
        PyObject_HEAD_INIT(NULL)                                \
        0,                         /* ob_size*/                 \
        (name_str),                /* tp_name*/                 \
        sizeof(name ## _Object),   /* tp_basicsize */           \
        0,                         /* tp_itemsize */            \
        0,                         /* tp_dealloc */             \
        0,                         /* tp_print */               \
        0,                         /* tp_getattr */             \
        0,                         /* tp_setattr */             \
        0,                         /* tp_compare */             \
        0,                         /* tp_repr */                \
        0,                         /* tp_as_number */           \
        0,                         /* tp_as_sequence */         \
        0,                         /* tp_as_mapping */          \
        0,                         /* tp_hash */                \
        0,                         /* tp_call */                \
        0,                         /* tp_str */                 \
        0,                         /* tp_getattro */            \
        0,                         /* tp_setattro */            \
        0,                         /* tp_as_buffer */           \
        Py_TPFLAGS_DEFAULT,        /* tp_flags */               \
        (doc),                     /* tp_doc */                 \
    }

#define ENUM_TYPE(name,namestr,docstr)                          \
    typedef struct                                              \
    {                                                           \
        PyObject_HEAD                                           \
    } name ## _Object;                                          \
    SIMPLE_TYPE(name,namestr,docstr)

#define ADD_ENUM(typename,name,value)                           \
    do {                                                        \
        PyObject *o = PyInt_FromLong((value));                  \
        PyDict_SetItemString((&typename ## _Type)->             \
                             tp_dict,(name),o);                 \
        Py_DECREF(o);                                           \
    } while(0)

#define REGISTER_SIMPLE(m,py_name,name)                         \
    do {                                                        \
        name ## _Type.tp_new = PyType_GenericNew;               \
        if(PyType_Ready(&name ## _Type) < 0)                    \
            return;                                             \
        Py_INCREF(&name ## _Type);                              \
        PyModule_AddObject(m,py_name,                           \
                           (PyObject *)&name ## _Type);         \
    } while(0)

#endif
