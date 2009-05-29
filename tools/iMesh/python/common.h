#pragma once

#include <Python.h>

#define SIMPLE_TYPE(name,type,name_str,doc)		\
  static PyTypeObject type = {				\
    PyObject_HEAD_INIT(NULL)				\
    0,                         /*ob_size*/		\
    (name_str),                /*tp_name*/		\
    sizeof(name),              /*tp_basicsize*/		\
    0,                         /*tp_itemsize*/		\
    0,                         /*tp_dealloc*/		\
    0,                         /*tp_print*/		\
    0,                         /*tp_getattr*/		\
    0,                         /*tp_setattr*/		\
    0,                         /*tp_compare*/		\
    0,                         /*tp_repr*/		\
    0,                         /*tp_as_number*/		\
    0,                         /*tp_as_sequence*/	\
    0,                         /*tp_as_mapping*/	\
    0,                         /*tp_hash */		\
    0,                         /*tp_call*/		\
    0,                         /*tp_str*/		\
    0,                         /*tp_getattro*/		\
    0,                         /*tp_setattro*/		\
    0,                         /*tp_as_buffer*/		\
    Py_TPFLAGS_DEFAULT,        /*tp_flags*/		\
    (doc),                     /* tp_doc */		\
}

#define ENUM_TYPE(name,namestr,docstr)			\
  typedef struct					\
  {							\
    PyObject_HEAD					\
  } name ## _Object;					\
  SIMPLE_TYPE(name ## _Object,name ## _Type,		\
		     namestr,docstr)

#define ADD_ENUM(typename,name,value)                   \
  do {                                                  \
    PyObject *o = Py_BuildValue("i",(value));           \
    PyDict_SetItemString((&typename ## _Type)->         \
                         tp_dict,(name),o);             \
    Py_DECREF(o);                                       \
  } while(0)

#define REGISTER_SIMPLE(m,name)				\
  do {							\
    name ## _Type.tp_new = PyType_GenericNew;		\
    if(PyType_Ready(&name ## _Type) < 0)		\
      return;						\
    Py_INCREF(&name ## _Type);				\
    PyModule_AddObject(m,#name,				\
		       (PyObject *)&name ## _Type);	\
  } while(0)

#define REGISTER_SIMPLE_SUB(base,name)			\
  do {                                                  \
    name ## _Type.tp_new = PyType_GenericNew;           \
    if(PyType_Ready(&name ## _Type) < 0)                \
      return;                                           \
    Py_INCREF(&name ## _Type);                          \
    PyDict_SetItemString(base.tp_dict,#name,		\
			 (PyObject *)&name ## _Type);	\
  } while(0)
