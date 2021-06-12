#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <structmember.h>

typedef struct {
  PyObject_HEAD double head;
  double tail;
} QuadrupleObject;

static void Quadruple_dealloc(QuadrupleObject *self) {
  Py_TYPE(self)->tp_free((PyObject *)self);
}

static void two_add(double first, double second, double* head, double* tail) {
    double result_head = first + second;
    double second_virtual = result_head - first;
    double first_virtual = result_head - second_virtual;
    double first_tail = first - first_virtual;
    double second_tail = second - second_virtual;
    *head = result_head;
    *tail = first_tail + second_tail;
}

static PyObject *Quadruple_new(PyTypeObject *cls, PyObject *args,
                               PyObject *kwargs) {
  if (!_PyArg_NoKeywords("Quadruple", kwargs)) return NULL;
  double head = 0.0, tail = 0.0;
  if (!PyArg_ParseTuple(args, "|dd", &head, &tail)) return NULL;
  QuadrupleObject *self = (QuadrupleObject *)(cls->tp_alloc(cls, 0));
  if (self) {
    if (tail)
      two_add(head, tail, &head, &tail);
    self->head = head;
    self->tail = tail;
  }
  return (PyObject *)self;
}

static PyObject *Quadruple_repr(QuadrupleObject *self) {
  PyObject* head = PyFloat_FromDouble(self->head);
  PyObject* result;
  if (self->tail) {
    PyObject* tail = PyFloat_FromDouble(self->tail);
    result = PyUnicode_FromFormat("Quadruple(%R, %R)", head, tail);
    Py_XDECREF(tail);
  }
  else
    result = PyUnicode_FromFormat("Quadruple(%R)", head);
  Py_XDECREF(head);
  return result;
}

static PyTypeObject QuadrupleType = {
    PyVarObject_HEAD_INIT(NULL, 0).tp_basicsize = sizeof(QuadrupleObject),
    .tp_dealloc = (destructor)Quadruple_dealloc,
    .tp_doc =
        PyDoc_STR("Represents quadruple precision floating point number."),
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_itemsize = 0,
    .tp_name = "shewchuk.Quadruple",
    .tp_new = Quadruple_new,
    .tp_repr = (reprfunc)Quadruple_repr,
};

static PyModuleDef _shewchuk_module = {
    PyModuleDef_HEAD_INIT,
    .m_name = "shewchuk",
    .m_doc = PyDoc_STR("Robust floating point operations."),
    .m_size = -1,
};

PyMODINIT_FUNC PyInit__shewchuk(void) {
  PyObject *result;
  if (PyType_Ready(&QuadrupleType) < 0) return NULL;
  result = PyModule_Create(&_shewchuk_module);
  if (result == NULL) return NULL;
  Py_INCREF(&QuadrupleType);
  if (PyModule_AddObject(result, "Quadruple", (PyObject *)&QuadrupleType) < 0) {
    Py_DECREF(&QuadrupleType);
    Py_DECREF(result);
    return NULL;
  }
  return result;
}
