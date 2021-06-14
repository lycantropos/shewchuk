#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <structmember.h>

static void two_add(double left, double right, double *result_head,
                    double *result_tail) {
  double head = left + right;
  double right_virtual = head - left;
  double left_virtual = head - right_virtual;
  double left_tail = left - left_virtual;
  double right_tail = right - right_virtual;
  double tail = left_tail + right_tail;
  *result_head = head;
  *result_tail = tail;
}

static void fast_two_add(double left, double right, double *result_head,
                         double *result_tail) {
  double head = left + right;
  double right_virtual = head - left;
  double tail = right - right_virtual;
  *result_head = head;
  *result_tail = tail;
}

static void two_sub(double left, double right, double *result_head,
                    double *result_tail) {
  two_add(left, -right, result_head, result_tail);
}

static void two_one_add(double left_head, double left_tail, double right,
                        double *result_head, double *result_first_tail,
                        double *result_second_tail) {
  double mid_head;
  two_add(left_tail, right, &mid_head, result_second_tail);
  two_add(left_head, mid_head, result_head, result_first_tail);
}

static void two_one_sub(double left_head, double left_tail, double right,
                        double *head, double *first_tail, double *second_tail) {
  double mid_head;
  two_sub(left_tail, right, &mid_head, second_tail);
  two_add(left_head, mid_head, head, first_tail);
}

static void two_two_add(double left_head, double left_tail, double right_head,
                        double right_tail, double *head, double *first_tail,
                        double *second_tail, double *third_tail) {
  double mid_head, mid_tail;
  two_one_add(left_head, left_tail, right_tail, &mid_head, &mid_tail,
              third_tail);
  two_one_add(mid_head, mid_tail, right_head, head, first_tail, second_tail);
}

static void two_two_sub(double left_head, double left_tail, double right_head,
                        double right_tail, double *head, double *first_tail,
                        double *second_tail, double *third_tail) {
  double mid_head, mid_tail;
  two_one_sub(left_head, left_tail, right_tail, &mid_head, &mid_tail,
              third_tail);
  two_one_sub(mid_head, mid_tail, right_head, head, first_tail, second_tail);
}

size_t compress_components(size_t size, double *components) {
  size_t bottom = size - 1;
  double cursor = components[bottom];
  for (Py_ssize_t index = (Py_ssize_t)(bottom)-1; index >= 0; --index) {
    double head, tail;
    two_add(cursor, components[index], &head, &tail);
    if (!!tail) {
      components[bottom--] = head;
      cursor = tail;
    } else
      cursor = head;
  }
  size_t top = 0;
  for (size_t index = bottom + 1; index < size; ++index) {
    double head, tail;
    two_add(components[index], cursor, &head, &tail);
    if (!!tail) components[top++] = tail;
    cursor = head;
  }
  components[top] = cursor;
  return top + 1;
}

size_t add_components_eliminating_zeros(size_t left_size, double *left,
                                        size_t right_size, double *right,
                                        double *result) {
  size_t left_index = 0, right_index = 0;
  double left_component = left[left_index];
  double right_component = right[right_index];
  double cursor;
  if ((right_component > left_component) ==
      (right_component > -left_component)) {
    cursor = left_component;
    left_component = left[++left_index];
  } else {
    cursor = right_component;
    right_component = right[++right_index];
  }
  size_t result_size = 0;
  double head, tail;
  if ((left_index < left_size) && (right_index < right_size)) {
    if ((right_component > left_component) ==
        (right_component > -left_component)) {
      fast_two_add(left_component, cursor, &head, &tail);
      left_component = left[++left_index];
    } else {
      fast_two_add(right_component, cursor, &head, &tail);
      right_component = right[++right_index];
    }
    cursor = head;
    if (!!tail) result[result_size++] = tail;
    while ((left_index < left_size) && (right_index < right_size)) {
      if ((right_component > left_component) ==
          (right_component > -left_component)) {
        two_add(cursor, left_component, &head, &tail);
        left_component = left[++left_index];
      } else {
        two_add(cursor, right_component, &head, &tail);
        right_component = right[++right_index];
      }
      cursor = head;
      if (!!tail) result[result_size++] = tail;
    }
  }
  while (left_index < left_size) {
    two_add(cursor, left_component, &head, &tail);
    left_component = left[++left_index];
    cursor = head;
    if (!!tail) result[result_size++] = tail;
  }
  while (right_index < right_size) {
    two_add(cursor, right_component, &head, &tail);
    right_component = right[++right_index];
    cursor = head;
    if (!!tail) result[result_size++] = tail;
  }
  if (!!cursor || !result_size) result[result_size++] = cursor;
  return result_size;
}

typedef struct {
  PyObject_HEAD size_t size;
  double *components;
} ExpansionObject;

typedef struct {
  PyObject_HEAD double head;
  double tail;
} QuadrupleObject;

static ExpansionObject *construct_Expansion(PyTypeObject *cls,
                                            double *components, size_t size) {
  ExpansionObject *result = (ExpansionObject *)(cls->tp_alloc(cls, 0));
  if (result) {
    result->components = components;
    result->size = size;
  }
  return result;
}

static PyTypeObject ExpansionType;

static ExpansionObject *Expansions_add(ExpansionObject *self,
                                       ExpansionObject *other) {
  double *result_components =
      PyMem_RawCalloc(self->size + other->size, sizeof(double));
  if (!result_components) return PyErr_NoMemory();
  size_t result_size = add_components_eliminating_zeros(
      self->size, self->components, other->size, other->components,
      result_components);
  result_components =
      PyMem_RawRealloc(result_components, result_size * sizeof(double));
  if (!result_components) return PyErr_NoMemory();
  return construct_Expansion(Py_TYPE(self), result_components, result_size);
}

static PyObject *Expansion_add(PyObject *self, PyObject *other) {
  if (PyObject_TypeCheck(self, &ExpansionType)) {
    if (PyObject_TypeCheck(other, &ExpansionType))
      return (PyObject *)Expansions_add((ExpansionObject *)self,
                                        (ExpansionObject *)other);
  }
}

static void Expansion_dealloc(ExpansionObject *self) {
  PyMem_RawFree(self->components);
  Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject *Expansion_new(PyTypeObject *cls, PyObject *args,
                               PyObject *kwargs) {
  if (!_PyArg_NoKeywords("Expansion", kwargs)) return NULL;
  double *components;
  Py_ssize_t raw_size = PyTuple_Size(args);
  if (raw_size < 0) return NULL;
  size_t size = (size_t)raw_size;
  if (size) {
    components = (double *)PyMem_RawCalloc(size, sizeof(double));
    if (!components) return PyErr_NoMemory();
    for (size_t index = 0; index < size; ++index) {
      PyObject *item = PyTuple_GET_ITEM(args, index);
      if (!item) {
        PyMem_RawFree(components);
        return NULL;
      }
      double component = components[index] = PyFloat_AsDouble(item);
      if (component == -1.0 && PyErr_Occurred()) {
        PyMem_RawFree(components);
        return NULL;
      }
    }
    size = compress_components(size, components);
    components = PyMem_RawRealloc(components, size * sizeof(double));
    if (!components) return PyErr_NoMemory();
  } else {
    components = (double *)PyMem_RawMalloc(sizeof(double));
    if (!components) return PyErr_NoMemory();
    components[0] = 0.0;
    size = 1;
  }
  return (PyObject *)construct_Expansion(cls, components, size);
}

static PyObject *Expansion_repr(ExpansionObject *self) {
  PyObject *result;
  if (self->size > 1) {
    PyObject *components_reprs = PyTuple_New(self->size);
    if (!components_reprs) return NULL;
    for (size_t index = 0; index < self->size; ++index) {
      PyObject *item = PyFloat_FromDouble(self->components[index]);
      if (!item) {
        Py_DECREF(components_reprs);
        return NULL;
      }
      PyTuple_SET_ITEM(components_reprs, index, PyObject_Repr(item));
      Py_DECREF(item);
    }
    PyObject *separator = PyUnicode_FromString(", ");
    if (!separator) {
      Py_DECREF(components_reprs);
      return NULL;
    }
    PyObject *joined_components_reprs =
        PyUnicode_Join(separator, components_reprs);
    Py_DECREF(separator);
    Py_DECREF(components_reprs);
    if (!joined_components_reprs) return NULL;
    result = PyUnicode_FromFormat("Expansion(%U)", joined_components_reprs);
    Py_DECREF(joined_components_reprs);
  } else {
    PyObject *head = PyFloat_FromDouble(self->components[0]);
    result = PyUnicode_FromFormat("Expansion(%R)", head);
    Py_DECREF(head);
  }
  return result;
}

static void Quadruple_dealloc(QuadrupleObject *self) {
  Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject *Quadruple_new(PyTypeObject *cls, PyObject *args,
                               PyObject *kwargs) {
  if (!_PyArg_NoKeywords("Quadruple", kwargs)) return NULL;
  double head = 0.0, tail = 0.0;
  if (!PyArg_ParseTuple(args, "|dd", &head, &tail)) return NULL;
  QuadrupleObject *self = (QuadrupleObject *)(cls->tp_alloc(cls, 0));
  if (self) {
    if (tail) two_add(head, tail, &head, &tail);
    self->head = head;
    self->tail = tail;
  }
  return (PyObject *)self;
}

static PyObject *Quadruple_repr(QuadrupleObject *self) {
  PyObject *head = PyFloat_FromDouble(self->head);
  PyObject *result;
  if (self->tail) {
    PyObject *tail = PyFloat_FromDouble(self->tail);
    result = PyUnicode_FromFormat("Quadruple(%R, %R)", head, tail);
    Py_XDECREF(tail);
  } else
    result = PyUnicode_FromFormat("Quadruple(%R)", head);
  Py_XDECREF(head);
  return result;
}

static PyNumberMethods Expansion_as_number = {
    .nb_add = Expansion_add,
};

static PyTypeObject ExpansionType = {
    PyVarObject_HEAD_INIT(NULL, 0).tp_as_number = &Expansion_as_number,
    .tp_basicsize = sizeof(ExpansionObject),
    .tp_dealloc = (destructor)Expansion_dealloc,
    .tp_doc = PyDoc_STR("Represents floating point number expansion."),
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_itemsize = 0,
    .tp_name = "shewchuk.Expansion",
    .tp_new = Expansion_new,
    .tp_repr = (reprfunc)Expansion_repr,
};

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
  if (PyType_Ready(&ExpansionType) < 0) return NULL;
  if (PyType_Ready(&QuadrupleType) < 0) return NULL;
  result = PyModule_Create(&_shewchuk_module);
  if (result == NULL) return NULL;
  Py_INCREF(&ExpansionType);
  if (PyModule_AddObject(result, "Expansion", (PyObject *)&ExpansionType) < 0) {
    Py_DECREF(&ExpansionType);
    Py_DECREF(result);
    return NULL;
  }
  Py_INCREF(&QuadrupleType);
  if (PyModule_AddObject(result, "Quadruple", (PyObject *)&QuadrupleType) < 0) {
    Py_DECREF(&QuadrupleType);
    Py_DECREF(result);
    return NULL;
  }
  return result;
}
