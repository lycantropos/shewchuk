#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <float.h>
#include <structmember.h>

static int are_components_equal(size_t left_size, double *left,
                                size_t right_size, double *right) {
  if (left_size != right_size) return 0;
  for (size_t offset = 1; offset <= left_size; ++offset)
    if (left[left_size - offset] != right[left_size - offset]) return 0;
  return 1;
}

static int are_components_lesser_than(size_t left_size, double *left,
                                      size_t right_size, double *right) {
  size_t min_size = left_size < right_size ? left_size : right_size;
  for (size_t offset = 1; offset <= min_size; ++offset)
    if (left[left_size - offset] < right[right_size - offset])
      return 1;
    else if (left[left_size - offset] > right[right_size - offset])
      return 0;
  return left_size != right_size &&
         (left_size < right_size ? right[right_size - left_size - 1] > 0.0
                                 : left[left_size - right_size - 1] < 0.0);
}

static void fast_two_add(double left, double right, double *result_head,
                         double *result_tail) {
  double head = left + right;
  double right_virtual = head - left;
  double tail = right - right_virtual;
  *result_head = head;
  *result_tail = tail;
}

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

static void split(double value, double *result_high, double *result_low) {
  static const double splitter = (1 << (size_t)((DBL_MANT_DIG + 1) / 2)) + 1;
  double base = splitter * value;
  double high = base - (base - value);
  double low = value - high;
  *result_high = high;
  *result_low = low;
}

static void two_subtract(double left, double right, double *result_head,
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

static void two_one_subtract(double left_head, double left_tail, double right,
                             double *head, double *first_tail,
                             double *second_tail) {
  double mid_head;
  two_subtract(left_tail, right, &mid_head, second_tail);
  two_add(left_head, mid_head, head, first_tail);
}

static void two_product_presplit(double left, double right, double right_high,
                                 double right_low, double *result_head,
                                 double *result_tail) {
  double head = left * right;
  double left_high, left_low;
  split(left, &left_high, &left_low);
  double first_error = head - left_high * right_high;
  double second_error = first_error - left_low * right_high;
  double third_error = second_error - left_high * right_low;
  double tail = left_low * right_low - third_error;
  *result_head = head;
  *result_tail = tail;
}

static void two_two_add(double left_head, double left_tail, double right_head,
                        double right_tail, double *head, double *first_tail,
                        double *second_tail, double *third_tail) {
  double mid_head, mid_tail;
  two_one_add(left_head, left_tail, right_tail, &mid_head, &mid_tail,
              third_tail);
  two_one_add(mid_head, mid_tail, right_head, head, first_tail, second_tail);
}

static void two_two_subtract(double left_head, double left_tail,
                             double right_head, double right_tail, double *head,
                             double *first_tail, double *second_tail,
                             double *third_tail) {
  double mid_head, mid_tail;
  two_one_subtract(left_head, left_tail, right_tail, &mid_head, &mid_tail,
                   third_tail);
  two_one_subtract(mid_head, mid_tail, right_head, head, first_tail,
                   second_tail);
}

static size_t compress_components_single(size_t size, double *components) {
  size_t bottom = size - 1;
  double accumulator = components[bottom];
  for (Py_ssize_t index = (Py_ssize_t)(bottom)-1; index >= 0; --index) {
    double head, tail;
    two_add(accumulator, components[index], &head, &tail);
    if (!!tail) {
      components[bottom--] = head;
      accumulator = tail;
    } else
      accumulator = head;
  }
  size_t top = 0;
  for (size_t index = bottom + 1; index < size; ++index) {
    double head, tail;
    two_add(components[index], accumulator, &head, &tail);
    if (!!tail) components[top++] = tail;
    accumulator = head;
  }
  if (!!accumulator || !top) components[top++] = accumulator;
  return top;
}

static size_t compress_components(size_t size, double *components) {
  const size_t original_size = size;
  double *next_components = PyMem_RawCalloc(size, sizeof(double));
  for (size_t index = 0; index < size; ++index)
    next_components[index] = components[index];
  for (size_t step = 0; step < original_size; ++step) {
    const size_t next_size = compress_components_single(size, next_components);
    if (are_components_equal(next_size, next_components, size, components))
      break;
    size = next_size;
    for (size_t index = 0; index < size; ++index)
      components[index] = next_components[index];
  }
  PyMem_RawFree(next_components);
  return size;
}

static size_t add_double_eliminating_zeros(size_t left_size, double *left,
                                           double right, double *result) {
  size_t result_size = 0;
  double accumulator = right;
  for (size_t index = 0; index < left_size; index++) {
    double tail;
    two_add(accumulator, left[index], &accumulator, &tail);
    if (!!tail) result[result_size++] = tail;
  }
  if (!!accumulator || !result_size) result[result_size++] = accumulator;
  return result_size;
}

static size_t scale_components(size_t size, double *components, double scalar,
                               double *result) {
  double scalar_high, scalar_low;
  split(scalar, &scalar_high, &scalar_low);
  double accumulator, tail;
  two_product_presplit(components[0], scalar, scalar_high, scalar_low,
                       &accumulator, &tail);
  size_t result_size = 0;
  if (!!tail) result[result_size++] = tail;
  for (size_t index = 1; index < size; ++index) {
    double product, product_tail;
    two_product_presplit(components[index], scalar, scalar_high, scalar_low,
                         &product, &product_tail);
    double interim;
    two_add(accumulator, product_tail, &interim, &tail);
    if (!!tail) result[result_size++] = tail;
    fast_two_add(product, interim, &accumulator, &tail);
    if (!!tail) result[result_size++] = tail;
  }
  if (!!accumulator || !result_size) result[result_size++] = accumulator;
  return result_size;
}

static size_t subtract_double_eliminating_zeros(size_t minuend_size,
                                                double *minuend,
                                                double subtrahend,
                                                double *result) {
  return add_double_eliminating_zeros(minuend_size, minuend, -subtrahend,
                                      result);
}

static size_t subtract_from_double_eliminating_zeros(double minuend,
                                                     size_t subtrahend_size,
                                                     double *subtrahend,
                                                     double *result) {
  size_t result_size = 0;
  double accumulator = minuend;
  for (size_t index = 0; index < subtrahend_size; index++) {
    double tail;
    two_add(accumulator, -subtrahend[index], &accumulator, &tail);
    if (!!tail) result[result_size++] = tail;
  }
  if (!!accumulator || !result_size) result[result_size++] = accumulator;
  return result_size;
}

static size_t add_components_eliminating_zeros(size_t left_size, double *left,
                                               size_t right_size, double *right,
                                               double *result) {
  size_t left_index = 0, right_index = 0;
  double left_component = left[left_index];
  double right_component = right[right_index];
  double accumulator;
  if ((right_component > left_component) ==
      (right_component > -left_component)) {
    accumulator = left_component;
    left_component = left[++left_index];
  } else {
    accumulator = right_component;
    right_component = right[++right_index];
  }
  size_t result_size = 0;
  double head, tail;
  if ((left_index < left_size) && (right_index < right_size)) {
    if ((right_component > left_component) ==
        (right_component > -left_component)) {
      fast_two_add(left_component, accumulator, &head, &tail);
      left_component = left[++left_index];
    } else {
      fast_two_add(right_component, accumulator, &head, &tail);
      right_component = right[++right_index];
    }
    accumulator = head;
    if (!!tail) result[result_size++] = tail;
    while ((left_index < left_size) && (right_index < right_size)) {
      if ((right_component > left_component) ==
          (right_component > -left_component)) {
        two_add(accumulator, left_component, &head, &tail);
        left_component = left[++left_index];
      } else {
        two_add(accumulator, right_component, &head, &tail);
        right_component = right[++right_index];
      }
      accumulator = head;
      if (!!tail) result[result_size++] = tail;
    }
  }
  while (left_index < left_size) {
    two_add(accumulator, left_component, &head, &tail);
    left_component = left[++left_index];
    accumulator = head;
    if (!!tail) result[result_size++] = tail;
  }
  while (right_index < right_size) {
    two_add(accumulator, right_component, &head, &tail);
    right_component = right[++right_index];
    accumulator = head;
    if (!!tail) result[result_size++] = tail;
  }
  if (!!accumulator || !result_size) result[result_size++] = accumulator;
  return result_size;
}

static size_t subtract_components_eliminating_zeros(size_t minuend_size,
                                                    double *minuend,
                                                    size_t subtrahend_size,
                                                    double *subtrahend,
                                                    double *result) {
  size_t minuend_index = 0, subtrahend_index = 0;
  double minuend_component = minuend[minuend_index];
  double subtrahend_component = -subtrahend[subtrahend_index];
  double accumulator;
  if ((subtrahend_component > minuend_component) ==
      (subtrahend_component > -minuend_component)) {
    accumulator = minuend_component;
    minuend_component = minuend[++minuend_index];
  } else {
    accumulator = subtrahend_component;
    subtrahend_component = -subtrahend[++subtrahend_index];
  }
  size_t result_size = 0;
  double head, tail;
  if ((minuend_index < minuend_size) && (subtrahend_index < subtrahend_size)) {
    if ((subtrahend_component > minuend_component) ==
        (subtrahend_component > -minuend_component)) {
      fast_two_add(minuend_component, accumulator, &head, &tail);
      minuend_component = minuend[++minuend_index];
    } else {
      fast_two_add(subtrahend_component, accumulator, &head, &tail);
      subtrahend_component = -subtrahend[++subtrahend_index];
    }
    accumulator = head;
    if (!!tail) result[result_size++] = tail;
    while ((minuend_index < minuend_size) &&
           (subtrahend_index < subtrahend_size)) {
      if ((subtrahend_component > minuend_component) ==
          (subtrahend_component > -minuend_component)) {
        two_add(accumulator, minuend_component, &head, &tail);
        minuend_component = minuend[++minuend_index];
      } else {
        two_add(accumulator, subtrahend_component, &head, &tail);
        subtrahend_component = -subtrahend[++subtrahend_index];
      }
      accumulator = head;
      if (!!tail) result[result_size++] = tail;
    }
  }
  while (minuend_index < minuend_size) {
    two_add(accumulator, minuend_component, &head, &tail);
    minuend_component = minuend[++minuend_index];
    accumulator = head;
    if (!!tail) result[result_size++] = tail;
  }
  while (subtrahend_index < subtrahend_size) {
    two_add(accumulator, subtrahend_component, &head, &tail);
    subtrahend_component = -subtrahend[++subtrahend_index];
    accumulator = head;
    if (!!tail) result[result_size++] = tail;
  }
  if (!!accumulator || !result_size) result[result_size++] = accumulator;
  return result_size;
}

int is_PyObject_convertible_to_Float(PyObject *self) {
  return !!Py_TYPE(self)->tp_as_number &&
         !!Py_TYPE(self)->tp_as_number->nb_float;
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
  if (!result_components) return (ExpansionObject *)PyErr_NoMemory();
  size_t result_size = add_components_eliminating_zeros(
      self->size, self->components, other->size, other->components,
      result_components);
  result_components =
      PyMem_RawRealloc(result_components, result_size * sizeof(double));
  if (!result_components) return (ExpansionObject *)PyErr_NoMemory();
  return construct_Expansion(Py_TYPE(self), result_components, result_size);
}

static ExpansionObject *Expansion_double_add(ExpansionObject *self,
                                             double other) {
  double *result_components = PyMem_RawCalloc(self->size + 1, sizeof(double));
  if (!result_components) return (ExpansionObject *)PyErr_NoMemory();
  size_t result_size = add_double_eliminating_zeros(
      self->size, self->components, other, result_components);
  result_components =
      PyMem_RawRealloc(result_components, result_size * sizeof(double));
  if (!result_components) return (ExpansionObject *)PyErr_NoMemory();
  return construct_Expansion(Py_TYPE(self), result_components, result_size);
}

static PyObject *Expansion_add(PyObject *self, PyObject *other) {
  if (PyObject_TypeCheck(self, &ExpansionType)) {
    if (PyObject_TypeCheck(other, &ExpansionType))
      return (PyObject *)Expansions_add((ExpansionObject *)self,
                                        (ExpansionObject *)other);
    else if (PyFloat_Check(other))
      return (PyObject *)Expansion_double_add((ExpansionObject *)self,
                                              PyFloat_AS_DOUBLE(other));
    else if (is_PyObject_convertible_to_Float(other)) {
      double other_value = PyFloat_AsDouble(other);
      return other_value == -1.0 && PyErr_Occurred()
                 ? NULL
                 : (PyObject *)Expansion_double_add((ExpansionObject *)self,
                                                    other_value);
    }
  } else if (PyFloat_Check(self))
    return (PyObject *)Expansion_double_add((ExpansionObject *)other,
                                            PyFloat_AS_DOUBLE(self));
  else if (!!Py_TYPE(self)->tp_as_number &&
           !!Py_TYPE(self)->tp_as_number->nb_float) {
    double value = PyFloat_AsDouble(self);
    return value == -1.0 && PyErr_Occurred()
               ? NULL
               : (PyObject *)Expansion_double_add((ExpansionObject *)other,
                                                  value);
  }
  Py_RETURN_NOTIMPLEMENTED;
}

static int Expansion_bool(ExpansionObject *self) {
  return !!self->components[self->size - 1];
}

static void Expansion_dealloc(ExpansionObject *self) {
  PyMem_RawFree(self->components);
  Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject *Expansion_float(ExpansionObject *self) {
  double result = self->components[0];
  for (size_t index = 1; index < self->size; ++index)
    result += self->components[index];
  return PyFloat_FromDouble(result);
}

static ExpansionObject *Expansion_double_multiply(ExpansionObject *self,
                                                  double other) {
  double *result_components =
      PyMem_RawCalloc(2 * self->size + 1, sizeof(double));
  if (!result_components) return (ExpansionObject *)PyErr_NoMemory();
  size_t result_size =
      scale_components(self->size, self->components, other, result_components);
  result_components =
      PyMem_RawRealloc(result_components, result_size * sizeof(double));
  if (!result_components) return (ExpansionObject *)PyErr_NoMemory();
  return construct_Expansion(Py_TYPE(self), result_components, result_size);
}

static PyObject *Expansion_multiply(PyObject *self, PyObject *other) {
  if (PyObject_TypeCheck(self, &ExpansionType)) {
    if (PyFloat_Check(other))
      return (PyObject *)Expansion_double_multiply((ExpansionObject *)self,
                                                   PyFloat_AS_DOUBLE(other));
    else if (is_PyObject_convertible_to_Float(other)) {
      double other_value = PyFloat_AsDouble(other);
      return other_value == -1.0 && PyErr_Occurred()
                 ? NULL
                 : (PyObject *)Expansion_double_multiply(
                       (ExpansionObject *)self, other_value);
    }
  } else if (PyFloat_Check(self))
    return (PyObject *)Expansion_double_multiply((ExpansionObject *)other,
                                                 PyFloat_AS_DOUBLE(self));
  else if (is_PyObject_convertible_to_Float(self)) {
    double value = PyFloat_AsDouble(self);
    return value == -1.0 && PyErr_Occurred()
               ? NULL
               : (PyObject *)Expansion_double_multiply((ExpansionObject *)other,
                                                       value);
  }
  Py_RETURN_NOTIMPLEMENTED;
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

static ExpansionObject *Expansion_negative(ExpansionObject *self) {
  double *result_components = PyMem_RawCalloc(self->size, sizeof(double));
  for (size_t index = 0; index < self->size; ++index)
    result_components[index] = -self->components[index];
  return construct_Expansion(Py_TYPE(self), result_components, self->size);
}

static ExpansionObject *Expansion_positive(ExpansionObject *self) {
  Py_INCREF(self);
  return self;
}

static ExpansionObject *Expansion_absolute(ExpansionObject *self) {
  return self->components[self->size - 1] < 0.0 ? Expansion_negative(self)
                                                : Expansion_positive(self);
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

static ExpansionObject *Expansions_subtract(ExpansionObject *self,
                                            ExpansionObject *other) {
  double *result_components =
      PyMem_RawCalloc(self->size + other->size, sizeof(double));
  if (!result_components) return (ExpansionObject *)PyErr_NoMemory();
  size_t result_size = subtract_components_eliminating_zeros(
      self->size, self->components, other->size, other->components,
      result_components);
  result_components =
      PyMem_RawRealloc(result_components, result_size * sizeof(double));
  if (!result_components) return (ExpansionObject *)PyErr_NoMemory();
  return construct_Expansion(Py_TYPE(self), result_components, result_size);
}

static ExpansionObject *Expansion_double_subtract(ExpansionObject *self,
                                                  double other) {
  double *result_components = PyMem_RawCalloc(self->size + 1, sizeof(double));
  if (!result_components) return (ExpansionObject *)PyErr_NoMemory();
  size_t result_size = subtract_double_eliminating_zeros(
      self->size, self->components, other, result_components);
  result_components =
      PyMem_RawRealloc(result_components, result_size * sizeof(double));
  if (!result_components) return (ExpansionObject *)PyErr_NoMemory();
  return construct_Expansion(Py_TYPE(self), result_components, result_size);
}

static ExpansionObject *double_Expansion_subtract(double self,
                                                  ExpansionObject *other) {
  double *result_components = PyMem_RawCalloc(other->size + 1, sizeof(double));
  if (!result_components) return (ExpansionObject *)PyErr_NoMemory();
  size_t result_size = subtract_from_double_eliminating_zeros(
      self, other->size, other->components, result_components);
  result_components =
      PyMem_RawRealloc(result_components, result_size * sizeof(double));
  if (!result_components) return (ExpansionObject *)PyErr_NoMemory();
  return construct_Expansion(Py_TYPE(other), result_components, result_size);
}

static PyObject *Expansion_subtract(PyObject *self, PyObject *other) {
  if (PyObject_TypeCheck(self, &ExpansionType)) {
    if (PyObject_TypeCheck(other, &ExpansionType))
      return (PyObject *)Expansions_subtract((ExpansionObject *)self,
                                             (ExpansionObject *)other);
    else if (PyFloat_Check(other))
      return (PyObject *)Expansion_double_subtract((ExpansionObject *)self,
                                                   PyFloat_AS_DOUBLE(other));
    else if (!!Py_TYPE(other)->tp_as_number &&
             !!Py_TYPE(other)->tp_as_number->nb_float) {
      double other_value = PyFloat_AsDouble(other);
      return other_value == -1.0 && PyErr_Occurred()
                 ? NULL
                 : (PyObject *)Expansion_double_subtract(
                       (ExpansionObject *)self, other_value);
    }
  } else if (PyFloat_Check(self))
    return (PyObject *)double_Expansion_subtract(PyFloat_AS_DOUBLE(self),
                                                 (ExpansionObject *)other);
  else if (!!Py_TYPE(self)->tp_as_number &&
           !!Py_TYPE(self)->tp_as_number->nb_float) {
    double value = PyFloat_AsDouble(self);
    return value == -1.0 && PyErr_Occurred()
               ? NULL
               : (PyObject *)double_Expansion_subtract(
                     value, (ExpansionObject *)other);
  }
  Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *Expansions_richcompare(ExpansionObject *self,
                                        ExpansionObject *other, int op) {
  switch (op) {
    case Py_EQ:
      return PyBool_FromLong(are_components_equal(
          self->size, self->components, other->size, other->components));
    case Py_GE:
      return PyBool_FromLong(!are_components_lesser_than(
          self->size, self->components, other->size, other->components));
    case Py_GT:
      return PyBool_FromLong(are_components_lesser_than(
          other->size, other->components, self->size, self->components));
    case Py_LE:
      return PyBool_FromLong(!are_components_lesser_than(
          other->size, other->components, self->size, self->components));
    case Py_LT:
      return PyBool_FromLong(are_components_lesser_than(
          self->size, self->components, other->size, other->components));
    case Py_NE:
      return PyBool_FromLong(!are_components_equal(
          self->size, self->components, other->size, other->components));
    default:
      Py_RETURN_NOTIMPLEMENTED;
  }
}

static PyObject *Expansion_double_richcompare(ExpansionObject *self,
                                              double other, int op) {
  switch (op) {
    case Py_EQ:
      return PyBool_FromLong(self->size == 1 && self->components[0] == other);
    case Py_GE:
      return PyBool_FromLong(
          self->components[self->size - 1] > other ||
          (self->components[self->size - 1] == other &&
           (self->size == 1 || self->components[self->size - 2] > 0.0)));
    case Py_GT:
      return PyBool_FromLong(
          self->components[self->size - 1] > other ||
          (self->components[self->size - 1] == other &&
           (self->size > 1 && self->components[self->size - 2] > 0.0)));
    case Py_LE:
      return PyBool_FromLong(
          self->components[self->size - 1] < other ||
          (self->components[self->size - 1] == other &&
           (self->size == 1 || self->components[self->size - 2] < 0.0)));
    case Py_LT:
      return PyBool_FromLong(
          self->components[self->size - 1] < other ||
          (self->components[self->size - 1] == other &&
           (self->size > 1 && self->components[self->size - 2] < 0.0)));
    case Py_NE:
      return PyBool_FromLong(self->size > 1 || self->components[0] != other);
    default:
      Py_RETURN_NOTIMPLEMENTED;
  }
}

static PyObject *Expansion_PyObject_richcompare(ExpansionObject *self,
                                                PyObject *other, int op) {
  switch (op) {
    case Py_EQ: {
      if (self->size > 1) Py_RETURN_FALSE;
      PyObject *head = PyFloat_FromDouble(self->components[0]);
      PyObject *result = PyObject_RichCompare(head, other, op);
      Py_DECREF(head);
      return result;
    }
    case Py_GE: {
      PyObject *head = PyFloat_FromDouble(self->components[self->size - 1]);
      PyObject *result = PyObject_RichCompare(head, other, Py_GT);
      if (!result || result == Py_True) return result;
      Py_DECREF(result);
      result = PyObject_RichCompare(head, other, Py_EQ);
      if (!result || result == Py_False) return result;
      Py_DECREF(result);
      return PyBool_FromLong(self->size == 1 ||
                             self->components[self->size - 2] > 0.0);
    }
    case Py_GT: {
      PyObject *head = PyFloat_FromDouble(self->components[self->size - 1]);
      PyObject *result = PyObject_RichCompare(head, other, Py_GT);
      if (!result || result == Py_True) return result;
      Py_DECREF(result);
      result = PyObject_RichCompare(head, other, Py_EQ);
      if (!result || result == Py_False) return result;
      Py_DECREF(result);
      return PyBool_FromLong(self->size > 1 &&
                             self->components[self->size - 2] > 0.0);
    }
    case Py_LE: {
      PyObject *head = PyFloat_FromDouble(self->components[self->size - 1]);
      PyObject *result = PyObject_RichCompare(head, other, Py_LT);
      if (!result || result == Py_True) return result;
      Py_DECREF(result);
      result = PyObject_RichCompare(head, other, Py_EQ);
      if (!result || result == Py_False) return result;
      Py_DECREF(result);
      return PyBool_FromLong(self->size == 1 ||
                             self->components[self->size - 2] < 0.0);
    }
    case Py_LT: {
      PyObject *head = PyFloat_FromDouble(self->components[self->size - 1]);
      PyObject *result = PyObject_RichCompare(head, other, Py_LT);
      if (!result || result == Py_True) return result;
      Py_DECREF(result);
      result = PyObject_RichCompare(head, other, Py_EQ);
      if (!result || result == Py_False) return result;
      Py_DECREF(result);
      return PyBool_FromLong(self->size > 1 &&
                             self->components[self->size - 2] < 0.0);
    }
    case Py_NE: {
      if (self->size > 1) Py_RETURN_TRUE;
      PyObject *head = PyFloat_FromDouble(self->components[0]);
      PyObject *result = PyObject_RichCompare(head, other, op);
      Py_DECREF(head);
      return result;
    }
    default:
      Py_RETURN_NOTIMPLEMENTED;
  }
}

static PyObject *Expansion_richcompare(ExpansionObject *self, PyObject *other,
                                       int op) {
  if (PyObject_TypeCheck(other, &ExpansionType))
    return Expansions_richcompare(self, (ExpansionObject *)other, op);
  else if (PyFloat_Check(other))
    return Expansion_double_richcompare(self, PyFloat_AS_DOUBLE(other), op);
  else
    return Expansion_PyObject_richcompare(self, other, op);
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
    .nb_absolute = (unaryfunc)Expansion_absolute,
    .nb_add = Expansion_add,
    .nb_bool = (inquiry)Expansion_bool,
    .nb_float = (unaryfunc)Expansion_float,
    .nb_multiply = Expansion_multiply,
    .nb_negative = (unaryfunc)Expansion_negative,
    .nb_positive = (unaryfunc)Expansion_positive,
    .nb_subtract = Expansion_subtract,
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
    .tp_richcompare = (richcmpfunc)Expansion_richcompare,
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
