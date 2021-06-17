#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <float.h>
#include <math.h>
#include <structmember.h>
#define EPSILON (DBL_EPSILON / 2.0)

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

static void two_multiply(double left, double right, double *result_head,
                         double *result_tail) {
  double head = left * right;
  double left_high, left_low;
  split(left, &left_high, &left_low);
  double right_high, right_low;
  split(right, &right_high, &right_low);
  double first_error = head - left_high * right_high;
  double second_error = first_error - left_low * right_high;
  double third_error = second_error - left_high * right_low;
  double tail = left_low * right_low - third_error;
  *result_head = head;
  *result_tail = tail;
}

static void two_multiply_presplit(double left, double right, double right_high,
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

static double two_subtract_tail(double left, double right, double head) {
  double right_virtual = left - head;
  double left_virtual = head + right_virtual;
  double right_error = right_virtual - right;
  double left_error = left - left_virtual;
  return left_error + right_error;
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

static void copy_components(double *source, size_t size, double *destination) {
  for (size_t index = 0; index < size; ++index)
    destination[index] = source[index];
}

static size_t compress_components(size_t size, double *components) {
  const size_t original_size = size;
  double *next_components = PyMem_RawCalloc(size, sizeof(double));
  copy_components(components, size, next_components);
  for (size_t step = 0; step < original_size; ++step) {
    const size_t next_size = compress_components_single(size, next_components);
    if (are_components_equal(next_size, next_components, size, components))
      break;
    size = next_size;
    copy_components(next_components, size, components);
  }
  PyMem_RawFree(next_components);
  return size;
}

double sum_components(size_t size, double *components) {
  double result = components[0];
  for (size_t index = 1; index < size; ++index) result += components[index];
  return result;
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
  two_multiply_presplit(components[0], scalar, scalar_high, scalar_low,
                        &accumulator, &tail);
  size_t result_size = 0;
  if (!!tail) result[result_size++] = tail;
  for (size_t index = 1; index < size; ++index) {
    double product, product_tail;
    two_multiply_presplit(components[index], scalar, scalar_high, scalar_low,
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

size_t adaptive_vectors_cross_product_impl(
    double first_start_x, double first_start_y, double first_end_x,
    double first_end_y, double second_start_x, double second_start_y,
    double second_end_x, double second_end_y, double upper_bound,
    double *result) {
  double minuend_x = first_end_x - first_start_x;
  double minuend_y = first_end_y - first_start_y;
  double subtrahend_x = second_end_x - second_start_x;
  double subtrahend_y = second_end_y - second_start_y;
  double minuend, minuend_tail;
  two_multiply(minuend_x, subtrahend_y, &minuend, &minuend_tail);
  double subtrahend, subtrahend_tail;
  two_multiply(minuend_y, subtrahend_x, &subtrahend, &subtrahend_tail);
  double first_components[4];
  two_two_subtract(minuend, minuend_tail, subtrahend, subtrahend_tail,
                   &first_components[3], &first_components[2],
                   &first_components[1], &first_components[0]);
  double estimation = sum_components(4, first_components);
  static const double first_upper_bound_coefficient =
      (2.0 + 12.0 * EPSILON) * EPSILON;
  double threshold = first_upper_bound_coefficient * upper_bound;
  if ((estimation >= threshold) || (-estimation >= threshold)) {
    copy_components(first_components, 4, result);
    return 4;
  }
  double minuend_x_tail =
      two_subtract_tail(first_end_x, first_start_x, minuend_x);
  double subtrahend_x_tail =
      two_subtract_tail(second_end_x, second_start_x, subtrahend_x);
  double minuend_y_tail =
      two_subtract_tail(first_end_y, first_start_y, minuend_y);
  double subtrahend_y_tail =
      two_subtract_tail(second_end_y, second_start_y, subtrahend_y);
  if (!minuend_x_tail && !minuend_y_tail && !subtrahend_x_tail &&
      !subtrahend_y_tail) {
    copy_components(first_components, 4, result);
    return 4;
  }
  static const double second_upper_bound_coefficient =
      (9.0 + 64.0 * EPSILON) * EPSILON * EPSILON;
  static const double estimation_coefficient = (3.0 + 8.0 * EPSILON) * EPSILON;
  threshold = second_upper_bound_coefficient * upper_bound +
              estimation_coefficient * fabs(estimation);
  double extra =
      (minuend_x * subtrahend_y_tail + subtrahend_y * minuend_x_tail) -
      (minuend_y * subtrahend_x_tail + subtrahend_x * minuend_y_tail);
  estimation += extra;
  if ((estimation >= threshold) || (-estimation >= threshold)) {
    size_t result_size =
        add_double_eliminating_zeros(4, first_components, extra, result);
    return result_size;
  }
  double minuend_x_subtrahend_y_head, minuend_x_subtrahend_y_tail;
  two_multiply(minuend_x_tail, subtrahend_y, &minuend_x_subtrahend_y_head,
               &minuend_x_subtrahend_y_tail);
  double minuend_y_subtrahend_x_head, minuend_y_subtrahend_x_tail;
  two_multiply(minuend_y_tail, subtrahend_x, &minuend_y_subtrahend_x_head,
               &minuend_y_subtrahend_x_tail);
  double extra_components[4];
  two_two_subtract(minuend_x_subtrahend_y_head, minuend_x_subtrahend_y_tail,
                   minuend_y_subtrahend_x_head, minuend_y_subtrahend_x_tail,
                   &extra_components[3], &extra_components[2],
                   &extra_components[1], &extra_components[0]);
  double second_components[8];
  size_t second_components_size = add_components_eliminating_zeros(
      4, first_components, 4, extra_components, second_components);
  two_multiply(minuend_x, subtrahend_y_tail, &minuend_x_subtrahend_y_head,
               &minuend_x_subtrahend_y_tail);
  two_multiply(minuend_y, subtrahend_x_tail, &minuend_y_subtrahend_x_head,
               &minuend_y_subtrahend_x_tail);
  two_two_subtract(minuend_x_subtrahend_y_head, minuend_x_subtrahend_y_tail,
                   minuend_y_subtrahend_x_head, minuend_y_subtrahend_x_tail,
                   &extra_components[3], &extra_components[2],
                   &extra_components[1], &extra_components[0]);
  double third_components[12];
  size_t third_components_size = add_components_eliminating_zeros(
      second_components_size, second_components, 4, extra_components,
      third_components);
  two_multiply(minuend_x_tail, subtrahend_y_tail, &minuend_x_subtrahend_y_head,
               &minuend_x_subtrahend_y_tail);
  two_multiply(minuend_y_tail, subtrahend_x_tail, &minuend_y_subtrahend_x_head,
               &minuend_y_subtrahend_x_tail);
  two_two_subtract(minuend_x_subtrahend_y_head, minuend_x_subtrahend_y_tail,
                   minuend_y_subtrahend_x_head, minuend_y_subtrahend_x_tail,
                   &extra_components[3], &extra_components[2],
                   &extra_components[1], &extra_components[0]);
  return add_components_eliminating_zeros(
      third_components_size, third_components, 4, extra_components, result);
}

double vectors_cross_product_impl(double first_start_x, double first_start_y,
                                  double first_end_x, double first_end_y,
                                  double second_start_x, double second_start_y,
                                  double second_end_x, double second_end_y,
                                  double *result) {
  double minuend =
      (first_end_x - first_start_x) * (second_end_y - second_start_y);
  double subtrahend =
      (first_end_y - first_start_y) * (second_end_x - second_start_x);
  double estimation = minuend - subtrahend;
  double upper_bound;
  if (minuend > 0.0) {
    if (subtrahend <= 0.0) {
      result[0] = estimation;
      return 1;
    } else
      upper_bound = minuend + subtrahend;
  } else if (minuend < 0.0) {
    if (subtrahend >= 0.0) {
      result[0] = estimation;
      return 1;
    } else
      upper_bound = -minuend - subtrahend;
  } else {
    result[0] = estimation;
    return 1;
  }
  static const double upper_bound_coefficient =
      (3.0 + 16.0 * EPSILON) * EPSILON;
  double threshold = upper_bound_coefficient * upper_bound;
  if ((estimation >= threshold) || (-estimation >= threshold)) {
    result[0] = estimation;
    return 1;
  }
  return adaptive_vectors_cross_product_impl(
      first_start_x, first_start_y, first_end_x, first_end_y, second_start_x,
      second_start_y, second_end_x, second_end_y, upper_bound, result);
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
  return construct_Expansion(&ExpansionType, result_components, result_size);
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
  return construct_Expansion(&ExpansionType, result_components, result_size);
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
  else if (is_PyObject_convertible_to_Float(self)) {
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
  return PyFloat_FromDouble(sum_components(self->size, self->components));
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
  return construct_Expansion(&ExpansionType, result_components, result_size);
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
  return construct_Expansion(&ExpansionType, result_components, self->size);
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
  return construct_Expansion(&ExpansionType, result_components, result_size);
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
  return construct_Expansion(&ExpansionType, result_components, result_size);
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
  else if (is_PyObject_convertible_to_Float(self)) {
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

static PyObject *vectors_cross_product(PyObject *Py_UNUSED(self),
                                       PyObject *args) {
  double first_start_x, first_start_y, first_end_x, first_end_y, second_start_x,
      second_start_y, second_end_x, second_end_y;
  if (!PyArg_ParseTuple(args, "dddddddd", &first_start_x, &first_start_y,
                        &first_end_x, &first_end_y, &second_start_x,
                        &second_start_y, &second_end_x, &second_end_y))
    return NULL;
  double components[16];
  size_t result_size = vectors_cross_product_impl(
      first_start_x, first_start_y, first_end_x, first_end_y, second_start_x,
      second_start_y, second_end_x, second_end_y, components);
  size_t offset = 0;
  for (; offset < result_size - 1 && !components[offset]; ++offset)
    ;
  result_size -= offset;
  double *result_components = PyMem_RawCalloc(result_size, sizeof(double));
  if (!result_components) return PyErr_NoMemory();
  copy_components(&components[offset], result_size, result_components);
  return (PyObject *)construct_Expansion(&ExpansionType, result_components,
                                         result_size);
}

static PyMethodDef _shewchuk_methods[] = {
    {"vectors_cross_product", (PyCFunction)vectors_cross_product, METH_VARARGS,
     PyDoc_STR("Computes cross product of two vectors given their endpoints "
               "coordinates.")},
    {NULL, NULL},
};

static PyModuleDef _shewchuk_module = {
    PyModuleDef_HEAD_INIT,
    .m_doc = PyDoc_STR("Robust floating point operations."),
    .m_methods = _shewchuk_methods,
    .m_name = "shewchuk",
    .m_size = -1,
};

PyMODINIT_FUNC PyInit__shewchuk(void) {
  PyObject *result;
  if (PyType_Ready(&ExpansionType) < 0) return NULL;
  result = PyModule_Create(&_shewchuk_module);
  if (result == NULL) return NULL;
  Py_INCREF(&ExpansionType);
  if (PyModule_AddObject(result, "Expansion", (PyObject *)&ExpansionType) < 0) {
    Py_DECREF(&ExpansionType);
    Py_DECREF(result);
    return NULL;
  }
  return result;
}
