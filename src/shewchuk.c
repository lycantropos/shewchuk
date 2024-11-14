#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <float.h>
#include <pymath.h>
#include <stddef.h>
#include <stdio.h>
#include <structmember.h>
#define EPSILON (DBL_EPSILON / 2.0)
#define PY39_OR_MORE PY_VERSION_HEX >= 0x03090000

static int to_sign(double value) {
  return value > 0.0 ? 1 : (value == 0.0 ? 0 : -1);
}

static const size_t BIT_LENGTHS_TABLE[32] = {0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4,
                                             4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5,
                                             5, 5, 5, 5, 5, 5, 5, 5, 5, 5};

size_t bit_length(const size_t value) {
  size_t result = 0;
  size_t step = value;
  while (step >= 32) {
    result += 6;
    step >>= 6;
  }
  result += BIT_LENGTHS_TABLE[step];
  return result;
}

size_t ceil_log2(const size_t value) {
  return bit_length(value) + !(value & (value - 1));
}

static void swap(double **first, double **second) {
  double *temp = *first;
  *first = *second;
  *second = temp;
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
  static const double splitter = 1 + (1 << (size_t)((DBL_MANT_DIG + 1) / 2));
  double base = splitter * value;
  double high = base - (base - value);
  double low = value - high;
  *result_high = high;
  *result_low = low;
}

static void square(double value, double *result_head, double *result_tail) {
  double head = value * value;
  double value_high, value_low;
  split(value, &value_high, &value_low);
  double first_error = head - value_high * value_high;
  double second_error = first_error - (value_high + value_high) * value_low;
  double tail = value_low * value_low - second_error;
  *result_head = head;
  *result_tail = tail;
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

static size_t add_double_in_place(size_t size, double *components, double value,
                                  double *result) {
  size_t result_size = 0;
  double accumulator = value;
  for (size_t index = 0; index < size; index++) {
    double tail;
    two_add(accumulator, components[index], &accumulator, &tail);
    if (tail != 0.0) result[result_size++] = tail;
  }
  if (accumulator != 0.0 || result_size == 0)
    result[result_size++] = accumulator;
  return result_size;
}

static int add_double(size_t size, double *components, double value,
                      size_t *result_size, double **result) {
  *result = (double *)PyMem_Malloc((size + 1) * sizeof(double));
  if (*result == NULL) {
    PyErr_NoMemory();
    return -1;
  }
  *result_size = add_double_in_place(size, components, value, *result);
  if (!PyMem_Resize(*result, double, *result_size)) {
    PyErr_NoMemory();
    return -1;
  }
  return 0;
}

static int are_components_equal(const size_t left_size,
                                const double *const left,
                                const size_t right_size,
                                const double *const right) {
  if (left_size != right_size) return 0;
  for (size_t index = 0; index < left_size; ++index)
    if (left[index] != right[index]) return 0;
  return 1;
}

static int are_components_lesser_than(const size_t left_size,
                                      const double *const left,
                                      const size_t right_size,
                                      const double *const right) {
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

static int are_components_equal_to_double(const size_t size,
                                          double *const components,
                                          double value) {
  return size == 1 && components[0] == value;
}

static int are_components_lesser_than_double(const size_t size,
                                             double *const components,
                                             double value) {
  return components[size - 1] < value ||
         (size > 1 && components[size - 1] == value &&
          components[size - 2] < 0.0);
}

static int components_to_fractions(const size_t size, double *const components,
                                   size_t *result_size,
                                   double **result_components) {
  *result_components = (double *)PyMem_Malloc(size * sizeof(double));
  if (*result_components == NULL) {
    PyErr_NoMemory();
    return -1;
  }
  double _;
  size_t index;
  for (index = 0; index < size; ++index)
    (*result_components)[index] = modf(components[index], &_);
  for (; index > 1 && (*result_components)[index - 1] == 0.0; --index) {
  }
  *result_size = index;
  if (!PyMem_Resize(*result_components, double, *result_size)) {
    PyErr_NoMemory();
    return -1;
  }
  return 0;
}

static double components_to_accumulated_fraction(
    const size_t size, const double *const components) {
  double _;
  double result = 0.0;
  size_t index;
  for (index = 0; index < size; ++index) {
    double component_fraction = modf(components[index], &_);
    if (component_fraction == 0.0) break;
    result += component_fraction;
  }
  assert(result == 0.0 && index == 0 ||
         result == modf(components[index - 1], &_));
  return result;
}

static int do_components_have_fraction(const double *const components) {
  double _;
  return modf(components[0], &_) != 0.0;
}

static int is_double_lesser_than_components(double value, const size_t size,
                                            double *const components) {
  return value < components[size - 1] ||
         (size > 1 && value == components[size - 1] &&
          components[size - 2] > 0.0);
}

static size_t compress_components_single(size_t size, double *components) {
  size_t bottom = size - 1;
  double accumulator = components[bottom];
  for (Py_ssize_t index = (Py_ssize_t)(bottom)-1; index >= 0; --index) {
    double head, tail;
    two_add(accumulator, components[index], &head, &tail);
    if (tail != 0.0) {
      components[bottom--] = head;
      accumulator = tail;
    } else
      accumulator = head;
  }
  size_t top = 0;
  for (size_t index = bottom + 1; index < size; ++index) {
    double head, tail;
    two_add(components[index], accumulator, &head, &tail);
    if (tail != 0.0) components[top++] = tail;
    accumulator = head;
  }
  if (accumulator != 0.0 || top == 0) components[top++] = accumulator;
  return top;
}

static void copy_components(size_t size, double *source, double *destination) {
  for (size_t index = 0; index < size; ++index)
    destination[index] = source[index];
}

static size_t compress_components(size_t size, double *components) {
  const size_t original_size = size;
  double *next_components = (double *)PyMem_Malloc(size * sizeof(double));
  if (next_components == NULL) {
    PyErr_NoMemory();
    return 0;
  }
  copy_components(size, components, next_components);
  for (size_t step = 0; step < original_size; ++step) {
    const size_t next_size = compress_components_single(size, next_components);
    if (are_components_equal(next_size, next_components, size, components))
      break;
    size = next_size;
    copy_components(size, next_components, components);
  }
  PyMem_Free(next_components);
  return size;
}

static size_t negate_components(size_t size, double *components,
                                double *result_components) {
  for (size_t index = 0; index < size; ++index)
    result_components[index] = -components[index];
  return size;
}

double sum_components(size_t size, double *components) {
  double result = components[0];
  for (size_t index = 1; index < size; ++index) result += components[index];
  return result;
}

static size_t scale_components_in_place(size_t size, double *components,
                                        double scalar, double *result) {
  double scalar_high, scalar_low;
  split(scalar, &scalar_high, &scalar_low);
  double accumulator, tail;
  two_multiply_presplit(components[0], scalar, scalar_high, scalar_low,
                        &accumulator, &tail);
  size_t result_size = 0;
  if (tail != 0.0) result[result_size++] = tail;
  for (size_t index = 1; index < size; ++index) {
    double product, product_tail;
    two_multiply_presplit(components[index], scalar, scalar_high, scalar_low,
                          &product, &product_tail);
    double interim;
    two_add(accumulator, product_tail, &interim, &tail);
    if (tail != 0.0) result[result_size++] = tail;
    fast_two_add(product, interim, &accumulator, &tail);
    if (tail != 0.0) result[result_size++] = tail;
  }
  if (accumulator != 0.0 || result_size == 0)
    result[result_size++] = accumulator;
  return result_size;
}

static int scale_components(size_t size, double *components, double scalar,
                            size_t *result_size, double **result) {
  *result = (double *)PyMem_Malloc(2 * size * sizeof(double));
  if (*result == NULL) {
    PyErr_NoMemory();
    return -1;
  }
  *result_size = scale_components_in_place(size, components, scalar, *result);
  if (!PyMem_Resize(*result, double, *result_size)) {
    PyErr_NoMemory();
    return -1;
  }
  return 0;
}

static size_t subtract_double_in_place(size_t minuend_size, double *minuend,
                                       double subtrahend, double *result) {
  return add_double_in_place(minuend_size, minuend, -subtrahend, result);
}

static int subtract_double(size_t minuend_size, double *minuend,
                           double subtrahend, size_t *result_size,
                           double **result) {
  *result = (double *)PyMem_Malloc((minuend_size + 1) * sizeof(double));
  if (*result == NULL) {
    PyErr_NoMemory();
    return -1;
  }
  *result_size =
      subtract_double_in_place(minuend_size, minuend, subtrahend, *result);
  if (!PyMem_Resize(*result, double, *result_size)) {
    PyErr_NoMemory();
    return -1;
  }
  return 0;
}

static size_t subtract_from_double_in_place(double minuend,
                                            size_t subtrahend_size,
                                            double *subtrahend,
                                            double *result) {
  size_t result_size = 0;
  double accumulator = minuend;
  for (size_t index = 0; index < subtrahend_size; index++) {
    double tail;
    two_add(accumulator, -subtrahend[index], &accumulator, &tail);
    if (tail != 0.0) result[result_size++] = tail;
  }
  if (accumulator != 0.0 || result_size == 0)
    result[result_size++] = accumulator;
  return result_size;
}

static int subtract_from_double(double minuend, size_t subtrahend_size,
                                double *subtrahend, size_t *result_size,
                                double **result) {
  *result = (double *)PyMem_Malloc((subtrahend_size + 1) * sizeof(double));
  if (*result == NULL) {
    PyErr_NoMemory();
    return -1;
  }
  *result_size = subtract_from_double_in_place(minuend, subtrahend_size,
                                               subtrahend, *result);
  if (!PyMem_Resize(*result, double, *result_size)) {
    PyErr_NoMemory();
    return -1;
  }
  return 0;
}

static size_t add_components_in_place(const size_t left_size,
                                      const double *const left,
                                      const size_t right_size,
                                      const double *const right,
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
    if (tail != 0.0) result[result_size++] = tail;
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
      if (tail != 0.0) result[result_size++] = tail;
    }
  }
  for (; left_index < left_size; ++left_index) {
    two_add(accumulator, left[left_index], &head, &tail);
    accumulator = head;
    if (tail != 0.0) result[result_size++] = tail;
  }
  for (; right_index < right_size; ++right_index) {
    two_add(accumulator, right[right_index], &head, &tail);
    accumulator = head;
    if (tail != 0.0) result[result_size++] = tail;
  }
  if (accumulator != 0.0 || result_size == 0)
    result[result_size++] = accumulator;
  return result_size;
}

static int add_components(const size_t left_size, const double *const left,
                          const size_t right_size, const double *const right,
                          size_t *result_size, double **result) {
  *result = (double *)PyMem_Malloc((left_size + right_size) * sizeof(double));
  if (*result == NULL) {
    PyErr_NoMemory();
    return -1;
  }
  *result_size =
      add_components_in_place(left_size, left, right_size, right, *result);
  if (!PyMem_Resize(*result, double, *result_size)) {
    PyErr_NoMemory();
    return -1;
  }
  return 0;
}

static size_t multiply_components_in_place(size_t left_size, double *left,
                                           size_t right_size, double *right,
                                           double *result) {
  double *next_components =
      (double *)PyMem_Malloc(2 * left_size * (right_size - 1) * sizeof(double));
  if (next_components == NULL) {
    PyErr_NoMemory();
    return 0;
  }
  double *step = (double *)PyMem_Malloc(2 * left_size * sizeof(double));
  if (step == NULL) {
    PyMem_Free(next_components);
    PyErr_NoMemory();
    return 0;
  }
  size_t result_size =
      scale_components_in_place(left_size, left, right[0], result);
  for (size_t index = 1; index < right_size; ++index) {
    size_t step_size =
        scale_components_in_place(left_size, left, right[index], step);
    copy_components(result_size, result, next_components);
    result_size = add_components_in_place(result_size, next_components,
                                          step_size, step, result);
  }
  PyMem_Free(next_components);
  PyMem_Free(step);
  return result_size;
}

static int multiply_components(size_t left_size, double *left,
                               size_t right_size, double *right,
                               size_t *result_size, double **result) {
  *result = (double *)PyMem_Malloc(2 * left_size * right_size * sizeof(double));
  if (*result == NULL) {
    PyErr_NoMemory();
    return -1;
  }
  *result_size = left_size < right_size
                     ? multiply_components_in_place(right_size, right,
                                                    left_size, left, *result)
                     : multiply_components_in_place(left_size, left, right_size,
                                                    right, *result);
  if (!PyMem_Resize(*result, double, *result_size)) {
    PyErr_NoMemory();
    return -1;
  }
  return 0;
}

static size_t subtract_components_in_place(size_t minuend_size, double *minuend,
                                           size_t subtrahend_size,
                                           double *subtrahend, double *result) {
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
    if (tail != 0.0) result[result_size++] = tail;
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
      if (tail != 0.0) result[result_size++] = tail;
    }
  }
  while (minuend_index < minuend_size) {
    two_add(accumulator, minuend_component, &head, &tail);
    minuend_component = minuend[++minuend_index];
    accumulator = head;
    if (tail != 0.0) result[result_size++] = tail;
  }
  while (subtrahend_index < subtrahend_size) {
    two_add(accumulator, subtrahend_component, &head, &tail);
    subtrahend_component = -subtrahend[++subtrahend_index];
    accumulator = head;
    if (tail != 0.0) result[result_size++] = tail;
  }
  if (accumulator != 0.0 || result_size == 0)
    result[result_size++] = accumulator;
  return result_size;
}

static int subtract_components(size_t minuend_size, double *minuend,
                               size_t subtrahend_size, double *subtrahend,
                               size_t *result_size, double **result) {
  *result =
      (double *)PyMem_Malloc((minuend_size + subtrahend_size) * sizeof(double));
  if (*result == NULL) {
    PyErr_NoMemory();
    return -1;
  }
  *result_size = subtract_components_in_place(
      minuend_size, minuend, subtrahend_size, subtrahend, *result);
  if (!PyMem_Resize(*result, double, *result_size)) {
    PyErr_NoMemory();
    return -1;
  }
  return 0;
}

static int invert_components(const size_t size, double *const components,
                             size_t *result_size, double **result_components) {
  double first_approximation = 1.0 / components[size - 1];
  double first_approximation_high, first_approximation_low;
  split(first_approximation, &first_approximation_high,
        &first_approximation_low);
  if (!isfinite(first_approximation_high)) {
    PyObject *components_object = PyList_New((Py_ssize_t)(size));
    if (components_object == NULL) return -1;
    for (size_t index = 0; index < size; ++index) {
      PyObject *component_object = PyFloat_FromDouble(components[index]);
      if (component_object == NULL) {
        Py_DECREF(components_object);
        return -1;
      }
      PyList_SET_ITEM(components_object, (Py_ssize_t)(index), component_object);
    }
    PyErr_Format(PyExc_OverflowError,
                 "Components %R are not finitely invertible.",
                 components_object);
    Py_DECREF(components_object);
    return -1;
  }
  size_t iterations_count = 6 + ceil_log2(size);
  size_t max_result_size = (iterations_count - 1) * iterations_count *
                           (1 + 2 * size * (2 * iterations_count - 1) / 3);
  double *step_components =
      (double *)PyMem_Malloc(max_result_size * sizeof(double));
  if (step_components == NULL) return -1;
  step_components[0] = first_approximation;
  size_t step_size = 1;
  double *negated_components = (double *)PyMem_Malloc(size * sizeof(double));
  if (negated_components == NULL) {
    PyMem_Free(step_components);
    return -1;
  }
  size_t negated_size = negate_components(size, components, negated_components);
  double *tmp_components =
      (double *)PyMem_Malloc(max_result_size * sizeof(double));
  if (tmp_components == NULL) {
    PyMem_Free(negated_components);
    PyMem_Free(step_components);
    return -1;
  }
  double *other_tmp_components =
      (double *)PyMem_Malloc(max_result_size * sizeof(double));
  if (other_tmp_components == NULL) {
    PyMem_Free(tmp_components);
    PyMem_Free(negated_components);
    PyMem_Free(step_components);
    return -1;
  }
  for (size_t _ = 0; _ < iterations_count; ++_) {
    size_t tmp_size =
        multiply_components_in_place(step_size, step_components, negated_size,
                                     negated_components, tmp_components);
    if (tmp_size == 0) {
      PyMem_Free(other_tmp_components);
      PyMem_Free(tmp_components);
      PyMem_Free(negated_components);
      PyMem_Free(step_components);
      return -1;
    }
    size_t other_tmp_size = add_double_in_place(tmp_size, tmp_components, 2.0,
                                                other_tmp_components);
    step_size = multiply_components_in_place(other_tmp_size,
                                             other_tmp_components, step_size,
                                             step_components, tmp_components);
    if (step_size == 0) {
      PyMem_Free(other_tmp_components);
      PyMem_Free(tmp_components);
      PyMem_Free(negated_components);
      PyMem_Free(step_components);
      return -1;
    }
    swap(&step_components, &tmp_components);
  }
  PyMem_Free(other_tmp_components);
  PyMem_Free(tmp_components);
  PyMem_Free(negated_components);
  step_size = compress_components(step_size, step_components);
  if (step_size == 0) return 0;
  if (!PyMem_Resize(step_components, double, step_size)) {
    PyErr_NoMemory();
    return 0;
  }
  *result_size = step_size;
  *result_components = step_components;
  return 0;
}

static int divide_components(size_t dividend_size, double *dividend,
                             size_t divisor_size, double *divisor,
                             size_t *result_size, double **result) {
  double *divisor_reciprocal = NULL;
  size_t divisor_reciprocal_size = 0;
  if (invert_components(divisor_size, divisor, &divisor_reciprocal_size,
                        &divisor_reciprocal) < 0)
    return -1;
  if (multiply_components(divisor_reciprocal_size, divisor_reciprocal,
                          dividend_size, dividend, result_size, result) < 0)
    return -1;
  PyMem_Free(divisor_reciprocal);
  if (*result_size == 0) {
    PyMem_Free(*result);
    return -1;
  }
  *result_size = compress_components(*result_size, *result);
  if (!PyMem_Resize(*result, double, *result_size)) return -1;
  return 0;
}

static int is_py_long_odd(PyObject *value) {
  PyObject *one = PyLong_FromLong(1);
  if (one == NULL) return -1;
  PyObject *first_bit = PyNumber_And(value, one);
  Py_DECREF(one);
  if (first_bit == NULL) return -1;
  const int result = PyObject_IsTrue(first_bit);
  Py_DECREF(first_bit);
  return result;
}

static int is_py_long_one(PyObject *value) {
  PyObject *one = PyLong_FromLong(1);
  if (one == NULL) return -1;
  int result = PyObject_RichCompareBool(value, one, Py_EQ);
  Py_DECREF(one);
  return result;
}

static PyObject *components_to_py_long(const size_t size,
                                       const double *const components) {
  PyObject *result = PyLong_FromDouble(components[size - 1]);
  if (result == NULL) return NULL;
  for (size_t offset = 2; offset <= size; ++offset) {
    PyObject *component_integer = PyLong_FromDouble(components[size - offset]);
    if (component_integer == NULL) {
      Py_DECREF(result);
      return NULL;
    }
    if (PyObject_Not(component_integer)) break;
    PyObject *tmp = result;
    result = PyNumber_InPlaceAdd(result, component_integer);
    Py_DECREF(tmp);
    Py_DECREF(component_integer);
    if (result == NULL) return NULL;
  }
  return result;
}

static int py_long_to_components(PyObject *value, size_t *size,
                                 double **components) {
  if (PyObject_Not(value)) {
    *components = (double *)PyMem_Malloc(sizeof(double));
    if (*components == NULL) {
      PyErr_NoMemory();
      return -1;
    }
    *size = 1;
    *components[0] = 0.0;
    return 0;
  }
  PyObject *rest = PyNumber_Long(value);
  if (rest == NULL) return -1;
  double component = PyFloat_AsDouble(rest);
  if (component == -1.0 && PyErr_Occurred()) {
    Py_DECREF(rest);
    return -1;
  }
  assert(component >= 1.0 || component <= -1.0);
  int exponent;
  frexp(component, &exponent);
  size_t max_components_count = 1 + ((size_t)exponent - 1) / DBL_MANT_DIG;
  double *reversed_components =
      (double *)PyMem_Malloc(max_components_count * sizeof(double));
  if (reversed_components == NULL) {
    PyErr_NoMemory();
    return -1;
  }
  size_t result_size = 0;
  for (;;) {
    reversed_components[result_size++] = component;
    assert(result_size <= max_components_count);
    PyObject *component_integer = PyLong_FromDouble(component);
    PyObject *tmp = rest;
    rest = PyNumber_InPlaceSubtract(rest, component_integer);
    Py_DECREF(tmp);
    if (rest == NULL) {
      PyMem_Free(reversed_components);
      return -1;
    }
    int is_zero = PyObject_Not(rest);
    if (is_zero < 0) {
      PyMem_Free(reversed_components);
      Py_DECREF(rest);
      return -1;
    }
    if (is_zero) break;
    component = PyFloat_AsDouble(rest);
    if (component == -1.0 && PyErr_Occurred()) {
      PyMem_Free(reversed_components);
      Py_DECREF(rest);
      return -1;
    }
  }
  Py_DECREF(rest);
  *components = (double *)PyMem_Malloc(result_size * sizeof(double));
  if (*components == NULL) {
    PyMem_Free(reversed_components);
    PyErr_NoMemory();
    return -1;
  }
  *size = result_size;
  for (size_t index = 0; index < result_size; ++index) {
    (*components)[result_size - 1 - index] = reversed_components[index];
  }
  PyMem_Free(reversed_components);
  return 0;
}

static int Rational_to_components(PyObject *value, size_t *size,
                                  double **components) {
  PyObject *denominator = PyObject_GetAttrString(value, "denominator");
  if (denominator == NULL) return -1;
  PyObject *numerator = PyObject_GetAttrString(value, "numerator");
  if (numerator == NULL) {
    Py_DECREF(denominator);
    return -1;
  }
  assert(PyLong_Check(denominator));
  assert(PyLong_Check(numerator));
  double *numerator_components;
  size_t numerator_size;
  if (py_long_to_components(numerator, &numerator_size, &numerator_components) <
      0) {
    Py_DECREF(numerator);
    Py_DECREF(denominator);
    return -1;
  }
  Py_DECREF(numerator);
  int is_denominator_one = is_py_long_one(denominator);
  if (is_denominator_one < 0) {
    PyMem_Free(numerator_components);
    Py_DECREF(denominator);
    return -1;
  } else if (is_denominator_one) {
    *components = numerator_components;
    *size = numerator_size;
    return 0;
  }
  double *denominator_components;
  size_t denominator_size;
  if (py_long_to_components(denominator, &denominator_size,
                            &denominator_components) < 0) {
    PyMem_Free(numerator_components);
    Py_DECREF(denominator);
    return -1;
  }
  Py_DECREF(denominator);
  if (divide_components(numerator_size, numerator_components, denominator_size,
                        denominator_components, size, components) < 0) {
    PyMem_Free(numerator_components);
    PyMem_Free(denominator_components);
    return -1;
  }
  PyMem_Free(numerator_components);
  PyMem_Free(denominator_components);
  return 0;
}

static int are_components_equal_to_py_long(const size_t size,
                                           const double *const components,
                                           PyObject *value) {
  if (do_components_have_fraction(components)) return 0;
  PyObject *components_integer = components_to_py_long(size, components);
  int result = PyObject_RichCompareBool(components_integer, value, Py_EQ);
  Py_DECREF(components_integer);
  return result;
}

static int are_components_lesser_than_py_long(const size_t size,
                                              const double *const components,
                                              PyObject *integral) {
  PyObject *components_integer = components_to_py_long(size, components);
  int components_integer_is_lesser_than_integral =
      PyObject_RichCompareBool(components_integer, integral, Py_LT);
  return components_integer_is_lesser_than_integral < 0
             ? components_integer_is_lesser_than_integral
             : (components_integer_is_lesser_than_integral ||
                (PyObject_RichCompareBool(components_integer, integral,
                                          Py_EQ) &&
                 components_to_accumulated_fraction(size, components) < 0.0));
}

static int is_py_long_lesser_than_components(PyObject *integral,
                                             const size_t size,
                                             const double *const components) {
  PyObject *components_integer = components_to_py_long(size, components);
  int components_integer_is_greater_than_integral =
      PyObject_RichCompareBool(components_integer, integral, Py_GT);
  return components_integer_is_greater_than_integral < 0
             ? components_integer_is_greater_than_integral
             : (components_integer_is_greater_than_integral ||
                (PyObject_RichCompareBool(components_integer, integral,
                                          Py_EQ) &&
                 components_to_accumulated_fraction(size, components) > 0.0));
}

static int are_components_equal_to_Rational(const size_t size,
                                            const double *const components,
                                            PyObject *value) {
  double *rational_components;
  size_t rational_size;
  if (Rational_to_components(value, &rational_size, &rational_components) < 0)
    return -1;
  int result = are_components_equal(size, components, rational_size,
                                    rational_components);
  PyMem_Free(rational_components);
  return result;
}

static int are_components_lesser_than_Rational(const size_t size,
                                               const double *const components,
                                               PyObject *value) {
  double *rational_components;
  size_t rational_size;
  if (Rational_to_components(value, &rational_size, &rational_components) < 0)
    return -1;
  int result = are_components_lesser_than(size, components, rational_size,
                                          rational_components);
  PyMem_Free(rational_components);
  return result;
}

static int is_Rational_lesser_than_components(PyObject *value,
                                              const size_t size,
                                              const double *const components) {
  double *rational_components;
  size_t rational_size;
  if (Rational_to_components(value, &rational_size, &rational_components) < 0)
    return -1;
  int result = are_components_lesser_than(rational_size, rational_components,
                                          size, components);
  PyMem_Free(rational_components);
  return result;
}

static void cross_product(double first_dx, double first_dy, double second_dx,
                          double second_dy, double *head, double *first_tail,
                          double *second_tail, double *third_tail) {
  double first_dx_second_dy_head, first_dx_second_dy_tail;
  two_multiply(first_dx, second_dy, &first_dx_second_dy_head,
               &first_dx_second_dy_tail);
  double second_dx_first_dy_head, second_dx_first_dy_tail;
  two_multiply(second_dx, first_dy, &second_dx_first_dy_head,
               &second_dx_first_dy_tail);
  two_two_subtract(first_dx_second_dy_head, first_dx_second_dy_tail,
                   second_dx_first_dy_head, second_dx_first_dy_tail, head,
                   first_tail, second_tail, third_tail);
}

static void squared_length(double dx, double dy, double *head,
                           double *first_tail, double *second_tail,
                           double *third_tail) {
  double dx_squared_head, dx_squared_tail;
  square(dx, &dx_squared_head, &dx_squared_tail);
  double dy_squared_head, dy_squared_tail;
  square(dy, &dy_squared_head, &dy_squared_tail);
  two_two_add(dx_squared_head, dx_squared_tail, dy_squared_head,
              dy_squared_tail, head, first_tail, second_tail, third_tail);
}

static size_t scale_by_squared_length(size_t size, double *components,
                                      double dx, double dy, double *result) {
  double dx_components[8], dx_squared_components[16], dy_components[8],
      dy_squared_components[16];
  size_t dx_components_size =
      scale_components_in_place(size, components, dx, dx_components);
  size_t dx_squared_components_size = scale_components_in_place(
      dx_components_size, dx_components, dx, dx_squared_components);
  size_t dy_components_size =
      scale_components_in_place(size, components, dy, dy_components);
  size_t dy_squared_components_size = scale_components_in_place(
      dy_components_size, dy_components, dy, dy_squared_components);
  return add_components_in_place(
      dx_squared_components_size, dx_squared_components,
      dy_squared_components_size, dy_squared_components, result);
}

static size_t add_extras(
    size_t final_size, double **final_components,
    double **accumulated_components, double first_dx, double first_dx_tail,
    double first_dy, double first_dy_tail, double second_dx,
    double second_dx_tail, double second_dy, double second_dy_tail,
    double third_dx, double third_dx_tail, double third_dy,
    double third_dy_tail, double second_third_cross_product[4],
    double second_squared_length[4], double third_squared_length[4]) {
  double first_buffer_16[16], second_buffer_16[16], third_buffer_16[16];
  double first_buffer_32[32], second_buffer_32[32], buffer_48[48];
  size_t first_buffer_16_limit, second_buffer_16_limit, third_buffer_16_limit;
  size_t first_buffer_32_limit, second_buffer_32_limit, buffer_48_limit;
  size_t first_dx_tail_second_third_cross_product_size = 0;
  double first_dx_tail_second_third_cross_product[8];
  if (first_dx_tail != 0.0) {
    first_dx_tail_second_third_cross_product_size =
        scale_components_in_place(4, second_third_cross_product, first_dx_tail,
                                  first_dx_tail_second_third_cross_product);
    first_buffer_16_limit =
        scale_components_in_place(first_dx_tail_second_third_cross_product_size,
                                  first_dx_tail_second_third_cross_product,
                                  2.0 * first_dx, first_buffer_16);
    double first_dx_tail_third_squared_length[8];
    size_t first_dx_tail_third_squared_length_size =
        scale_components_in_place(4, third_squared_length, first_dx_tail,
                                  first_dx_tail_third_squared_length);
    second_buffer_16_limit = scale_components_in_place(
        first_dx_tail_third_squared_length_size,
        first_dx_tail_third_squared_length, second_dy, second_buffer_16);
    double first_dx_tail_second_squared_length[8];
    size_t first_dx_tail_second_squared_length_size =
        scale_components_in_place(4, second_squared_length, first_dx_tail,
                                  first_dx_tail_second_squared_length);
    third_buffer_16_limit = scale_components_in_place(
        first_dx_tail_second_squared_length_size,
        first_dx_tail_second_squared_length, -third_dy, third_buffer_16);
    first_buffer_32_limit = add_components_in_place(
        first_buffer_16_limit, first_buffer_16, second_buffer_16_limit,
        second_buffer_16, first_buffer_32);
    buffer_48_limit = add_components_in_place(
        third_buffer_16_limit, third_buffer_16, first_buffer_32_limit,
        first_buffer_32, buffer_48);
    final_size =
        add_components_in_place(final_size, *final_components, buffer_48_limit,
                                buffer_48, *accumulated_components);
    swap(final_components, accumulated_components);
  }
  size_t first_dy_tail_second_third_cross_product_size = 0;
  double first_dy_tail_second_third_cross_product[8];
  if (first_dy_tail != 0.0) {
    first_dy_tail_second_third_cross_product_size =
        scale_components_in_place(4, second_third_cross_product, first_dy_tail,
                                  first_dy_tail_second_third_cross_product);
    first_buffer_16_limit =
        scale_components_in_place(first_dy_tail_second_third_cross_product_size,
                                  first_dy_tail_second_third_cross_product,
                                  2.0 * first_dy, first_buffer_16);
    double first_dy_tail_second_squared_length[8];
    size_t first_dy_tail_second_squared_length_size =
        scale_components_in_place(4, second_squared_length, first_dy_tail,
                                  first_dy_tail_second_squared_length);
    second_buffer_16_limit = scale_components_in_place(
        first_dy_tail_second_squared_length_size,
        first_dy_tail_second_squared_length, third_dx, second_buffer_16);
    double first_dy_tail_third_squared_length[8];
    size_t first_dy_tail_third_squared_length_size =
        scale_components_in_place(4, third_squared_length, first_dy_tail,
                                  first_dy_tail_third_squared_length);
    third_buffer_16_limit = scale_components_in_place(
        first_dy_tail_third_squared_length_size,
        first_dy_tail_third_squared_length, -second_dx, third_buffer_16);
    first_buffer_32_limit = add_components_in_place(
        first_buffer_16_limit, first_buffer_16, second_buffer_16_limit,
        second_buffer_16, first_buffer_32);
    buffer_48_limit = add_components_in_place(
        third_buffer_16_limit, third_buffer_16, first_buffer_32_limit,
        first_buffer_32, buffer_48);
    final_size =
        add_components_in_place(final_size, *final_components, buffer_48_limit,
                                buffer_48, *accumulated_components);
    swap(final_components, accumulated_components);
  }
  double dx_tail_dy_head_head, dx_head_dy_tail_head;
  double dx_tail_dy_head_tail, dx_head_dy_tail_tail;
  double buffer_8[8], buffer_64[64];
  size_t buffer_8_limit, buffer_64_limit;
  double first_buffer_4[4], second_buffer_4[4];
  if (first_dx_tail != 0.0 || first_dy_tail != 0.0) {
    size_t second_third_cross_product_bodies_size,
        second_third_cross_product_tails_size;
    double second_third_cross_product_bodies[8],
        second_third_cross_product_tails[4];
    if (second_dx_tail != 0.0 || second_dy_tail != 0.0 ||
        third_dx_tail != 0.0 || third_dy_tail != 0.0) {
      two_multiply(second_dx_tail, third_dy, &dx_tail_dy_head_head,
                   &dx_tail_dy_head_tail);
      two_multiply(second_dx, third_dy_tail, &dx_head_dy_tail_head,
                   &dx_head_dy_tail_tail);
      two_two_add(dx_tail_dy_head_head, dx_tail_dy_head_tail,
                  dx_head_dy_tail_head, dx_head_dy_tail_tail,
                  &first_buffer_4[3], &first_buffer_4[2], &first_buffer_4[1],
                  &first_buffer_4[0]);
      two_multiply(third_dx_tail, -second_dy, &dx_tail_dy_head_head,
                   &dx_tail_dy_head_tail);
      two_multiply(third_dx, -second_dy_tail, &dx_head_dy_tail_head,
                   &dx_head_dy_tail_tail);
      two_two_add(dx_tail_dy_head_head, dx_tail_dy_head_tail,
                  dx_head_dy_tail_head, dx_head_dy_tail_tail,
                  &second_buffer_4[3], &second_buffer_4[2], &second_buffer_4[1],
                  &second_buffer_4[0]);
      second_third_cross_product_bodies_size =
          add_components_in_place(4, first_buffer_4, 4, second_buffer_4,
                                  second_third_cross_product_bodies);
      two_multiply(second_dx_tail, third_dy_tail, &dx_tail_dy_head_head,
                   &dx_tail_dy_head_tail);
      two_multiply(third_dx_tail, second_dy_tail, &dx_head_dy_tail_head,
                   &dx_head_dy_tail_tail);
      two_two_subtract(dx_tail_dy_head_head, dx_tail_dy_head_tail,
                       dx_head_dy_tail_head, dx_head_dy_tail_tail,
                       &second_third_cross_product_tails[3],
                       &second_third_cross_product_tails[2],
                       &second_third_cross_product_tails[1],
                       &second_third_cross_product_tails[0]);
      second_third_cross_product_tails_size = 4;
    } else {
      second_third_cross_product_bodies[0] = 0.0;
      second_third_cross_product_bodies_size = 1;
      second_third_cross_product_tails[0] = 0.0;
      second_third_cross_product_tails_size = 1;
    }
    if (first_dx_tail != 0.0) {
      first_buffer_16_limit = scale_components_in_place(
          first_dx_tail_second_third_cross_product_size,
          first_dx_tail_second_third_cross_product, first_dx_tail,
          first_buffer_16);
      double first_dx_tail_second_third_cross_product_bodies[16];
      size_t first_dx_tail_second_third_cross_product_bodies_size =
          scale_components_in_place(
              second_third_cross_product_bodies_size,
              second_third_cross_product_bodies, first_dx_tail,
              first_dx_tail_second_third_cross_product_bodies);
      first_buffer_32_limit = scale_components_in_place(
          first_dx_tail_second_third_cross_product_bodies_size,
          first_dx_tail_second_third_cross_product_bodies, 2.0 * first_dx,
          first_buffer_32);
      buffer_48_limit = add_components_in_place(
          first_buffer_16_limit, first_buffer_16, first_buffer_32_limit,
          first_buffer_32, buffer_48);
      final_size = add_components_in_place(final_size, *final_components,
                                           buffer_48_limit, buffer_48,
                                           *accumulated_components);
      swap(final_components, accumulated_components);
      if (second_dy_tail != 0.0) {
        buffer_8_limit = scale_components_in_place(4, third_squared_length,
                                                   first_dx_tail, buffer_8);
        first_buffer_16_limit = scale_components_in_place(
            buffer_8_limit, buffer_8, second_dy_tail, first_buffer_16);
        final_size = add_components_in_place(
            final_size, *final_components, first_buffer_16_limit,
            first_buffer_16, *accumulated_components);
        swap(final_components, accumulated_components);
      }
      if (third_dy_tail != 0.0) {
        buffer_8_limit = scale_components_in_place(4, second_squared_length,
                                                   -first_dx_tail, buffer_8);
        first_buffer_16_limit = scale_components_in_place(
            buffer_8_limit, buffer_8, third_dy_tail, first_buffer_16);
        final_size = add_components_in_place(
            final_size, *final_components, first_buffer_16_limit,
            first_buffer_16, *accumulated_components);
        swap(final_components, accumulated_components);
      }
      first_buffer_32_limit = scale_components_in_place(
          first_dx_tail_second_third_cross_product_bodies_size,
          first_dx_tail_second_third_cross_product_bodies, first_dx_tail,
          first_buffer_32);
      double first_dx_tail_second_third_cross_product_tails[8];
      size_t first_dx_tail_second_third_cross_product_tails_size =
          scale_components_in_place(
              second_third_cross_product_tails_size,
              second_third_cross_product_tails, first_dx_tail,
              first_dx_tail_second_third_cross_product_tails);
      first_buffer_16_limit = scale_components_in_place(
          first_dx_tail_second_third_cross_product_tails_size,
          first_dx_tail_second_third_cross_product_tails, 2.0 * first_dx,
          first_buffer_16);
      second_buffer_16_limit = scale_components_in_place(
          first_dx_tail_second_third_cross_product_tails_size,
          first_dx_tail_second_third_cross_product_tails, first_dx_tail,
          second_buffer_16);
      second_buffer_32_limit = add_components_in_place(
          first_buffer_16_limit, first_buffer_16, second_buffer_16_limit,
          second_buffer_16, second_buffer_32);
      buffer_64_limit = add_components_in_place(
          first_buffer_32_limit, first_buffer_32, second_buffer_32_limit,
          second_buffer_32, buffer_64);
      final_size = add_components_in_place(final_size, *final_components,
                                           buffer_64_limit, buffer_64,
                                           *accumulated_components);
      swap(final_components, accumulated_components);
    }
    if (first_dy_tail != 0.0) {
      first_buffer_16_limit = scale_components_in_place(
          first_dy_tail_second_third_cross_product_size,
          first_dy_tail_second_third_cross_product, first_dy_tail,
          first_buffer_16);
      double first_dy_tail_second_third_cross_product_bodies[16];
      size_t first_dy_tail_second_third_cross_product_bodies_size =
          scale_components_in_place(
              second_third_cross_product_bodies_size,
              second_third_cross_product_bodies, first_dy_tail,
              first_dy_tail_second_third_cross_product_bodies);
      first_buffer_32_limit = scale_components_in_place(
          first_dy_tail_second_third_cross_product_bodies_size,
          first_dy_tail_second_third_cross_product_bodies, 2.0 * first_dy,
          first_buffer_32);
      buffer_48_limit = add_components_in_place(
          first_buffer_16_limit, first_buffer_16, first_buffer_32_limit,
          first_buffer_32, buffer_48);
      final_size = add_components_in_place(final_size, *final_components,
                                           buffer_48_limit, buffer_48,
                                           *accumulated_components);
      swap(final_components, accumulated_components);
      first_buffer_32_limit = scale_components_in_place(
          first_dy_tail_second_third_cross_product_bodies_size,
          first_dy_tail_second_third_cross_product_bodies, first_dy_tail,
          first_buffer_32);
      double first_dy_tail_second_third_cross_product_tails[8];
      size_t first_dy_tail_second_third_cross_product_tails_size =
          scale_components_in_place(
              second_third_cross_product_tails_size,
              second_third_cross_product_tails, first_dy_tail,
              first_dy_tail_second_third_cross_product_tails);
      first_buffer_16_limit = scale_components_in_place(
          first_dy_tail_second_third_cross_product_tails_size,
          first_dy_tail_second_third_cross_product_tails, 2.0 * first_dy,
          first_buffer_16);
      second_buffer_16_limit = scale_components_in_place(
          first_dy_tail_second_third_cross_product_tails_size,
          first_dy_tail_second_third_cross_product_tails, first_dy_tail,
          second_buffer_16);
      second_buffer_32_limit = add_components_in_place(
          first_buffer_16_limit, first_buffer_16, second_buffer_16_limit,
          second_buffer_16, second_buffer_32);
      buffer_64_limit = add_components_in_place(
          first_buffer_32_limit, first_buffer_32, second_buffer_32_limit,
          second_buffer_32, buffer_64);
      final_size = add_components_in_place(final_size, *final_components,
                                           buffer_64_limit, buffer_64,
                                           *accumulated_components);
      swap(final_components, accumulated_components);
    }
  }
  return final_size;
}

double adaptive_incircle_determinant_estimation(double point_x, double point_y,
                                                double first_x, double first_y,
                                                double second_x,
                                                double second_y, double third_x,
                                                double third_y,
                                                double upper_bound) {
  double first_dx = first_x - point_x;
  double second_dx = second_x - point_x;
  double third_dx = third_x - point_x;
  double first_dy = first_y - point_y;
  double second_dy = second_y - point_y;
  double third_dy = third_y - point_y;
  double second_third_cross_product[4];
  cross_product(second_dx, second_dy, third_dx, third_dy,
                &second_third_cross_product[3], &second_third_cross_product[2],
                &second_third_cross_product[1], &second_third_cross_product[0]);
  double first_components[32];
  size_t first_components_size = scale_by_squared_length(
      4, second_third_cross_product, first_dx, first_dy, first_components);
  double third_first_cross_product[4];
  cross_product(third_dx, third_dy, first_dx, first_dy,
                &third_first_cross_product[3], &third_first_cross_product[2],
                &third_first_cross_product[1], &third_first_cross_product[0]);
  double second_components[32];
  size_t second_components_size = scale_by_squared_length(
      4, third_first_cross_product, second_dx, second_dy, second_components);
  double first_second_cross_product[4];
  cross_product(first_dx, first_dy, second_dx, second_dy,
                &first_second_cross_product[3], &first_second_cross_product[2],
                &first_second_cross_product[1], &first_second_cross_product[0]);
  double third_components[32];
  size_t third_components_size = scale_by_squared_length(
      4, first_second_cross_product, third_dx, third_dy, third_components);
  double first_second_sum_components[64];
  size_t first_second_sum_size = add_components_in_place(
      first_components_size, first_components, second_components_size,
      second_components, first_second_sum_components);
  double first_buffer[1152];
  size_t final_size = add_components_in_place(
      first_second_sum_size, first_second_sum_components, third_components_size,
      third_components, first_buffer);
  double result = sum_components(final_size, first_buffer);
  static const double first_upper_bound_coefficient =
      (4.0 + 48.0 * EPSILON) * EPSILON;
  double threshold = first_upper_bound_coefficient * upper_bound;
  if ((result >= threshold) || (-result >= threshold)) return result;
  double first_dx_tail = two_subtract_tail(first_x, point_x, first_dx);
  double first_dy_tail = two_subtract_tail(first_y, point_y, first_dy);
  double second_dx_tail = two_subtract_tail(second_x, point_x, second_dx);
  double second_dy_tail = two_subtract_tail(second_y, point_y, second_dy);
  double third_dx_tail = two_subtract_tail(third_x, point_x, third_dx);
  double third_dy_tail = two_subtract_tail(third_y, point_y, third_dy);
  if (first_dx_tail == 0.0 && second_dx_tail == 0.0 && third_dx_tail == 0.0 &&
      first_dy_tail == 0.0 && second_dy_tail == 0.0 && third_dy_tail == 0.0)
    return result;
  static const double second_upper_bound_coefficient =
      (44.0 + 576.0 * EPSILON) * EPSILON * EPSILON;
  static const double result_coefficient = (3.0 + 8.0 * EPSILON) * EPSILON;
  threshold = second_upper_bound_coefficient * upper_bound +
              result_coefficient * fabs(result);
  result += ((first_dx * first_dx + first_dy * first_dy) *
                 ((second_dx * third_dy_tail + third_dy * second_dx_tail) -
                  (second_dy * third_dx_tail + third_dx * second_dy_tail)) +
             2.0 * (first_dx * first_dx_tail + first_dy * first_dy_tail) *
                 (second_dx * third_dy - second_dy * third_dx)) +
            ((second_dx * second_dx + second_dy * second_dy) *
                 ((third_dx * first_dy_tail + first_dy * third_dx_tail) -
                  (third_dy * first_dx_tail + first_dx * third_dy_tail)) +
             2.0 * (second_dx * second_dx_tail + second_dy * second_dy_tail) *
                 (third_dx * first_dy - third_dy * first_dx)) +
            ((third_dx * third_dx + third_dy * third_dy) *
                 ((first_dx * second_dy_tail + second_dy * first_dx_tail) -
                  (first_dy * second_dx_tail + second_dx * first_dy_tail)) +
             2.0 * (third_dx * third_dx_tail + third_dy * third_dy_tail) *
                 (first_dx * second_dy - first_dy * second_dx));
  if ((result >= threshold) || (-result >= threshold)) return result;
  double first_squared_length[4] = {0};
  if (second_dx_tail != 0.0 || second_dy_tail != 0.0 || third_dx_tail != 0.0 ||
      third_dy_tail != 0.0)
    squared_length(first_dx, first_dy, &first_squared_length[3],
                   &first_squared_length[2], &first_squared_length[1],
                   &first_squared_length[0]);
  double second_squared_length[4] = {0};
  if (third_dx_tail != 0.0 || third_dy_tail != 0.0 || first_dx_tail != 0.0 ||
      first_dy_tail != 0.0)
    squared_length(second_dx, second_dy, &second_squared_length[3],
                   &second_squared_length[2], &second_squared_length[1],
                   &second_squared_length[0]);
  double third_squared_length[4] = {0};
  if (first_dx_tail != 0.0 || first_dy_tail != 0.0 || second_dx_tail != 0.0 ||
      second_dy_tail != 0.0)
    squared_length(third_dx, third_dy, &third_squared_length[3],
                   &third_squared_length[2], &third_squared_length[1],
                   &third_squared_length[0]);
  double *final_components = first_buffer;
  double second_buffer_1152[1152];
  double *accumulated_components = second_buffer_1152;
  final_size = add_extras(final_size, &final_components,
                          &accumulated_components, first_dx, first_dx_tail,
                          first_dy, first_dy_tail, second_dx, second_dx_tail,
                          second_dy, second_dy_tail, third_dx, third_dx_tail,
                          third_dy, third_dy_tail, second_third_cross_product,
                          second_squared_length, third_squared_length);
  final_size = add_extras(
      final_size, &final_components, &accumulated_components, second_dx,
      second_dx_tail, second_dy, second_dy_tail, third_dx, third_dx_tail,
      third_dy, third_dy_tail, first_dx, first_dx_tail, first_dy, first_dy_tail,
      third_first_cross_product, third_squared_length, first_squared_length);
  final_size = add_extras(
      final_size, &final_components, &accumulated_components, third_dx,
      third_dx_tail, third_dy, third_dy_tail, first_dx, first_dx_tail, first_dy,
      first_dy_tail, second_dx, second_dx_tail, second_dy, second_dy_tail,
      first_second_cross_product, first_squared_length, second_squared_length);
  return final_components[final_size - 1];
}

double incircle_determinant_estimation(double point_x, double point_y,
                                       double first_x, double first_y,
                                       double second_x, double second_y,
                                       double third_x, double third_y) {
  double first_dx = first_x - point_x;
  double second_dx = second_x - point_x;
  double third_dx = third_x - point_x;
  double first_dy = first_y - point_y;
  double second_dy = second_y - point_y;
  double third_dy = third_y - point_y;
  double second_dx_third_dy = second_dx * third_dy;
  double third_dx_second_dy = third_dx * second_dy;
  double first_squared_distance = first_dx * first_dx + first_dy * first_dy;
  double third_dx_first_dy = third_dx * first_dy;
  double first_dx_third_dy = first_dx * third_dy;
  double second_squared_distance =
      second_dx * second_dx + second_dy * second_dy;
  double first_dx_second_dy = first_dx * second_dy;
  double second_dx_first_dy = second_dx * first_dy;
  double third_squared_distance = third_dx * third_dx + third_dy * third_dy;
  double result =
      first_squared_distance * (second_dx_third_dy - third_dx_second_dy) +
      second_squared_distance * (third_dx_first_dy - first_dx_third_dy) +
      third_squared_distance * (first_dx_second_dy - second_dx_first_dy);
  double upper_bound = (fabs(second_dx_third_dy) + fabs(third_dx_second_dy)) *
                           first_squared_distance +
                       (fabs(third_dx_first_dy) + fabs(first_dx_third_dy)) *
                           second_squared_distance +
                       (fabs(first_dx_second_dy) + fabs(second_dx_first_dy)) *
                           third_squared_distance;
  static const double upper_bound_coefficient =
      (10.0 + 96.0 * EPSILON) * EPSILON;
  double threshold = upper_bound_coefficient * upper_bound;
  if ((result > threshold) || (-result > threshold)) return result;
  return adaptive_incircle_determinant_estimation(
      point_x, point_y, first_x, first_y, second_x, second_y, third_x, third_y,
      upper_bound);
}

double adaptive_vectors_cross_product_estimation(
    double first_start_x, double first_start_y, double first_end_x,
    double first_end_y, double second_start_x, double second_start_y,
    double second_end_x, double second_end_y, double upper_bound) {
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
  double result = sum_components(4, first_components);
  static const double first_upper_bound_coefficient =
      (2.0 + 12.0 * EPSILON) * EPSILON;
  double threshold = first_upper_bound_coefficient * upper_bound;
  if ((result >= threshold) || (-result >= threshold)) return result;
  double minuend_x_tail =
      two_subtract_tail(first_end_x, first_start_x, minuend_x);
  double subtrahend_x_tail =
      two_subtract_tail(second_end_x, second_start_x, subtrahend_x);
  double minuend_y_tail =
      two_subtract_tail(first_end_y, first_start_y, minuend_y);
  double subtrahend_y_tail =
      two_subtract_tail(second_end_y, second_start_y, subtrahend_y);
  if (minuend_x_tail == 0.0 && minuend_y_tail == 0.0 &&
      subtrahend_x_tail == 0.0 && subtrahend_y_tail == 0.0)
    return result;
  static const double second_upper_bound_coefficient =
      (9.0 + 64.0 * EPSILON) * EPSILON * EPSILON;
  static const double result_coefficient = (3.0 + 8.0 * EPSILON) * EPSILON;
  threshold = second_upper_bound_coefficient * upper_bound +
              result_coefficient * fabs(result);
  result += (minuend_x * subtrahend_y_tail + subtrahend_y * minuend_x_tail) -
            (minuend_y * subtrahend_x_tail + subtrahend_x * minuend_y_tail);
  if ((result >= threshold) || (-result >= threshold)) return result;
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
  size_t second_components_size = add_components_in_place(
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
  size_t third_components_size =
      add_components_in_place(second_components_size, second_components, 4,
                              extra_components, third_components);
  two_multiply(minuend_x_tail, subtrahend_y_tail, &minuend_x_subtrahend_y_head,
               &minuend_x_subtrahend_y_tail);
  two_multiply(minuend_y_tail, subtrahend_x_tail, &minuend_y_subtrahend_x_head,
               &minuend_y_subtrahend_x_tail);
  two_two_subtract(minuend_x_subtrahend_y_head, minuend_x_subtrahend_y_tail,
                   minuend_y_subtrahend_x_head, minuend_y_subtrahend_x_tail,
                   &extra_components[3], &extra_components[2],
                   &extra_components[1], &extra_components[0]);
  double final_components[16];
  size_t final_components_size =
      add_components_in_place(third_components_size, third_components, 4,
                              extra_components, final_components);
  return final_components[final_components_size - 1];
}

double vectors_cross_product_estimation(
    double first_start_x, double first_start_y, double first_end_x,
    double first_end_y, double second_start_x, double second_start_y,
    double second_end_x, double second_end_y) {
  double minuend =
      (first_end_x - first_start_x) * (second_end_y - second_start_y);
  double subtrahend =
      (first_end_y - first_start_y) * (second_end_x - second_start_x);
  double result = minuend - subtrahend;
  double upper_bound;
  if (minuend > 0.0) {
    if (subtrahend <= 0.0)
      return result;
    else
      upper_bound = minuend + subtrahend;
  } else if (minuend < 0.0) {
    if (subtrahend >= 0.0)
      return result;
    else
      upper_bound = -minuend - subtrahend;
  } else
    return result;
  static const double upper_bound_coefficient =
      (3.0 + 16.0 * EPSILON) * EPSILON;
  double threshold = upper_bound_coefficient * upper_bound;
  if ((result >= threshold) || (-result >= threshold)) return result;
  return adaptive_vectors_cross_product_estimation(
      first_start_x, first_start_y, first_end_x, first_end_y, second_start_x,
      second_start_y, second_end_x, second_end_y, upper_bound);
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
    size_t result_size = compress_components_single(4, first_components);
    copy_components(result_size, first_components, result);
    return result_size;
  }
  double minuend_x_tail =
      two_subtract_tail(first_end_x, first_start_x, minuend_x);
  double subtrahend_x_tail =
      two_subtract_tail(second_end_x, second_start_x, subtrahend_x);
  double minuend_y_tail =
      two_subtract_tail(first_end_y, first_start_y, minuend_y);
  double subtrahend_y_tail =
      two_subtract_tail(second_end_y, second_start_y, subtrahend_y);
  if (minuend_x_tail == 0.0 && minuend_y_tail == 0.0 &&
      subtrahend_x_tail == 0.0 && subtrahend_y_tail == 0.0) {
    size_t result_size = compress_components_single(4, first_components);
    copy_components(result_size, first_components, result);
    return result_size;
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
  if ((estimation >= threshold) || (-estimation >= threshold))
    return add_double_in_place(4, first_components, extra, result);
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
  size_t second_components_size = add_components_in_place(
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
  size_t third_components_size =
      add_components_in_place(second_components_size, second_components, 4,
                              extra_components, third_components);
  two_multiply(minuend_x_tail, subtrahend_y_tail, &minuend_x_subtrahend_y_head,
               &minuend_x_subtrahend_y_tail);
  two_multiply(minuend_y_tail, subtrahend_x_tail, &minuend_y_subtrahend_x_head,
               &minuend_y_subtrahend_x_tail);
  two_two_subtract(minuend_x_subtrahend_y_head, minuend_x_subtrahend_y_tail,
                   minuend_y_subtrahend_x_head, minuend_y_subtrahend_x_tail,
                   &extra_components[3], &extra_components[2],
                   &extra_components[1], &extra_components[0]);
  return add_components_in_place(third_components_size, third_components, 4,
                                 extra_components, result);
}

size_t vectors_cross_product_impl(double first_start_x, double first_start_y,
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

static PyObject *round_method_name = NULL;
static PyObject *zero = NULL;
static PyObject *Rational = NULL;
static PyObject *Real = NULL;

typedef struct {
  PyObject_HEAD size_t size;
  double *components;
} ExpansionObject;

static ExpansionObject *construct_expansion(PyTypeObject *cls, size_t size,
                                            double *components) {
  for (size_t index = 0; index < size; ++index)
    if (!isfinite(components[index])) {
      PyObject *component = PyFloat_FromDouble(components[index]);
      if (component) {
        PyErr_Format(PyExc_ValueError,
                     "Components should be finite, but found: %R.", component);
        Py_DECREF(component);
      }
      PyMem_Free(components);
      return NULL;
    }
  ExpansionObject *result = (ExpansionObject *)(cls->tp_alloc(cls, 0));
  if (result != NULL) {
    result->components = components;
    result->size = size;
  } else
    PyMem_Free(components);
  return result;
}

static PyTypeObject ExpansionType;

static ExpansionObject *Expansions_add(ExpansionObject *self,
                                       ExpansionObject *other) {
  double *result_components =
      (double *)PyMem_Malloc((self->size + other->size) * sizeof(double));
  size_t result_size;
  if (add_components(self->size, self->components, other->size,
                     other->components, &result_size, &result_components) < 0)
    return NULL;
  result_size = compress_components(result_size, result_components);
  if (result_size == 0) {
    PyMem_Free(result_components);
    return NULL;
  }
  if (!PyMem_Resize(result_components, double, result_size))
    return (ExpansionObject *)PyErr_NoMemory();
  return construct_expansion(&ExpansionType, result_size, result_components);
}

static ExpansionObject *expansion_double_add(ExpansionObject *self,
                                             double other) {
  double *result_components;
  size_t result_size;
  if (add_double(self->size, self->components, other, &result_size,
                 &result_components) < 0)
    return NULL;
  result_size = compress_components(result_size, result_components);
  if (!PyMem_Resize(result_components, double, result_size))
    return (ExpansionObject *)PyErr_NoMemory();
  return construct_expansion(&ExpansionType, result_size, result_components);
}

static PyObject *expansion_PyObject_add(ExpansionObject *self,
                                        PyObject *other) {
  if (PyFloat_Check(other))
    return (PyObject *)expansion_double_add((ExpansionObject *)self,
                                            PyFloat_AS_DOUBLE(other));
  else if (PyLong_Check(other)) {
    double *other_components;
    size_t other_size;
    if (py_long_to_components(other, &other_size, &other_components) < 0)
      return NULL;
    double *result_components;
    size_t result_size;
    if (add_components(self->size, self->components, other_size,
                       other_components, &result_size, &result_components) < 0)
      return NULL;
    result_size = compress_components(result_size, result_components);
    if (!PyMem_Resize(result_components, double, result_size))
      return PyErr_NoMemory();
    return (PyObject *)construct_expansion(&ExpansionType, result_size,
                                           result_components);
  } else if (PyObject_IsInstance(other, Rational)) {
    double *other_components;
    size_t other_size;
    if (Rational_to_components(other, &other_size, &other_components) < 0)
      return NULL;
    double *result_components;
    size_t result_size;
    if (add_components(self->size, self->components, other_size,
                       other_components, &result_size, &result_components) < 0)
      return NULL;
    result_size = compress_components(result_size, result_components);
    if (!PyMem_Resize(result_components, double, result_size))
      return PyErr_NoMemory();
    return (PyObject *)construct_expansion(&ExpansionType, result_size,
                                           result_components);
  }
  Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *expansion_add(PyObject *self, PyObject *other) {
  if (PyObject_TypeCheck(self, &ExpansionType))
    return PyObject_TypeCheck(other, &ExpansionType)
               ? (PyObject *)Expansions_add((ExpansionObject *)self,
                                            (ExpansionObject *)other)
               : expansion_PyObject_add((ExpansionObject *)self, other);
  else
    return expansion_PyObject_add((ExpansionObject *)other, self);
}

static int expansion_bool(ExpansionObject *self) {
  return self->components[self->size - 1] != 0.0;
}

static void expansion_dealloc(ExpansionObject *self) {
  PyMem_Free(self->components);
  Py_TYPE(self)->tp_free((PyObject *)self);
}

static double expansion_double(ExpansionObject *self) {
  assert(sum_components(self->size, self->components) ==
         self->components[self->size - 1]);
  return self->components[self->size - 1];
}

static PyObject *expansion_ceil(ExpansionObject *self,
                                PyObject *Py_UNUSED(args)) {
  PyObject *result = components_to_py_long(self->size, self->components);
  if (result == NULL) return NULL;
  double fraction =
      components_to_accumulated_fraction(self->size, self->components);
  assert(fabs(fraction) < 1.0);
  PyObject *fraction_ceil = PyLong_FromLong((long)ceil(fraction));
  if (fraction_ceil == NULL) {
    Py_DECREF(result);
    return NULL;
  }
  PyObject *tmp = result;
  result = PyNumber_InPlaceAdd(result, fraction_ceil);
  Py_DECREF(tmp);
  Py_DECREF(fraction_ceil);
  return result;
}

static PyObject *expansion_float(ExpansionObject *self) {
  return PyFloat_FromDouble(expansion_double(self));
}

static PyObject *expansion_floor(ExpansionObject *self,
                                 PyObject *Py_UNUSED(args)) {
  PyObject *result = components_to_py_long(self->size, self->components);
  if (result == NULL) return NULL;
  double fraction =
      components_to_accumulated_fraction(self->size, self->components);
  assert(fabs(fraction) < 1.0);
  PyObject *fraction_floor = PyLong_FromLong((long)floor(fraction));
  if (fraction_floor == NULL) {
    Py_DECREF(result);
    return NULL;
  }
  PyObject *tmp = result;
  result = PyNumber_InPlaceAdd(result, fraction_floor);
  Py_DECREF(tmp);
  Py_DECREF(fraction_floor);
  return result;
}

static PyObject *expansion_getnewargs(ExpansionObject *self,
                                      PyObject *Py_UNUSED(args)) {
  PyObject *result = PyTuple_New((Py_ssize_t)(self->size));
  if (result == NULL) return NULL;
  for (size_t index = 0; index < self->size; ++index) {
    PyObject *component = PyFloat_FromDouble(self->components[index]);
    if (component == NULL) {
      Py_DECREF(component);
      return NULL;
    }
    PyTuple_SET_ITEM(result, (Py_ssize_t)(index), component);
  }
  return result;
}

static Py_hash_t expansion_hash(ExpansionObject *self) {
  PyObject *components = PyTuple_New((Py_ssize_t)(self->size));
  if (components == NULL) return -1;
  for (size_t index = 0; index < self->size; ++index)
    PyTuple_SET_ITEM(components, (Py_ssize_t)(index),
                     PyFloat_FromDouble(self->components[index]));
  Py_hash_t result = PyObject_Hash(components);
  Py_DECREF(components);
  return result;
}

static ExpansionObject *Expansions_multiply(ExpansionObject *self,
                                            ExpansionObject *other) {
  double *result_components;
  size_t result_size;
  if (multiply_components(self->size, self->components, other->size,
                          other->components, &result_size,
                          &result_components) < 0)
    return NULL;
  result_size = compress_components(result_size, result_components);
  if (result_size == 0) {
    PyMem_Free(result_components);
    return NULL;
  }
  if (!PyMem_Resize(result_components, double, result_size))
    return (ExpansionObject *)PyErr_NoMemory();
  return construct_expansion(&ExpansionType, result_size, result_components);
}

static ExpansionObject *expansion_double_multiply(ExpansionObject *self,
                                                  double other) {
  double *result_components;
  size_t result_size;
  if (scale_components(self->size, self->components, other, &result_size,
                       &result_components) < 0)
    return NULL;
  result_size = compress_components(result_size, result_components);
  if (result_size == 0) {
    PyMem_Free(result_components);
    return NULL;
  }
  if (!PyMem_Resize(result_components, double, result_size))
    return (ExpansionObject *)PyErr_NoMemory();
  return construct_expansion(&ExpansionType, result_size, result_components);
}

static PyObject *expansion_PyObject_multiply(ExpansionObject *self,
                                             PyObject *other) {
  if (PyFloat_Check(other))
    return (PyObject *)expansion_double_multiply(self,
                                                 PyFloat_AS_DOUBLE(other));
  else if (PyLong_Check(other)) {
    double *other_components;
    size_t other_size;
    if (py_long_to_components(other, &other_size, &other_components) < 0)
      return NULL;
    double *result_components;
    size_t result_size;
    if (multiply_components(self->size, self->components, other_size,
                            other_components, &result_size,
                            &result_components) < 0)
      return NULL;
    result_size = compress_components(result_size, result_components);
    if (!PyMem_Resize(result_components, double, result_size))
      return PyErr_NoMemory();
    return (PyObject *)construct_expansion(&ExpansionType, result_size,
                                           result_components);
  } else if (PyObject_IsInstance(other, Rational)) {
    double *other_components;
    size_t other_size;
    if (Rational_to_components(other, &other_size, &other_components) < 0)
      return NULL;
    double *result_components;
    size_t result_size;
    if (multiply_components(self->size, self->components, other_size,
                            other_components, &result_size,
                            &result_components) < 0)
      return NULL;
    result_size = compress_components(result_size, result_components);
    if (!PyMem_Resize(result_components, double, result_size))
      return PyErr_NoMemory();
    return (PyObject *)construct_expansion(&ExpansionType, result_size,
                                           result_components);
  }
  Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *expansion_multiply(PyObject *self, PyObject *other) {
  if (PyObject_TypeCheck(self, &ExpansionType))
    return PyObject_TypeCheck(other, &ExpansionType)
               ? (PyObject *)Expansions_multiply((ExpansionObject *)self,
                                                 (ExpansionObject *)other)

               : expansion_PyObject_multiply((ExpansionObject *)self, other);
  else
    return expansion_PyObject_multiply((ExpansionObject *)other, self);
}

static int is_unit_py_object(PyObject *self) {
  PyObject *tmp = PyLong_FromLong(1);
  int result = PyObject_RichCompareBool(self, tmp, Py_EQ);
  Py_DECREF(tmp);
  return result;
}

static int normalize_fraction_components(PyObject **result_numerator,
                                         PyObject **result_denominator) {
  PyObject *gcd = _PyLong_GCD(*result_numerator, *result_denominator);
  if (gcd == NULL) return -1;
  int is_gcd_unit = is_unit_py_object(gcd);
  if (is_gcd_unit < 0) {
    Py_DECREF(gcd);
    return -1;
  } else if (!is_gcd_unit) {
    PyObject *numerator = PyNumber_FloorDivide(*result_numerator, gcd);
    if (numerator == NULL) {
      Py_DECREF(gcd);
      return -1;
    }
    PyObject *denominator = PyNumber_FloorDivide(*result_denominator, gcd);
    if (denominator == NULL) {
      Py_DECREF(numerator);
      Py_DECREF(gcd);
      return -1;
    }
    PyObject *tmp = *result_numerator;
    *result_numerator = numerator;
    Py_DECREF(tmp);
    tmp = *result_denominator;
    *result_denominator = denominator;
    Py_DECREF(tmp);
  }
  Py_DECREF(gcd);
  return 0;
}

static int fractions_components_add(PyObject *numerator, PyObject *denominator,
                                    PyObject *other_numerator,
                                    PyObject *other_denominator,
                                    PyObject **result_numerator,
                                    PyObject **result_denominator) {
  PyObject *first_result_numerator_component =
      PyNumber_Multiply(numerator, other_denominator);
  if (first_result_numerator_component == NULL) return -1;
  PyObject *second_result_numerator_component =
      PyNumber_Multiply(other_numerator, denominator);
  if (second_result_numerator_component == NULL) {
    Py_DECREF(first_result_numerator_component);
    return -1;
  }
  *result_numerator = PyNumber_Add(first_result_numerator_component,
                                   second_result_numerator_component);
  Py_DECREF(second_result_numerator_component);
  Py_DECREF(first_result_numerator_component);
  if (*result_numerator == NULL) return -1;
  *result_denominator = PyNumber_Multiply(denominator, other_denominator);
  if (*result_denominator == NULL) {
    Py_DECREF(*result_numerator);
    return -1;
  }
  if (normalize_fraction_components(result_numerator, result_denominator)) {
    Py_DECREF(*result_denominator);
    Py_DECREF(*result_numerator);
    return -1;
  }
  return 0;
}

static PyObject *double_as_integer_ratio(double value) {
  PyObject *as_integer_ratio_method_name =
      PyUnicode_FromString("as_integer_ratio");
  if (as_integer_ratio_method_name == NULL) return NULL;
  PyObject *tmp = PyFloat_FromDouble(value);
  PyObject *result =
#if PY39_OR_MORE
      PyObject_CallMethodNoArgs(tmp, as_integer_ratio_method_name);
#else
      PyObject_CallMethodObjArgs(tmp, as_integer_ratio_method_name, NULL)
#endif
  ;
  Py_DECREF(tmp);
  return result;
}

static PyObject *expansion_as_integer_ratio(ExpansionObject *self,
                                            PyObject *Py_UNUSED(args)) {
  PyObject *result = double_as_integer_ratio(self->components[0]);
  if (self->size == 1 || result == NULL) {
    return result;
  }
  PyObject *result_numerator = PyTuple_GET_ITEM(result, 0);
  PyObject *result_denominator = PyTuple_GET_ITEM(result, 1);
  Py_INCREF(result_numerator);
  Py_INCREF(result_denominator);
  Py_DECREF(result);
  for (size_t index = 1; index < self->size; ++index) {
    PyObject *step = double_as_integer_ratio(self->components[index]);
    if (step == NULL) {
      return NULL;
    }
    PyObject *step_numerator = PyTuple_GET_ITEM(step, 0);
    PyObject *step_denominator = PyTuple_GET_ITEM(step, 1);
    PyObject *next_result_numerator, *next_result_denominator;
    if (fractions_components_add(result_numerator, result_denominator,
                                 step_numerator, step_denominator,
                                 &next_result_numerator,
                                 &next_result_denominator) < 0) {
      Py_DECREF(step);
      Py_DECREF(result_denominator);
      Py_DECREF(result_numerator);
      return NULL;
    }
    Py_DECREF(step);
    Py_DECREF(result_denominator);
    Py_DECREF(result_numerator);
    result_numerator = next_result_numerator;
    result_denominator = next_result_denominator;
  }
  return PyTuple_Pack(2, result_numerator, result_denominator);
}

int are_kwargs_passed(PyObject *kwargs) {
  return kwargs != NULL &&
         (!PyDict_CheckExact(kwargs) || PyDict_GET_SIZE(kwargs) != 0);
}

static PyObject *expansion_new(PyTypeObject *cls, PyObject *args,
                               PyObject *kwargs) {
  if (are_kwargs_passed(kwargs)) {
    PyErr_Format(PyExc_TypeError, "Expansion() takes no keyword arguments");
    return NULL;
  }
  double *components;
  Py_ssize_t raw_size = PyTuple_Size(args);
  if (raw_size < 0) return NULL;
  size_t size = (size_t)raw_size;
  if (size == 1) {
    PyObject *argument = PyTuple_GET_ITEM(args, 0);
    if (PyObject_TypeCheck(argument, &ExpansionType)) {
      ExpansionObject *expansion_argument = (ExpansionObject *)argument;
      components =
          (double *)PyMem_Malloc(expansion_argument->size * sizeof(double));
      if (components == NULL) return NULL;
      copy_components(expansion_argument->size, expansion_argument->components,
                      components);
      size = expansion_argument->size;
    } else if (PyFloat_Check(argument)) {
      components = (double *)PyMem_Malloc(sizeof(double));
      if (components == NULL) return PyErr_NoMemory();
      components[0] = PyFloat_AS_DOUBLE(argument);
      size = 1;
    } else if (PyLong_Check(argument)) {
      if (py_long_to_components(argument, &size, &components) < 0) return NULL;
    } else if (PyObject_IsInstance(argument, Rational)) {
      if (Rational_to_components(argument, &size, &components) < 0) return NULL;
    } else {
      PyErr_Format(PyExc_TypeError,
                   "Argument should be of type %R, `Rational` or `float`, but "
                   "found: %R.",
                   &ExpansionType, Py_TYPE(argument));
      return NULL;
    }
  } else if (size != 0) {
    components = (double *)PyMem_Malloc(size * sizeof(double));
    if (components == NULL) return PyErr_NoMemory();
    for (size_t index = 0; index < size; ++index) {
      PyObject *item = PyTuple_GET_ITEM(args, index);
      if (item == NULL) {
        PyMem_Free(components);
        return NULL;
      }
      if (!PyFloat_Check(item)) {
        PyErr_Format(PyExc_TypeError,
                     "Components should be of type `float`, but found: %R.",
                     Py_TYPE(item));
        PyMem_Free(components);
        return NULL;
      }
      components[index] = PyFloat_AS_DOUBLE(item);
    }
    size = compress_components(size, components);
    if (size == 0) return NULL;
    if (!PyMem_Resize(components, double, size)) return PyErr_NoMemory();
  } else {
    components = (double *)PyMem_Malloc(sizeof(double));
    if (components == NULL) return PyErr_NoMemory();
    components[0] = 0.0;
    size = 1;
  }
  return (PyObject *)construct_expansion(cls, size, components);
}

static ExpansionObject *expansion_negative(ExpansionObject *self) {
  double *result_components =
      (double *)PyMem_Malloc(self->size * sizeof(double));
  size_t result_size =
      negate_components(self->size, self->components, result_components);
  return construct_expansion(&ExpansionType, result_size, result_components);
}

static ExpansionObject *expansion_positive(ExpansionObject *self) {
  Py_INCREF(self);
  return self;
}

static ExpansionObject *expansion_absolute(ExpansionObject *self) {
  return self->components[self->size - 1] < 0.0 ? expansion_negative(self)
                                                : expansion_positive(self);
}

static PyObject *expansion_repr(ExpansionObject *self) {
  PyObject *result;
  if (self->size > 1) {
    PyObject *components_reprs = PyTuple_New((Py_ssize_t)(self->size));
    if (components_reprs == NULL) return NULL;
    for (size_t index = 0; index < self->size; ++index) {
      PyObject *item = PyFloat_FromDouble(self->components[index]);
      if (item == NULL) {
        Py_DECREF(components_reprs);
        return NULL;
      }
      PyTuple_SET_ITEM(components_reprs, (Py_ssize_t)(index),
                       PyObject_Repr(item));
      Py_DECREF(item);
    }
    PyObject *separator = PyUnicode_FromString(", ");
    if (separator == NULL) {
      Py_DECREF(components_reprs);
      return NULL;
    }
    PyObject *joined_components_reprs =
        PyUnicode_Join(separator, components_reprs);
    Py_DECREF(separator);
    Py_DECREF(components_reprs);
    if (joined_components_reprs == NULL) return NULL;
    result = PyUnicode_FromFormat("Expansion(%U)", joined_components_reprs);
    Py_DECREF(joined_components_reprs);
  } else {
    PyObject *head = PyFloat_FromDouble(self->components[0]);
    result = PyUnicode_FromFormat("Expansion(%R)", head);
    Py_DECREF(head);
  }
  return result;
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

static PyObject *expansion_py_long_richcompare(ExpansionObject *self,
                                               PyObject *other, int op) {
  switch (op) {
    case Py_EQ:
      return PyBool_FromLong(
          are_components_equal_to_py_long(self->size, self->components, other));
    case Py_GE:
      return PyBool_FromLong(!are_components_lesser_than_py_long(
          self->size, self->components, other));
    case Py_GT:
      return PyBool_FromLong(is_py_long_lesser_than_components(
          other, self->size, self->components));
    case Py_LE:
      return PyBool_FromLong(!is_py_long_lesser_than_components(
          other, self->size, self->components));
    case Py_LT:
      return PyBool_FromLong(are_components_lesser_than_py_long(
          self->size, self->components, other));
    case Py_NE:
      return PyBool_FromLong(!are_components_equal_to_py_long(
          self->size, self->components, other));
    default:
      Py_RETURN_NOTIMPLEMENTED;
  }
}

static PyObject *expansion_Rational_richcompare(ExpansionObject *self,
                                                PyObject *other, int op) {
  switch (op) {
    case Py_EQ:
      return PyBool_FromLong(are_components_equal_to_Rational(
          self->size, self->components, other));
    case Py_GE:
      return PyBool_FromLong(!are_components_lesser_than_Rational(
          self->size, self->components, other));
    case Py_GT:
      return PyBool_FromLong(is_Rational_lesser_than_components(
          other, self->size, self->components));
    case Py_LE:
      return PyBool_FromLong(!is_Rational_lesser_than_components(
          other, self->size, self->components));
    case Py_LT:
      return PyBool_FromLong(are_components_lesser_than_Rational(
          self->size, self->components, other));
    case Py_NE:
      return PyBool_FromLong(!are_components_equal_to_Rational(
          self->size, self->components, other));
    default:
      Py_RETURN_NOTIMPLEMENTED;
  }
}

static PyObject *expansion_double_richcompare(ExpansionObject *self,
                                              double other, int op) {
  switch (op) {
    case Py_EQ:
      return PyBool_FromLong(
          are_components_equal_to_double(self->size, self->components, other));
    case Py_GE:
      return PyBool_FromLong(!are_components_lesser_than_double(
          self->size, self->components, other));
    case Py_GT:
      return PyBool_FromLong(is_double_lesser_than_components(
          other, self->size, self->components));
    case Py_LE:
      return PyBool_FromLong(!is_double_lesser_than_components(
          other, self->size, self->components));
    case Py_LT:
      return PyBool_FromLong(are_components_lesser_than_double(
          self->size, self->components, other));
    case Py_NE:
      return PyBool_FromLong(
          !are_components_equal_to_double(self->size, self->components, other));
    default:
      Py_RETURN_NOTIMPLEMENTED;
  }
}

static PyObject *expansion_richcompare(ExpansionObject *self, PyObject *other,
                                       int op) {
  if (PyObject_TypeCheck(other, &ExpansionType))
    return Expansions_richcompare(self, (ExpansionObject *)other, op);
  else if (PyFloat_Check(other))
    return expansion_double_richcompare(self, PyFloat_AS_DOUBLE(other), op);
  else if (PyLong_Check(other))
    return expansion_py_long_richcompare(self, other, op);
  else if (PyObject_IsInstance(other, Rational))
    return expansion_Rational_richcompare(self, other, op);
  else
    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *expansion_round_plain(ExpansionObject *self) {
  PyObject *result = components_to_py_long(self->size, self->components);
  if (result == NULL) return NULL;
  size_t fractions_size;
  double *fractions_components;
  if (components_to_fractions(self->size, self->components, &fractions_size,
                              &fractions_components) < 0) {
    Py_DECREF(result);
    return NULL;
  }
  if (are_components_equal_to_double(fractions_size, fractions_components,
                                     0.5) ||
      are_components_equal_to_double(fractions_size, fractions_components,
                                     -0.5)) {
    double fraction =
        components_to_accumulated_fraction(self->size, self->components);
    if ((fraction > 0.0 && Py_SIZE(result) < 0) ||
        (fraction < 0.0 && Py_SIZE(result) > 0)) {
      PyObject *sign = PyLong_FromLong(Py_SIZE(result) > 0 ? 1 : -1);
      if (sign == NULL) {
        Py_DECREF(result);
        return NULL;
      }
      PyObject *tmp = result;
      result = PyNumber_InPlaceSubtract(result, sign);
      Py_DECREF(tmp);
      Py_DECREF(sign);
    }
    const int is_truncation_odd = is_py_long_odd(result);
    if (is_truncation_odd < 0) {
      Py_DECREF(result);
      return NULL;
    }
    if (is_truncation_odd) {
      PyObject *sign = PyLong_FromLong(Py_SIZE(result) > 0 ? 1 : -1);
      if (sign == NULL) {
        Py_DECREF(result);
        return NULL;
      }
      PyObject *tmp = result;
      result = PyNumber_InPlaceAdd(result, sign);
      Py_DECREF(tmp);
      Py_DECREF(sign);
    }
  } else if (are_components_lesser_than_double(fractions_size,
                                               fractions_components, -0.5) ||
             is_double_lesser_than_components(0.5, fractions_size,
                                              fractions_components)) {
    PyObject *sign =
        PyLong_FromLong(is_double_lesser_than_components(0.0, fractions_size,
                                                         fractions_components)
                            ? 1
                            : -1);
    if (sign == NULL) {
      Py_DECREF(result);
      return NULL;
    }
    PyObject *tmp = result;
    result = PyNumber_InPlaceAdd(result, sign);
    Py_DECREF(tmp);
    Py_DECREF(sign);
  }
  return result;
}

static PyObject *expansion_round(ExpansionObject *self, PyObject *args) {
  PyObject *precision = NULL;
  if (!PyArg_ParseTuple(args, "|O", &precision)) return NULL;
  if (precision == NULL) return expansion_round_plain(self);
  const size_t size = self->size;
  double *const components = self->components;
  size_t result_size = size;
  double *result_components =
      (double *)PyMem_Malloc(result_size * sizeof(double));
  if (result_components == NULL) {
    Py_DECREF(precision);
    return PyErr_NoMemory();
  }
  for (size_t index = 0; index < self->size; ++index) {
    PyObject *component = PyFloat_FromDouble(components[index]);
    if (component == NULL) {
      PyMem_Free(result_components);
      return NULL;
    }
    PyObject *rounded_component =
#if PY39_OR_MORE
        PyObject_CallMethodOneArg(component, round_method_name, precision)
#else
        PyObject_CallMethodObjArgs(component, round_method_name, precision,
                                   NULL)
#endif
        ;
    if (rounded_component == NULL) {
      PyMem_Free(result_components);
      return NULL;
    }
    result_components[index] = PyFloat_AS_DOUBLE(rounded_component);
    Py_DECREF(rounded_component);
  }
  result_size = compress_components(result_size, result_components);
  if (result_size == 0) {
    PyMem_Free(result_components);
    return NULL;
  }
  if (!PyMem_Resize(result_components, double, result_size))
    return PyErr_NoMemory();
  return (PyObject *)construct_expansion(&ExpansionType, result_size,
                                         result_components);
}

static ExpansionObject *Expansions_subtract(ExpansionObject *self,
                                            ExpansionObject *other) {
  double *result_components;
  size_t result_size;
  if (subtract_components(self->size, self->components, other->size,
                          other->components, &result_size,
                          &result_components) < 0)
    return NULL;
  result_size = compress_components(result_size, result_components);
  if (!PyMem_Resize(result_components, double, result_size))
    return (ExpansionObject *)PyErr_NoMemory();
  return construct_expansion(&ExpansionType, result_size, result_components);
}

static ExpansionObject *expansion_py_long_subtract(ExpansionObject *self,
                                                   PyObject *other) {
  double *other_components;
  size_t other_size;
  if (py_long_to_components(other, &other_size, &other_components) < 0)
    return NULL;
  size_t result_size;
  double *result_components;
  if (subtract_components(self->size, self->components, other_size,
                          other_components, &result_size,
                          &result_components) < 0) {
    PyMem_Free(other_components);
    return NULL;
  }
  PyMem_Free(other_components);
  result_size = compress_components(result_size, result_components);
  if (!PyMem_Resize(result_components, double, result_size))
    return (ExpansionObject *)PyErr_NoMemory();
  return construct_expansion(&ExpansionType, result_size, result_components);
}

static ExpansionObject *py_long_expansion_subtract(PyObject *self,
                                                   ExpansionObject *other) {
  double *components;
  size_t size;
  if (py_long_to_components(self, &size, &components) < 0) return NULL;
  size_t result_size;
  double *result_components;
  if (subtract_components(size, components, other->size, other->components,
                          &result_size, &result_components) < 0) {
    PyMem_Free(components);
    return NULL;
  }
  PyMem_Free(components);
  result_size = compress_components(result_size, result_components);
  if (!PyMem_Resize(result_components, double, result_size))
    return (ExpansionObject *)PyErr_NoMemory();
  return construct_expansion(&ExpansionType, result_size, result_components);
}

static ExpansionObject *expansion_Rational_subtract(ExpansionObject *self,
                                                    PyObject *other) {
  double *other_components;
  size_t other_size;
  if (Rational_to_components(other, &other_size, &other_components) < 0)
    return NULL;
  size_t result_size;
  double *result_components;
  if (subtract_components(self->size, self->components, other_size,
                          other_components, &result_size,
                          &result_components) < 0) {
    PyMem_Free(other_components);
    return NULL;
  }
  PyMem_Free(other_components);
  result_size = compress_components(result_size, result_components);
  if (!PyMem_Resize(result_components, double, result_size))
    return (ExpansionObject *)PyErr_NoMemory();
  return construct_expansion(&ExpansionType, result_size, result_components);
}

static ExpansionObject *Rational_expansion_subtract(PyObject *self,
                                                    ExpansionObject *other) {
  double *components;
  size_t size;
  if (Rational_to_components(self, &size, &components) < 0) return NULL;
  size_t result_size;
  double *result_components;
  if (subtract_components(size, components, other->size, other->components,
                          &result_size, &result_components) < 0) {
    PyMem_Free(components);
    return NULL;
  }
  PyMem_Free(components);
  result_size = compress_components(result_size, result_components);
  if (!PyMem_Resize(result_components, double, result_size))
    return (ExpansionObject *)PyErr_NoMemory();
  return construct_expansion(&ExpansionType, result_size, result_components);
}

static ExpansionObject *expansion_double_subtract(ExpansionObject *self,
                                                  double other) {
  double *result_components;
  size_t result_size = 0;
  if (subtract_double(self->size, self->components, other, &result_size,
                      &result_components) < 0)
    return NULL;
  result_size = compress_components(result_size, result_components);
  if (!PyMem_Resize(result_components, double, result_size))
    return (ExpansionObject *)PyErr_NoMemory();
  return construct_expansion(&ExpansionType, result_size, result_components);
}

static ExpansionObject *double_expansion_subtract(double self,
                                                  ExpansionObject *other) {
  double *result_components;
  size_t result_size;
  if (subtract_from_double(self, other->size, other->components, &result_size,
                           &result_components) < 0)
    return NULL;
  result_size = compress_components(result_size, result_components);
  if (!PyMem_Resize(result_components, double, result_size))
    return (ExpansionObject *)PyErr_NoMemory();
  return construct_expansion(&ExpansionType, result_size, result_components);
}

static PyObject *expansion_subtract(PyObject *self, PyObject *other) {
  if (PyObject_TypeCheck(self, &ExpansionType)) {
    if (PyObject_TypeCheck(other, &ExpansionType))
      return (PyObject *)Expansions_subtract((ExpansionObject *)self,
                                             (ExpansionObject *)other);
    else if (PyFloat_Check(other))
      return (PyObject *)expansion_double_subtract((ExpansionObject *)self,
                                                   PyFloat_AS_DOUBLE(other));
    else if (PyLong_Check(other))
      return (PyObject *)expansion_py_long_subtract((ExpansionObject *)self,
                                                    other);
    else if (PyObject_IsInstance(other, Rational))
      return (PyObject *)expansion_Rational_subtract((ExpansionObject *)self,
                                                     other);
  } else if (PyFloat_Check(self))
    return (PyObject *)double_expansion_subtract(PyFloat_AS_DOUBLE(self),
                                                 (ExpansionObject *)other);
  else if (PyLong_Check(self))
    return (PyObject *)py_long_expansion_subtract(self,
                                                  (ExpansionObject *)other);
  else if (PyObject_IsInstance(self, Rational))
    return (PyObject *)Rational_expansion_subtract(self,
                                                   (ExpansionObject *)other);
  Py_RETURN_NOTIMPLEMENTED;
}

static ExpansionObject *Expansions_true_divide(ExpansionObject *self,
                                               ExpansionObject *other) {
  if (!expansion_bool(other)) {
    PyErr_Format(PyExc_ZeroDivisionError, "Divisor is zero.");
    return NULL;
  }
  double *result_components;
  size_t result_size;
  if (divide_components(self->size, self->components, other->size,
                        other->components, &result_size,
                        &result_components) < 0)
    return NULL;
  return construct_expansion(&ExpansionType, result_size, result_components);
}

static ExpansionObject *expansion_py_long_true_divide(ExpansionObject *self,
                                                      PyObject *other) {
  if (PyObject_Not(other)) {
    PyErr_Format(PyExc_ZeroDivisionError, "Divisor is zero.");
    return NULL;
  }
  double *other_components;
  size_t other_size;
  if (py_long_to_components(other, &other_size, &other_components) < 0)
    return NULL;
  size_t result_size;
  double *result_components;
  if (divide_components(self->size, self->components, other_size,
                        other_components, &result_size,
                        &result_components) < 0) {
    PyMem_Free(other_components);
    return NULL;
  }
  PyMem_Free(other_components);
  return construct_expansion(&ExpansionType, result_size, result_components);
}

static ExpansionObject *expansion_Rational_true_divide(ExpansionObject *self,
                                                       PyObject *other) {
  if (PyObject_Not(other)) {
    PyErr_Format(PyExc_ZeroDivisionError, "Divisor is zero.");
    return NULL;
  }
  double *other_components;
  size_t other_size;
  if (Rational_to_components(other, &other_size, &other_components) < 0)
    return NULL;
  size_t result_size;
  double *result_components;
  if (divide_components(self->size, self->components, other_size,
                        other_components, &result_size,
                        &result_components) < 0) {
    PyMem_Free(other_components);
    return NULL;
  }
  PyMem_Free(other_components);
  return construct_expansion(&ExpansionType, result_size, result_components);
}

static ExpansionObject *py_long_expansion_true_divide(PyObject *self,
                                                      ExpansionObject *other) {
  if (!expansion_bool(other)) {
    PyErr_Format(PyExc_ZeroDivisionError, "Divisor is zero.");
    return NULL;
  }
  double *components;
  size_t size;
  if (py_long_to_components(self, &size, &components) < 0) return NULL;
  size_t result_size;
  double *result_components;
  if (divide_components(size, components, other->size, other->components,
                        &result_size, &result_components) < 0) {
    PyMem_Free(components);
    return NULL;
  }
  PyMem_Free(components);
  return construct_expansion(&ExpansionType, result_size, result_components);
}

static ExpansionObject *Rational_expansion_true_divide(PyObject *self,
                                                       ExpansionObject *other) {
  if (!expansion_bool(other)) {
    PyErr_Format(PyExc_ZeroDivisionError, "Divisor is zero.");
    return NULL;
  }
  double *components;
  size_t size;
  if (Rational_to_components(self, &size, &components) < 0) return NULL;
  size_t result_size;
  double *result_components;
  if (divide_components(size, components, other->size, other->components,
                        &result_size, &result_components) < 0) {
    PyMem_Free(components);
    return NULL;
  }
  PyMem_Free(components);
  return construct_expansion(&ExpansionType, result_size, result_components);
}

static ExpansionObject *expansion_double_true_divide(ExpansionObject *self,
                                                     double other) {
  if (other == 0.0) {
    PyErr_Format(PyExc_ZeroDivisionError, "Divisor is zero.");
    return NULL;
  }
  double *other_components = (double *)PyMem_Malloc(sizeof(double));
  if (other_components == NULL) return (ExpansionObject *)PyErr_NoMemory();
  other_components[0] = other;
  double *result_components;
  size_t result_size;
  if (divide_components(self->size, self->components, 1, other_components,
                        &result_size, &result_components) < 0) {
    PyMem_Free(other_components);
    return NULL;
  }
  PyMem_Free(other_components);
  return construct_expansion(&ExpansionType, result_size, result_components);
}

static ExpansionObject *double_expansion_true_divide(double self,
                                                     ExpansionObject *other) {
  if (!expansion_bool(other)) {
    PyErr_Format(PyExc_ZeroDivisionError, "Divisor is zero.");
    return NULL;
  }
  double *components = (double *)PyMem_Malloc(sizeof(double));
  if (components == NULL) return (ExpansionObject *)PyErr_NoMemory();
  components[0] = self;
  double *result_components;
  size_t result_size;
  if (divide_components(1, components, other->size, other->components,
                        &result_size, &result_components) < 0) {
    PyMem_Free(components);
    return NULL;
  }
  PyMem_Free(components);
  return construct_expansion(&ExpansionType, result_size, result_components);
}

static PyObject *expansion_true_divide(PyObject *self, PyObject *other) {
  if (PyObject_TypeCheck(self, &ExpansionType)) {
    if (PyObject_TypeCheck(other, &ExpansionType))
      return (PyObject *)Expansions_true_divide((ExpansionObject *)self,
                                                (ExpansionObject *)other);
    else if (PyFloat_Check(other))
      return (PyObject *)expansion_double_true_divide((ExpansionObject *)self,
                                                      PyFloat_AS_DOUBLE(other));
    else if (PyLong_Check(other))
      return (PyObject *)expansion_py_long_true_divide((ExpansionObject *)self,
                                                       other);
    else if (PyObject_IsInstance(other, Rational))
      return (PyObject *)expansion_Rational_true_divide((ExpansionObject *)self,
                                                        other);
  } else if (PyFloat_Check(self))
    return (PyObject *)double_expansion_true_divide(PyFloat_AS_DOUBLE(self),
                                                    (ExpansionObject *)other);
  else if (PyLong_Check(self))
    return (PyObject *)py_long_expansion_true_divide(self,
                                                     (ExpansionObject *)other);
  else if (PyObject_IsInstance(self, Rational))
    return (PyObject *)Rational_expansion_true_divide(self,
                                                      (ExpansionObject *)other);
  Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *expansion_trunc(ExpansionObject *self,
                                 PyObject *Py_UNUSED(args)) {
  PyObject *result = components_to_py_long(self->size, self->components);
  if (result == NULL) return NULL;
  double fraction =
      components_to_accumulated_fraction(self->size, self->components);
  if (fraction < 0.0) {
    int is_whole_part_positive = PyObject_RichCompareBool(result, zero, Py_GT);
    if (is_whole_part_positive < 0) {
      Py_DECREF(result);
      return NULL;
    } else if (is_whole_part_positive > 0) {
      PyObject *sign = PyLong_FromLong(1);
      result = PyNumber_InPlaceSubtract(result, sign);
      Py_DECREF(sign);
    }
  } else if (fraction > 0.0) {
    int is_whole_part_negative = PyObject_RichCompareBool(result, zero, Py_LT);
    if (is_whole_part_negative < 0) {
      Py_DECREF(result);
      return NULL;
    } else if (is_whole_part_negative > 0) {
      PyObject *sign = PyLong_FromLong(-1);
      result = PyNumber_InPlaceSubtract(result, sign);
      Py_DECREF(sign);
    }
  }
  return result;
}

static PyNumberMethods expansion_as_number = {
    .nb_absolute = (unaryfunc)expansion_absolute,
    .nb_add = expansion_add,
    .nb_bool = (inquiry)expansion_bool,
    .nb_float = (unaryfunc)expansion_float,
    .nb_multiply = expansion_multiply,
    .nb_negative = (unaryfunc)expansion_negative,
    .nb_positive = (unaryfunc)expansion_positive,
    .nb_subtract = expansion_subtract,
    .nb_true_divide = expansion_true_divide,
};

PyObject *expansion_getimag(ExpansionObject *Py_UNUSED(self),
                            void *Py_UNUSED(closure)) {
  return PyLong_FromLong(0);
}

PyObject *expansion_getreal(ExpansionObject *self, void *Py_UNUSED(closure)) {
  return (PyObject *)expansion_positive(self);
}

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#endif
static PyGetSetDef expansion_getset[] = {
    {"real", (getter)expansion_getreal, (setter)NULL,
     "The real part of the expansion.", NULL},
    {"imag", (getter)expansion_getimag, (setter)NULL,
     "The imaginary part of the expansion.", NULL},
    {NULL} /* sentinel */
};
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#endif
static PyMethodDef expansion_methods[] = {
    {"as_integer_ratio", (PyCFunction)expansion_as_integer_ratio, METH_NOARGS,
     NULL},
    {"__ceil__", (PyCFunction)expansion_ceil, METH_NOARGS, NULL},
    {"__floor__", (PyCFunction)expansion_floor, METH_NOARGS, NULL},
    {"__getnewargs__", (PyCFunction)expansion_getnewargs, METH_NOARGS, NULL},
    {"__round__", (PyCFunction)expansion_round, METH_VARARGS, NULL},
    {"__trunc__", (PyCFunction)expansion_trunc, METH_NOARGS, NULL},
    {NULL, NULL} /* sentinel */
};
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

static PyTypeObject ExpansionType = {
    PyVarObject_HEAD_INIT(NULL, 0).tp_as_number = &expansion_as_number,
    .tp_basicsize = sizeof(ExpansionObject),
    .tp_dealloc = (destructor)expansion_dealloc,
    .tp_doc = PyDoc_STR("Represents floating point number expansion."),
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_getset = expansion_getset,
    .tp_hash = (hashfunc)expansion_hash,
    .tp_itemsize = 0,
    .tp_methods = expansion_methods,
    .tp_name = "shewchuk.Expansion",
    .tp_new = expansion_new,
    .tp_repr = (reprfunc)expansion_repr,
    .tp_richcompare = (richcmpfunc)expansion_richcompare,
};

static PyObject *incircle_test(PyObject *Py_UNUSED(self), PyObject *args) {
  double point_x, point_y, first_x, first_y, second_x, second_y, third_x,
      third_y;
  if (!PyArg_ParseTuple(args, "dddddddd", &point_x, &point_y, &first_x,
                        &first_y, &second_x, &second_y, &third_x, &third_y))
    return NULL;
  double estimation = incircle_determinant_estimation(
      point_x, point_y, first_x, first_y, second_x, second_y, third_x, third_y);
  return PyLong_FromLong(to_sign(estimation));
}

static PyObject *kind(PyObject *Py_UNUSED(self), PyObject *args) {
  double vertex_x, vertex_y, first_ray_point_x, first_ray_point_y,
      second_ray_point_x, second_ray_point_y;
  if (!PyArg_ParseTuple(args, "dddddd", &vertex_x, &vertex_y,
                        &first_ray_point_x, &first_ray_point_y,
                        &second_ray_point_x, &second_ray_point_y))
    return NULL;
  return PyLong_FromLong(to_sign(vectors_cross_product_estimation(
      vertex_x, vertex_y, first_ray_point_x, first_ray_point_y, -vertex_y,
      vertex_x, -second_ray_point_y, second_ray_point_x)));
}

static PyObject *orientation(PyObject *Py_UNUSED(self), PyObject *args) {
  double start_x, start_y, end_x, end_y, point_x, point_y;
  if (!PyArg_ParseTuple(args, "dddddd", &start_x, &start_y, &end_x, &end_y,
                        &point_x, &point_y))
    return NULL;
  return PyLong_FromLong(to_sign(vectors_cross_product_estimation(
      start_x, start_y, end_x, end_y, start_x, start_y, point_x, point_y)));
}

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
  double *result_components =
      (double *)PyMem_Malloc(result_size * sizeof(double));
  if (result_components == NULL) return PyErr_NoMemory();
  copy_components(result_size, components, result_components);
  return (PyObject *)construct_expansion(&ExpansionType, result_size,
                                         result_components);
}

static PyObject *vectors_dot_product(PyObject *Py_UNUSED(self),
                                     PyObject *args) {
  double first_start_x, first_start_y, first_end_x, first_end_y, second_start_x,
      second_start_y, second_end_x, second_end_y;
  if (!PyArg_ParseTuple(args, "dddddddd", &first_start_x, &first_start_y,
                        &first_end_x, &first_end_y, &second_start_x,
                        &second_start_y, &second_end_x, &second_end_y))
    return NULL;
  double components[16];
  size_t result_size = vectors_cross_product_impl(
      first_start_x, first_start_y, first_end_x, first_end_y, -second_start_y,
      second_start_x, -second_end_y, second_end_x, components);
  double *result_components =
      (double *)PyMem_Malloc(result_size * sizeof(double));
  if (result_components == NULL) return PyErr_NoMemory();
  copy_components(result_size, components, result_components);
  return (PyObject *)construct_expansion(&ExpansionType, result_size,
                                         result_components);
}

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#endif
static PyMethodDef _cshewchuk_methods[] = {
    {"incircle_test", incircle_test, METH_VARARGS,
     PyDoc_STR("incircle_test(point_x, point_y, first_x, first_y, second_x, "
               "second_y, third_x, third_y, /)\n--\n\n"
               "Computes location  of point relative to a circle formed by "
               "three others given their coordinates.")},
    {"kind", kind, METH_VARARGS,
     PyDoc_STR("kind(vertex_x, vertex_y, first_ray_point_x, first_ray_point_y, "
               "second_ray_point_x, second_ray_point_y, /)\n--\n\n"
               "Computes kind of angle given its endpoints coordinates.")},
    {"orientation", orientation, METH_VARARGS,
     PyDoc_STR("orientation(start_x, start_y, end_x, end_y, point_x, point_y, "
               "/)\n--\n\n"
               "Computes orientation of point relative to segment given their "
               "coordinates.")},
    {"vectors_cross_product", vectors_cross_product, METH_VARARGS,
     PyDoc_STR("vectors_cross_product(first_start_x, first_start_y, "
               "first_end_x, first_end_y, second_start_x, second_start_y, "
               "second_end_x, second_end_y, /)\n--\n\n"
               "Computes cross product of two vectors given their endpoints "
               "coordinates.")},
    {"vectors_dot_product", vectors_dot_product, METH_VARARGS,
     PyDoc_STR("vectors_dot_product(first_start_x, first_start_y, first_end_x, "
               "first_end_y, second_start_x, second_start_y, second_end_x, "
               "second_end_y, /)\n--\n\n"
               "Computes dot product of two vectors given their endpoints "
               "coordinates.")},
    {NULL, NULL}, /* sentinel */
};
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

static PyModuleDef _cshewchuk_module = {
    PyModuleDef_HEAD_INIT,
    .m_doc = PyDoc_STR("Robust floating point operations."),
    .m_methods = _cshewchuk_methods,
    .m_name = "shewchuk",
    .m_size = -1,
};

static int load_number_interfaces() {
  PyObject *numbers_module = PyImport_ImportModule("numbers");
  if (numbers_module == NULL) return -1;
  Rational = PyObject_GetAttrString(numbers_module, "Rational");
  if (Rational == NULL) {
    Py_DECREF(numbers_module);
    return -1;
  }
  Real = PyObject_GetAttrString(numbers_module, "Real");
  Py_DECREF(numbers_module);
  if (Real == NULL) {
    Py_DECREF(Rational);
    return -1;
  }
  return 0;
}

static int mark_as_real(PyObject *python_type) {
  PyObject *register_method_name = PyUnicode_FromString("register");
  if (register_method_name == NULL) return -1;
  PyObject *tmp =
#if PY39_OR_MORE
      PyObject_CallMethodOneArg(Real, register_method_name, python_type);
#else
      PyObject_CallMethodObjArgs(Real, register_method_name, python_type, NULL)
#endif
  ;
  Py_DECREF(register_method_name);
  if (tmp == NULL) return -1;
  Py_DECREF(tmp);
  return 0;
}

PyMODINIT_FUNC PyInit__cshewchuk(void) {
  PyObject *result;
  if (PyType_Ready(&ExpansionType) < 0) return NULL;
  result = PyModule_Create(&_cshewchuk_module);
  if (result == NULL) return NULL;
  Py_INCREF(&ExpansionType);
  if (PyModule_AddObject(result, "Expansion", (PyObject *)&ExpansionType) < 0) {
    Py_DECREF(&ExpansionType);
    Py_DECREF(result);
    return NULL;
  }
  round_method_name = PyUnicode_InternFromString("__round__");
  if (round_method_name == NULL) {
    Py_DECREF(result);
    return NULL;
  }
  zero = PyLong_FromLong(0);
  if (zero == NULL) {
    Py_DECREF(round_method_name);
    Py_DECREF(result);
    return NULL;
  }
  if (load_number_interfaces() < 0) {
    Py_DECREF(zero);
    Py_DECREF(round_method_name);
    Py_DECREF(result);
    return NULL;
  }
  if (mark_as_real((PyObject *)&ExpansionType) < 0) {
    Py_DECREF(zero);
    Py_DECREF(round_method_name);
    Py_DECREF(Rational);
    Py_DECREF(Real);
    Py_DECREF(result);
    return NULL;
  }
  return result;
}
