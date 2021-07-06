#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <float.h>
#include <pymath.h>
#include <stddef.h>
#include <stdio.h>
#include <structmember.h>
#define EPSILON (DBL_EPSILON / 2.0)
#define PY39_OR_MORE PY_VERSION_HEX >= 0x03090000

static int to_sign(double value) { return value > 0.0 ? 1 : (!value ? 0 : -1); }

static double double_remainder(double dividend, double divisor) {
  double result = fmod(dividend, divisor);
  // ensure the remainder has the same sign as the denominator
  if (result) {
    if ((divisor < 0) != (result < 0)) result += divisor;
  } else
    result = copysign(0.0, divisor);
  return result;
}

static void swap(double **first, double **second) {
  double *temp = *first;
  *first = *second;
  *second = temp;
}

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
  double *next_components = PyMem_Calloc(size, sizeof(double));
  copy_components(components, size, next_components);
  for (size_t step = 0; step < original_size; ++step) {
    const size_t next_size = compress_components_single(size, next_components);
    if (are_components_equal(next_size, next_components, size, components))
      break;
    size = next_size;
    copy_components(next_components, size, components);
  }
  PyMem_Free(next_components);
  return size;
}

double sum_components(size_t size, double *components) {
  double result = components[0];
  for (size_t index = 1; index < size; ++index) result += components[index];
  return result;
}

static size_t add_double_eliminating_zeros(size_t size, double *components,
                                           double value, double *result) {
  size_t result_size = 0;
  double accumulator = value;
  for (size_t index = 0; index < size; index++) {
    double tail;
    two_add(accumulator, components[index], &accumulator, &tail);
    if (!!tail) result[result_size++] = tail;
  }
  if (!!accumulator || !result_size) result[result_size++] = accumulator;
  return result_size;
}

static size_t scale_components_eliminating_zeros(size_t size,
                                                 double *components,
                                                 double scalar,
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

static size_t multiply_components_eliminating_zeros(size_t left_size,
                                                    double *left,
                                                    size_t right_size,
                                                    double *right, double *step,
                                                    double *result) {
  size_t result_size =
      scale_components_eliminating_zeros(left_size, left, right[0], result);
  for (size_t index = 1; index < right_size; ++index) {
    size_t step_size =
        scale_components_eliminating_zeros(left_size, left, right[index], step);
    result_size = add_components_eliminating_zeros(result_size, result,
                                                   step_size, step, result);
  }
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
      scale_components_eliminating_zeros(size, components, dx, dx_components);
  size_t dx_squared_components_size = scale_components_eliminating_zeros(
      dx_components_size, dx_components, dx, dx_squared_components);
  size_t dy_components_size =
      scale_components_eliminating_zeros(size, components, dy, dy_components);
  size_t dy_squared_components_size = scale_components_eliminating_zeros(
      dy_components_size, dy_components, dy, dy_squared_components);
  return add_components_eliminating_zeros(
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
  size_t first_dx_tail_second_third_cross_product_size;
  double first_dx_tail_second_third_cross_product[8];
  if (!!first_dx_tail) {
    first_dx_tail_second_third_cross_product_size =
        scale_components_eliminating_zeros(
            4, second_third_cross_product, first_dx_tail,
            first_dx_tail_second_third_cross_product);
    first_buffer_16_limit = scale_components_eliminating_zeros(
        first_dx_tail_second_third_cross_product_size,
        first_dx_tail_second_third_cross_product, 2.0 * first_dx,
        first_buffer_16);
    double first_dx_tail_third_squared_length[8];
    size_t first_dx_tail_third_squared_length_size =
        scale_components_eliminating_zeros(4, third_squared_length,
                                           first_dx_tail,
                                           first_dx_tail_third_squared_length);
    second_buffer_16_limit = scale_components_eliminating_zeros(
        first_dx_tail_third_squared_length_size,
        first_dx_tail_third_squared_length, second_dy, second_buffer_16);
    double first_dx_tail_second_squared_length[8];
    size_t first_dx_tail_second_squared_length_size =
        scale_components_eliminating_zeros(4, second_squared_length,
                                           first_dx_tail,
                                           first_dx_tail_second_squared_length);
    third_buffer_16_limit = scale_components_eliminating_zeros(
        first_dx_tail_second_squared_length_size,
        first_dx_tail_second_squared_length, -third_dy, third_buffer_16);
    first_buffer_32_limit = add_components_eliminating_zeros(
        first_buffer_16_limit, first_buffer_16, second_buffer_16_limit,
        second_buffer_16, first_buffer_32);
    buffer_48_limit = add_components_eliminating_zeros(
        third_buffer_16_limit, third_buffer_16, first_buffer_32_limit,
        first_buffer_32, buffer_48);
    final_size = add_components_eliminating_zeros(final_size, *final_components,
                                                  buffer_48_limit, buffer_48,
                                                  *accumulated_components);
    swap(final_components, accumulated_components);
  }
  size_t first_dy_tail_second_third_cross_product_size;
  double first_dy_tail_second_third_cross_product[8];
  if (!!first_dy_tail) {
    first_dy_tail_second_third_cross_product_size =
        scale_components_eliminating_zeros(
            4, second_third_cross_product, first_dy_tail,
            first_dy_tail_second_third_cross_product);
    first_buffer_16_limit = scale_components_eliminating_zeros(
        first_dy_tail_second_third_cross_product_size,
        first_dy_tail_second_third_cross_product, 2.0 * first_dy,
        first_buffer_16);
    double first_dy_tail_second_squared_length[8];
    size_t first_dy_tail_second_squared_length_size =
        scale_components_eliminating_zeros(4, second_squared_length,
                                           first_dy_tail,
                                           first_dy_tail_second_squared_length);
    second_buffer_16_limit = scale_components_eliminating_zeros(
        first_dy_tail_second_squared_length_size,
        first_dy_tail_second_squared_length, third_dx, second_buffer_16);
    double first_dy_tail_third_squared_length[8];
    size_t first_dy_tail_third_squared_length_size =
        scale_components_eliminating_zeros(4, third_squared_length,
                                           first_dy_tail,
                                           first_dy_tail_third_squared_length);
    third_buffer_16_limit = scale_components_eliminating_zeros(
        first_dy_tail_third_squared_length_size,
        first_dy_tail_third_squared_length, -second_dx, third_buffer_16);
    first_buffer_32_limit = add_components_eliminating_zeros(
        first_buffer_16_limit, first_buffer_16, second_buffer_16_limit,
        second_buffer_16, first_buffer_32);
    buffer_48_limit = add_components_eliminating_zeros(
        third_buffer_16_limit, third_buffer_16, first_buffer_32_limit,
        first_buffer_32, buffer_48);
    final_size = add_components_eliminating_zeros(final_size, *final_components,
                                                  buffer_48_limit, buffer_48,
                                                  *accumulated_components);
    swap(final_components, accumulated_components);
  }
  double dx_tail_dy_head_head, dx_head_dy_tail_head;
  double dx_tail_dy_head_tail, dx_head_dy_tail_tail;
  double buffer_8[8], buffer_64[64];
  size_t buffer_8_limit, buffer_64_limit;
  double first_buffer_4[4], second_buffer_4[4];
  if (!!first_dx_tail || !!first_dy_tail) {
    size_t second_third_cross_product_bodies_size,
        second_third_cross_product_tails_size;
    double second_third_cross_product_bodies[8],
        second_third_cross_product_tails[4];
    if (!!second_dx_tail || !!second_dy_tail || !!third_dx_tail ||
        !!third_dy_tail) {
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
      second_third_cross_product_bodies_size = add_components_eliminating_zeros(
          4, first_buffer_4, 4, second_buffer_4,
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
    if (!!first_dx_tail) {
      first_buffer_16_limit = scale_components_eliminating_zeros(
          first_dx_tail_second_third_cross_product_size,
          first_dx_tail_second_third_cross_product, first_dx_tail,
          first_buffer_16);
      double first_dx_tail_second_third_cross_product_bodies[16];
      size_t first_dx_tail_second_third_cross_product_bodies_size =
          scale_components_eliminating_zeros(
              second_third_cross_product_bodies_size,
              second_third_cross_product_bodies, first_dx_tail,
              first_dx_tail_second_third_cross_product_bodies);
      first_buffer_32_limit = scale_components_eliminating_zeros(
          first_dx_tail_second_third_cross_product_bodies_size,
          first_dx_tail_second_third_cross_product_bodies, 2.0 * first_dx,
          first_buffer_32);
      buffer_48_limit = add_components_eliminating_zeros(
          first_buffer_16_limit, first_buffer_16, first_buffer_32_limit,
          first_buffer_32, buffer_48);
      final_size = add_components_eliminating_zeros(
          final_size, *final_components, buffer_48_limit, buffer_48,
          *accumulated_components);
      swap(final_components, accumulated_components);
      if (!!second_dy_tail) {
        buffer_8_limit = scale_components_eliminating_zeros(
            4, third_squared_length, first_dx_tail, buffer_8);
        first_buffer_16_limit = scale_components_eliminating_zeros(
            buffer_8_limit, buffer_8, second_dy_tail, first_buffer_16);
        final_size = add_components_eliminating_zeros(
            final_size, *final_components, first_buffer_16_limit,
            first_buffer_16, *accumulated_components);
        swap(final_components, accumulated_components);
      }
      if (!!third_dy_tail) {
        buffer_8_limit = scale_components_eliminating_zeros(
            4, second_squared_length, -first_dx_tail, buffer_8);
        first_buffer_16_limit = scale_components_eliminating_zeros(
            buffer_8_limit, buffer_8, third_dy_tail, first_buffer_16);
        final_size = add_components_eliminating_zeros(
            final_size, *final_components, first_buffer_16_limit,
            first_buffer_16, *accumulated_components);
        swap(final_components, accumulated_components);
      }
      first_buffer_32_limit = scale_components_eliminating_zeros(
          first_dx_tail_second_third_cross_product_bodies_size,
          first_dx_tail_second_third_cross_product_bodies, first_dx_tail,
          first_buffer_32);
      double first_dx_tail_second_third_cross_product_tails[8];
      size_t first_dx_tail_second_third_cross_product_tails_size =
          scale_components_eliminating_zeros(
              second_third_cross_product_tails_size,
              second_third_cross_product_tails, first_dx_tail,
              first_dx_tail_second_third_cross_product_tails);
      first_buffer_16_limit = scale_components_eliminating_zeros(
          first_dx_tail_second_third_cross_product_tails_size,
          first_dx_tail_second_third_cross_product_tails, 2.0 * first_dx,
          first_buffer_16);
      second_buffer_16_limit = scale_components_eliminating_zeros(
          first_dx_tail_second_third_cross_product_tails_size,
          first_dx_tail_second_third_cross_product_tails, first_dx_tail,
          second_buffer_16);
      second_buffer_32_limit = add_components_eliminating_zeros(
          first_buffer_16_limit, first_buffer_16, second_buffer_16_limit,
          second_buffer_16, second_buffer_32);
      buffer_64_limit = add_components_eliminating_zeros(
          first_buffer_32_limit, first_buffer_32, second_buffer_32_limit,
          second_buffer_32, buffer_64);
      final_size = add_components_eliminating_zeros(
          final_size, *final_components, buffer_64_limit, buffer_64,
          *accumulated_components);
      swap(final_components, accumulated_components);
    }
    if (!!first_dy_tail) {
      first_buffer_16_limit = scale_components_eliminating_zeros(
          first_dy_tail_second_third_cross_product_size,
          first_dy_tail_second_third_cross_product, first_dy_tail,
          first_buffer_16);
      double first_dy_tail_second_third_cross_product_bodies[16];
      size_t first_dy_tail_second_third_cross_product_bodies_size =
          scale_components_eliminating_zeros(
              second_third_cross_product_bodies_size,
              second_third_cross_product_bodies, first_dy_tail,
              first_dy_tail_second_third_cross_product_bodies);
      first_buffer_32_limit = scale_components_eliminating_zeros(
          first_dy_tail_second_third_cross_product_bodies_size,
          first_dy_tail_second_third_cross_product_bodies, 2.0 * first_dy,
          first_buffer_32);
      buffer_48_limit = add_components_eliminating_zeros(
          first_buffer_16_limit, first_buffer_16, first_buffer_32_limit,
          first_buffer_32, buffer_48);
      final_size = add_components_eliminating_zeros(
          final_size, *final_components, buffer_48_limit, buffer_48,
          *accumulated_components);
      swap(final_components, accumulated_components);
      first_buffer_32_limit = scale_components_eliminating_zeros(
          first_dy_tail_second_third_cross_product_bodies_size,
          first_dy_tail_second_third_cross_product_bodies, first_dy_tail,
          first_buffer_32);
      double first_dy_tail_second_third_cross_product_tails[8];
      size_t first_dy_tail_second_third_cross_product_tails_size =
          scale_components_eliminating_zeros(
              second_third_cross_product_tails_size,
              second_third_cross_product_tails, first_dy_tail,
              first_dy_tail_second_third_cross_product_tails);
      first_buffer_16_limit = scale_components_eliminating_zeros(
          first_dy_tail_second_third_cross_product_tails_size,
          first_dy_tail_second_third_cross_product_tails, 2.0 * first_dy,
          first_buffer_16);
      second_buffer_16_limit = scale_components_eliminating_zeros(
          first_dy_tail_second_third_cross_product_tails_size,
          first_dy_tail_second_third_cross_product_tails, first_dy_tail,
          second_buffer_16);
      second_buffer_32_limit = add_components_eliminating_zeros(
          first_buffer_16_limit, first_buffer_16, second_buffer_16_limit,
          second_buffer_16, second_buffer_32);
      buffer_64_limit = add_components_eliminating_zeros(
          first_buffer_32_limit, first_buffer_32, second_buffer_32_limit,
          second_buffer_32, buffer_64);
      final_size = add_components_eliminating_zeros(
          final_size, *final_components, buffer_64_limit, buffer_64,
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
  size_t first_second_sum_size = add_components_eliminating_zeros(
      first_components_size, first_components, second_components_size,
      second_components, first_second_sum_components);
  double first_buffer[1152];
  size_t final_size = add_components_eliminating_zeros(
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
  if (!first_dx_tail && !second_dx_tail && !third_dx_tail && !first_dy_tail &&
      !second_dy_tail && !third_dy_tail)
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
  double first_squared_length[4];
  if (!!second_dx_tail || !!second_dy_tail || !!third_dx_tail ||
      !!third_dy_tail)
    squared_length(first_dx, first_dy, &first_squared_length[3],
                   &first_squared_length[2], &first_squared_length[1],
                   &first_squared_length[0]);
  double second_squared_length[4];
  if (!!third_dx_tail || !!third_dy_tail || !!first_dx_tail || !!first_dy_tail)
    squared_length(second_dx, second_dy, &second_squared_length[3],
                   &second_squared_length[2], &second_squared_length[1],
                   &second_squared_length[0]);
  double third_squared_length[4];
  if (!!first_dx_tail || !!first_dy_tail || !!second_dx_tail ||
      !!second_dy_tail)
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
  if (!minuend_x_tail && !minuend_y_tail && !subtrahend_x_tail &&
      !subtrahend_y_tail)
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
  double final_components[16];
  size_t final_components_size =
      add_components_eliminating_zeros(third_components_size, third_components,
                                       4, extra_components, final_components);
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
    size_t offset = 0, result_size = 4;
    for (; offset < result_size - 1 && !first_components[offset]; ++offset)
      ;
    result_size -= offset;
    copy_components(&first_components[offset], result_size, result);
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
  if (!minuend_x_tail && !minuend_y_tail && !subtrahend_x_tail &&
      !subtrahend_y_tail) {
    size_t offset = 0, result_size = 4;
    for (; offset < result_size - 1 && !first_components[offset]; ++offset)
      ;
    result_size -= offset;
    copy_components(&first_components[offset], result_size, result);
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
    return add_double_eliminating_zeros(4, first_components, extra, result);
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

static PyObject *PyObject_round = NULL;
static PyObject *Real = NULL;

typedef struct {
  PyObject_HEAD size_t size;
  double *components;
} ExpansionObject;

static ExpansionObject *construct_Expansion(PyTypeObject *cls,
                                            double *components, size_t size) {
  ExpansionObject *result = (ExpansionObject *)(cls->tp_alloc(cls, 0));
  if (result) {
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
      PyMem_Calloc(self->size + other->size, sizeof(double));
  if (!result_components) return (ExpansionObject *)PyErr_NoMemory();
  size_t result_size = add_components_eliminating_zeros(
      self->size, self->components, other->size, other->components,
      result_components);
  result_size = compress_components(result_size, result_components);
  if (!PyMem_Resize(result_components, double, result_size))
    return (ExpansionObject *)PyErr_NoMemory();
  return construct_Expansion(&ExpansionType, result_components, result_size);
}

static ExpansionObject *Expansion_double_add(ExpansionObject *self,
                                             double other) {
  double *result_components = PyMem_Calloc(self->size + 1, sizeof(double));
  if (!result_components) return (ExpansionObject *)PyErr_NoMemory();
  size_t result_size = add_double_eliminating_zeros(
      self->size, self->components, other, result_components);
  result_size = compress_components(result_size, result_components);
  if (!PyMem_Resize(result_components, double, result_size))
    return (ExpansionObject *)PyErr_NoMemory();
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
    else if (PyObject_IsInstance(other, Real)) {
      double other_value = PyFloat_AsDouble(other);
      return other_value == -1.0 && PyErr_Occurred()
                 ? NULL
                 : (PyObject *)Expansion_double_add((ExpansionObject *)self,
                                                    other_value);
    }
  } else if (PyFloat_Check(self))
    return (PyObject *)Expansion_double_add((ExpansionObject *)other,
                                            PyFloat_AS_DOUBLE(self));
  else if (PyObject_IsInstance(self, Real)) {
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
  PyMem_Free(self->components);
  Py_TYPE(self)->tp_free((PyObject *)self);
}

static double Expansion_double(ExpansionObject *self) {
  return sum_components(self->size, self->components);
}

static PyObject *Expansion_ceil(ExpansionObject *self,
                                PyObject *Py_UNUSED(args)) {
  double value = Expansion_double(self);
  PyObject *result = PyLong_FromDouble(value);
  if (value > 0.0) {
    double integer_part;
    PyObject *increment = PyLong_FromLong(!!modf(value, &integer_part));
    if (!increment) {
      Py_DECREF(result);
      return NULL;
    }
    PyObject *tmp = result;
    result = PyNumber_Add(result, increment);
    Py_DECREF(tmp);
    Py_DECREF(increment);
  }
  return result;
}

static PyObject *Expansion_float(ExpansionObject *self) {
  return PyFloat_FromDouble(Expansion_double(self));
}

static PyObject *Expansion_floor(ExpansionObject *self,
                                 PyObject *Py_UNUSED(args)) {
  double value = Expansion_double(self);
  PyObject *result = PyLong_FromDouble(value);
  if (value < 0.0) {
    double integer_part;
    PyObject *increment = PyLong_FromLong(!!modf(value, &integer_part));
    if (!increment) {
      Py_DECREF(result);
      return NULL;
    }
    PyObject *tmp = result;
    result = PyNumber_Add(result, increment);
    Py_DECREF(tmp);
    Py_DECREF(increment);
  }
  return result;
}

static PyObject *Expansion_floor_divide(PyObject *self, PyObject *other) {
  if (PyObject_TypeCheck(self, &ExpansionType)) {
    PyObject *self_float = Expansion_float((ExpansionObject *)self);
    if (!self_float) return NULL;
    PyObject *result = PyNumber_FloorDivide(self_float, other);
    Py_DECREF(self_float);
    return result;
  } else if (PyObject_IsInstance(self, Real)) {
    PyObject *other_float = Expansion_float((ExpansionObject *)other);
    if (!other_float) return NULL;
    PyObject *result = PyNumber_FloorDivide(self, other_float);
    Py_DECREF(other_float);
    return result;
  }
  Py_RETURN_NOTIMPLEMENTED;
}

static Py_hash_t Expansion_hash(ExpansionObject *self) {
  PyObject *self_float = Expansion_float(self);
  Py_hash_t result = PyObject_Hash(self_float);
  Py_DECREF(self_float);
  return result;
}

static ExpansionObject *Expansions_multiply(ExpansionObject *self,
                                            ExpansionObject *other) {
  if (self->size < other->size) {
    ExpansionObject *tmp = self;
    self = other;
    other = tmp;
  }
  double *result_components =
      PyMem_Calloc(2 * self->size * other->size, sizeof(double));
  if (!result_components) return (ExpansionObject *)PyErr_NoMemory();
  double *step_components = PyMem_Calloc(2 * self->size, sizeof(double));
  if (!step_components) {
    PyMem_Free(result_components);
    return (ExpansionObject *)PyErr_NoMemory();
  }
  size_t result_size = multiply_components_eliminating_zeros(
      self->size, self->components, other->size, other->components,
      step_components, result_components);
  PyMem_Free(step_components);
  result_size = compress_components(result_size, result_components);
  if (!PyMem_Resize(result_components, double, result_size))
    return (ExpansionObject *)PyErr_NoMemory();
  return construct_Expansion(&ExpansionType, result_components, result_size);
}

static ExpansionObject *Expansion_double_multiply(ExpansionObject *self,
                                                  double other) {
  double *result_components = PyMem_Calloc(2 * self->size, sizeof(double));
  if (!result_components) return (ExpansionObject *)PyErr_NoMemory();
  size_t result_size = scale_components_eliminating_zeros(
      self->size, self->components, other, result_components);
  result_size = compress_components(result_size, result_components);
  if (!PyMem_Resize(result_components, double, result_size))
    return (ExpansionObject *)PyErr_NoMemory();
  return construct_Expansion(&ExpansionType, result_components, result_size);
}

static PyObject *Expansion_multiply(PyObject *self, PyObject *other) {
  if (PyObject_TypeCheck(self, &ExpansionType)) {
    if (PyFloat_Check(other))
      return (PyObject *)Expansion_double_multiply((ExpansionObject *)self,
                                                   PyFloat_AS_DOUBLE(other));
    else if (PyObject_TypeCheck(other, &ExpansionType))
      return (PyObject *)Expansions_multiply((ExpansionObject *)self,
                                             (ExpansionObject *)other);
    else if (PyObject_IsInstance(other, Real)) {
      double other_value = PyFloat_AsDouble(other);
      return other_value == -1.0 && PyErr_Occurred()
                 ? NULL
                 : (PyObject *)Expansion_double_multiply(
                       (ExpansionObject *)self, other_value);
    }
  } else if (PyFloat_Check(self))
    return (PyObject *)Expansion_double_multiply((ExpansionObject *)other,
                                                 PyFloat_AS_DOUBLE(self));
  else if (PyObject_IsInstance(self, Real)) {
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
  if (size == 1) {
    PyObject *argument = PyTuple_GET_ITEM(args, 0);
    if (PyObject_TypeCheck(argument, &ExpansionType)) {
      ExpansionObject *expansion_argument = (ExpansionObject *)argument;
      components = PyMem_Calloc(expansion_argument->size, sizeof(double));
      if (!components) return NULL;
      copy_components(expansion_argument->components, expansion_argument->size,
                      components);
      size = expansion_argument->size;
    } else {
      components = (double *)PyMem_Malloc(sizeof(double));
      if (!components) return PyErr_NoMemory();
      double value = PyFloat_AsDouble(argument);
      if (value == 1.0 && PyErr_Occurred()) return NULL;
      components[0] = value;
      size = 1;
    }
  } else if (size) {
    components = (double *)PyMem_Calloc(size, sizeof(double));
    if (!components) return PyErr_NoMemory();
    for (size_t index = 0; index < size; ++index) {
      PyObject *item = PyTuple_GET_ITEM(args, index);
      if (!item) {
        PyMem_Free(components);
        return NULL;
      }
      double component = components[index] = PyFloat_AsDouble(item);
      if (component == -1.0 && PyErr_Occurred()) {
        PyMem_Free(components);
        return NULL;
      }
    }
    size = compress_components(size, components);
    if (!PyMem_Resize(components, double, size)) return PyErr_NoMemory();
  } else {
    components = (double *)PyMem_Malloc(sizeof(double));
    if (!components) return PyErr_NoMemory();
    components[0] = 0.0;
    size = 1;
  }
  return (PyObject *)construct_Expansion(cls, components, size);
}

static ExpansionObject *Expansion_negative(ExpansionObject *self) {
  double *result_components = PyMem_Calloc(self->size, sizeof(double));
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

static ExpansionObject *Expansion_double_remainder(ExpansionObject *self,
                                                   double other) {
  if (!other) {
    PyErr_SetString(PyExc_ZeroDivisionError, "Zero divisor.");
    return NULL;
  }
  double *result_components = PyMem_Calloc(self->size, sizeof(double));
  size_t result_size;
  result_components[0] = double_remainder(self->components[0], other);
  result_size = 1;
  for (size_t index = 1; index < self->size; ++index)
    result_size = add_double_eliminating_zeros(
        result_size, result_components,
        double_remainder(self->components[index], other), result_components);
  result_components[result_size - 1] =
      double_remainder(result_components[result_size - 1], other);
  result_size = compress_components_single(result_size, result_components);
  return construct_Expansion(&ExpansionType, result_components, result_size);
}

static PyObject *Expansion_remainder(PyObject *self, PyObject *other) {
  if (PyObject_TypeCheck(self, &ExpansionType)) {
    if (PyObject_IsInstance(other, Real)) {
      double other_value = PyFloat_AsDouble(other);
      return other_value == -1.0 && PyErr_Occurred()
                 ? NULL
                 : (PyObject *)Expansion_double_remainder(
                       (ExpansionObject *)self, other_value);
    }
  } else if (PyObject_IsInstance(self, Real)) {
    PyObject *other_float = Expansion_float((ExpansionObject *)other);
    if (!other_float) return NULL;
    PyObject *result = PyNumber_Remainder(self, other_float);
    Py_DECREF(other_float);
    return result;
  }
  Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *Expansion_power(PyObject *self, PyObject *exponent,
                                 PyObject *modulo) {
  if (PyObject_TypeCheck(self, &ExpansionType)) {
    PyObject *self_float = Expansion_float((ExpansionObject *)self);
    PyObject *result = PyNumber_Power(self_float, exponent, modulo);
    Py_DECREF(self_float);
    return result;
  } else if (PyObject_TypeCheck(exponent, &ExpansionType)) {
    PyObject *exponent_float = Expansion_float((ExpansionObject *)self);
    PyObject *result = PyNumber_Power(self, exponent_float, modulo);
    Py_DECREF(exponent_float);
    return result;
  }
  Py_RETURN_NOTIMPLEMENTED;
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

static PyObject *Expansion_PyObject_richcompare(ExpansionObject *self,
                                                PyObject *other, int op) {
  PyObject *self_float = Expansion_float(self);
  if (!self_float) return NULL;
  PyObject *result = PyObject_RichCompare(self_float, other, op);
  Py_DECREF(self_float);
  return result;
}

static PyObject *Expansion_richcompare(ExpansionObject *self, PyObject *other,
                                       int op) {
  if (PyObject_TypeCheck(other, &ExpansionType))
    return Expansions_richcompare(self, (ExpansionObject *)other, op);
  else if (PyObject_IsInstance(other, Real))
    return Expansion_PyObject_richcompare(self, other, op);
  else
    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *Expansion_round(ExpansionObject *self, PyObject *args) {
  PyObject *self_float = Expansion_float(self);
  if (!self_float) return NULL;
  Py_ssize_t size = PyTuple_Size(args);
  if (size < 0) {
    Py_DECREF(self_float);
    return NULL;
  }
  PyObject *full_args = PyTuple_New(size + 1);
  if (!full_args) {
    Py_DECREF(self_float);
    return NULL;
  }
  PyTuple_SET_ITEM(full_args, 0, self_float);
  for (Py_ssize_t index = 0; index < size; ++index) {
    PyObject *argument = PyTuple_GET_ITEM(args, index);
    Py_INCREF(argument);
    PyTuple_SET_ITEM(full_args, index + 1, argument);
  }
  PyObject *result = PyObject_CallObject(PyObject_round, full_args);
  Py_DECREF(full_args);
  return result;
}

static ExpansionObject *Expansions_subtract(ExpansionObject *self,
                                            ExpansionObject *other) {
  double *result_components =
      PyMem_Calloc(self->size + other->size, sizeof(double));
  if (!result_components) return (ExpansionObject *)PyErr_NoMemory();
  size_t result_size = subtract_components_eliminating_zeros(
      self->size, self->components, other->size, other->components,
      result_components);
  result_size = compress_components(result_size, result_components);
  if (!PyMem_Resize(result_components, double, result_size))
    return (ExpansionObject *)PyErr_NoMemory();
  return construct_Expansion(&ExpansionType, result_components, result_size);
}

static ExpansionObject *Expansion_double_subtract(ExpansionObject *self,
                                                  double other) {
  double *result_components = PyMem_Calloc(self->size + 1, sizeof(double));
  if (!result_components) return (ExpansionObject *)PyErr_NoMemory();
  size_t result_size = subtract_double_eliminating_zeros(
      self->size, self->components, other, result_components);
  result_size = compress_components(result_size, result_components);
  if (!PyMem_Resize(result_components, double, result_size))
    return (ExpansionObject *)PyErr_NoMemory();
  return construct_Expansion(&ExpansionType, result_components, result_size);
}

static ExpansionObject *double_Expansion_subtract(double self,
                                                  ExpansionObject *other) {
  double *result_components = PyMem_Calloc(other->size + 1, sizeof(double));
  if (!result_components) return (ExpansionObject *)PyErr_NoMemory();
  size_t result_size = subtract_from_double_eliminating_zeros(
      self, other->size, other->components, result_components);
  result_size = compress_components(result_size, result_components);
  if (!PyMem_Resize(result_components, double, result_size))
    return (ExpansionObject *)PyErr_NoMemory();
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
    else if (PyObject_IsInstance(other, Real)) {
      double other_value = PyFloat_AsDouble(other);
      return other_value == -1.0 && PyErr_Occurred()
                 ? NULL
                 : (PyObject *)Expansion_double_subtract(
                       (ExpansionObject *)self, other_value);
    }
  } else if (PyFloat_Check(self))
    return (PyObject *)double_Expansion_subtract(PyFloat_AS_DOUBLE(self),
                                                 (ExpansionObject *)other);
  else if (PyObject_IsInstance(self, Real)) {
    double value = PyFloat_AsDouble(self);
    return value == -1.0 && PyErr_Occurred()
               ? NULL
               : (PyObject *)double_Expansion_subtract(
                     value, (ExpansionObject *)other);
  }
  Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *PyObject_Expansion_true_divide(PyObject *self,
                                                ExpansionObject *other) {
  if (!Expansion_bool(other)) {
    PyErr_Format(PyExc_ZeroDivisionError, "Divisor is zero.");
    return NULL;
  }
  PyObject *other_float = Expansion_float(other);
  if (!other_float) return NULL;
  PyObject *result = PyNumber_TrueDivide(self, other_float);
  Py_DECREF(other_float);
  return result;
}

static PyObject *Expansion_true_divide(PyObject *self, PyObject *other) {
  if (PyObject_TypeCheck(self, &ExpansionType)) {
    if (PyFloat_Check(other)) {
      double other_value = PyFloat_AS_DOUBLE(other);
      if (!other_value) {
        PyErr_Format(PyExc_ZeroDivisionError, "Divisor is zero.");
        return NULL;
      }
      return (PyObject *)Expansion_double_multiply((ExpansionObject *)self,
                                                   1.0 / other_value);
    } else if (PyObject_TypeCheck(other, &ExpansionType) ||
               PyObject_IsInstance(other, Real)) {
      double other_value = PyFloat_AsDouble(other);
      if (other_value == -1.0 && PyErr_Occurred())
        return NULL;
      else if (!other_value) {
        PyErr_Format(PyExc_ZeroDivisionError, "Divisor is zero.");
        return NULL;
      }
      return (PyObject *)Expansion_double_multiply((ExpansionObject *)self,
                                                   1.0 / other_value);
    }
  } else if (PyFloat_Check(self) || PyObject_IsInstance(self, Real))
    return PyObject_Expansion_true_divide(self, (ExpansionObject *)other);
  Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *Expansion_trunc(ExpansionObject *self,
                                 PyObject *Py_UNUSED(args)) {
  return PyLong_FromDouble(Expansion_double(self));
}

static PyNumberMethods Expansion_as_number = {
    .nb_absolute = (unaryfunc)Expansion_absolute,
    .nb_add = Expansion_add,
    .nb_bool = (inquiry)Expansion_bool,
    .nb_float = (unaryfunc)Expansion_float,
    .nb_floor_divide = Expansion_floor_divide,
    .nb_multiply = Expansion_multiply,
    .nb_negative = (unaryfunc)Expansion_negative,
    .nb_positive = (unaryfunc)Expansion_positive,
    .nb_power = Expansion_power,
    .nb_remainder = Expansion_remainder,
    .nb_subtract = Expansion_subtract,
    .nb_true_divide = Expansion_true_divide,
};

PyObject *Expansion_getimag(ExpansionObject *self, void *closure) {
  return PyLong_FromLong(0);
}

PyObject *Expansion_getreal(ExpansionObject *self, void *closure) {
  return (PyObject *)Expansion_positive(self);
}

static PyGetSetDef Expansion_getset[] = {
    {"real", (getter)Expansion_getreal, (setter)NULL,
     "The real part of the expansion.", NULL},
    {"imag", (getter)Expansion_getimag, (setter)NULL,
     "The imaginary part of the expansion.", NULL},
    {NULL} /* Sentinel */
};

static PyMethodDef Expansion_methods[] = {
    {"__ceil__", (PyCFunction)Expansion_ceil, METH_NOARGS, NULL},
    {"__floor__", (PyCFunction)Expansion_floor, METH_NOARGS, NULL},
    {"__round__", (PyCFunction)Expansion_round, METH_VARARGS, NULL},
    {"__trunc__", (PyCFunction)Expansion_trunc, METH_NOARGS, NULL},
    {NULL, NULL} /* sentinel */
};

static PyTypeObject ExpansionType = {
    PyVarObject_HEAD_INIT(NULL, 0).tp_as_number = &Expansion_as_number,
    .tp_basicsize = sizeof(ExpansionObject),
    .tp_dealloc = (destructor)Expansion_dealloc,
    .tp_doc = PyDoc_STR("Represents floating point number expansion."),
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_getset = Expansion_getset,
    .tp_hash = (hashfunc)Expansion_hash,
    .tp_itemsize = 0,
    .tp_methods = Expansion_methods,
    .tp_name = "shewchuk.Expansion",
    .tp_new = Expansion_new,
    .tp_repr = (reprfunc)Expansion_repr,
    .tp_richcompare = (richcmpfunc)Expansion_richcompare,
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
  double *result_components = PyMem_Calloc(result_size, sizeof(double));
  if (!result_components) return PyErr_NoMemory();
  copy_components(components, result_size, result_components);
  return (PyObject *)construct_Expansion(&ExpansionType, result_components,
                                         result_size);
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
  double *result_components = PyMem_Calloc(result_size, sizeof(double));
  if (!result_components) return PyErr_NoMemory();
  copy_components(components, result_size, result_components);
  return (PyObject *)construct_Expansion(&ExpansionType, result_components,
                                         result_size);
}

static PyMethodDef _shewchuk_methods[] = {
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
    {NULL, NULL},
};

static PyModuleDef _shewchuk_module = {
    PyModuleDef_HEAD_INIT,
    .m_doc = PyDoc_STR("Robust floating point operations."),
    .m_methods = _shewchuk_methods,
    .m_name = "shewchuk",
    .m_size = -1,
};

static int load_PyObject_round() {
  PyObject *builtins_module = PyImport_ImportModule("builtins");
  if (!builtins_module) return -1;
  PyObject_round = PyObject_GetAttrString(builtins_module, "round");
  Py_DECREF(builtins_module);
  return !PyObject_round ? -1 : 0;
}

static int load_real() {
  PyObject *numbers_module = PyImport_ImportModule("numbers");
  if (!numbers_module) return -1;
  Real = PyObject_GetAttrString(numbers_module, "Real");
  Py_DECREF(numbers_module);
  return !Real ? -1 : 0;
}

static int mark_as_real(PyObject *python_type) {
  PyObject *register_method_name = PyUnicode_FromString("register");
  if (!register_method_name) return -1;
  PyObject *tmp =
#if PY39_OR_MORE
      PyObject_CallMethodOneArg(Real, register_method_name, python_type);
#else
      PyObject_CallMethodObjArgs(Real, register_method_name, python_type, NULL)
#endif
  ;
  Py_DECREF(register_method_name);
  if (!tmp) return -1;
  Py_DECREF(tmp);
  return 0;
}

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
  if (load_PyObject_round() < 0) {
    Py_DECREF(result);
    return NULL;
  }
  if (load_real() < 0) {
    Py_DECREF(PyObject_round);
    Py_DECREF(result);
    return NULL;
  }
  if (mark_as_real((PyObject *)&ExpansionType) < 0) {
    Py_DECREF(PyObject_round);
    Py_DECREF(Real);
    Py_DECREF(result);
    return NULL;
  }
  return result;
}
