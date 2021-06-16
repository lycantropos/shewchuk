"""Robust floating point operations."""

__version__ = '0.0.0'

try:
    from _shewchuk import (Expansion,
                           vectors_cross_product)
except ImportError:
    from sys import float_info as _float_info
    from itertools import repeat as _repeat
    from numbers import Real as _Real
    from typing import (Sequence as _Sequence,
                        Tuple as _Tuple,
                        Union as _Union)


    class Expansion:
        __slots__ = '_components',

        def __new__(cls, *components: float,
                    _compress: bool = True) -> 'Expansion':
            self = super().__new__(cls)
            components = [float(component) for component in components]
            if _compress and len(components) > 1:
                components = _compress_components(components)
            elif not components:
                components = [0.]
            self._components = components
            return self

        def __abs__(self) -> 'Expansion':
            return +self if self._components[-1] > 0. else -self

        def __add__(self, other: _Union[_Real, 'Expansion']) -> 'Expansion':
            return (Expansion(*_add_components_eliminating_zeros(
                    self._components, other._components))
                    if isinstance(other, Expansion)
                    else self.__radd__(other))

        def __bool__(self) -> bool:
            return bool(self._components[-1])

        def __eq__(self, other: _Union[_Real, 'Expansion']) -> bool:
            return (self._components == other._components
                    if isinstance(other, Expansion)
                    else (len(self._components) == 1
                          and self._components[0] == other
                          if isinstance(other, _Real)
                          else NotImplemented))

        def __float__(self) -> float:
            return sum(self._components)

        def __ge__(self, other: _Union[_Real, 'Expansion']) -> bool:
            return (not _are_components_lesser_than(self._components,
                                                    other._components)
                    if isinstance(other, Expansion)
                    else (self._components[-1] > other
                          or self._components[-1] == other
                          and (len(self._components) == 1
                               or self._components[-2] > 0.)
                          if isinstance(other, _Real)
                          else NotImplemented))

        def __gt__(self, other: _Union[_Real, 'Expansion']) -> bool:
            return (_are_components_lesser_than(other._components,
                                                self._components)
                    if isinstance(other, Expansion)
                    else (self._components[-1] > other
                          or (self._components[-1] == other
                              and len(self._components) > 1
                              and self._components[-2] > 0.)
                          if isinstance(other, _Real)
                          else NotImplemented))

        def __le__(self, other: _Union[_Real, 'Expansion']) -> bool:
            return (not _are_components_lesser_than(other._components,
                                                    self._components)
                    if isinstance(other, Expansion)
                    else (self._components[-1] < other
                          or self._components[-1] == other
                          and (len(self._components) == 1
                               or self._components[-2] < 0.)
                          if isinstance(other, _Real)
                          else NotImplemented))

        def __lt__(self, other: _Union[_Real, 'Expansion']) -> bool:
            return (_are_components_lesser_than(self._components,
                                                other._components)
                    if isinstance(other, Expansion)
                    else (self._components[-1] < other
                          or (self._components[-1] == other
                              and len(self._components) > 1
                              and self._components[-2] < 0.)
                          if isinstance(other, _Real)
                          else NotImplemented))

        def __rmul__(self, other: _Real) -> 'Expansion':
            return (Expansion(*_scale_components(self._components,
                                                 float(other)))
                    if isinstance(other, _Real)
                    else NotImplemented)

        __mul__ = __rmul__

        def __neg__(self) -> 'Expansion':
            return Expansion(*[-component for component in self._components])

        def __pos__(self) -> 'Expansion':
            return self

        def __radd__(self, other: _Real) -> 'Expansion':
            return (Expansion(*_add_float_eliminating_zeros(self._components,
                                                            float(other)))
                    if isinstance(other, _Real)
                    else NotImplemented)

        def __repr__(self) -> str:
            return (type(self).__qualname__
                    + '({})'.format(', '.join(map(str, self._components))))

        def __rsub__(self, other: _Union[_Real, 'Expansion']) -> 'Expansion':
            return (Expansion(*_subtract_from_double_eliminating_zeros(
                    float(other), self._components))
                    if isinstance(other, _Real)
                    else NotImplemented)

        def __sub__(self, other: _Union[_Real, 'Expansion']) -> 'Expansion':
            return (Expansion(*_subtract_components_eliminating_zeros(
                    self._components, other._components))
                    if isinstance(other, Expansion)
                    else
                    (Expansion(*_subtract_double_eliminating_zeros(
                            self._components, float(other)))
                     if isinstance(other, _Real)
                     else NotImplemented))


    def _are_components_lesser_than(left: _Sequence[float],
                                    right: _Sequence[float]) -> bool:
        left_size, right_size = len(left), len(right)
        for offset in range(min(left_size, right_size)):
            if left[left_size - 1 - offset] < right[right_size - 1 - offset]:
                return True
            elif left[left_size - 1 - offset] > right[right_size - 1 - offset]:
                return False
        return (left_size != right_size
                and (right[right_size - left_size - 1] > 0.
                     if left_size < right_size
                     else left[left_size - right_size - 1] < 0.))


    _EPSILON = _float_info.epsilon / 2.0


    def _adaptive_vectors_cross_product(first_start_x: float,
                                        first_start_y: float,
                                        first_end_x: float,
                                        first_end_y: float,
                                        second_start_x: float,
                                        second_start_y: float,
                                        second_end_x: float,
                                        second_end_y: float,
                                        upper_bound: float,
                                        first_upper_bound_coefficient: float
                                        = (2.0 + 12.0 * _EPSILON) * _EPSILON,
                                        second_upper_bound_coefficient: float
                                        = ((9.0 + 64.0 * _EPSILON) * _EPSILON
                                           * _EPSILON),
                                        estimation_coefficient: float
                                        = (3.0 + 8.0 * _EPSILON) * _EPSILON
                                        ) -> _Sequence[float]:
        minuend_x = first_end_x - first_start_x
        minuend_y = first_end_y - first_start_y
        subtrahend_x = second_end_x - second_start_x
        subtrahend_y = second_end_y - second_start_y
        minuend_tail, minuend = _two_multiply(minuend_x, subtrahend_y)
        subtrahend_tail, subtrahend = _two_multiply(minuend_y, subtrahend_x)
        first_components = _two_two_subtract(minuend, minuend_tail, subtrahend,
                                             subtrahend_tail)
        estimation = sum(first_components)
        threshold = first_upper_bound_coefficient * upper_bound
        if (estimation >= threshold) or (-estimation >= threshold):
            return first_components
        minuend_x_tail = _two_subtract_tail(first_end_x, first_start_x,
                                            minuend_x)
        subtrahend_x_tail = _two_subtract_tail(second_end_x, second_start_x,
                                               subtrahend_x)
        minuend_y_tail = _two_subtract_tail(first_end_y, first_start_y,
                                            minuend_y)
        subtrahend_y_tail = _two_subtract_tail(second_end_y, second_start_y,
                                               subtrahend_y)
        if (not minuend_x_tail and not minuend_y_tail and not subtrahend_x_tail
                and not subtrahend_y_tail):
            return first_components
        threshold = (second_upper_bound_coefficient * upper_bound
                     + estimation_coefficient * abs(estimation))
        extra = ((minuend_x * subtrahend_y_tail
                  + subtrahend_y * minuend_x_tail)
                 - (minuend_y * subtrahend_x_tail
                    + subtrahend_x * minuend_y_tail))
        estimation += extra
        if (estimation >= threshold) or (-estimation >= threshold):
            return _add_float_eliminating_zeros(first_components, extra)
        minuend_x_subtrahend_y_tail, minuend_x_subtrahend_y = _two_multiply(
                minuend_x_tail, subtrahend_y)
        minuend_y_subtrahend_x_tail, minuend_y_subtrahend_x = _two_multiply(
                minuend_y_tail, subtrahend_x)
        extra_components = _two_two_subtract(
                minuend_x_subtrahend_y, minuend_x_subtrahend_y_tail,
                minuend_y_subtrahend_x, minuend_y_subtrahend_x_tail)
        second_components = _add_components_eliminating_zeros(
                first_components, extra_components)
        minuend_x_subtrahend_y_tail, minuend_x_subtrahend_y = _two_multiply(
                minuend_x, subtrahend_y_tail)
        minuend_y_subtrahend_x_tail, minuend_y_subtrahend_x = _two_multiply(
                minuend_y, subtrahend_x_tail)
        extra_components = _two_two_subtract(
                minuend_x_subtrahend_y, minuend_x_subtrahend_y_tail,
                minuend_y_subtrahend_x, minuend_y_subtrahend_x_tail)
        third_components = _add_components_eliminating_zeros(
                second_components, extra_components)
        minuend_x_subtrahend_y_tail, minuend_x_subtrahend_y = _two_multiply(
                minuend_x_tail, subtrahend_y_tail)
        minuend_y_subtrahend_x_tail, minuend_y_subtrahend_x = _two_multiply(
                minuend_y_tail, subtrahend_x_tail)
        extra_components = _two_two_subtract(
                minuend_x_subtrahend_y, minuend_x_subtrahend_y_tail,
                minuend_y_subtrahend_x, minuend_y_subtrahend_x_tail)
        return _add_components_eliminating_zeros(third_components,
                                                 extra_components)


    def vectors_cross_product(first_start_x: float,
                              first_start_y: float,
                              first_end_x: float,
                              first_end_y: float,
                              second_start_x: float,
                              second_start_y: float,
                              second_end_x: float,
                              second_end_y: float) -> Expansion:
        return Expansion(*_vectors_cross_product(first_start_x, first_start_y,
                                                 first_end_x, first_end_y,
                                                 second_start_x,
                                                 second_start_y, second_end_x,
                                                 second_end_y),
                         _compress=False)


    def _vectors_cross_product(first_start_x: float,
                               first_start_y: float,
                               first_end_x: float,
                               first_end_y: float,
                               second_start_x: float,
                               second_start_y: float,
                               second_end_x: float,
                               second_end_y: float,
                               upper_bound_coefficient: float
                               = (3.0 + 16.0 * _EPSILON) * _EPSILON
                               ) -> _Sequence[float]:
        minuend = ((first_end_x - first_start_x)
                   * (second_end_y - second_start_y))
        subtrahend = ((first_end_y - first_start_y)
                      * (second_end_x - second_start_x))
        estimation = minuend - subtrahend
        if minuend > 0.0:
            if subtrahend <= 0.0:
                return [estimation]
            else:
                upper_bound = minuend + subtrahend
        elif minuend < 0.0:
            if subtrahend >= 0.0:
                return [estimation]
            else:
                upper_bound = -minuend - subtrahend
        else:
            return [estimation]
        threshold = upper_bound_coefficient * upper_bound
        return ([estimation]
                if (estimation >= threshold) or (-estimation >= threshold)
                else
                _adaptive_vectors_cross_product(first_start_x, first_start_y,
                                                first_end_x, first_end_y,
                                                second_start_x, second_start_y,
                                                second_end_x, second_end_y,
                                                upper_bound))


    def _add_components_eliminating_zeros(left: _Sequence[float],
                                          right: _Sequence[float]
                                          ) -> _Sequence[float]:
        left_length, right_length = len(left), len(right)
        left_component, right_component = left[0], right[0]
        left_index = right_index = 0
        if ((right_component > left_component)
                is (right_component > -left_component)):
            accumulator = left_component
            left_index += 1
        else:
            accumulator = right_component
            right_index += 1
        result = []
        if (left_index < left_length) and (right_index < right_length):
            left_component, right_component = (left[left_index],
                                               right[right_index])
            if ((right_component > left_component)
                    is (right_component > -left_component)):
                tail, accumulator = _fast_two_add(left_component, accumulator)
                left_index += 1
            else:
                tail, accumulator = _fast_two_add(right_component, accumulator)
                right_index += 1
            if tail:
                result.append(tail)
            while (left_index < left_length) and (right_index < right_length):
                left_component, right_component = (left[left_index],
                                                   right[right_index])
                if ((right_component > left_component)
                        is (right_component > -left_component)):
                    tail, accumulator = _two_add(accumulator, left_component)
                    left_index += 1
                else:
                    tail, accumulator = _two_add(accumulator, right_component)
                    right_index += 1
                if tail:
                    result.append(tail)
        for left_index in range(left_index, left_length):
            left_component = left[left_index]
            tail, accumulator = _two_add(accumulator, left_component)
            if tail:
                result.append(tail)
        for right_index in range(right_index, right_length):
            right_component = right[right_index]
            tail, accumulator = _two_add(accumulator, right_component)
            if tail:
                result.append(tail)
        if accumulator or not result:
            result.append(accumulator)
        return result


    def _add_float_eliminating_zeros(left: _Sequence[float],
                                     right: float) -> _Sequence[float]:
        result = []
        accumulator = right
        for left_component in left:
            tail, accumulator = _two_add(accumulator, left_component)
            if tail:
                result.append(tail)
        if accumulator or not result:
            result.append(accumulator)
        return result


    def _ceil_divide(dividend: int, divisor: int) -> int:
        return -(-dividend // divisor)


    def _compress_components(components: _Sequence[float]) -> _Sequence[float]:
        for _ in _repeat(None, len(components)):
            next_components = _compress_components_single(components)
            if next_components == components:
                break
            components = next_components
        return components


    def _compress_components_single(components: _Sequence[float]
                                    ) -> _Sequence[float]:
        bottom = len(components) - 1
        cursor = components[bottom]
        result = [None] * len(components)
        for index in range(bottom - 1, -1, -1):
            tail, head = _two_add(cursor, components[index])
            if tail:
                result[bottom] = head
                bottom -= 1
                cursor = tail
            else:
                cursor = head
        top = 0
        for index in range(bottom + 1, len(result)):
            tail, head = _two_add(result[index], cursor)
            if tail:
                result[top] = tail
                top += 1
            cursor = head
        if cursor or not top:
            result[top] = cursor
            top += 1
        return result[:top]


    def _fast_two_add(left: float, right: float) -> _Tuple[float, float]:
        head = left + right
        right_virtual = head - left
        tail = right - right_virtual
        return tail, head


    def _scale_components(components: _Sequence[float],
                          scalar: float) -> _Sequence[float]:
        components_iterator = iter(components)
        scalar_low, scalar_high = _split(scalar)
        tail, accumulator = _two_multiply_presplit(next(components_iterator),
                                                   scalar, scalar_high,
                                                   scalar_low)
        result = []
        if tail:
            result.append(tail)
        for component in components_iterator:
            product_tail, product = _two_multiply_presplit(component, scalar,
                                                           scalar_high,
                                                           scalar_low)
            tail, interim = _two_add(accumulator, product_tail)
            if tail:
                result.append(tail)
            tail, accumulator = _fast_two_add(product, interim)
            if tail:
                result.append(tail)
        if accumulator or not result:
            result.append(accumulator)
        return result


    def _split(value: float,
               *,
               splitter: float
               = float(1 << _ceil_divide(_float_info.mant_dig, 2) + 1)
               ) -> _Tuple[float, float]:
        base = splitter * value
        high = base - (base - value)
        low = value - high
        return low, high


    def _subtract_components_eliminating_zeros(minuend: _Sequence[float],
                                               subtrahend: _Sequence[float]
                                               ) -> _Sequence[float]:
        minuend_length, subtrahend_length = len(minuend), len(subtrahend)
        minuend_component, subtrahend_component = minuend[0], -subtrahend[0]
        minuend_index = subtrahend_index = 0
        if ((subtrahend_component > minuend_component)
                is (subtrahend_component > -minuend_component)):
            accumulator = minuend_component
            minuend_index += 1
        else:
            accumulator = subtrahend_component
            subtrahend_index += 1
        result = []
        if ((minuend_index < minuend_length)
                and (subtrahend_index < subtrahend_length)):
            minuend_component, subtrahend_component = (minuend[minuend_index],
                                                       -subtrahend[
                                                           subtrahend_index])
            if ((subtrahend_component > minuend_component)
                    is (subtrahend_component > -minuend_component)):
                tail, accumulator = _fast_two_add(minuend_component,
                                                  accumulator)
                minuend_index += 1
            else:
                tail, accumulator = _fast_two_add(subtrahend_component,
                                                  accumulator)
                subtrahend_index += 1
            if tail:
                result.append(tail)
            while ((minuend_index < minuend_length)
                   and (subtrahend_index < subtrahend_length)):
                minuend_component, subtrahend_component = (
                    minuend[minuend_index], -subtrahend[subtrahend_index])
                if ((subtrahend_component > minuend_component)
                        is (subtrahend_component > -minuend_component)):
                    tail, accumulator = _two_add(accumulator,
                                                 minuend_component)
                    minuend_index += 1
                else:
                    tail, accumulator = _two_add(accumulator,
                                                 subtrahend_component)
                    subtrahend_index += 1
                if tail:
                    result.append(tail)
        for minuend_index in range(minuend_index, minuend_length):
            minuend_component = minuend[minuend_index]
            tail, accumulator = _two_add(accumulator, minuend_component)
            if tail:
                result.append(tail)
        for subtrahend_index in range(subtrahend_index, subtrahend_length):
            subtrahend_component = -subtrahend[subtrahend_index]
            tail, accumulator = _two_add(accumulator, subtrahend_component)
            if tail:
                result.append(tail)
        if accumulator or not result:
            result.append(accumulator)
        return result


    def _subtract_double_eliminating_zeros(minuend: _Sequence[float],
                                           subtrahend: float
                                           ) -> _Sequence[float]:
        return _add_float_eliminating_zeros(minuend, -subtrahend)


    def _subtract_from_double_eliminating_zeros(minuend: float,
                                                subtrahend: _Sequence[float]
                                                ) -> _Sequence[float]:
        result = []
        accumulator = minuend
        for subtrahend_component in subtrahend:
            tail, accumulator = _two_add(accumulator, -subtrahend_component)
            if tail:
                result.append(tail)
        if accumulator or not result:
            result.append(accumulator)
        return result


    def _two_add(left: float, right: float) -> _Tuple[float, float]:
        head = left + right
        right_virtual = head - left
        left_virtual = head - right_virtual
        right_tail = right - right_virtual
        left_tail = left - left_virtual
        tail = left_tail + right_tail
        return tail, head


    def _two_subtract(left: float, right: float) -> _Tuple[float, float]:
        return _two_add(left, -right)


    def _two_one_subtract(left_tail: float,
                          left_head: float,
                          right: float) -> _Tuple[float, float, float]:
        second_tail, mid_head = _two_subtract(left_tail, right)
        first_tail, head = _two_subtract(left_head, mid_head)
        return second_tail, first_tail, head


    def _two_multiply(left: float, right: float) -> _Tuple[float, float]:
        head = left * right
        left_low, left_high = _split(left)
        right_low, right_high = _split(right)
        first_error = head - left_high * right_high
        second_error = first_error - left_low * right_high
        third_error = second_error - left_high * right_low
        tail = left_low * right_low - third_error
        return tail, head


    def _two_multiply_presplit(left: float,
                               right: float,
                               right_high: float,
                               right_low: float) -> _Tuple[float, float]:
        head = left * right
        left_low, left_high = _split(left)
        first_error = head - left_high * right_high
        second_error = first_error - left_low * right_high
        third_error = second_error - left_high * right_low
        tail = left_low * right_low - third_error
        return tail, head


    def _two_subtract_tail(left: float, right: float, head: float) -> float:
        right_virtual = left - head
        left_virtual = head + right_virtual
        right_tail = right_virtual - right
        left_tail = left - left_virtual
        return left_tail + right_tail


    def _two_two_subtract(left_head: float,
                          left_tail: float,
                          right_head: float,
                          right_tail: float
                          ) -> _Tuple[float, float, float, float]:
        third_tail, mid_tail, mid_head = _two_one_subtract(
                left_tail, left_head, right_tail)
        second_tail, first_tail, head = _two_one_subtract(mid_tail, mid_head,
                                                          right_head)
        return third_tail, second_tail, first_tail, head
