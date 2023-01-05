from __future__ import annotations

import typing as _t
from functools import reduce as _reduce
from itertools import (dropwhile as _dropwhile,
                       repeat as _repeat)
from math import (ceil as _ceil,
                  floor as _floor,
                  isfinite as _isfinite,
                  modf as _modf)
from numbers import (Rational as _Rational,
                     Real as _Real)
from operator import not_ as _not
from sys import float_info as _float_info

_Number = _t.Union[_Rational, float, int]


@_Real.register
class Expansion:
    """Represents floating point number expansion."""

    @property
    def real(self) -> Expansion:
        """The imaginary part of the expansion."""
        return self

    @property
    def imag(self) -> _Number:
        """The real part of the expansion."""
        return 0

    _components: _t.Sequence[float]

    __slots__ = '_components',

    def __new__(cls,
                _argument: _t.Union[Expansion, _Number] = 0.0,
                *args: float,
                _compress: bool = True) -> Expansion:
        self = super().__new__(cls)
        components: _t.Sequence[float]
        if args:
            try:
                invalid_component = (
                    _argument
                    if not isinstance(_argument, float)
                    else next(argument
                              for argument in args
                              if not isinstance(argument, float))
                )
            except StopIteration:
                pass
            else:
                raise TypeError(f'Components should be of type {float!r}, '
                                f'but found: {type(invalid_component)!r}.')
            assert isinstance(_argument, float), _argument
            components = [_argument, *args]
            if _compress:
                components = _compress_components(components)
        elif isinstance(_argument, Expansion):
            components = _argument._components
        elif isinstance(_argument, float):
            components = [_argument]
        elif isinstance(_argument, int):
            components = _int_to_components(_argument)
        elif isinstance(_argument, _Rational):
            components = _rational_to_components(_argument)
        else:
            raise TypeError('Argument should be of type '
                            f'{Expansion!r}, {_Rational!r}, '
                            f'{int!r} or {float!r}, '
                            f'but found: {type(_argument)!r}.')
        try:
            invalid_component = next(component
                                     for component in components
                                     if not _isfinite(component))
        except StopIteration:
            pass
        else:
            raise ValueError('Components should be finite, '
                             'but found: {!r}.'
                             .format(invalid_component))
        self._components = tuple(components)
        return self

    def __abs__(self) -> Expansion:
        return +self if self._components[-1] > 0. else -self

    @_t.overload
    def __add__(self, other: _t.Union[Expansion, _Number]) -> Expansion:
        ...

    @_t.overload
    def __add__(self, other: _t.Any) -> _t.Any:
        ...

    def __add__(
            self, other: _t.Union[Expansion, _Number, _t.Any]
    ) -> _t.Union[Expansion, _t.Any]:
        return (Expansion(*_add_components(self._components,
                                           other._components))
                if isinstance(other, Expansion)
                else self.__radd__(other))

    def __bool__(self) -> bool:
        return bool(self._components[-1])

    def __ceil__(self) -> int:
        return (_components_to_integer(self._components)
                + _ceil(_components_to_accumulated_fraction(
                        self._components)))

    @_t.overload
    def __eq__(self, other: _t.Union[Expansion, _Number]) -> bool:
        ...

    @_t.overload
    def __eq__(self, other: _t.Any) -> _t.Any:
        ...

    def __eq__(self,
               other: _t.Union[Expansion, _Number]) -> _t.Union[_t.Any, bool]:
        return (self._components == other._components
                if isinstance(other, Expansion)
                else
                (_are_components_equal_to_float(self._components, other)
                 if isinstance(other, float)
                 else
                 (_are_components_equal_to_int(self._components, other)
                  if isinstance(other, int)
                  else (_are_components_equal_to_rational(self._components,
                                                          other)
                        if isinstance(other, _Rational)
                        else NotImplemented))))

    def __float__(self) -> float:
        assert sum(self._components) == self._components[-1], self
        return self._components[-1]

    def __floor__(self) -> int:
        return (_components_to_integer(self._components)
                + _floor(_components_to_accumulated_fraction(
                        self._components)))

    @_t.overload
    def __ge__(self, other: _t.Union[Expansion, _Number]) -> bool:
        ...

    @_t.overload
    def __ge__(self, other: _t.Any) -> _t.Any:
        ...

    def __ge__(self,
               other: _t.Union[Expansion, _Number]) -> _t.Union[_t.Any, bool]:
        return (not _are_components_lesser_than(self._components,
                                                other._components)
                if isinstance(other, Expansion)
                else
                (not _are_components_lesser_than_float(self._components,
                                                       other)
                 if isinstance(other, float)
                 else
                 (not _are_components_lesser_than_int(self._components,
                                                      other)
                  if isinstance(other, int)
                  else
                  (not _are_components_lesser_than_rational(
                          self._components, other)
                   if isinstance(other, _Rational)
                   else NotImplemented))))

    @_t.overload
    def __gt__(self, other: _t.Union[Expansion, _Number]) -> bool:
        ...

    @_t.overload
    def __gt__(self, other: _t.Any) -> _t.Any:
        ...

    def __gt__(self,
               other: _t.Union[Expansion, _Number]) -> _t.Union[_t.Any, bool]:
        return (_are_components_lesser_than(other._components,
                                            self._components)
                if isinstance(other, Expansion)
                else
                (_is_float_lesser_than_components(other, self._components)
                 if isinstance(other, float)
                 else
                 (_is_int_lesser_than_components(other, self._components)
                  if isinstance(other, int)
                  else
                  (_is_rational_lesser_than_components(other,
                                                       self._components)
                   if isinstance(other, _Rational)
                   else NotImplemented))))

    def __hash__(self) -> int:
        return hash(self._components)

    @_t.overload
    def __le__(self, other: _t.Union[Expansion, _Number]) -> bool:
        ...

    @_t.overload
    def __le__(self, other: _t.Any) -> _t.Any:
        ...

    def __le__(self,
               other: _t.Union[Expansion, _Number]) -> _t.Union[_t.Any, bool]:
        return (not _are_components_lesser_than(other._components,
                                                self._components)
                if isinstance(other, Expansion)
                else
                (not _is_float_lesser_than_components(other,
                                                      self._components)
                 if isinstance(other, float)
                 else
                 (not _is_int_lesser_than_components(other,
                                                     self._components)
                  if isinstance(other, int)
                  else
                  (not _is_rational_lesser_than_components(
                          other, self._components)
                   if isinstance(other, _Rational)
                   else NotImplemented))))

    @_t.overload
    def __lt__(self, other: _t.Union[Expansion, _Number]) -> bool:
        ...

    @_t.overload
    def __lt__(self, other: _t.Any) -> _t.Any:
        ...

    def __lt__(self,
               other: _t.Union[Expansion, _Number]) -> _t.Union[_t.Any, bool]:
        return (_are_components_lesser_than(self._components,
                                            other._components)
                if isinstance(other, Expansion)
                else
                (_are_components_lesser_than_float(self._components, other)
                 if isinstance(other, float)
                 else
                 (_are_components_lesser_than_int(self._components,
                                                  other)
                  if isinstance(other, int)
                  else
                  (_are_components_lesser_than_rational(self._components,
                                                        other)
                   if isinstance(other, _Rational)
                   else NotImplemented))))

    @_t.overload
    def __mul__(self, other: _t.Union[Expansion, _Number]) -> Expansion:
        ...

    @_t.overload
    def __mul__(self, other: _t.Any) -> _t.Any:
        ...

    def __mul__(self, other: _t.Any) -> _t.Any:
        return (Expansion(*_multiply_components(self._components,
                                                other._components))
                if isinstance(other, Expansion)
                else self.__rmul__(other))

    def __neg__(self) -> Expansion:
        return Expansion(*_negate_components(self._components),
                         _compress=False)

    def __pos__(self) -> Expansion:
        return self

    @_t.overload
    def __radd__(self, other: _t.Union[Expansion, _Number]) -> Expansion:
        ...

    @_t.overload
    def __radd__(self, other: _t.Any) -> _t.Any:
        ...

    def __radd__(
            self, other: _t.Union[Expansion, _Number, _t.Any]
    ) -> _t.Union[Expansion, _t.Any]:
        return (Expansion(*_add_float(self._components, other))
                if isinstance(other, float)
                else
                (Expansion(*_add_components(self._components,
                                            _int_to_components(other)))
                 if isinstance(other, int)
                 else
                 (Expansion(*_add_components(
                         self._components, _rational_to_components(other)))
                  if isinstance(other, _Rational)
                  else NotImplemented)))

    def __repr__(self) -> str:
        return (type(self).__qualname__
                + '({})'.format(', '.join(map(str, self._components))))

    @_t.overload
    def __rmul__(self, other: _Number) -> Expansion:
        ...

    @_t.overload
    def __rmul__(self, other: _t.Any) -> _t.Any:
        ...

    def __rmul__(
            self, other: _t.Union[_t.Any, _Number]
    ) -> _t.Union[_t.Any, Expansion]:
        return (Expansion(*_scale_components(self._components, other))
                if isinstance(other, float)
                else
                (Expansion(*_multiply_components(
                        self._components, _int_to_components(other)))
                 if isinstance(other, int)
                 else
                 (Expansion(*_multiply_components(
                         self._components, _rational_to_components(other)))
                  if isinstance(other, _Rational)
                  else NotImplemented)))

    @_t.overload
    def __round__(self, precision: None = ...) -> int:
        ...

    @_t.overload
    def __round__(self, precision: int) -> Expansion:
        ...

    def __round__(
            self, precision: _t.Optional[int] = None
    ) -> _t.Union[Expansion, int]:
        if precision is None:
            result = _components_to_integer(self._components)
            fractions = _components_to_fractions(self._components)
            fraction_sign = (
                -1
                if _are_components_lesser_than_float(fractions, 0)
                else 1
            )
            if _are_components_equal_to_float(fractions,
                                              0.5 * fraction_sign):
                sign = _to_sign(result)
                if sign != fraction_sign:
                    result -= sign
                if result & 1:
                    result += sign
            elif (_is_float_lesser_than_components(0.5, fractions)
                  or _are_components_lesser_than_float(fractions, -0.5)):
                result += fraction_sign
            return result
        else:
            return Expansion(*[round(component, precision)
                               for component in self._components])

    @_t.overload
    def __rsub__(self, other: _Number) -> Expansion:
        ...

    @_t.overload
    def __rsub__(self, other: _t.Any) -> _t.Any:
        ...

    def __rsub__(
            self, other: _t.Union[_t.Any, _Number]
    ) -> _t.Union[_t.Any, Expansion]:
        return (Expansion(*_subtract_from_double(other, self._components))
                if isinstance(other, float)
                else
                (Expansion(*_subtract_components(_int_to_components(other),
                                                 self._components))
                 if isinstance(other, int)
                 else
                 (Expansion(*_subtract_components(
                         _rational_to_components(other), self._components))
                  if isinstance(other, _Rational)
                  else NotImplemented)))

    @_t.overload
    def __rtruediv__(self, other: _Number) -> Expansion:
        ...

    @_t.overload
    def __rtruediv__(self, other: _t.Any) -> _t.Any:
        ...

    def __rtruediv__(
            self, other: _t.Union[_t.Any, _Number]
    ) -> _t.Union[_t.Any, Expansion]:
        return (Expansion(*_divide_components([other], self._components))
                if isinstance(other, float)
                else
                (Expansion(*_divide_components(_int_to_components(other),
                                               self._components))
                 if isinstance(other, int)
                 else
                 (Expansion(*_divide_components(
                         _rational_to_components(other), self._components))
                  if isinstance(other, _Rational)
                  else NotImplemented)))

    @_t.overload
    def __sub__(self, other: _t.Union[Expansion, _Number]) -> Expansion:
        ...

    @_t.overload
    def __sub__(self, other: _t.Any) -> _t.Any:
        ...

    def __sub__(self, other: _t.Union[Expansion, _Number]) -> Expansion:
        return (Expansion(*_subtract_components(self._components,
                                                other._components))
                if isinstance(other, Expansion)
                else
                (Expansion(*_subtract_double(self._components, other))
                 if isinstance(other, float)
                 else
                 (Expansion(*_subtract_components(
                         self._components, _int_to_components(other)))
                  if isinstance(other, int)
                  else (Expansion(*_subtract_components(
                         self._components, _rational_to_components(other)))
                        if isinstance(other, _Rational)
                        else NotImplemented))))

    @_t.overload
    def __truediv__(self, other: _t.Union[Expansion, _Number]) -> Expansion:
        ...

    @_t.overload
    def __truediv__(self, other: _t.Any) -> _t.Any:
        ...

    def __truediv__(self, other: _t.Union[Expansion, _Number]) -> Expansion:
        return (Expansion(*_divide_components(self._components,
                                              other._components))
                if isinstance(other, Expansion)
                else
                (Expansion(*_divide_components(self._components, [other]))
                 if isinstance(other, float)
                 else
                 (Expansion(*_divide_components(self._components,
                                                _int_to_components(other)))
                  if isinstance(other, int)
                  else
                  (Expansion(*_divide_components(
                          self._components, _rational_to_components(other)
                  ))
                   if isinstance(other, _Rational)
                   else NotImplemented))))

    def __trunc__(self) -> int:
        integer = _components_to_integer(self._components)
        integer_sign = _to_sign(integer)
        fraction_sign = _to_sign(_components_to_accumulated_fraction(
                self._components
        ))
        return (integer - integer_sign
                if (integer_sign and fraction_sign
                    and fraction_sign != integer_sign)
                else integer)


def _rational_to_components(value: _Rational) -> _t.List[float]:
    return (_int_to_components(value.numerator)
            if value.denominator == 1
            else _divide_components(_int_to_components(value.numerator),
                                    _int_to_components(value.denominator)))


def _are_components_equal_to_float(components: _t.Sequence[float],
                                   value: float) -> bool:
    return len(components) == 1 and components[0] == value


def _are_components_equal_to_int(components: _t.Sequence[float], value: int
                                 ) -> bool:
    return (_to_fraction(components[0]) == 0.
            and _components_to_integer(components) == value)


def _are_components_equal_to_rational(components: _t.Sequence[float],
                                      value: _Rational) -> bool:
    return components == tuple(_rational_to_components(value))


def _are_components_lesser_than_float(components: _t.Sequence[float],
                                      value: float) -> bool:
    return components[-1] < value or (len(components) > 1
                                      and components[-1] == value
                                      and components[-2] < 0.)


def _are_components_lesser_than_int(components: _t.Sequence[float],
                                    value: int) -> bool:
    components_integer = _components_to_integer(components)
    return (components_integer < value
            or (components_integer == value
                and _components_to_accumulated_fraction(components) < 0.))


def _are_components_lesser_than_rational(components: _t.Sequence[float],
                                         value: _Rational) -> bool:
    return _are_components_lesser_than(components,
                                       _rational_to_components(value))


def _divide_components(dividend: _t.Sequence[float],
                       divisor: _t.Sequence[float]) -> _t.List[float]:
    return _multiply_components(dividend, _invert_components(divisor))


def _int_to_components(value: int) -> _t.List[float]:
    if not value:
        return [0.]
    result = []
    while value:
        component = float(value)
        result.append(component)
        value -= int(component)
    return result[::-1]


def _invert_components(components: _t.Sequence[float]) -> _t.Sequence[float]:
    # based on Newton's method
    # https://en.wikipedia.org/wiki/Newton%27s_method
    # for f(x) = 1 / x
    first_approximation = 1. / components[-1]
    _, high = _split(first_approximation)
    if not _isfinite(high):
        raise OverflowError(f'Components {components} '
                            'are not finitely invertible.')
    result: _t.Sequence[float] = [first_approximation]
    negated_components = _negate_components(components)
    for _ in _repeat(None, 6 + _ceil_log2(len(components))):
        result = _multiply_components(
                _add_float(_multiply_components(result, negated_components),
                           2.),
                result
        )
    return result


def _ceil_log2(number: int) -> int:
    return number.bit_length() - (not (number & (number - 1)))


def _is_float_lesser_than_components(value: float,
                                     components: _t.Sequence[float]) -> bool:
    return value < components[-1] or (len(components) > 1
                                      and components[-1] == value
                                      and components[-2] > 0.)


def _is_int_lesser_than_components(value: int, components: _t.Sequence[float]
                                   ) -> bool:
    components_integer = _components_to_integer(components);
    return (value < components_integer
            or (value == components_integer
                and _components_to_accumulated_fraction(components) > 0.))


def _is_rational_lesser_than_components(value: _Rational,
                                        components: _t.Sequence[float]
                                        ) -> bool:
    return _are_components_lesser_than(_rational_to_components(value),
                                       components)


def _components_to_accumulated_fraction(components: _t.Sequence[float]
                                        ) -> float:
    result = 0.
    for component in components:
        component_fraction = _to_fraction(component)
        if not component_fraction:
            break
        result += component_fraction
    assert abs(result) < 1.0, components
    return result


def _components_to_fractions(components: _t.Sequence[float]
                             ) -> _t.Sequence[float]:
    result = []
    for component in components:
        component_fraction = _to_fraction(component)
        if not component_fraction:
            break
        result.append(component_fraction)
    return result or [0.]


def _components_to_integer(components: _t.Sequence[float]) -> int:
    result = 0
    for component in reversed(components):
        component_integer = int(component)
        if not component_integer:
            break
        result += component_integer
    return result


def _negate_components(components: _t.Sequence[float]) -> _t.Sequence[float]:
    return [-component for component in components]


def incircle_test(_point_x: float,
                  _point_y: float,
                  _first_x: float,
                  _first_y: float,
                  _second_x: float,
                  _second_y: float,
                  _third_x: float,
                  _third_y: float) -> int:
    """
    Computes location of point relative to a circle formed by three others
    given their coordinates.
    """
    return _to_sign(_incircle_determinant_estimation(
            _point_x, _point_y, _first_x, _first_y, _second_x, _second_y,
            _third_x, _third_y
    ))


def kind(_vertex_x: float,
         _vertex_y: float,
         _first_ray_point_x: float,
         _first_ray_point_y: float,
         _second_ray_point_x: float,
         _second_ray_point_y: float) -> int:
    """Computes kind of angle given its endpoints coordinates."""
    return _to_sign(_vectors_cross_product_estimation(
            _vertex_x, _vertex_y, _first_ray_point_x, _first_ray_point_y,
            -_vertex_y, _vertex_x, -_second_ray_point_y, _second_ray_point_x
    ))


def orientation(_start_x: float,
                _start_y: float,
                _end_x: float,
                _end_y: float,
                _point_x: float,
                _point_y: float) -> int:
    """
    Computes orientation of point relative to segment
    given their coordinates.
    """
    return _to_sign(_vectors_cross_product_estimation(_start_x, _start_y,
                                                      _end_x, _end_y,
                                                      _start_x, _start_y,
                                                      _point_x, _point_y))


def vectors_cross_product(_first_start_x: float,
                          _first_start_y: float,
                          _first_end_x: float,
                          _first_end_y: float,
                          _second_start_x: float,
                          _second_start_y: float,
                          _second_end_x: float,
                          _second_end_y: float) -> Expansion:
    """
    Computes cross product of two vectors
    given their endpoints coordinates.
    """
    return Expansion(
            *_vectors_cross_product(_first_start_x, _first_start_y,
                                    _first_end_x, _first_end_y,
                                    _second_start_x, _second_start_y,
                                    _second_end_x, _second_end_y),
            _compress=False
    )


def vectors_dot_product(_first_start_x: float,
                        _first_start_y: float,
                        _first_end_x: float,
                        _first_end_y: float,
                        _second_start_x: float,
                        _second_start_y: float,
                        _second_end_x: float,
                        _second_end_y: float) -> Expansion:
    """
    Computes dot product of two vectors given their endpoints coordinates.
    """
    return Expansion(
            *_vectors_cross_product(_first_start_x, _first_start_y,
                                    _first_end_x, _first_end_y,
                                    -_second_start_y, _second_start_x,
                                    -_second_end_y, _second_end_x),
            _compress=False
    )


_EPSILON = _float_info.epsilon / 2.


def _are_components_lesser_than(left: _t.Sequence[float],
                                right: _t.Sequence[float]) -> bool:
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


def _add_components(left: _t.Sequence[float],
                    right: _t.Sequence[float]) -> _t.List[float]:
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
        left_component, right_component = left[left_index], right[right_index]
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


def _add_float(components: _t.Sequence[float], value: float
               ) -> _t.Sequence[float]:
    result = []
    accumulator = value
    for component in components:
        tail, accumulator = _two_add(accumulator, component)
        if tail:
            result.append(tail)
    if accumulator or not result:
        result.append(accumulator)
    return result


def _compress_components(components: _t.List[float]) -> _t.List[float]:
    for _ in _repeat(None, len(components)):
        next_components = _compress_components_single(components)
        if next_components == components:
            break
        components = next_components
    return components


def _compress_components_single(
        components: _t.Sequence[float]
) -> _t.List[float]:
    bottom = len(components) - 1
    cursor = components[bottom]
    result = [0.0] * len(components)
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


def _fast_two_add(left: float, right: float) -> _t.Tuple[float, float]:
    head = left + right
    right_virtual = head - left
    tail = right - right_virtual
    return tail, head


def _multiply_components(left: _t.Sequence[float],
                         right: _t.Sequence[float]) -> _t.List[float]:
    return _reduce(_add_components,
                   [_scale_components(left, right_component)
                    for right_component in right])


def _scale_components(components: _t.Sequence[float],
                      scalar: float) -> _t.List[float]:
    components_iterator = iter(components)
    scalar_low, scalar_high = _split(scalar)
    tail, accumulator = _two_multiply_presplit(next(components_iterator),
                                               scalar, scalar_low, scalar_high)
    result = []
    if tail:
        result.append(tail)
    for component in components_iterator:
        product_tail, product = _two_multiply_presplit(component, scalar,
                                                       scalar_low, scalar_high)
        tail, interim = _two_add(accumulator, product_tail)
        if tail:
            result.append(tail)
        tail, accumulator = _fast_two_add(product, interim)
        if tail:
            result.append(tail)
    if accumulator or not result:
        result.append(accumulator)
    return result


def _split(
        value: float,
        *,
        splitter: float = float(1 + (1 << ((_float_info.mant_dig + 1) // 2)))
) -> _t.Tuple[float, float]:
    base = splitter * value
    high = base - (base - value)
    low = value - high
    return low, high


def _square(value: float) -> _t.Tuple[float, float]:
    head = value * value
    value_low, value_high = _split(value)
    first_error = head - value_high * value_high
    second_error = first_error - (value_high + value_high) * value_low
    tail = value_low * value_low - second_error
    return tail, head


def _subtract_components(minuend: _t.Sequence[float],
                         subtrahend: _t.Sequence[float]) -> _t.Sequence[float]:
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


def _subtract_double(minuend: _t.Sequence[float], subtrahend: float
                     ) -> _t.Sequence[float]:
    return _add_float(minuend, -subtrahend)


def _subtract_from_double(minuend: float, subtrahend: _t.Sequence[float]
                          ) -> _t.Sequence[float]:
    result = []
    accumulator = minuend
    for subtrahend_component in subtrahend:
        tail, accumulator = _two_add(accumulator, -subtrahend_component)
        if tail:
            result.append(tail)
    if accumulator or not result:
        result.append(accumulator)
    return result


def _to_sign(value: float) -> int:
    return 1 if value > 0. else (0 if not value else -1)


def _two_add(left: float, right: float) -> _t.Tuple[float, float]:
    head = left + right
    right_virtual = head - left
    left_virtual = head - right_virtual
    right_tail = right - right_virtual
    left_tail = left - left_virtual
    tail = left_tail + right_tail
    return tail, head


def _two_multiply(left: float, right: float) -> _t.Tuple[float, float]:
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
                           right_low: float,
                           right_high: float) -> _t.Tuple[float, float]:
    head = left * right
    left_low, left_high = _split(left)
    first_error = head - left_high * right_high
    second_error = first_error - left_low * right_high
    third_error = second_error - left_high * right_low
    tail = left_low * right_low - third_error
    return tail, head


def _two_one_add(left_tail: float,
                 left_head: float,
                 right: float) -> _t.Tuple[float, float, float]:
    second_tail, mid_head = _two_add(left_tail, right)
    first_tail, head = _two_add(left_head, mid_head)
    return second_tail, first_tail, head


def _two_one_subtract(left_tail: float,
                      left_head: float,
                      right: float) -> _t.Tuple[float, float, float]:
    second_tail, mid_head = _two_subtract(left_tail, right)
    first_tail, head = _two_add(left_head, mid_head)
    return second_tail, first_tail, head


def _two_subtract(left: float, right: float) -> _t.Tuple[float, float]:
    head = left - right
    return _two_subtract_tail(left, right, head), head


def _two_subtract_tail(left: float, right: float, head: float) -> float:
    right_virtual = left - head
    left_virtual = head + right_virtual
    right_tail = right_virtual - right
    left_tail = left - left_virtual
    return left_tail + right_tail


def _two_two_add(left_tail: float,
                 left_head: float,
                 right_tail: float,
                 right_head: float) -> _t.Tuple[float, float, float, float]:
    third_tail, mid_tail, mid_head = _two_one_add(left_tail, left_head,
                                                  right_tail)
    second_tail, first_tail, head = _two_one_add(mid_tail, mid_head,
                                                 right_head)
    return third_tail, second_tail, first_tail, head


def _two_two_subtract(
        left_tail: float,
        left_head: float,
        right_tail: float,
        right_head: float
) -> _t.Tuple[float, float, float, float]:
    third_tail, mid_tail, mid_head = _two_one_subtract(left_tail, left_head,
                                                       right_tail)
    second_tail, first_tail, head = _two_one_subtract(mid_tail, mid_head,
                                                      right_head)
    return third_tail, second_tail, first_tail, head


def _cross_product(first_dx: float,
                   first_dy: float,
                   second_dx: float,
                   second_dy: float) -> _t.Tuple[float, float, float, float]:
    first_dx_second_dy_tail, first_dx_second_dy_head = _two_multiply(first_dx,
                                                                     second_dy)
    second_dx_first_dy_tail, second_dx_first_dy_head = _two_multiply(second_dx,
                                                                     first_dy)
    return _two_two_subtract(first_dx_second_dy_tail, first_dx_second_dy_head,
                             second_dx_first_dy_tail, second_dx_first_dy_head)


def _scale_by_squared_length(components: _t.Sequence[float],
                             dx: float,
                             dy: float) -> _t.Sequence[float]:
    dx_components = _scale_components(components, dx)
    dx_squared_components = _scale_components(dx_components, dx)
    dy_components = _scale_components(components, dy)
    dy_squared_components = _scale_components(dy_components, dy)
    return _add_components(dx_squared_components, dy_squared_components)


def _squared_length(dx: float,
                    dy: float) -> _t.Tuple[float, float, float, float]:
    dx_squared_tail, dx_squared_head = _square(dx)
    dy_squared_tail, dy_squared_head = _square(dy)
    return _two_two_add(dx_squared_tail, dx_squared_head, dy_squared_tail,
                        dy_squared_head)


def _to_fraction(value: float) -> float:
    result, _ = _modf(value)
    return result


def _add_extras(
        final_components: _t.Sequence[float],
        first_dx: float,
        first_dx_tail: float,
        first_dy: float,
        first_dy_tail: float,
        second_dx: float,
        second_dx_tail: float,
        second_dy: float,
        second_dy_tail: float,
        third_dx: float,
        third_dx_tail: float,
        third_dy: float,
        third_dy_tail: float,
        second_third_cross_product: _t.Tuple[float, float, float, float],
        second_squared_length: _t.Tuple[float, float, float, float],
        third_squared_length: _t.Tuple[float, float, float, float]
) -> _t.Sequence[float]:
    if first_dx_tail:
        first_dx_tail_second_third_cross_product = _scale_components(
                second_third_cross_product, first_dx_tail
        )
        first_buffer_16 = _scale_components(
                first_dx_tail_second_third_cross_product, 2. * first_dx
        )
        first_dx_tail_third_squared_length = _scale_components(
                third_squared_length, first_dx_tail
        )
        second_buffer_16 = _scale_components(
                first_dx_tail_third_squared_length, second_dy)
        first_dx_tail_second_squared_length = _scale_components(
                second_squared_length, first_dx_tail
        )
        third_buffer_16 = _scale_components(
                first_dx_tail_second_squared_length, -third_dy
        )
        first_buffer_32 = _add_components(first_buffer_16,
                                          second_buffer_16)
        buffer_48 = _add_components(third_buffer_16, first_buffer_32)
        final_components = _add_components(final_components, buffer_48)
    first_dy_tail_second_third_cross_product: _t.Sequence[float]
    if first_dy_tail:
        first_dy_tail_second_third_cross_product = _scale_components(
                second_third_cross_product, first_dy_tail
        )
        first_buffer_16 = _scale_components(
                first_dy_tail_second_third_cross_product, 2. * first_dy
        )
        first_dy_tail_second_squared_length = _scale_components(
                second_squared_length, first_dy_tail
        )
        second_buffer_16 = _scale_components(
                first_dy_tail_second_squared_length, third_dx
        )
        first_dy_tail_third_squared_length = _scale_components(
                third_squared_length, first_dy_tail
        )
        third_buffer_16 = _scale_components(
                first_dy_tail_third_squared_length, -second_dx
        )
        first_buffer_32 = _add_components(first_buffer_16,
                                          second_buffer_16)
        buffer_48 = _add_components(third_buffer_16, first_buffer_32)
        final_components = _add_components(final_components, buffer_48)
    if first_dx_tail or first_dy_tail:
        second_third_cross_product_bodies: _t.Sequence[float]
        second_third_cross_product_tails: _t.Sequence[float]
        if (second_dx_tail or second_dy_tail or third_dx_tail
                or third_dy_tail):
            dx_tail_dy_head_tail, dx_tail_dy_head_head = _two_multiply(
                    second_dx_tail, third_dy
            )
            dx_head_dy_tail_tail, dx_head_dy_tail_head = _two_multiply(
                    second_dx, third_dy_tail
            )
            first_buffer_4 = _two_two_add(dx_tail_dy_head_tail,
                                          dx_tail_dy_head_head,
                                          dx_head_dy_tail_tail,
                                          dx_head_dy_tail_head)
            dx_tail_dy_head_tail, dx_tail_dy_head_head = _two_multiply(
                    third_dx_tail, -second_dy
            )
            dx_head_dy_tail_tail, dx_head_dy_tail_head = _two_multiply(
                    third_dx, -second_dy_tail
            )
            second_buffer_4 = _two_two_add(dx_tail_dy_head_tail,
                                           dx_tail_dy_head_head,
                                           dx_head_dy_tail_tail,
                                           dx_head_dy_tail_head)
            second_third_cross_product_bodies = _add_components(
                    first_buffer_4, second_buffer_4
            )
            dx_tail_dy_head_tail, dx_tail_dy_head_head = _two_multiply(
                    second_dx_tail, third_dy_tail
            )
            dx_head_dy_tail_tail, dx_head_dy_tail_head = _two_multiply(
                    third_dx_tail, second_dy_tail
            )
            second_third_cross_product_tails = _two_two_subtract(
                    dx_tail_dy_head_tail, dx_tail_dy_head_head,
                    dx_head_dy_tail_tail, dx_head_dy_tail_head
            )
        else:
            second_third_cross_product_bodies = [0.]
            second_third_cross_product_tails = [0.]
        if first_dx_tail:
            first_buffer_16 = _scale_components(
                    first_dx_tail_second_third_cross_product, first_dx_tail
            )
            first_dx_tail_second_third_cross_product_bodies = (
                _scale_components(second_third_cross_product_bodies,
                                  first_dx_tail)
            )
            first_buffer_32 = _scale_components(
                    first_dx_tail_second_third_cross_product_bodies,
                    2. * first_dx
            )
            buffer_48 = _add_components(first_buffer_16, first_buffer_32)
            final_components = _add_components(buffer_48, final_components)
            if second_dy_tail:
                buffer_8 = _scale_components(third_squared_length,
                                             first_dx_tail)
                first_buffer_16 = _scale_components(buffer_8,
                                                    second_dy_tail)
                final_components = _add_components(final_components,
                                                   first_buffer_16)
            if third_dy_tail:
                buffer_8 = _scale_components(second_squared_length,
                                             -first_dx_tail)
                first_buffer_16 = _scale_components(buffer_8,
                                                    third_dy_tail)
                final_components = _add_components(final_components,
                                                   first_buffer_16)
            first_buffer_32 = _scale_components(
                    first_dx_tail_second_third_cross_product_bodies,
                    first_dx_tail
            )
            first_dx_tail_second_third_cross_product_tails = (
                _scale_components(second_third_cross_product_tails,
                                  first_dx_tail)
            )
            first_buffer_16 = _scale_components(
                    first_dx_tail_second_third_cross_product_tails,
                    2. * first_dx
            )
            second_buffer_16 = _scale_components(
                    first_dx_tail_second_third_cross_product_tails,
                    first_dx_tail
            )
            second_buffer_32 = _add_components(first_buffer_16,
                                               second_buffer_16)
            buffer_64 = _add_components(first_buffer_32, second_buffer_32)
            final_components = _add_components(final_components, buffer_64)
        if first_dy_tail:
            first_buffer_16 = _scale_components(
                    first_dy_tail_second_third_cross_product, first_dy_tail
            )
            first_dy_tail_second_third_cross_product_bodies = (
                _scale_components(second_third_cross_product_bodies,
                                  first_dy_tail)
            )
            first_buffer_32 = _scale_components(
                    first_dy_tail_second_third_cross_product_bodies,
                    2. * first_dy
            )
            buffer_48 = _add_components(first_buffer_16, first_buffer_32)
            final_components = _add_components(final_components, buffer_48)
            first_buffer_32 = _scale_components(
                    first_dy_tail_second_third_cross_product_bodies,
                    first_dy_tail
            )
            first_dy_tail_second_third_cross_product_tails = (
                _scale_components(second_third_cross_product_tails,
                                  first_dy_tail)
            )
            first_buffer_16 = _scale_components(
                    first_dy_tail_second_third_cross_product_tails,
                    2. * first_dy
            )
            second_buffer_16 = _scale_components(
                    first_dy_tail_second_third_cross_product_tails,
                    first_dy_tail
            )
            second_buffer_32 = _add_components(first_buffer_16,
                                               second_buffer_16)
            buffer_64 = _add_components(first_buffer_32, second_buffer_32)
            final_components = _add_components(final_components, buffer_64)
    return final_components


def _adaptive_incircle_determinant_estimation(
        point_x: float,
        point_y: float,
        first_x: float,
        first_y: float,
        second_x: float,
        second_y: float,
        third_x: float,
        third_y: float,
        upper_bound: float,
        second_upper_bound_coefficient: float
        = (44. + 576. * _EPSILON) * _EPSILON * _EPSILON,
        result_coefficient: float = (3. + 8. * _EPSILON) * _EPSILON
) -> float:
    first_dx = first_x - point_x
    second_dx = second_x - point_x
    third_dx = third_x - point_x
    first_dy = first_y - point_y
    second_dy = second_y - point_y
    third_dy = third_y - point_y
    first_second_cross_product = _cross_product(first_dx, first_dy,
                                                second_dx, second_dy)
    second_third_cross_product = _cross_product(second_dx, second_dy,
                                                third_dx, third_dy)
    third_first_cross_product = _cross_product(third_dx, third_dy,
                                               first_dx, first_dy)
    first_components = _scale_by_squared_length(second_third_cross_product,
                                                first_dx, first_dy)
    second_components = _scale_by_squared_length(third_first_cross_product,
                                                 second_dx, second_dy)
    third_components = _scale_by_squared_length(first_second_cross_product,
                                                third_dx, third_dy)
    first_second_sum_components = _add_components(first_components,
                                                  second_components)
    first_buffer = _add_components(first_second_sum_components,
                                   third_components)
    result = sum(first_buffer)
    first_upper_bound_coefficient = (4. + 48. * _EPSILON) * _EPSILON
    threshold = first_upper_bound_coefficient * upper_bound
    if (result >= threshold) or (-result >= threshold):
        return result
    first_dx_tail = _two_subtract_tail(first_x, point_x, first_dx)
    first_dy_tail = _two_subtract_tail(first_y, point_y, first_dy)
    second_dx_tail = _two_subtract_tail(second_x, point_x, second_dx)
    second_dy_tail = _two_subtract_tail(second_y, point_y, second_dy)
    third_dx_tail = _two_subtract_tail(third_x, point_x, third_dx)
    third_dy_tail = _two_subtract_tail(third_y, point_y, third_dy)
    if (not first_dx_tail and not second_dx_tail and not third_dx_tail
            and not first_dy_tail and not second_dy_tail
            and not third_dy_tail):
        return result
    threshold = (second_upper_bound_coefficient * upper_bound
                 + result_coefficient * abs(result))
    result += (((first_dx * first_dx + first_dy * first_dy)
                * ((second_dx * third_dy_tail + third_dy * second_dx_tail)
                   - (second_dy * third_dx_tail
                      + third_dx * second_dy_tail))
                + 2. * (first_dx * first_dx_tail
                        + first_dy * first_dy_tail)
                * (second_dx * third_dy - second_dy * third_dx))
               + ((second_dx * second_dx + second_dy * second_dy)
                  * ((third_dx * first_dy_tail + first_dy * third_dx_tail)
                     - (third_dy * first_dx_tail
                        + first_dx * third_dy_tail))
                  + 2. * (second_dx * second_dx_tail
                          + second_dy * second_dy_tail)
                  * (third_dx * first_dy - third_dy * first_dx))
               + ((third_dx * third_dx + third_dy * third_dy)
                  * ((first_dx * second_dy_tail
                      + second_dy * first_dx_tail)
                     - (first_dy * second_dx_tail
                        + second_dx * first_dy_tail))
                  + 2. * (third_dx * third_dx_tail
                          + third_dy * third_dy_tail)
                  * (first_dx * second_dy - first_dy * second_dx)))
    if (result >= threshold) or (-result >= threshold):
        return result
    first_squared_length = (_squared_length(first_dx, first_dy)
                            if (second_dx_tail or second_dy_tail
                                or third_dx_tail or third_dy_tail)
                            else (0, 0, 0, 0))
    second_squared_length = (_squared_length(second_dx, second_dy)
                             if (third_dx_tail or third_dy_tail
                                 or first_dx_tail or first_dy_tail)
                             else (0, 0, 0, 0))
    third_squared_length = (_squared_length(third_dx, third_dy)
                            if (first_dx_tail or first_dy_tail
                                or second_dx_tail or second_dy_tail)
                            else (0, 0, 0, 0))
    final_components = _add_extras(
            first_buffer, first_dx, first_dx_tail, first_dy, first_dy_tail,
            second_dx, second_dx_tail, second_dy, second_dy_tail, third_dx,
            third_dx_tail, third_dy, third_dy_tail, second_third_cross_product,
            second_squared_length, third_squared_length
    )
    final_components = _add_extras(
            final_components, second_dx, second_dx_tail, second_dy,
            second_dy_tail, third_dx, third_dx_tail, third_dy, third_dy_tail,
            first_dx, first_dx_tail, first_dy, first_dy_tail,
            third_first_cross_product, third_squared_length,
            first_squared_length
    )
    final_components = _add_extras(
            final_components, third_dx, third_dx_tail, third_dy, third_dy_tail,
            first_dx, first_dx_tail, first_dy, first_dy_tail, second_dx,
            second_dx_tail, second_dy, second_dy_tail,
            first_second_cross_product, first_squared_length,
            second_squared_length
    )
    return final_components[-1]


def _incircle_determinant_estimation(
        point_x: float,
        point_y: float,
        first_x: float,
        first_y: float,
        second_x: float,
        second_y: float,
        third_x: float,
        third_y: float,
        upper_bound_coefficient: float = (10. + 96. * _EPSILON) * _EPSILON
) -> float:
    first_dx = first_x - point_x
    second_dx = second_x - point_x
    third_dx = third_x - point_x
    first_dy = first_y - point_y
    second_dy = second_y - point_y
    third_dy = third_y - point_y
    second_dx_third_dy = second_dx * third_dy
    third_dx_second_dy = third_dx * second_dy
    first_squared_distance = first_dx * first_dx + first_dy * first_dy
    third_dx_first_dy = third_dx * first_dy
    first_dx_third_dy = first_dx * third_dy
    second_squared_distance = second_dx * second_dx + second_dy * second_dy
    first_dx_second_dy = first_dx * second_dy
    second_dx_first_dy = second_dx * first_dy
    third_squared_distance = third_dx * third_dx + third_dy * third_dy
    result = (first_squared_distance * (second_dx_third_dy
                                        - third_dx_second_dy) +
              second_squared_distance * (third_dx_first_dy
                                         - first_dx_third_dy) +
              third_squared_distance * (first_dx_second_dy
                                        - second_dx_first_dy))
    upper_bound = ((abs(second_dx_third_dy) + abs(third_dx_second_dy))
                   * first_squared_distance
                   + (abs(third_dx_first_dy) + abs(first_dx_third_dy))
                   * second_squared_distance
                   + (abs(first_dx_second_dy) + abs(second_dx_first_dy))
                   * third_squared_distance)
    threshold = upper_bound_coefficient * upper_bound
    return (
        result
        if (result > threshold) or (-result > threshold)
        else _adaptive_incircle_determinant_estimation(
                point_x, point_y, first_x, first_y, second_x, second_y,
                third_x, third_y, upper_bound
        )
    )


def _vectors_cross_product_estimation(
        first_start_x: float,
        first_start_y: float,
        first_end_x: float,
        first_end_y: float,
        second_start_x: float,
        second_start_y: float,
        second_end_x: float,
        second_end_y: float,
        upper_bound_coefficient: float = (3. + 16. * _EPSILON) * _EPSILON
) -> float:
    minuend = ((first_end_x - first_start_x)
               * (second_end_y - second_start_y))
    subtrahend = ((first_end_y - first_start_y)
                  * (second_end_x - second_start_x))
    result = minuend - subtrahend
    if minuend > 0.:
        if subtrahend <= 0.:
            return result
        else:
            upper_bound = minuend + subtrahend
    elif minuend < 0.:
        if subtrahend >= 0.:
            return result
        else:
            upper_bound = -minuend - subtrahend
    else:
        return result
    threshold = upper_bound_coefficient * upper_bound
    return (result
            if (result >= threshold) or (-result >= threshold)
            else
            _adaptive_vectors_cross_product_estimation(
                    first_start_x, first_start_y, first_end_x, first_end_y,
                    second_start_x, second_start_y, second_end_x,
                    second_end_y, upper_bound))


def _adaptive_vectors_cross_product_estimation(
        first_start_x: float,
        first_start_y: float,
        first_end_x: float,
        first_end_y: float,
        second_start_x: float,
        second_start_y: float,
        second_end_x: float,
        second_end_y: float,
        upper_bound: float,
        first_upper_bound_coefficient: float
        = (2. + 12. * _EPSILON) * _EPSILON,
        second_upper_bound_coefficient: float
        = (9. + 64. * _EPSILON) * _EPSILON * _EPSILON,
        estimation_coefficient: float = (3. + 8. * _EPSILON) * _EPSILON
) -> float:
    minuend_x = first_end_x - first_start_x
    minuend_y = first_end_y - first_start_y
    subtrahend_x = second_end_x - second_start_x
    subtrahend_y = second_end_y - second_start_y
    minuend_tail, minuend = _two_multiply(minuend_x, subtrahend_y)
    subtrahend_tail, subtrahend = _two_multiply(minuend_y, subtrahend_x)
    first_components = _two_two_subtract(minuend_tail, minuend,
                                         subtrahend_tail, subtrahend)
    result = sum(first_components)
    threshold = first_upper_bound_coefficient * upper_bound
    if (result >= threshold) or (-result >= threshold):
        return result
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
        return result
    threshold = (second_upper_bound_coefficient * upper_bound
                 + estimation_coefficient * abs(result))
    result += ((minuend_x * subtrahend_y_tail
                + subtrahend_y * minuend_x_tail)
               - (minuend_y * subtrahend_x_tail
                  + subtrahend_x * minuend_y_tail))
    if (result >= threshold) or (-result >= threshold):
        return result
    minuend_x_subtrahend_y_tail, minuend_x_subtrahend_y = _two_multiply(
            minuend_x_tail, subtrahend_y)
    minuend_y_subtrahend_x_tail, minuend_y_subtrahend_x = _two_multiply(
            minuend_y_tail, subtrahend_x)
    extra_components = _two_two_subtract(minuend_x_subtrahend_y_tail,
                                         minuend_x_subtrahend_y,
                                         minuend_y_subtrahend_x_tail,
                                         minuend_y_subtrahend_x)
    second_components = _add_components(first_components, extra_components)
    minuend_x_subtrahend_y_tail, minuend_x_subtrahend_y = _two_multiply(
            minuend_x, subtrahend_y_tail)
    minuend_y_subtrahend_x_tail, minuend_y_subtrahend_x = _two_multiply(
            minuend_y, subtrahend_x_tail)
    extra_components = _two_two_subtract(minuend_x_subtrahend_y_tail,
                                         minuend_x_subtrahend_y,
                                         minuend_y_subtrahend_x_tail,
                                         minuend_y_subtrahend_x)
    third_components = _add_components(second_components, extra_components)
    minuend_x_subtrahend_y_tail, minuend_x_subtrahend_y = _two_multiply(
            minuend_x_tail, subtrahend_y_tail)
    minuend_y_subtrahend_x_tail, minuend_y_subtrahend_x = _two_multiply(
            minuend_y_tail, subtrahend_x_tail)
    extra_components = _two_two_subtract(minuend_x_subtrahend_y_tail,
                                         minuend_x_subtrahend_y,
                                         minuend_y_subtrahend_x_tail,
                                         minuend_y_subtrahend_x)
    return _add_components(third_components, extra_components)[-1]


def _vectors_cross_product(first_start_x: float,
                           first_start_y: float,
                           first_end_x: float,
                           first_end_y: float,
                           second_start_x: float,
                           second_start_y: float,
                           second_end_x: float,
                           second_end_y: float,
                           upper_bound_coefficient: float
                           = (3. + 16. * _EPSILON) * _EPSILON
                           ) -> _t.Sequence[float]:
    minuend = ((first_end_x - first_start_x)
               * (second_end_y - second_start_y))
    subtrahend = ((first_end_y - first_start_y)
                  * (second_end_x - second_start_x))
    estimation = minuend - subtrahend
    if minuend > 0.:
        if subtrahend <= 0.:
            return [estimation]
        else:
            upper_bound = minuend + subtrahend
    elif minuend < 0.:
        if subtrahend >= 0.:
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
                                    = (2. + 12. * _EPSILON) * _EPSILON,
                                    second_upper_bound_coefficient: float
                                    = ((9. + 64. * _EPSILON) * _EPSILON
                                       * _EPSILON),
                                    estimation_coefficient: float
                                    = (3. + 8. * _EPSILON) * _EPSILON
                                    ) -> _t.Sequence[float]:
    minuend_x = first_end_x - first_start_x
    minuend_y = first_end_y - first_start_y
    subtrahend_x = second_end_x - second_start_x
    subtrahend_y = second_end_y - second_start_y
    minuend_tail, minuend = _two_multiply(minuend_x, subtrahend_y)
    subtrahend_tail, subtrahend = _two_multiply(minuend_y, subtrahend_x)
    first_components = _two_two_subtract(minuend_tail, minuend,
                                         subtrahend_tail, subtrahend)
    estimation = sum(first_components)
    threshold = first_upper_bound_coefficient * upper_bound
    if (estimation >= threshold) or (-estimation >= threshold):
        return _compress_components_single(first_components)
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
        return _compress_components_single(first_components)
    threshold = (second_upper_bound_coefficient * upper_bound
                 + estimation_coefficient * abs(estimation))
    extra = ((minuend_x * subtrahend_y_tail
              + subtrahend_y * minuend_x_tail)
             - (minuend_y * subtrahend_x_tail
                + subtrahend_x * minuend_y_tail))
    estimation += extra
    if (estimation >= threshold) or (-estimation >= threshold):
        return _add_float(first_components, extra)
    minuend_x_subtrahend_y_tail, minuend_x_subtrahend_y = _two_multiply(
            minuend_x_tail, subtrahend_y)
    minuend_y_subtrahend_x_tail, minuend_y_subtrahend_x = _two_multiply(
            minuend_y_tail, subtrahend_x)
    extra_components = _two_two_subtract(minuend_x_subtrahend_y_tail,
                                         minuend_x_subtrahend_y,
                                         minuend_y_subtrahend_x_tail,
                                         minuend_y_subtrahend_x)
    second_components = _add_components(first_components, extra_components)
    minuend_x_subtrahend_y_tail, minuend_x_subtrahend_y = _two_multiply(
            minuend_x, subtrahend_y_tail)
    minuend_y_subtrahend_x_tail, minuend_y_subtrahend_x = _two_multiply(
            minuend_y, subtrahend_x_tail)
    extra_components = _two_two_subtract(minuend_x_subtrahend_y_tail,
                                         minuend_x_subtrahend_y,
                                         minuend_y_subtrahend_x_tail,
                                         minuend_y_subtrahend_x)
    third_components = _add_components(second_components, extra_components)
    minuend_x_subtrahend_y_tail, minuend_x_subtrahend_y = _two_multiply(
            minuend_x_tail, subtrahend_y_tail)
    minuend_y_subtrahend_x_tail, minuend_y_subtrahend_x = _two_multiply(
            minuend_y_tail, subtrahend_x_tail)
    extra_components = _two_two_subtract(minuend_x_subtrahend_y_tail,
                                         minuend_x_subtrahend_y,
                                         minuend_y_subtrahend_x_tail,
                                         minuend_y_subtrahend_x)
    return _add_components(third_components, extra_components)
