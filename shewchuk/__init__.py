"""Robust floating point operations."""

__version__ = '0.0.0'

try:
    from _shewchuk import (Expansion,
                           Quadruple)
except ImportError:
    from numbers import Real as _Real
    from typing import (Sequence as _Sequence,
                        Tuple as _Tuple,
                        Union as _Union)


    class Expansion:
        __slots__ = '_components',

        def __new__(cls, *components: float) -> 'Expansion':
            self = super().__new__(cls)
            if len(components) > 1:
                components = tuple(_compress_components(components))
            self._components = components or (0.,)
            return self

        def __abs__(self) -> 'Expansion':
            return +self if self._components[-1] > 0. else -self

        def __add__(self, other: _Union[_Real, 'Expansion']) -> 'Expansion':
            return (Expansion(*_add_components_eliminating_zeros(
                    self._components, other._components))
                    if isinstance(other, Expansion)
                    else self.__radd__(other))

        def __neg__(self) -> 'Expansion':
            return Expansion(*[-component for component in self._components])

        def __pos__(self) -> 'Expansion':
            return self

        def __radd__(self, other: _Real) -> 'Expansion':
            return (Expansion(*_add_component_eliminating_zeros(
                    self._components, float(other)))
                    if isinstance(other, _Real)
                    else NotImplemented)

        def __repr__(self) -> str:
            return (type(self).__qualname__
                    + '({})'.format(', '.join(map(str, self._components))))


    class Quadruple:
        __slots__ = '_head', '_tail'

        def __new__(cls, _head: float = 0., _tail: float = 0.) -> 'Quadruple':
            self = super().__new__(cls)
            self._head, self._tail = _two_add(_head, _tail)
            return self

        def __repr__(self) -> str:
            return (type(self).__qualname__
                    + (('({}, {})'
                        if self._tail
                        else '({})').format(self._head, self._tail)))


    def _add_component_eliminating_zeros(left: _Sequence[float],
                                         right: float) -> _Sequence[float]:
        result = []
        accumulator = right
        for left_component in left:
            accumulator, tail = _two_add(accumulator, left_component)
            if tail:
                result.append(tail)
        if accumulator or not result:
            result.append(accumulator)
        return result


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
            left_component, right_component = left[left_index], right[
                right_index]
            if ((right_component > left_component)
                    is (right_component > -left_component)):
                accumulator, tail = _fast_two_add(left_component, accumulator)
                left_index += 1
            else:
                accumulator, tail = _fast_two_add(right_component, accumulator)
                right_index += 1
            if tail:
                result.append(tail)
            while (left_index < left_length) and (right_index < right_length):
                left_component, right_component = (left[left_index],
                                                   right[right_index])
                if ((right_component > left_component)
                        is (right_component > -left_component)):
                    accumulator, tail = _two_add(accumulator, left_component)
                    left_index += 1
                else:
                    accumulator, tail = _two_add(accumulator, right_component)
                    right_index += 1
                if tail:
                    result.append(tail)
        for left_index in range(left_index, left_length):
            left_component = left[left_index]
            accumulator, tail = _two_add(accumulator, left_component)
            if tail:
                result.append(tail)
        for right_index in range(right_index, right_length):
            right_component = right[right_index]
            accumulator, tail = _two_add(accumulator, right_component)
            if tail:
                result.append(tail)
        if accumulator or not result:
            result.append(accumulator)
        return result


    def _compress_components(components: _Sequence[float]) -> _Sequence[float]:
        bottom = len(components) - 1
        cursor = components[bottom]
        result = [None] * len(components)
        for index in range(bottom - 1, -1, -1):
            head, tail = _two_add(cursor, components[index])
            if tail:
                result[bottom] = head
                bottom -= 1
                cursor = tail
            else:
                cursor = head
        top = 0
        for index in range(bottom + 1, len(result)):
            head, tail = _two_add(result[index], cursor)
            if tail:
                result[top] = tail
                top += 1
            cursor = head
        result[top] = cursor
        return result[:top + 1]


    def _fast_two_add(left: float, right: float) -> _Tuple[float, float]:
        head = left + right
        right_virtual = head - left
        tail = right - right_virtual
        return head, tail


    def _two_add(left: float, right: float) -> _Tuple[float, float]:
        head = left + right
        right_virtual = head - left
        left_virtual = head - right_virtual
        right_tail = right - right_virtual
        left_tail = left - left_virtual
        tail = left_tail + right_tail
        return head, tail
