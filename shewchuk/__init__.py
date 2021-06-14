"""Robust floating point operations."""

__version__ = '0.0.0'

try:
    from _shewchuk import (Expansion,
                           Quadruple)
except ImportError:
    from typing import (Sequence as _Sequence,
                        Tuple as _Tuple)


    class Expansion:
        __slots__ = '_components',

        def __new__(cls, *components: float) -> 'Expansion':
            self = super().__new__(cls)
            if len(components) > 1:
                components = tuple(_compress_components(components))
            self._components = components or (0.,)
            return self

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


    def _two_add(left: float, right: float) -> _Tuple[float, float]:
        head = left + right
        right_virtual = head - left
        left_virtual = head - right_virtual
        right_tail = right - right_virtual
        left_tail = left - left_virtual
        tail = left_tail + right_tail
        return head, tail
