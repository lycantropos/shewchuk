"""Robust floating point operations."""

__version__ = '0.0.0'

try:
    from _shewchuk import Quadruple
except ImportError:
    from typing import Tuple


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


    def _two_add(left: float, right: float) -> Tuple[float, float]:
        head = left + right
        right_virtual = head - left
        left_virtual = head - right_virtual
        right_tail = right - right_virtual
        left_tail = left - left_virtual
        tail = left_tail + right_tail
        return head, tail
