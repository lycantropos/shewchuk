"""Robust floating point operations."""

from __future__ import annotations

from typing import (
    Any as _Any,
    TYPE_CHECKING,
    final as _final,
    overload as _overload,
)

__version__ = '6.10.0'

if TYPE_CHECKING:
    from typing_extensions import Self as _Self

    from ._shewchuk import Number

    @_final
    class Expansion:
        @property
        def real(self) -> _Self: ...

        @property
        def imag(self) -> Number: ...

        def as_integer_ratio(self) -> tuple[int, int]: ...

        def __new__(
            cls, _argument: _Self | Number = 0.0, *args: float
        ) -> _Self: ...

        def __abs__(self) -> _Self: ...

        @_overload
        def __add__(self, other: _Self | Number) -> _Self: ...

        @_overload
        def __add__(self, other: _Any) -> _Any: ...

        def __add__(self, other: _Any) -> _Any: ...

        def __bool__(self) -> bool: ...

        def __ceil__(self) -> int: ...

        @_overload
        def __eq__(self, other: _Self | Number) -> bool: ...

        @_overload
        def __eq__(self, other: _Any) -> _Any: ...

        def __eq__(self, other: _Any) -> _Any: ...

        def __float__(self) -> float: ...

        def __floor__(self) -> int: ...

        @_overload
        def __ge__(self, other: _Self | Number) -> bool: ...

        @_overload
        def __ge__(self, other: _Any) -> _Any: ...

        def __ge__(self, other: _Any) -> _Any: ...

        @_overload
        def __gt__(self, other: _Self | Number) -> bool: ...

        @_overload
        def __gt__(self, other: _Any) -> _Any: ...

        def __gt__(self, other: _Any) -> _Any: ...

        def __hash__(self) -> int: ...

        @_overload
        def __le__(self, other: _Self | Number) -> bool: ...

        @_overload
        def __le__(self, other: _Any) -> _Any: ...

        def __le__(self, other: _Self | Number) -> _Any | bool: ...

        @_overload
        def __lt__(self, other: _Self | Number) -> bool: ...

        @_overload
        def __lt__(self, other: _Any) -> _Any: ...

        def __lt__(self, other: _Self | Number) -> _Any | bool: ...

        @_overload
        def __mul__(self, other: _Self | Number) -> _Self: ...

        @_overload
        def __mul__(self, other: _Any) -> _Any: ...

        def __mul__(self, other: _Any) -> _Any: ...

        def __neg__(self) -> _Self: ...

        def __pos__(self) -> _Self: ...

        def __radd__(self, other: _Any) -> _Any: ...

        def __repr__(self) -> str: ...

        def __rmul__(self, other: _Any) -> _Any: ...

        @_overload
        def __round__(self, precision: None = ...) -> int: ...

        @_overload
        def __round__(self, precision: int) -> _Self: ...

        def __round__(self, precision: int | None = None) -> _Self | int: ...

        def __rsub__(self, other: _Any) -> _Any: ...

        def __rtruediv__(self, other: _Any) -> _Any: ...

        @_overload
        def __sub__(self, other: _Self | Number) -> _Self: ...

        @_overload
        def __sub__(self, other: _Any) -> _Any: ...

        def __sub__(self, other: _Any) -> _Any: ...

        @_overload
        def __truediv__(self, other: _Self | Number) -> _Self: ...

        @_overload
        def __truediv__(self, other: _Any) -> _Any: ...

        def __truediv__(self, other: _Any) -> _Any: ...

        def __trunc__(self) -> int: ...

    def incircle_test(
        point_x: float,
        point_y: float,
        first_x: float,
        first_y: float,
        second_x: float,
        second_y: float,
        third_x: float,
        third_y: float,
        /,
    ) -> int: ...

    def kind(
        vertex_x: float,
        vertex_y: float,
        first_ray_point_x: float,
        first_ray_point_y: float,
        second_ray_point_x: float,
        second_ray_point_y: float,
        /,
    ) -> int: ...

    def orientation(
        start_x: float,
        start_y: float,
        end_x: float,
        end_y: float,
        point_x: float,
        point_y: float,
        /,
    ) -> int: ...

    def vectors_cross_product(
        first_start_x: float,
        first_start_y: float,
        first_end_x: float,
        first_end_y: float,
        second_start_x: float,
        second_start_y: float,
        second_end_x: float,
        second_end_y: float,
        /,
    ) -> Expansion: ...

    def vectors_dot_product(
        first_start_x: float,
        first_start_y: float,
        first_end_x: float,
        first_end_y: float,
        second_start_x: float,
        second_start_y: float,
        second_end_x: float,
        second_end_y: float,
        /,
    ) -> Expansion: ...

else:
    try:
        from . import _cshewchuk
    except ImportError:
        from ._shewchuk import (
            Expansion,
            incircle_test,
            kind,
            orientation,
            vectors_cross_product,
            vectors_dot_product,
        )
    else:
        Expansion = _cshewchuk.Expansion
        incircle_test = _cshewchuk.incircle_test
        kind = _cshewchuk.kind
        orientation = _cshewchuk.orientation
        vectors_cross_product = _cshewchuk.vectors_cross_product
        vectors_dot_product = _cshewchuk.vectors_dot_product
