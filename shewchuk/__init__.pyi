from numbers import (Integral as _Integral,
                     Real as _Real)
from typing import (Any as _Any,
                    Union as _Union,
                    overload as _overload)

__version__: str = ...


class Expansion:
    @property
    def real(self) -> _Real:
        ...

    @property
    def imag(self) -> _Real:
        ...

    @_overload
    def __new__(cls, value: _Union['Expansion', _Integral, float]
                ) -> 'Expansion':
        ...

    @_overload
    def __new__(cls, *args: float) -> 'Expansion':
        ...

    def __abs__(self) -> 'Expansion':
        ...

    def __add__(self, other: _Union['Expansion', _Integral, float]
                ) -> 'Expansion':
        ...

    def __bool__(self) -> bool:
        ...

    def __ceil__(self) -> int:
        ...

    @_overload
    def __eq__(self, other: _Union['Expansion', _Integral, float]) -> bool:
        ...

    @_overload
    def __eq__(self, other: _Any) -> _Any:
        ...

    def __float__(self) -> float:
        ...

    def __floor__(self) -> int:
        ...

    def __ge__(self, other: _Union['Expansion', _Integral, float]) -> bool:
        ...

    def __gt__(self, other: _Union['Expansion', _Integral, float]) -> bool:
        ...

    def __hash__(self) -> int:
        ...

    def __le__(self, other: _Union['Expansion', _Integral, float]) -> bool:
        ...

    def __lt__(self, other: _Union['Expansion', _Integral, float]) -> bool:
        ...

    def __mul__(self, other: _Union['Expansion', _Integral, float]
                ) -> 'Expansion':
        ...

    def __neg__(self) -> 'Expansion':
        ...

    def __pos__(self) -> 'Expansion':
        ...

    def __radd__(self, other: _Union[_Integral, float]) -> 'Expansion':
        ...

    def __rmul__(self, other: _Union[_Integral, float]) -> 'Expansion':
        ...

    @_overload
    def __round__(self, precision: None = ...) -> int:
        ...

    @_overload
    def __round__(self, precision: int) -> 'Expansion':
        ...

    def __rsub__(self, other: _Union[_Integral, float]) -> 'Expansion':
        ...

    def __rtruediv__(self, other: _Union[_Integral, float]) -> 'Expansion':
        ...

    def __sub__(self, other: _Union['Expansion', _Integral, float]
                ) -> 'Expansion':
        ...

    def __truediv__(self, other: _Union['Expansion', _Integral, float]
                    ) -> 'Expansion':
        ...

    def __trunc__(self) -> int:
        ...


def incircle_test(_point_x: float,
                  _point_y: float,
                  _first_x: float,
                  _first_y: float,
                  _second_x: float,
                  _second_y: float,
                  _third_x: float,
                  _third_y: float) -> int:
    ...


def kind(_vertex_x: float,
         _vertex_y: float,
         _first_ray_point_x: float,
         _first_ray_point_y: float,
         _second_ray_point_x: float,
         _second_ray_point_y: float) -> int:
    ...


def orientation(_start_x: float,
                _start_y: float,
                _end_x: float,
                _end_y: float,
                _point_x: float,
                _point_y: float) -> int:
    ...


def vectors_cross_product(_first_start_x: float,
                          _first_start_y: float,
                          _first_end_x: float,
                          _first_end_y: float,
                          _second_start_x: float,
                          _second_start_y: float,
                          _second_end_x: float,
                          _second_end_y: float) -> Expansion:
    ...


def vectors_dot_product(_first_start_x: float,
                        _first_start_y: float,
                        _first_end_x: float,
                        _first_end_y: float,
                        _second_start_x: float,
                        _second_start_y: float,
                        _second_end_x: float,
                        _second_end_y: float) -> Expansion:
    ...
