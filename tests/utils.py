from __future__ import annotations

import pickle
import platform
from collections.abc import Iterable
from fractions import Fraction
from functools import partial
from numbers import Rational
from typing import Any, Callable, TypeVar, Union

import pytest
from hypothesis.strategies import SearchStrategy as Strategy

import shewchuk
from shewchuk import Expansion

LeftOperand = Union[Rational, float]
RightOperand = Union[Expansion, Rational, float]
Strategy = Strategy
Domain = TypeVar('Domain')
Range = TypeVar('Range')

MAX_VALUE = 10**50


def apply(function: Callable[..., Range], args: tuple[Domain, ...]) -> Range:
    return function(*args)


def equivalence(
    left_statement: bool,  # noqa: FBT001
    right_statement: bool,  # noqa: FBT001
    /,
) -> bool:
    return left_statement is right_statement


def exact_incircle_test(
    point_x: float,
    point_y: float,
    first_x: float,
    first_y: float,
    second_x: float,
    second_y: float,
    third_x: float,
    third_y: float,
) -> int:
    exact_point_x, exact_point_y = Fraction(point_x), Fraction(point_y)
    first_dx, first_dy = (
        Fraction(first_x) - exact_point_x,
        Fraction(first_y) - exact_point_y,
    )
    second_dx, second_dy = (
        Fraction(second_x) - exact_point_x,
        Fraction(second_y) - exact_point_y,
    )
    third_dx, third_dy = (
        Fraction(third_x) - exact_point_x,
        Fraction(third_y) - exact_point_y,
    )
    return to_sign(
        (
            (first_dx * first_dx + first_dy * first_dy)
            * (second_dx * third_dy - second_dy * third_dx)
        )
        - (
            (second_dx * second_dx + second_dy * second_dy)
            * (first_dx * third_dy - first_dy * third_dx)
        )
        + (
            (third_dx * third_dx + third_dy * third_dy)
            * (first_dx * second_dy - first_dy * second_dx)
        )
    )


def exact_kind(
    vertex_x: float,
    vertex_y: float,
    first_ray_point_x: float,
    first_ray_point_y: float,
    second_ray_point_x: float,
    second_ray_point_y: float,
) -> int:
    return to_sign(
        (
            (Fraction(first_ray_point_x) - Fraction(vertex_x))
            * (Fraction(second_ray_point_x) - Fraction(vertex_x))
        )
        + (
            (Fraction(first_ray_point_y) - Fraction(vertex_y))
            * (Fraction(second_ray_point_y) - Fraction(vertex_y))
        )
    )


def exact_orientation(
    start_x: float,
    start_y: float,
    end_x: float,
    end_y: float,
    point_x: float,
    point_y: float,
) -> int:
    return to_sign(
        (Fraction(end_x) - Fraction(start_x))
        * (Fraction(point_y) - Fraction(start_y))
        - (Fraction(end_y) - Fraction(start_y))
        * (Fraction(point_x) - Fraction(start_x))
    )


def implication(
    antecedent: bool,  # noqa: FBT001
    consequent: bool,  # noqa: FBT001
    /,
) -> bool:
    return not antecedent or consequent


def is_expansion_valid(expansion: Expansion) -> bool:
    result = eval(repr(expansion), vars(shewchuk)) == expansion
    assert isinstance(result, bool), result
    return result


def pack(
    function: Callable[..., Range],
) -> Callable[[Iterable[Domain]], Range]:
    return partial(apply, function)


def pickle_round_trip(value: Any) -> Any:
    return pickle.loads(pickle.dumps(value))


skip_reference_counter_test = pytest.mark.skipif(
    platform.python_implementation() == 'PyPy',
    reason='PyPy garbage collection is not based on reference counting.',
)


def to_sign(value: Any) -> int:
    return 1 if value > 0 else (0 if not value else -1)
