import platform
from fractions import Fraction
from functools import partial
from typing import (Callable,
                    Iterable,
                    Tuple,
                    TypeVar)

import pytest
from hypothesis.strategies import SearchStrategy as Strategy

Strategy = Strategy
Domain = TypeVar('Domain')
Range = TypeVar('Range')


def equivalence(left_statement: bool, right_statement: bool) -> bool:
    return left_statement is right_statement


def implication(antecedent: bool, consequent: bool) -> bool:
    return not antecedent or consequent


def pack(function: Callable[..., Range]
         ) -> Callable[[Iterable[Domain]], Range]:
    return partial(apply, function)


def apply(function: Callable[..., Range], args: Tuple[Domain, ...]) -> Range:
    return function(*args)


skip_reference_counter_test = pytest.mark.skipif(
        platform.python_implementation() == 'PyPy',
        reason='PyPy\'s garbage collection '
               'is not based on reference counting.')


def to_sign(value: Domain) -> int:
    return 1 if value > 0 else (0 if not value else -1)


def exact_kind(vertex_x: float,
               vertex_y: float,
               first_ray_point_x: float,
               first_ray_point_y: float,
               second_ray_point_x: float,
               second_ray_point_y: float) -> int:
    return to_sign((Fraction(first_ray_point_x) - Fraction(vertex_x))
                   * (Fraction(second_ray_point_x) - Fraction(vertex_x))
                   + (Fraction(first_ray_point_y) - Fraction(vertex_y))
                   * (Fraction(second_ray_point_y) - Fraction(vertex_y)))


def exact_orientation(start_x: float,
                      start_y: float,
                      end_x: float,
                      end_y: float,
                      point_x: float,
                      point_y: float) -> int:
    return to_sign((Fraction(end_x) - Fraction(start_x))
                   * (Fraction(point_y) - Fraction(start_y))
                   - (Fraction(end_y) - Fraction(start_y))
                   * (Fraction(point_x) - Fraction(start_x)))
