from typing import Tuple

from hypothesis import given

from shewchuk import incircle_test
from . import strategies
from ..utils import exact_incircle_test


@given(strategies.floats_octuplets)
def test_basic(octuplet: Tuple[float, float, float, float, float, float, float,
                               float]) -> None:
    (first_x, first_y, second_x, second_y, third_x, third_y, fourth_x,
     fourth_y) = octuplet

    result = incircle_test(first_x, first_y, second_x, second_y, third_x,
                           third_y, fourth_x, fourth_y)

    assert isinstance(result, int)
    assert result in (-1, 0, 1)


@given(strategies.floats_octuplets)
def test_alternatives(octuplet: Tuple[float, float, float, float, float, float,
                                      float, float]) -> None:
    (first_x, first_y, second_x, second_y, third_x, third_y, fourth_x,
     fourth_y) = octuplet

    result = incircle_test(first_x, first_y, second_x, second_y, third_x,
                           third_y, fourth_x, fourth_y)

    assert result == exact_incircle_test(first_x, first_y, second_x, second_y,
                                         third_x, third_y, fourth_x, fourth_y)
