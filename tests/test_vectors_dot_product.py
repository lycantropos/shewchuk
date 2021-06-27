from typing import Tuple

from hypothesis import given

from shewchuk import (Expansion,
                      vectors_dot_product)
from . import strategies


@given(strategies.floats_octuplets)
def test_basic(octuplet: Tuple[float, float, float, float, float, float, float,
                               float]) -> None:
    (first_start_x, first_start_y, first_end_x, first_end_y, second_start_x,
     second_start_y, second_end_x, second_end_y) = octuplet

    result = vectors_dot_product(first_start_x, first_start_y, first_end_x,
                                 first_end_y, second_start_x, second_start_y,
                                 second_end_x, second_end_y)

    assert isinstance(result, Expansion)


@given(strategies.floats_quadruplets)
def test_perpendicular_endpoints(quadruplet: Tuple[float, float, float, float]
                                 ) -> None:
    first_start_x, first_start_y, first_end_x, first_end_y = quadruplet

    assert not vectors_dot_product(first_start_x, first_start_y, first_end_x,
                                   first_end_y, -first_start_y, first_start_x,
                                   -first_end_y, first_end_x)


@given(strategies.floats_octuplets)
def test_segments_permutation(octuplet: Tuple[float, float, float, float,
                                              float, float, float, float]
                              ) -> None:
    (first_start_x, first_start_y, first_end_x, first_end_y, second_start_x,
     second_start_y, second_end_x, second_end_y) = octuplet

    result = vectors_dot_product(first_start_x, first_start_y, first_end_x,
                                 first_end_y, second_start_x, second_start_y,
                                 second_end_x, second_end_y)

    assert result == vectors_dot_product(second_start_x, second_start_y,
                                         second_end_x, second_end_y,
                                         first_start_x, first_start_y,
                                         first_end_x, first_end_y)
