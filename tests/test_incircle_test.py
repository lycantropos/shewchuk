from hypothesis import given

from shewchuk import incircle_test
from tests.utils import exact_incircle_test

from . import strategies


@given(strategies.floats_octuplets)
def test_basic(
    octuplet: tuple[float, float, float, float, float, float, float, float],
) -> None:
    (
        point_x,
        point_y,
        first_x,
        first_y,
        second_x,
        second_y,
        third_x,
        third_y,
    ) = octuplet

    result = incircle_test(
        point_x,
        point_y,
        first_x,
        first_y,
        second_x,
        second_y,
        third_x,
        third_y,
    )

    assert isinstance(result, int)
    assert result in (-1, 0, 1)


@given(strategies.floats_octuplets)
def test_alternatives(
    octuplet: tuple[float, float, float, float, float, float, float, float],
) -> None:
    (
        point_x,
        point_y,
        first_x,
        first_y,
        second_x,
        second_y,
        third_x,
        third_y,
    ) = octuplet

    result = incircle_test(
        point_x,
        point_y,
        first_x,
        first_y,
        second_x,
        second_y,
        third_x,
        third_y,
    )

    assert result == exact_incircle_test(
        point_x,
        point_y,
        first_x,
        first_y,
        second_x,
        second_y,
        third_x,
        third_y,
    )
