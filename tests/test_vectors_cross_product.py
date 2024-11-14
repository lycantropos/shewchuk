from hypothesis import given

from shewchuk import Expansion, vectors_cross_product

from . import strategies


@given(strategies.floats_octuplets)
def test_basic(
    octuplet: tuple[float, float, float, float, float, float, float, float],
) -> None:
    (
        first_start_x,
        first_start_y,
        first_end_x,
        first_end_y,
        second_start_x,
        second_start_y,
        second_end_x,
        second_end_y,
    ) = octuplet

    result = vectors_cross_product(
        first_start_x,
        first_start_y,
        first_end_x,
        first_end_y,
        second_start_x,
        second_start_y,
        second_end_x,
        second_end_y,
    )

    assert isinstance(result, Expansion)


@given(strategies.floats_quadruplets)
def test_same_endpoints(quadruplet: tuple[float, float, float, float]) -> None:
    first_start_x, first_start_y, first_end_x, first_end_y = quadruplet

    assert not vectors_cross_product(
        first_start_x,
        first_start_y,
        first_end_x,
        first_end_y,
        first_start_x,
        first_start_y,
        first_end_x,
        first_end_y,
    )


@given(strategies.floats_octuplets)
def test_segments_permutation(
    octuplet: tuple[float, float, float, float, float, float, float, float],
) -> None:
    (
        first_start_x,
        first_start_y,
        first_end_x,
        first_end_y,
        second_start_x,
        second_start_y,
        second_end_x,
        second_end_y,
    ) = octuplet

    result = vectors_cross_product(
        first_start_x,
        first_start_y,
        first_end_x,
        first_end_y,
        second_start_x,
        second_start_y,
        second_end_x,
        second_end_y,
    )

    assert result == -vectors_cross_product(
        second_start_x,
        second_start_y,
        second_end_x,
        second_end_y,
        first_start_x,
        first_start_y,
        first_end_x,
        first_end_y,
    )
