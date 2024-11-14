from hypothesis import given

from shewchuk import orientation, vectors_cross_product
from tests.utils import exact_orientation, to_sign

from . import strategies


@given(strategies.floats_sextuplets)
def test_basic(
    sextuplet: tuple[float, float, float, float, float, float],
) -> None:
    start_x, start_y, end_x, end_y, point_x, point_y = sextuplet

    result = orientation(start_x, start_y, end_x, end_y, point_x, point_y)

    assert isinstance(result, int)
    assert result in (-1, 0, 1)


@given(strategies.floats_quadruplets)
def test_endpoints(quadruplet: tuple[float, float, float, float]) -> None:
    start_x, start_y, end_x, end_y = quadruplet

    assert not orientation(start_x, start_y, end_x, end_y, start_x, start_y)
    assert not orientation(start_x, start_y, end_x, end_y, end_x, end_y)


@given(strategies.floats_sextuplets)
def test_endpoints_permutation(
    sextuplet: tuple[float, float, float, float, float, float],
) -> None:
    start_x, start_y, end_x, end_y, point_x, point_y = sextuplet

    result = orientation(start_x, start_y, end_x, end_y, point_x, point_y)

    assert result == -orientation(
        end_x, end_y, start_x, start_y, point_x, point_y
    )


@given(strategies.floats_sextuplets)
def test_alternatives(
    sextuplet: tuple[float, float, float, float, float, float],
) -> None:
    start_x, start_y, end_x, end_y, point_x, point_y = sextuplet

    result = orientation(start_x, start_y, end_x, end_y, point_x, point_y)

    assert result == to_sign(
        vectors_cross_product(
            start_x, start_y, end_x, end_y, start_x, start_y, point_x, point_y
        )
    )
    assert result == exact_orientation(
        start_x, start_y, end_x, end_y, point_x, point_y
    )
