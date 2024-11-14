from hypothesis import given

from shewchuk import kind, vectors_dot_product
from tests.utils import exact_kind, to_sign

from . import strategies


@given(strategies.floats_sextuplets)
def test_basic(
    sextuplet: tuple[float, float, float, float, float, float],
) -> None:
    (
        vertex_x,
        vertex_y,
        first_ray_second_ray_point_x,
        first_ray_second_ray_point_y,
        second_ray_point_x,
        second_ray_point_y,
    ) = sextuplet

    result = kind(
        vertex_x,
        vertex_y,
        first_ray_second_ray_point_x,
        first_ray_second_ray_point_y,
        second_ray_point_x,
        second_ray_point_y,
    )

    assert isinstance(result, int)
    assert result in (-1, 0, 1)


@given(strategies.floats_quadruplets)
def test_endpoints(quadruplet: tuple[float, float, float, float]) -> None:
    (
        vertex_x,
        vertex_y,
        first_ray_second_ray_point_x,
        first_ray_second_ray_point_y,
    ) = quadruplet

    assert not kind(
        vertex_x,
        vertex_y,
        first_ray_second_ray_point_x,
        first_ray_second_ray_point_y,
        vertex_x,
        vertex_y,
    )
    assert kind(
        vertex_x,
        vertex_y,
        first_ray_second_ray_point_x,
        first_ray_second_ray_point_y,
        first_ray_second_ray_point_x,
        first_ray_second_ray_point_y,
    ) == (
        vertex_x != first_ray_second_ray_point_x
        or vertex_y != first_ray_second_ray_point_y
    )


@given(strategies.floats_sextuplets)
def test_endpoints_permutation(
    sextuplet: tuple[float, float, float, float, float, float],
) -> None:
    (
        vertex_x,
        vertex_y,
        first_ray_second_ray_point_x,
        first_ray_second_ray_point_y,
        second_ray_point_x,
        second_ray_point_y,
    ) = sextuplet

    result = kind(
        vertex_x,
        vertex_y,
        first_ray_second_ray_point_x,
        first_ray_second_ray_point_y,
        second_ray_point_x,
        second_ray_point_y,
    )

    assert result == kind(
        vertex_x,
        vertex_y,
        second_ray_point_x,
        second_ray_point_y,
        first_ray_second_ray_point_x,
        first_ray_second_ray_point_y,
    )


@given(strategies.floats_sextuplets)
def test_alternatives(
    sextuplet: tuple[float, float, float, float, float, float],
) -> None:
    (
        vertex_x,
        vertex_y,
        first_ray_second_ray_point_x,
        first_ray_second_ray_point_y,
        second_ray_point_x,
        second_ray_point_y,
    ) = sextuplet

    result = kind(
        vertex_x,
        vertex_y,
        first_ray_second_ray_point_x,
        first_ray_second_ray_point_y,
        second_ray_point_x,
        second_ray_point_y,
    )

    assert result == to_sign(
        vectors_dot_product(
            vertex_x,
            vertex_y,
            first_ray_second_ray_point_x,
            first_ray_second_ray_point_y,
            vertex_x,
            vertex_y,
            second_ray_point_x,
            second_ray_point_y,
        )
    )
    assert result == exact_kind(
        vertex_x,
        vertex_y,
        first_ray_second_ray_point_x,
        first_ray_second_ray_point_y,
        second_ray_point_x,
        second_ray_point_y,
    )
