import sys

from hypothesis import strategies

from shewchuk import incircle_test

from tests.utils import MAX_VALUE, exact_incircle_test


def to_min_positive_float() -> float:
    result = 0.5
    while incircle_test(
        0.0, 0.0, result, 0.0, 0.0, result, -result, -result
    ) == exact_incircle_test(
        0.0, 0.0, result, 0.0, 0.0, result, -result, -result
    ):
        result /= 2.0
        assert result > sys.float_info.min, result
    return result * 2.0


MIN_POSITIVE_FLOAT = to_min_positive_float()
assert MIN_POSITIVE_FLOAT < 1.0, MIN_POSITIVE_FLOAT
finite_floats = (
    strategies.floats(-float(MAX_VALUE), -MIN_POSITIVE_FLOAT)
    | strategies.just(0.0)
    | strategies.floats(MIN_POSITIVE_FLOAT, float(MAX_VALUE))
)
floats_quadruplets = strategies.lists(
    finite_floats, min_size=4, max_size=4
).map(tuple)
floats_sextuplets = strategies.lists(
    finite_floats, min_size=6, max_size=6
).map(tuple)
floats_octuplets = strategies.lists(finite_floats, min_size=8, max_size=8).map(
    tuple
)
