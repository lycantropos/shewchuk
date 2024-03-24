from hypothesis import strategies

from tests.utils import MAX_VALUE

MIN_POSITIVE_FLOAT = 1.0 / MAX_VALUE
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
