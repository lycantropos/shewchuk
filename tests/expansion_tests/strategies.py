from fractions import Fraction

import math
from typing import Sequence

from hypothesis import strategies

from shewchuk import Expansion
from tests.strategies import finite_floats
from tests.utils import (MAX_VALUE,
                         pack)

precisions = strategies.none() | strategies.integers(-10, 10)
finite_floats = finite_floats
integers = strategies.integers(-MAX_VALUE, MAX_VALUE)
rationals = integers | strategies.fractions(-MAX_VALUE, MAX_VALUE)
reals = rationals | finite_floats


def is_floats_sequence_sum_finite(values: Sequence[float]) -> bool:
    return math.isfinite(2 * sum(map(abs, values)))


finite_floats_sequences = (strategies.lists(finite_floats)
                           .filter(is_floats_sequence_sum_finite))
floats = strategies.floats(allow_infinity=True,
                           allow_nan=True)


def is_floats_sequence_sum_non_finite(values: Sequence[float]) -> bool:
    return not math.isfinite(sum(values))


non_finite_floats_sequences = (strategies.lists(floats)
                               .filter(is_floats_sequence_sum_non_finite))
expansions = strategies.builds(pack(Expansion), finite_floats_sequences)
invalid_components = strategies.lists(rationals | expansions,
                                      min_size=2)
non_zero_expansions = expansions.filter(bool)
non_zero_reals = reals.filter(bool)
non_zero_reals_or_expansions = non_zero_reals | non_zero_expansions
reals_or_expansions = reals | expansions
zero_integers = strategies.builds(int)
zero_rationals = zero_integers | strategies.builds(Fraction)
zero_reals = zero_rationals | strategies.builds(float)
zero_expansions = strategies.builds(Expansion)
zero_reals_or_expansions = zero_reals | zero_expansions
ones = strategies.just(1) | strategies.just(1.0)
ones |= strategies.builds(Expansion, ones)
