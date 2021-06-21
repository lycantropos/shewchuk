import math
from fractions import Fraction
from typing import Sequence

from hypothesis import strategies

from shewchuk import Expansion
from tests.strategies import finite_floats
from tests.utils import pack

finite_floats = finite_floats
floats_quadruplets = strategies.lists(finite_floats,
                                      min_size=4,
                                      max_size=4)
floats_sextuplets = strategies.lists(finite_floats,
                                     min_size=6,
                                     max_size=6)
floats_octuplets = strategies.lists(finite_floats,
                                    min_size=8,
                                    max_size=8)
reals = strategies.integers() | strategies.fractions() | finite_floats


def is_floats_sequence_sum_finite(values: Sequence[float]) -> bool:
    return math.isfinite(2 * sum(map(abs, values)))


finite_floats_sequences = (strategies.lists(finite_floats)
                           .filter(is_floats_sequence_sum_finite))
expansions = strategies.builds(pack(Expansion), finite_floats_sequences)
non_zero_expansions = expansions.filter(bool)
non_zero_reals = reals.filter(bool)
non_zero_reals_or_expansions = non_zero_reals | non_zero_expansions
reals_or_expansions = reals | expansions
zero_reals = (strategies.builds(int) | strategies.builds(Fraction)
              | strategies.builds(float))
zero_expansions = strategies.builds(Expansion)
zero_reals_or_expansions = zero_reals | zero_expansions
