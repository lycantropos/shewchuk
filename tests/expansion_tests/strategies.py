import math
import sys
from fractions import Fraction
from typing import Sequence

from hypothesis import strategies

from shewchuk import Expansion
from tests.strategies import finite_floats
from tests.utils import pack

negative_integers = strategies.integers(max_value=-1)
ones = strategies.just(1)
small_positive_integers = strategies.integers(1, 5)
small_non_negative_integers = strategies.just(0) | small_positive_integers
precisions = strategies.none() | strategies.integers(-10, 10)
finite_floats = finite_floats
MAX_FLOAT_REPRESENTABLE_INTEGER = int(sys.float_info.max)
integers = strategies.integers(-MAX_FLOAT_REPRESENTABLE_INTEGER,
                               MAX_FLOAT_REPRESENTABLE_INTEGER)
reals = integers | finite_floats


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
