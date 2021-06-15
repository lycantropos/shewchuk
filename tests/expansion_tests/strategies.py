import math
from typing import Sequence

from hypothesis import strategies

from shewchuk import Expansion
from tests.strategies import finite_floats
from tests.utils import pack

finite_floats = finite_floats


def is_floats_sequence_sum_finite(values: Sequence[float]) -> bool:
    return math.isfinite(2 * sum(map(abs, values)))


finite_floats_sequences = (strategies.lists(finite_floats)
                           .filter(is_floats_sequence_sum_finite))
expansions = strategies.builds(pack(Expansion), finite_floats_sequences)
reals = strategies.integers() | strategies.fractions() | finite_floats
reals_or_expansions = reals | expansions
