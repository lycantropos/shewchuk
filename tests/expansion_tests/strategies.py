import math
from typing import (Sequence,
                    Tuple)

from hypothesis import strategies

from tests.strategies import finite_floats
from tests.utils import Strategy


def is_floats_sequence_sum_finite(values: Sequence[float]) -> bool:
    return math.isfinite(sum(map(abs, values)))


finite_floats_sequences = (strategies.lists(finite_floats)
                           .filter(is_floats_sequence_sum_finite))


def to_finite_floats_sequences_with_permutations(
        values: Sequence[float]) -> Strategy[Tuple[Sequence[float],
                                                   Sequence[float]]]:
    return strategies.tuples(strategies.just(values),
                             strategies.permutations(values))


finite_floats_sequences_with_permutations = finite_floats_sequences.flatmap(
        to_finite_floats_sequences_with_permutations)
