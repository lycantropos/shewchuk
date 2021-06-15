from typing import (Sequence,
                    Tuple)

from hypothesis import given

from shewchuk import Expansion
from . import strategies


@given(strategies.finite_floats_sequences)
def test_basic(components: Sequence[float]) -> None:
    result = Expansion(*components)

    assert isinstance(result, Expansion)


@given(strategies.finite_floats_sequences_with_permutations)
def test_permutation(components_with_permutation
                     : Tuple[Sequence[float], Sequence[float]]) -> None:
    components, permuted_components = components_with_permutation

    assert Expansion(*components) == Expansion(*permuted_components)
