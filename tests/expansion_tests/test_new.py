from itertools import permutations
from typing import Sequence

from hypothesis import given

from shewchuk import Expansion
from . import strategies
from ..utils import permute


@given(strategies.finite_floats_sequences)
def test_basic(components: Sequence[float]) -> None:
    result = Expansion(*components)

    assert isinstance(result, Expansion)


def test_no_argument() -> None:
    result = Expansion()

    assert not result


@given(strategies.finite_floats)
def test_float_argument(value: float) -> None:
    result = Expansion(value)

    assert result == value


@given(strategies.finite_floats_sequences)
def test_determinism(components: Sequence[float]) -> None:
    result = Expansion(*components)

    assert result == Expansion(*components)


@given(strategies.finite_floats_sequences, strategies.indices)
def test_permutations(components: Sequence[float], index: int) -> None:
    result = Expansion(*components)

    assert result == Expansion(*permute(components, index))
