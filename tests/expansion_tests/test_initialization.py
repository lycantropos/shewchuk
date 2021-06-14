from typing import Sequence

from hypothesis import given

from shewchuk import Expansion
from . import strategies


@given(strategies.finite_floats_sequences)
def test_basic(components: Sequence[float]) -> None:
    result = Expansion(*components)

    assert isinstance(result, Expansion)
