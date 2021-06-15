from numbers import Real

from hypothesis import given

from shewchuk import Expansion
from . import strategies


@given(strategies.reals, strategies.expansions)
def test_connection_with_add(first: Real, second: Expansion) -> None:
    result = first + second

    assert result == second + first
