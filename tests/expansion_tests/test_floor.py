import math

from hypothesis import given

from shewchuk import Expansion
from . import strategies


@given(strategies.expansions)
def test_basic(expansion: Expansion) -> None:
    result = math.floor(expansion)

    assert isinstance(result, int)


@given(strategies.expansions)
def test_value(expansion: Expansion) -> None:
    result = math.floor(expansion)

    assert result <= expansion
    assert result > expansion - 1
