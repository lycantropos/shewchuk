import math

from hypothesis import given

from shewchuk import Expansion

from . import strategies


@given(strategies.expansions)
def test_basic(expansion: Expansion) -> None:
    result = float(expansion)

    assert isinstance(result, float)
    assert not math.isnan(result)
    assert math.isfinite(result)


@given(strategies.expansions)
def test_round_trip(expansion: Expansion) -> None:
    result = float(expansion)

    assert float(Expansion(result)) == result
