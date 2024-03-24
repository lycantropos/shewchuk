from hypothesis import given

from shewchuk import Expansion

from . import strategies


@given(strategies.expansions)
def test_basic(expansion: Expansion) -> None:
    result = expansion.as_integer_ratio()

    assert isinstance(result, tuple)
    assert all(isinstance(element, int) for element in result)
    assert len(result) == 2


@given(strategies.expansions)
def test_components(expansion: Expansion) -> None:
    result_numerator, result_denominator = expansion.as_integer_ratio()

    assert result_denominator > 0
    assert result_denominator != 1 or result_numerator == int(expansion)
