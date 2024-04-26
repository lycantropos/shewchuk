import pytest
from hypothesis import given

from shewchuk import Expansion
from tests.utils import LeftOperand, is_expansion_valid

from . import strategies


@given(strategies.reals, strategies.non_zero_expansions)
def test_alternatives(first: LeftOperand, second: Expansion) -> None:
    result = first / second

    assert isinstance(result, Expansion)
    assert is_expansion_valid(result)


@given(strategies.reals, strategies.zero_expansions)
def test_zero_divisor(first: LeftOperand, second: Expansion) -> None:
    with pytest.raises(ZeroDivisionError):
        first / second
