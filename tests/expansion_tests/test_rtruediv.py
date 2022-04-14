import pytest
from hypothesis import given

from shewchuk import Expansion
from tests.utils import LeftOperand
from . import strategies


@given(strategies.reals, strategies.zero_expansions)
def test_zero_divisor(first: LeftOperand, second: Expansion) -> None:
    with pytest.raises(ZeroDivisionError):
        first / second
