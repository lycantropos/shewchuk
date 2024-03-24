from hypothesis import given

from shewchuk import Expansion

from tests.utils import LeftOperand
from . import strategies


@given(strategies.reals, strategies.expansions)
def test_alternatives(first: LeftOperand, second: Expansion) -> None:
    result = first - second

    assert result == first + (-second)
