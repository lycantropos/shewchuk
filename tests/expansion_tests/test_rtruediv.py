from numbers import Real

import pytest
from hypothesis import given

from shewchuk import Expansion
from . import strategies


@given(strategies.reals, strategies.zero_expansions)
def test_zero_divisor(first: Real, second: Expansion) -> None:
    with pytest.raises(ZeroDivisionError):
        first / second
