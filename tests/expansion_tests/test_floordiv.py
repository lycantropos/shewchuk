import math
import sys
from numbers import Real
from typing import Union

import pytest
from hypothesis import given

from shewchuk import Expansion
from tests.utils import (implication,
                         skip_reference_counter_test)
from . import strategies


@given(strategies.expansions, strategies.non_zero_reals_or_expansions)
def test_basic(dividend: Expansion, divisor: Union[Real, Expansion]) -> None:
    result = dividend // divisor

    assert isinstance(result, Expansion)
    assert math.isfinite(result)


@given(strategies.expansions, strategies.non_zero_reals_or_expansions)
def test_value(dividend: Expansion, divisor: Union[Real, Expansion]) -> None:
    result = dividend // divisor

    assert implication(not dividend, not result)


@given(strategies.expansions, strategies.non_zero_reals_or_expansions)
def test_connection_with_mod(dividend: Expansion,
                             divisor: Union[Real, Expansion]) -> None:
    result = dividend // divisor

    assert result * divisor + (dividend % divisor) == dividend


@skip_reference_counter_test
@given(strategies.expansions, strategies.non_zero_expansions)
def test_reference_counter(dividend: Expansion, divisor: Expansion) -> None:
    first_refcount_before = sys.getrefcount(dividend)
    second_refcount_before = sys.getrefcount(divisor)

    result = dividend // divisor

    first_refcount_after = sys.getrefcount(dividend)
    second_refcount_after = sys.getrefcount(divisor)
    assert first_refcount_after == first_refcount_before
    assert second_refcount_after == second_refcount_before


@given(strategies.expansions, strategies.zero_reals_or_expansions)
def test_zero_divisor(dividend: Expansion, divisor: Union[Real, Expansion]
                      ) -> None:
    with pytest.raises(ZeroDivisionError):
        dividend // divisor
