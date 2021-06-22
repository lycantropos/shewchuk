import sys
from numbers import Real
from typing import Union

import pytest
from hypothesis import given

from shewchuk import Expansion
from tests.utils import skip_reference_counter_test
from . import strategies


@given(strategies.expansions, strategies.non_zero_reals_or_expansions)
def test_basic(first: Expansion, second: Union[Real, Expansion]) -> None:
    result = first // second

    assert isinstance(result, Expansion)


@skip_reference_counter_test
@given(strategies.expansions, strategies.non_zero_expansions)
def test_reference_counter(first: Expansion, second: Expansion) -> None:
    first_refcount_before = sys.getrefcount(first)
    second_refcount_before = sys.getrefcount(second)

    result = first // second

    first_refcount_after = sys.getrefcount(first)
    second_refcount_after = sys.getrefcount(second)
    assert first_refcount_after == first_refcount_before
    assert second_refcount_after == second_refcount_before


@given(strategies.expansions, strategies.zero_reals_or_expansions)
def test_zero_divisor(first: Expansion,
                      second: Union[Real, Expansion]) -> None:
    with pytest.raises(ZeroDivisionError):
        first // second
