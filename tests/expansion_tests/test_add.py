import sys
from numbers import Real
from typing import Union

from hypothesis import given

from shewchuk import Expansion
from tests.utils import (is_expansion_valid,
                         skip_reference_counter_test)
from . import strategies


@given(strategies.expansions, strategies.reals_or_expansions)
def test_basic(first: Expansion, second: Union[Real, Expansion]) -> None:
    result = first + second

    assert isinstance(result, Expansion)
    assert is_expansion_valid(result)


@given(strategies.expansions, strategies.expansions)
def test_commutativity(first: Expansion, second: Expansion) -> None:
    assert first + second == second + first


@given(strategies.expansions, strategies.zero_reals_or_expansions)
def test_neutral_element(first: Expansion,
                         second: Union[Real, Expansion]) -> None:
    assert first + second == first == second + first


@skip_reference_counter_test
@given(strategies.expansions, strategies.expansions)
def test_reference_counter(first: Expansion, second: Expansion) -> None:
    first_refcount_before = sys.getrefcount(first)
    second_refcount_before = sys.getrefcount(second)

    result = first + second

    first_refcount_after = sys.getrefcount(first)
    second_refcount_after = sys.getrefcount(second)
    assert first_refcount_after == first_refcount_before
    assert second_refcount_after == second_refcount_before
