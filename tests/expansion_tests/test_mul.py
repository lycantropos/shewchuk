import sys
from numbers import Real
from typing import Union

from hypothesis import given

from shewchuk import Expansion
from tests.utils import skip_reference_counter_test
from . import strategies


@given(strategies.expansions, strategies.reals_or_expansions)
def test_basic(first: Expansion, second: Union[Real, Expansion]) -> None:
    result = first * second

    assert isinstance(result, Expansion)


@given(strategies.expansions, strategies.reals)
def test_real_argument(first: Expansion, second: Real) -> None:
    assert first * second == first * Expansion(second)


@given(strategies.expansions, strategies.zero_expansions)
def test_absorbing_element(first: Expansion, second: Expansion) -> None:
    assert first * second == second == second * first


@skip_reference_counter_test
@given(strategies.expansions, strategies.reals)
def test_reference_counter(first: Expansion, second: Real) -> None:
    first_refcount_before = sys.getrefcount(first)

    result = first * second

    first_refcount_after = sys.getrefcount(first)
    assert first_refcount_after == first_refcount_before
