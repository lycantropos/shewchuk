import sys

import pytest
from hypothesis import given

from shewchuk import Expansion
from tests.utils import skip_reference_counter_test
from . import strategies


@given(strategies.expansions, strategies.small_non_negative_integers)
def test_basic(base: Expansion, exponent: int) -> None:
    result = base ** exponent

    assert isinstance(result, float)


@given(strategies.zero_expansions, strategies.small_positive_integers)
def test_left_absorbing_element(base: Expansion, exponent: int) -> None:
    assert base ** exponent == base


@given(strategies.expansions, strategies.ones)
def test_right_neutral_element(base: Expansion, exponent: int) -> None:
    assert base ** exponent == base


@skip_reference_counter_test
@given(strategies.expansions, strategies.small_non_negative_integers)
def test_reference_counter(base: Expansion, exponent: int) -> None:
    base_refcount_before = sys.getrefcount(base)

    result = base ** exponent

    base_refcount_after = sys.getrefcount(base)
    assert base_refcount_after == base_refcount_before


@given(strategies.zero_expansions, strategies.negative_integers)
def test_zero_base(base: Expansion, exponent: int) -> None:
    with pytest.raises(ZeroDivisionError):
        base ** exponent
