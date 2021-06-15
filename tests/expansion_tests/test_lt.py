import sys
from numbers import Real

from hypothesis import given

from shewchuk import Expansion
from tests.utils import (equivalence,
                         implication,
                         skip_reference_counter_test)
from . import strategies


@given(strategies.expansions)
def test_irreflexivity(expansion: Expansion) -> None:
    assert not expansion < expansion


@given(strategies.expansions, strategies.reals)
def test_asymmetry(first: Expansion, second: Real) -> None:
    assert implication(first < second, not second < first)


@given(strategies.expansions, strategies.expansions, strategies.expansions)
def test_transitivity(first: Expansion,
                      second: Expansion,
                      third: Expansion) -> None:
    assert implication(first < second < third,
                       first < third)


@given(strategies.expansions, strategies.reals)
def test_equivalents(first: Expansion, second: Real) -> None:
    result = first < second

    assert equivalence(result, second > first)
    assert equivalence(result, second >= first != second)
    assert equivalence(result, first <= second != first)


@given(strategies.expansions, strategies.finite_floats)
def test_float_operand(first: Expansion, second: float) -> None:
    assert implication(first < second, float(first) < second)


@skip_reference_counter_test
@given(strategies.expansions, strategies.expansions)
def test_reference_counter(first: Expansion, second: Expansion) -> None:
    first_refcount_before = sys.getrefcount(first)
    second_refcount_before = sys.getrefcount(second)

    result = first < second

    first_refcount_after = sys.getrefcount(first)
    second_refcount_after = sys.getrefcount(second)
    assert first_refcount_after == first_refcount_before
    assert second_refcount_after == second_refcount_before
