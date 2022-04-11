import sys

from hypothesis import given

from shewchuk import Expansion
from tests.utils import (RightOperand,
                         equivalence,
                         implication,
                         skip_reference_counter_test)
from . import strategies


@given(strategies.expansions)
def test_reflexivity(expansion: Expansion) -> None:
    assert expansion >= expansion


@given(strategies.expansions, strategies.reals_or_expansions)
def test_antisymmetry(first: Expansion, second: RightOperand) -> None:
    assert equivalence(first >= second >= first, first == second)


@given(strategies.expansions, strategies.expansions, strategies.expansions)
def test_transitivity(first: Expansion, second: Expansion, third: Expansion
                      ) -> None:
    assert implication(first >= second >= third, first >= third)


@given(strategies.expansions, strategies.reals_or_expansions)
def test_equivalents(first: Expansion, second: RightOperand) -> None:
    result = first >= second

    assert equivalence(result, second <= first)
    assert equivalence(result, first > second or first == second)
    assert equivalence(result, second < first or first == second)


@skip_reference_counter_test
@given(strategies.expansions, strategies.expansions)
def test_reference_counter(first: Expansion, second: Expansion) -> None:
    first_refcount_before = sys.getrefcount(first)
    second_refcount_before = sys.getrefcount(second)

    result = first >= second

    first_refcount_after = sys.getrefcount(first)
    second_refcount_after = sys.getrefcount(second)
    assert first_refcount_after == first_refcount_before
    assert second_refcount_after == second_refcount_before
