import sys

from hypothesis import given

from shewchuk import Expansion

from tests.utils import (
    RightOperand,
    equivalence,
    is_expansion_valid,
    skip_reference_counter_test,
)
from . import strategies


@given(strategies.expansions, strategies.reals_or_expansions)
def test_basic(first: Expansion, second: RightOperand) -> None:
    result = first + second

    assert isinstance(result, Expansion)
    assert is_expansion_valid(result)


@given(strategies.expansions)
def test_diagonal(expansion: Expansion) -> None:
    assert not expansion - expansion


@given(strategies.expansions, strategies.expansions)
def test_commutative_case(first: Expansion, second: Expansion) -> None:
    assert equivalence(first - second == second - first, first == second)


@given(strategies.expansions, strategies.zero_reals_or_expansions)
def test_right_neutral_element(first: Expansion, second: RightOperand) -> None:
    assert first - second == first


@given(strategies.expansions, strategies.expansions)
def test_alternatives(first: Expansion, second: Expansion) -> None:
    result = first - second

    assert result == first + (-second)


@skip_reference_counter_test
@given(strategies.expansions, strategies.expansions)
def test_reference_counter(first: Expansion, second: Expansion) -> None:
    first_refcount_before = sys.getrefcount(first)
    second_refcount_before = sys.getrefcount(second)

    first - second

    first_refcount_after = sys.getrefcount(first)
    second_refcount_after = sys.getrefcount(second)
    assert first_refcount_after == first_refcount_before
    assert second_refcount_after == second_refcount_before
