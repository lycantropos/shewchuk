import sys

from hypothesis import given

from shewchuk import Expansion
from tests.utils import RightOperand, equivalence, skip_reference_counter_test

from . import strategies


@given(strategies.expansions)
def test_basic(expansion: Expansion) -> None:
    result = abs(expansion)

    assert isinstance(result, Expansion)


@given(strategies.expansions)
def test_idempotence(expansion: Expansion) -> None:
    result = abs(expansion)

    assert result == abs(result)


@given(strategies.expansions)
def test_positive_definiteness(expansion: Expansion) -> None:
    result = abs(expansion)

    assert equivalence(not result, not expansion)


@given(strategies.expansions)
def test_evenness(expansion: Expansion) -> None:
    result = abs(expansion)

    assert result == abs(-expansion)


@given(strategies.expansions, strategies.reals_or_expansions)
def test_triangle_inequality(first: Expansion, second: RightOperand) -> None:
    result = abs(first + second)

    assert result <= abs(first) + abs(second)


@skip_reference_counter_test
@given(strategies.expansions)
def test_reference_counter(expansion: Expansion) -> None:
    expansion_refcount_before = sys.getrefcount(expansion)

    result = abs(expansion)

    expansion_refcount_after = sys.getrefcount(expansion)
    assert expansion_refcount_after == (
        expansion_refcount_before + (result == expansion)
    )
