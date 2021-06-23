import sys

from hypothesis import given

from shewchuk import Expansion
from tests.utils import (implication,
                         skip_reference_counter_test)
from . import strategies


@given(strategies.expansions)
def test_basic(expansion: Expansion) -> None:
    result = hash(expansion)

    assert isinstance(result, int)


@given(strategies.expansions)
def test_determinism(expansion: Expansion) -> None:
    result = hash(expansion)

    assert result == hash(expansion)


@given(strategies.expansions, strategies.expansions)
def test_connection_with_equality(left: Expansion, right: Expansion) -> None:
    assert implication(left == right, hash(left) == hash(right))


@skip_reference_counter_test
@given(strategies.expansions)
def test_reference_counter(expansion: Expansion) -> None:
    expansion_refcount_before = sys.getrefcount(expansion)

    result = hash(expansion)

    expansion_refcount_after = sys.getrefcount(expansion)
    assert expansion_refcount_after == expansion_refcount_before
