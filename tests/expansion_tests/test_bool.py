import sys

from hypothesis import given

from shewchuk import Expansion
from tests.utils import equivalence, skip_reference_counter_test

from . import strategies


@given(strategies.expansions)
def test_properties(expansion: Expansion) -> None:
    assert equivalence(bool(expansion), bool(float(expansion)))


@skip_reference_counter_test
@given(strategies.expansions)
def test_reference_counter(expansion: Expansion) -> None:
    expansion_refcount_before = sys.getrefcount(expansion)

    bool(expansion)

    expansion_refcount_after = sys.getrefcount(expansion)
    assert expansion_refcount_after == expansion_refcount_before
