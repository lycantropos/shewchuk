import sys

from hypothesis import given

from shewchuk import Expansion
from tests.utils import skip_reference_counter_test

from . import strategies


@given(strategies.expansions)
def test_basic(expansion: Expansion) -> None:
    result = expansion.real

    assert isinstance(result, Expansion)


@given(strategies.expansions)
def test_identity(expansion: Expansion) -> None:
    result = expansion.real

    assert result == expansion


@skip_reference_counter_test
@given(strategies.expansions)
def test_reference_counter(expansion: Expansion) -> None:
    expansion_refcount_before = sys.getrefcount(expansion)

    _ = expansion.real

    expansion_refcount_after = sys.getrefcount(expansion)
    assert expansion_refcount_after == expansion_refcount_before + 1
