import sys

from hypothesis import given

from shewchuk import Expansion

from tests.utils import skip_reference_counter_test
from . import strategies


@given(strategies.expansions)
def test_basic(expansion: Expansion) -> None:
    result = expansion.imag

    assert isinstance(result, int)


@given(strategies.expansions)
def test_zeroness(expansion: Expansion) -> None:
    result = expansion.imag

    assert result == 0


@skip_reference_counter_test
@given(strategies.expansions)
def test_reference_counter(expansion: Expansion) -> None:
    expansion_refcount_before = sys.getrefcount(expansion)

    _ = expansion.imag

    expansion_refcount_after = sys.getrefcount(expansion)
    assert expansion_refcount_after == expansion_refcount_before
