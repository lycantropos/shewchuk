import sys

from hypothesis import given

import shewchuk
from shewchuk import Expansion
from tests.utils import skip_reference_counter_test

from . import strategies


@given(strategies.expansions)
def test_basic(expansion: Expansion) -> None:
    result = repr(expansion)

    assert result.startswith(Expansion.__qualname__)


@given(strategies.expansions)
def test_round_trip(expansion: Expansion) -> None:
    result = repr(expansion)

    assert eval(result, vars(shewchuk)) == expansion


@skip_reference_counter_test
@given(strategies.expansions)
def test_reference_counter(expansion: Expansion) -> None:
    expansion_refcount_before = sys.getrefcount(expansion)

    repr(expansion)

    expansion_refcount_after = sys.getrefcount(expansion)
    assert expansion_refcount_after == expansion_refcount_before
