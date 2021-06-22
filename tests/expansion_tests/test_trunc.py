import math

from hypothesis import given

from shewchuk import Expansion
from . import strategies
from ..utils import equivalence


@given(strategies.expansions)
def test_basic(expansion: Expansion) -> None:
    result = math.trunc(expansion)

    assert isinstance(result, int)


@given(strategies.expansions)
def test_value(expansion: Expansion) -> None:
    result = math.trunc(expansion)

    assert abs(result) <= abs(expansion)
    assert equivalence(result >= 0, expansion >= 0)
    assert equivalence(not result, not expansion)
