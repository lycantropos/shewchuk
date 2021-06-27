from typing import Optional

from hypothesis import given

from shewchuk import Expansion
from . import strategies


@given(strategies.expansions, strategies.precisions)
def test_basic(expansion: Expansion, precision: Optional[int]) -> None:
    result = round(expansion, precision)

    assert isinstance(result, int if precision is None else float)


@given(strategies.expansions, strategies.precisions)
def test_value(expansion: Expansion, precision: Optional[int]) -> None:
    result = round(expansion, precision)

    assert result == round(float(expansion), precision)
