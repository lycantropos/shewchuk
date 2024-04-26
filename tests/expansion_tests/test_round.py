from typing import Optional

from hypothesis import given

from shewchuk import Expansion
from tests.utils import to_sign

from . import strategies


@given(strategies.expansions, strategies.precisions)
def test_basic(expansion: Expansion, precision: Optional[int]) -> None:
    result = round(expansion, precision)

    assert isinstance(result, int if precision is None else Expansion)


@given(strategies.expansions)
def test_value(expansion: Expansion) -> None:
    result = round(expansion)

    truncated_expansion = int(expansion)
    distance = abs(truncated_expansion - expansion)
    assert (
        result == truncated_expansion + to_sign(expansion)
        if (distance > 0.5 or distance == 0.5 and truncated_expansion % 2 == 1)
        else result == truncated_expansion
    )
