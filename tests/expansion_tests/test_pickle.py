from hypothesis import given

from shewchuk import Expansion
from tests.utils import pickle_round_trip
from . import strategies


@given(strategies.expansions)
def test_round_trip(expansion: Expansion) -> None:
    result = pickle_round_trip(expansion)

    assert result == expansion
    assert result is not expansion
