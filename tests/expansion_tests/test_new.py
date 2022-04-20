from numbers import Rational
from typing import (Any,
                    Sequence)

import pytest
from hypothesis import given

from shewchuk import Expansion
from . import strategies


@given(strategies.finite_floats_sequences)
def test_basic(components: Sequence[float]) -> None:
    result = Expansion(*components)

    assert isinstance(result, Expansion)


def test_no_argument() -> None:
    result = Expansion()

    assert not result


@given(strategies.expansions)
def test_expansion_argument(value: Expansion) -> None:
    result = Expansion(value)

    assert result == value


@given(strategies.finite_floats)
def test_float_argument(value: float) -> None:
    result = Expansion(value)

    assert result == value


@given(strategies.rationals)
def test_rational_argument(value: Rational) -> None:
    result = Expansion(value)

    assert result == value


@given(strategies.finite_floats_sequences)
def test_determinism(components: Sequence[float]) -> None:
    result = Expansion(*components)

    assert result == Expansion(*components)


@given(strategies.invalid_components)
def test_invalid_components_types(components: Sequence[Any]) -> None:
    with pytest.raises(TypeError):
        Expansion(*components)


@given(strategies.non_finite_floats_sequences)
def test_invalid_components_values(components: Sequence[float]) -> None:
    with pytest.raises(ValueError):
        Expansion(*components)
