from __future__ import annotations

import math
from fractions import Fraction
from functools import reduce
from typing import Sequence

from hypothesis import strategies as _st

from shewchuk import Expansion
from tests import strategies as _strategies
from tests.utils import MAX_VALUE, pack

precisions = _st.none() | _st.integers(-10, 10)
finite_floats = _strategies.finite_floats
integers = _st.integers(-MAX_VALUE, MAX_VALUE)
rationals = integers | _st.fractions(-MAX_VALUE, MAX_VALUE)
reals = rationals | finite_floats


def is_floats_sequence_sum_finite(values: Sequence[float]) -> bool:
    return math.isfinite(2 * sum(map(abs, values)))


finite_floats_sequences = _st.lists(finite_floats).filter(
    is_floats_sequence_sum_finite
)
floats = _st.floats(allow_infinity=True, allow_nan=True)


def is_invalid_floats_sequence(values: Sequence[float]) -> bool:
    def to_next_component(component: float, value: float) -> float:
        tail, head = _two_add(component, value)
        return tail or head

    return not math.isfinite(reduce(to_next_component, reversed(values), 0.0))


def _are_non_overlapping(values: Sequence[float]) -> bool:
    return all(
        _do_not_overlap(value, values[next_value_index])
        for index, value in enumerate(values)
        for next_value_index in range(index + 1, len(values))
    )


def _do_not_overlap(first: float, second: float) -> bool:
    return _two_add(first, second) == (min(first, second), max(first, second))


def _two_add(left: float, right: float) -> tuple[float, float]:
    head = left + right
    right_virtual = head - left
    left_virtual = head - right_virtual
    right_tail = right - right_virtual
    left_tail = left - left_virtual
    tail = left_tail + right_tail
    return tail, head


invalid_floats_sequences = _st.lists(floats).filter(is_invalid_floats_sequence)
expansions = _st.builds(pack(Expansion), finite_floats_sequences)
invalid_components = _st.lists(rationals | expansions, min_size=2)
non_zero_expansions = expansions.filter(bool)
non_zero_reals = reals.filter(bool)
non_zero_reals_or_expansions = non_zero_reals | non_zero_expansions
reals_or_expansions = reals | expansions
zero_integers = _st.builds(int)
zero_rationals = zero_integers | _st.builds(Fraction)
zero_reals = zero_rationals | _st.builds(float)
zero_expansions = _st.builds(Expansion)
zero_reals_or_expansions = zero_reals | zero_expansions
ones: _st.SearchStrategy[Expansion | float | int] = _st.just(1) | _st.just(1.0)
ones |= _st.builds(Expansion, ones)
