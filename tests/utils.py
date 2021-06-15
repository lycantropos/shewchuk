import math
import platform
from functools import partial
from typing import (Callable,
                    Iterable,
                    Sequence,
                    Tuple,
                    TypeVar)

import pytest
from hypothesis.strategies import SearchStrategy as Strategy

Strategy = Strategy
Domain = TypeVar('Domain')
Range = TypeVar('Range')
Permutation = Sequence[int]


def apply(function: Callable[..., Range], args: Tuple[Domain, ...]) -> Range:
    return function(*args)


def equivalence(left_statement: bool, right_statement: bool) -> bool:
    return left_statement is right_statement


def implication(antecedent: bool, consequent: bool) -> bool:
    return not antecedent or consequent


def nth_permutation(index: int, size: int) -> Permutation:
    permutations_count = math.factorial(size)
    index %= permutations_count
    indices = list(range(size))
    result = []
    for rest_size in range(size, 0, -1):
        permutations_count //= rest_size
        step, index = divmod(index, permutations_count)
        result.append(indices.pop(step))
    return result


def pack(function: Callable[..., Range]
         ) -> Callable[[Iterable[Domain]], Range]:
    return partial(apply, function)


def permute(sequence: Sequence[Domain], index: int) -> Sequence[Domain]:
    return [sequence[index] for index in nth_permutation(index, len(sequence))]


skip_reference_counter_test = pytest.mark.skipif(
        platform.python_implementation() == 'PyPy',
        reason='PyPy\'s garbage collection '
               'is not based on reference counting.')
