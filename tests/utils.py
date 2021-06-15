import platform
from functools import partial
from typing import (Callable,
                    Iterable,
                    Tuple,
                    TypeVar)

import pytest
from hypothesis.strategies import SearchStrategy as Strategy

Strategy = Strategy
Domain = TypeVar('Domain')
Range = TypeVar('Range')


def equivalence(left_statement: bool, right_statement: bool) -> bool:
    return left_statement is right_statement


def implication(antecedent: bool, consequent: bool) -> bool:
    return not antecedent or consequent


def pack(function: Callable[..., Range]
         ) -> Callable[[Iterable[Domain]], Range]:
    return partial(apply, function)


def apply(function: Callable[..., Range], args: Tuple[Domain, ...]) -> Range:
    return function(*args)


skip_reference_counter_test = pytest.mark.skipif(
        platform.python_implementation() == 'PyPy',
        reason='PyPy\'s garbage collection '
               'is not based on reference counting.')
