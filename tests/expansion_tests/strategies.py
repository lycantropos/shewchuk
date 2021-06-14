from hypothesis import strategies

from tests.strategies import finite_floats

finite_floats_sequences = strategies.lists(finite_floats)
