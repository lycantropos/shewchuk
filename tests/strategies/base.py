from hypothesis import strategies

from tests.utils import MAX_VALUE

finite_floats = strategies.floats(-float(MAX_VALUE), float(MAX_VALUE),
                                  allow_infinity=False,
                                  allow_nan=False)
floats_quadruplets = strategies.lists(finite_floats,
                                      min_size=4,
                                      max_size=4)
floats_sextuplets = strategies.lists(finite_floats,
                                     min_size=6,
                                     max_size=6)
floats_octuplets = strategies.lists(finite_floats,
                                    min_size=8,
                                    max_size=8)
