from hypothesis import strategies

finite_floats = strategies.floats(allow_infinity=False,
                                  allow_nan=False)
