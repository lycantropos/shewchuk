from hypothesis import strategies

MAX_VALUE = 10 ** 50
finite_floats = strategies.floats(-float(MAX_VALUE), float(MAX_VALUE),
                                  allow_infinity=False,
                                  allow_nan=False)
