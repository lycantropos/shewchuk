"""Robust floating point operations."""

__version__ = '5.0.0'

try:
    from ._c import (Expansion,
                     incircle_test,
                     kind,
                     orientation,
                     vectors_cross_product,
                     vectors_dot_product)
except ImportError:
    from ._py import (Expansion,
                      incircle_test,
                      kind,
                      orientation,
                      vectors_cross_product,
                      vectors_dot_product)
