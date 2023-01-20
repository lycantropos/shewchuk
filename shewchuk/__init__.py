"""Robust floating point operations."""

__version__ = '6.6.0'

try:
    from ._cshewchuk import (Expansion,
                             incircle_test,
                             kind,
                             orientation,
                             vectors_cross_product,
                             vectors_dot_product)
except ImportError:
    from ._shewchuk import (Expansion,
                            incircle_test,
                            kind,
                            orientation,
                            vectors_cross_product,
                            vectors_dot_product)
