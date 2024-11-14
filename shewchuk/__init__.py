"""Robust floating point operations."""

from __future__ import annotations

from typing import TYPE_CHECKING

__version__ = '6.9.0'

if TYPE_CHECKING:
    from ._shewchuk import (
        Expansion,
        incircle_test,
        kind,
        orientation,
        vectors_cross_product,
        vectors_dot_product,
    )
else:
    try:
        from . import _cshewchuk
    except ImportError:
        from ._shewchuk import (
            Expansion,
            incircle_test,
            kind,
            orientation,
            vectors_cross_product,
            vectors_dot_product,
        )
    else:
        Expansion = _cshewchuk.Expansion
        incircle_test = _cshewchuk.incircle_test
        kind = _cshewchuk.kind
        orientation = _cshewchuk.orientation
        vectors_cross_product = _cshewchuk.vectors_cross_product
        vectors_dot_product = _cshewchuk.vectors_dot_product
