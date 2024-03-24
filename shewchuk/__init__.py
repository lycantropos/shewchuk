"""Robust floating point operations."""

from __future__ import annotations

import typing as _t

__version__ = '6.9.0'

if _t.TYPE_CHECKING:
    from . import _shewchuk

    Expansion = _shewchuk.Expansion
    incircle_test = _shewchuk.incircle_test
    kind = _shewchuk.kind
    orientation = _shewchuk.orientation
    vectors_cross_product = _shewchuk.vectors_cross_product
    vectors_dot_product = _shewchuk.vectors_dot_product
else:
    try:
        from . import _cshewchuk
    except ImportError:
        from . import _shewchuk

        Expansion = _shewchuk.Expansion
        incircle_test = _shewchuk.incircle_test
        kind = _shewchuk.kind
        orientation = _shewchuk.orientation
        vectors_cross_product = _shewchuk.vectors_cross_product
        vectors_dot_product = _shewchuk.vectors_dot_product
    else:
        Expansion = _cshewchuk.Expansion
        incircle_test = _cshewchuk.incircle_test
        kind = _cshewchuk.kind
        orientation = _cshewchuk.orientation
        vectors_cross_product = _cshewchuk.vectors_cross_product
        vectors_dot_product = _cshewchuk.vectors_dot_product
