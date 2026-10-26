"""Force-field adapter for the self-contained cg2aa pipeline.

The cg2aa pipeline uses bundled Multiscale OPLS implementation and package assets
under :mod:`multiscale.forcefield`.
"""

from __future__ import annotations

from typing import Any

from multiscale.forcefield.assign import ff as _internal_opls_ff


def ff(molecule: Any, **kwargs: Any):
    """Assign force-field parameters using bundled Multiscale OPLS assets."""
    return _internal_opls_ff(molecule, **kwargs)
