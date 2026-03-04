"""
two_fragment_free_energy
========================
Package for computing Gibbs free energies from Gaussian frequency calculations
with an optional two-fragment correction for transition states involving two
loosely-associated molecular fragments.

Public API
----------
    from two_fragment_free_energy import two_fragment_correction

    result = two_fragment_correction("path/to/file.log", frag_a="0-31")
    print(result["G_corrected"])   # Hartree
    print(result["correction_kcal"])
"""

import os as _os
import sys as _sys

# Ensure the package directory is on sys.path so the internal modules
# (api, parse_gaussian, thermochem, etc.) can resolve their own absolute
# imports regardless of how this package was found.
_pkg_dir = _os.path.dirname(_os.path.abspath(__file__))
if _pkg_dir not in _sys.path:
    _sys.path.insert(0, _pkg_dir)

from api import two_fragment_correction  # noqa: E402

del _pkg_dir

__all__ = ["two_fragment_correction"]
