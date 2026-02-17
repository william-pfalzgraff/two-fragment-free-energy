"""
Public library API for the Free Energy Calculator.

Usage
-----
    from api import two_fragment_correction

    result = two_fragment_correction(
        logfile="path/to/file.log",
        frag_a="0-31",
    )
    print(result["correction_kcal"])
"""

import warnings

from parse_gaussian import parse_log
from thermochem import (
    compute_thermochem, compute_two_fragment_thermochem,
    compute_fragment_rot_temps, kcal_to_hartree,
)
from main import parse_index_spec, validate_fragments, run_standard_thermochem
from output import run_validation_checks


def two_fragment_correction(
    logfile,
    frag_a=None,
    frag_b=None,
    sigma_a=1,
    sigma_b=1,
    validate=True,
):
    """Compute the two-fragment free energy correction for a Gaussian log file.

    Parameters
    ----------
    logfile : str
        Path to a Gaussian frequency log file.
    frag_a : str or list[int] or None
        Atom indices for fragment A. String specs like ``"0-31"`` or
        ``"0,1,2,3"`` are parsed automatically. At least one of
        *frag_a* / *frag_b* must be provided.
    frag_b : str or list[int] or None
        Atom indices for fragment B. If omitted, inferred as the
        complement of *frag_a*.
    sigma_a : int
        Rotational symmetry number for fragment A (default 1).
    sigma_b : int
        Rotational symmetry number for fragment B (default 1).
    validate : bool
        If True (default), run validation checks against Gaussian output
        and issue warnings on failure. If False, skip validation.

    Returns
    -------
    dict
        Keys: ``G_standard``, ``G_corrected``, ``correction_hartree``,
        ``correction_kcal``, ``E_elec``, ``standard``, ``two_fragment``.

    Raises
    ------
    ValueError
        If fragment indices are invalid (out of range, overlapping, or
        incomplete coverage).
    """
    if frag_a is None and frag_b is None:
        raise ValueError("At least one of frag_a or frag_b must be provided.")

    # --- Parse log file ---
    data = parse_log(logfile)

    # --- Standard thermochemistry (silent) ---
    result_std, real_freqs, E_elec, our_vib_temps = run_standard_thermochem(
        data, verbosity=0
    )

    # --- Validation ---
    if validate:
        tol = 5e-4
        _checks, failures = run_validation_checks(
            result_std, data, E_elec, our_vib_temps, tol
        )
        if failures:
            details = "; ".join(
                f"{label}: off by {diff:.2e} {unit}"
                for label, diff, unit in failures
            )
            warnings.warn(
                f"Validation against Gaussian output failed: {details}",
                stacklevel=2,
            )

    # --- Resolve fragments ---
    n_atoms = len(data.atom_numbers)
    all_indices = set(range(n_atoms))

    idx_a = _normalize_fragment(frag_a) if frag_a is not None else None
    idx_b = _normalize_fragment(frag_b) if frag_b is not None else None

    if idx_a is not None and idx_b is not None:
        pass  # both provided
    elif idx_a is not None:
        idx_b = sorted(all_indices - set(idx_a))
    else:
        idx_a = sorted(all_indices - set(idx_b))

    frag_errors = validate_fragments(idx_a, idx_b, n_atoms)
    if frag_errors:
        raise ValueError("\n".join(frag_errors))

    # --- Fragment rotational temperatures ---
    rot_a, moments_a, mass_a = compute_fragment_rot_temps(
        data.coordinates, data.atom_masses, idx_a
    )
    rot_b, moments_b, mass_b = compute_fragment_rot_temps(
        data.coordinates, data.atom_masses, idx_b
    )

    # --- Two-fragment thermochemistry ---
    result_2frag = compute_two_fragment_thermochem(
        freqs_all_real=real_freqs,
        mass_A=mass_a,
        mass_B=mass_b,
        rot_temps_A=rot_a,
        rot_temps_B=rot_b,
        sigma_A=sigma_a,
        sigma_B=sigma_b,
        T=data.temperature,
        P_atm=data.pressure,
        n_remove=6,
    )

    # --- Build result ---
    G_standard = E_elec + result_std['thermal_correction_gibbs']
    G_corrected = E_elec + result_2frag['thermal_correction_gibbs']
    correction_hartree = G_corrected - G_standard
    correction_kcal = correction_hartree / kcal_to_hartree

    return {
        "G_standard": G_standard,
        "G_corrected": G_corrected,
        "correction_hartree": correction_hartree,
        "correction_kcal": correction_kcal,
        "E_elec": E_elec,
        "standard": result_std,
        "two_fragment": result_2frag,
    }


def _normalize_fragment(frag):
    """Convert a fragment specification to a sorted list of ints."""
    if isinstance(frag, str):
        return parse_index_spec(frag)
    return sorted(set(int(i) for i in frag))
