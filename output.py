"""
Output formatting functions for the Free Energy Calculator.

All print/display logic lives here. Each function takes computed data in,
prints formatted output, and returns nothing (unless noted).
"""

import sys

from thermochem import kcal_to_hartree, freq_to_vib_temp


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def print_comparison_table(label, our_value, gauss_value, unit=''):
    """Print a single comparison line."""
    diff = our_value - gauss_value
    print(f"  {label:<45s} {our_value:>14.6f}  {gauss_value:>14.6f}  {diff:>+12.6f}  {unit}")


# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------

def print_parse_summary(data):
    """Print the post-parse summary block."""
    print(f"  Atoms: {len(data.atom_numbers)}")
    print(f"  Temperature: {data.temperature:.3f} K")
    print(f"  Pressure: {data.pressure:.5f} atm")
    print(f"  Molecular mass: {data.molecular_mass:.5f} amu")
    print(f"  SCF energy: {data.scf_energy:.8f} Hartree")
    print(f"  Frequencies: {len(data.frequencies)} total, {data.n_imaginary} imaginary")
    print(f"  Rotational temperatures: {data.rot_temperatures}")
    print(f"  Rotational symmetry number: {data.rot_symmetry_number}")
    print()


def print_vib_temp_comparison(our_vib_temps, gauss_vib_temps):
    """Print a table comparing computed vs Gaussian vibrational temperatures."""
    if not gauss_vib_temps:
        print("  (No vibrational temperatures found in Gaussian output)")
        print()
        return

    n = min(len(our_vib_temps), len(gauss_vib_temps))
    print("  Vibrational temperatures (K):")
    print(f"  {'Mode':<8s} {'Ours':>12s}  {'Gaussian':>12s}  {'Difference':>12s}")
    print(f"  {'-'*8} {'-'*12}  {'-'*12}  {'-'*12}")
    max_diff = 0.0
    for i in range(n):
        diff = our_vib_temps[i] - gauss_vib_temps[i]
        max_diff = max(max_diff, abs(diff))
        print(f"  {i+1:<8d} {our_vib_temps[i]:>12.2f}  {gauss_vib_temps[i]:>12.2f}  {diff:>+12.4f}")
    print()
    print(f"  Max vibrational temperature difference: {max_diff:.4f} K")
    if len(our_vib_temps) != len(gauss_vib_temps):
        print(f"  WARNING: mode count mismatch (ours: {len(our_vib_temps)}, Gaussian: {len(gauss_vib_temps)})")
    print()


def print_gaussian_comparison(result, data, E_elec, our_vib_temps):
    """Print the full E/Cv/S tables, corrections, final energies, ZPE, and vib temps."""
    print("=" * 100)
    print("VALIDATION: Comparison with Gaussian output")
    print("=" * 100)
    print()

    # --- E (Thermal) ---
    print("  E (Thermal) in kcal/mol:")
    print(f"  {'Component':<20s} {'Ours':>14s}  {'Gaussian':>14s}  {'Difference':>12s}")
    print(f"  {'-'*20} {'-'*14}  {'-'*14}  {'-'*12}")
    for comp, ours, gauss in [
        ('Electronic', result['E_elec'], data.e_electronic),
        ('Translational', result['E_trans'], data.e_translational),
        ('Rotational', result['E_rot'], data.e_rotational),
        ('Vibrational', result['E_vib'], data.e_vibrational),
        ('Total', result['E_total'], data.e_total),
    ]:
        diff = ours - gauss
        print(f"  {comp:<20s} {ours:>14.3f}  {gauss:>14.3f}  {diff:>+12.6f}")
    print()

    # --- Cv ---
    print("  Cv in cal/(mol·K):")
    print(f"  {'Component':<20s} {'Ours':>14s}  {'Gaussian':>14s}  {'Difference':>12s}")
    print(f"  {'-'*20} {'-'*14}  {'-'*14}  {'-'*12}")
    for comp, ours, gauss in [
        ('Electronic', result['Cv_elec'], data.cv_electronic),
        ('Translational', result['Cv_trans'], data.cv_translational),
        ('Rotational', result['Cv_rot'], data.cv_rotational),
        ('Vibrational', result['Cv_vib'], data.cv_vibrational),
        ('Total', result['Cv_total'], data.cv_total),
    ]:
        diff = ours - gauss
        print(f"  {comp:<20s} {ours:>14.3f}  {gauss:>14.3f}  {diff:>+12.6f}")
    print()

    # --- S ---
    print("  S in cal/(mol·K):")
    print(f"  {'Component':<20s} {'Ours':>14s}  {'Gaussian':>14s}  {'Difference':>12s}")
    print(f"  {'-'*20} {'-'*14}  {'-'*14}  {'-'*12}")
    for comp, ours, gauss in [
        ('Electronic', result['S_elec'], data.s_electronic),
        ('Translational', result['S_trans'], data.s_translational),
        ('Rotational', result['S_rot'], data.s_rotational),
        ('Vibrational', result['S_vib'], data.s_vibrational),
        ('Total', result['S_total'], data.s_total),
    ]:
        diff = ours - gauss
        print(f"  {comp:<20s} {ours:>14.3f}  {gauss:>14.3f}  {diff:>+12.6f}")
    print()

    # --- Thermal corrections ---
    print("  Thermal corrections (Hartree):")
    print(f"  {'Quantity':<45s} {'Ours':>14s}  {'Gaussian':>14s}  {'Difference':>12s}")
    print(f"  {'-'*45} {'-'*14}  {'-'*14}  {'-'*12}")
    for label, ours, gauss in [
        ('Zero-point correction', result['ZPE_hartree'], data.zpe_correction),
        ('Thermal correction to Energy', result['thermal_correction_energy'], data.thermal_correction_energy),
        ('Thermal correction to Enthalpy', result['thermal_correction_enthalpy'], data.thermal_correction_enthalpy),
        ('Thermal correction to Gibbs', result['thermal_correction_gibbs'], data.thermal_correction_gibbs),
    ]:
        diff = ours - gauss
        print(f"  {label:<45s} {ours:>14.6f}  {gauss:>14.6f}  {diff:>+12.6f}")
    print()

    # --- Final energies ---
    print("  Final energies (Hartree):")
    print(f"  {'Quantity':<45s} {'Ours':>14s}  {'Gaussian':>14s}  {'Difference':>12s}")
    print(f"  {'-'*45} {'-'*14}  {'-'*14}  {'-'*12}")
    for label, correction, gauss_sum in [
        ('E_elec + ZPE', result['ZPE_hartree'], data.sum_elec_zpe),
        ('E_elec + E_thermal', result['thermal_correction_energy'], data.sum_elec_thermal_e),
        ('E_elec + H', result['thermal_correction_enthalpy'], data.sum_elec_thermal_h),
        ('E_elec + G', result['thermal_correction_gibbs'], data.sum_elec_thermal_g),
    ]:
        ours = E_elec + correction
        diff = ours - gauss_sum
        print(f"  {label:<45s} {ours:>14.6f}  {gauss_sum:>14.6f}  {diff:>+12.6f}")
    print()

    # --- ZPE ---
    print(f"  ZPE: {result['ZPE_kcal']:.5f} kcal/mol  (Gaussian: {data.zpe_kcal_per_mol:.5f})")
    print()

    # --- Vibrational temperatures ---
    print_vib_temp_comparison(our_vib_temps, data.vib_temperatures)


def run_validation_checks(result, data, E_elec, our_vib_temps, tol):
    """Build validation checks and identify failures.

    Returns
    -------
    checks : list of (label, ours, ref, unit)
    failures : list of (label, diff, unit)  — subset that exceeded tolerance.
    """
    checks = [
        ('ZPE correction', result['ZPE_hartree'], data.zpe_correction, 'Hartree'),
        ('Thermal correction to Energy', result['thermal_correction_energy'], data.thermal_correction_energy, 'Hartree'),
        ('Thermal correction to Enthalpy', result['thermal_correction_enthalpy'], data.thermal_correction_enthalpy, 'Hartree'),
        ('Thermal correction to Gibbs', result['thermal_correction_gibbs'], data.thermal_correction_gibbs, 'Hartree'),
        ('E_elec + ZPE', E_elec + result['ZPE_hartree'], data.sum_elec_zpe, 'Hartree'),
        ('E_elec + G', E_elec + result['thermal_correction_gibbs'], data.sum_elec_thermal_g, 'Hartree'),
    ]

    # Vibrational temperature check
    vib_temp_tol = 0.05  # K
    if our_vib_temps and data.vib_temperatures:
        n = min(len(our_vib_temps), len(data.vib_temperatures))
        max_diff = max(abs(our_vib_temps[i] - data.vib_temperatures[i]) for i in range(n))
        checks.append(('Vibrational temperatures (max diff)', max_diff, 0.0, 'K'))

    failures = []
    for entry in checks:
        label, ours, ref, unit = entry
        diff = abs(ours - ref)
        entry_tol = vib_temp_tol if unit == 'K' else tol
        if diff >= entry_tol:
            failures.append((label, diff, unit))
    return checks, failures


def print_validation_test(checks, tol):
    """Print the PASS/FAIL table for ``--validate``."""
    vib_temp_tol = 0.05  # K
    print("=" * 70)
    print("VALIDATION TEST")
    print("=" * 70)
    for label, ours, ref, unit in checks:
        diff = abs(ours - ref)
        entry_tol = vib_temp_tol if unit == 'K' else tol
        status = "PASS" if diff < entry_tol else "FAIL"
        print(f"  {status}  {label:<40s}  diff = {diff:.2e} {unit}")
    print()


def print_two_fragment_results(result_std, result_2frag, E_elec, frag_info, verbose=False):
    """Print full two-fragment comparison output.

    Parameters
    ----------
    result_std : dict — standard thermochem result
    result_2frag : dict — two-fragment thermochem result
    E_elec : float — electronic energy in Hartree
    frag_info : dict with keys:
        idx_a, idx_b, mass_a, mass_b, rot_a, rot_b, moments_a, moments_b,
        n_remove, sigma_a, sigma_b
    """
    idx_a = frag_info['idx_a']
    idx_b = frag_info['idx_b']
    mass_a = frag_info['mass_a']
    mass_b = frag_info['mass_b']
    rot_a = frag_info['rot_a']
    rot_b = frag_info['rot_b']
    moments_a = frag_info['moments_a']
    moments_b = frag_info['moments_b']
    n_remove = frag_info['n_remove']

    print()
    print("=" * 100)
    print("TWO-FRAGMENT CORRECTED THERMOCHEMISTRY")
    print("=" * 100)
    print()

    print(f"  Fragment A: atoms {idx_a[0]}-{idx_a[-1]} ({len(idx_a)} atoms)")
    print(f"  Fragment B: atoms {idx_b[0]}-{idx_b[-1]} ({len(idx_b)} atoms)")
    print()

    print(f"  Fragment A: mass = {mass_a:.5f} amu")
    print(f"    Rotational temperatures: {rot_a[0]:.5f}  {rot_a[1]:.5f}  {rot_a[2]:.5f} K")
    print(f"    Principal moments (amu·A²): {moments_a[0]:.3f}  {moments_a[1]:.3f}  {moments_a[2]:.3f}")
    print()
    print(f"  Fragment B: mass = {mass_b:.5f} amu")
    print(f"    Rotational temperatures: {rot_b[0]:.5f}  {rot_b[1]:.5f}  {rot_b[2]:.5f} K")
    print(f"    Principal moments (amu·A²): {moments_b[0]:.3f}  {moments_b[1]:.3f}  {moments_b[2]:.3f}")
    print()

    print(f"  Removed {n_remove} lowest frequencies (cm⁻¹):")
    for i, f in enumerate(result_2frag['removed_freqs']):
        print(f"    {i+1}. {f:.4f}")
    print(f"  Remaining vibrational modes: {len(result_2frag['kept_freqs'])}")
    print()

    # --- E/Cv/S comparison tables (verbose only) ---
    result = result_std
    result2 = result_2frag

    if verbose:
        print("  E (Thermal) in kcal/mol:")
        print(f"  {'Component':<20s} {'Standard':>14s}  {'Two-Fragment':>14s}  {'Difference':>12s}")
        print(f"  {'-'*20} {'-'*14}  {'-'*14}  {'-'*12}")
        for comp in ['E_trans', 'E_rot', 'E_vib', 'E_elec', 'E_total']:
            label = comp.replace('E_', '').capitalize()
            v1, v2 = result[comp], result2[comp]
            print(f"  {label:<20s} {v1:>14.3f}  {v2:>14.3f}  {v2-v1:>+12.3f}")
        print()

        print("  Cv in cal/(mol·K):")
        print(f"  {'Component':<20s} {'Standard':>14s}  {'Two-Fragment':>14s}  {'Difference':>12s}")
        print(f"  {'-'*20} {'-'*14}  {'-'*14}  {'-'*12}")
        for comp in ['Cv_trans', 'Cv_rot', 'Cv_vib', 'Cv_elec', 'Cv_total']:
            label = comp.replace('Cv_', '').capitalize()
            v1, v2 = result[comp], result2[comp]
            print(f"  {label:<20s} {v1:>14.3f}  {v2:>14.3f}  {v2-v1:>+12.3f}")
        print()

        print("  S in cal/(mol·K):")
        print(f"  {'Component':<20s} {'Standard':>14s}  {'Two-Fragment':>14s}  {'Difference':>12s}")
        print(f"  {'-'*20} {'-'*14}  {'-'*14}  {'-'*12}")
        for comp in ['S_trans', 'S_rot', 'S_vib', 'S_elec', 'S_total']:
            label = comp.replace('S_', '').capitalize()
            v1, v2 = result[comp], result2[comp]
            print(f"  {label:<20s} {v1:>14.3f}  {v2:>14.3f}  {v2-v1:>+12.3f}")
        print()

    # --- Corrections ---
    print("  Thermal corrections (Hartree):")
    print(f"  {'Quantity':<45s} {'Standard':>14s}  {'Two-Fragment':>14s}  {'Difference':>12s}")
    print(f"  {'-'*45} {'-'*14}  {'-'*14}  {'-'*12}")
    for label, key in [
        ('Zero-point correction', 'ZPE_hartree'),
        ('Thermal correction to Energy', 'thermal_correction_energy'),
        ('Thermal correction to Enthalpy', 'thermal_correction_enthalpy'),
        ('Thermal correction to Gibbs', 'thermal_correction_gibbs'),
    ]:
        v1, v2 = result[key], result2[key]
        print(f"  {label:<45s} {v1:>14.6f}  {v2:>14.6f}  {v2-v1:>+12.6f}")
    print()

    # --- Final energies ---
    print("  Final energies (Hartree):")
    print(f"  {'Quantity':<45s} {'Standard':>14s}  {'Two-Fragment':>14s}  {'Difference':>12s}")
    print(f"  {'-'*45} {'-'*14}  {'-'*14}  {'-'*12}")
    for label, key in [
        ('E_elec + ZPE', 'ZPE_hartree'),
        ('E_elec + E_thermal', 'thermal_correction_energy'),
        ('E_elec + H', 'thermal_correction_enthalpy'),
        ('E_elec + G', 'thermal_correction_gibbs'),
    ]:
        v1 = E_elec + result[key]
        v2 = E_elec + result2[key]
        print(f"  {label:<45s} {v1:>14.6f}  {v2:>14.6f}  {v2-v1:>+12.6f}")
    print()

    # --- Key result ---
    G_standard = E_elec + result['thermal_correction_gibbs']
    G_corrected = E_elec + result2['thermal_correction_gibbs']
    print(f"  *** E_elec + G (standard):      {G_standard:.6f} Hartree")
    print(f"  *** E_elec + G (two-fragment):   {G_corrected:.6f} Hartree")
    print(f"  *** Correction:                  {(G_corrected - G_standard):.6f} Hartree")
    print(f"  ***                              {(G_corrected - G_standard) / kcal_to_hartree:.3f} kcal/mol")
    print()
