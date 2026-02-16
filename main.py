#!/usr/bin/env python3
"""
Free Energy Calculator — Main driver script.

Parses a Gaussian frequency log file, reproduces standard IGRRHO
thermochemistry, and optionally computes two-fragment corrected
free energies.

Usage
-----
Standard validation:
    python main.py <logfile>

Two-fragment correction:
    python main.py <logfile> --frag_a 0-4 --frag_b 5-83

Fragment indices are 0-based and can be specified as:
    --frag_a 0,1,2,3,4      (comma-separated)
    --frag_a 0-4             (range, inclusive)
    --frag_a 0-4,10,12-15   (mixed)
"""

import argparse
import sys

from parse_gaussian import parse_log
from thermochem import (
    compute_thermochem, compute_two_fragment_thermochem,
    freq_to_vib_temp, kcal_to_hartree
)
from inertia import compute_fragment_rot_temps


def parse_index_spec(spec):
    """Parse an index specification like '0-4,10,12-15' into a list of ints."""
    indices = []
    for part in spec.split(','):
        part = part.strip()
        if '-' in part:
            start, end = part.split('-', 1)
            indices.extend(range(int(start), int(end) + 1))
        else:
            indices.append(int(part))
    return sorted(set(indices))


def print_comparison_table(label, our_value, gauss_value, unit=''):
    """Print a single comparison line."""
    diff = our_value - gauss_value
    print(f"  {label:<45s} {our_value:>14.6f}  {gauss_value:>14.6f}  {diff:>+12.6f}  {unit}")


def main():
    parser = argparse.ArgumentParser(description='Free Energy Calculator')
    parser.add_argument('logfile', help='Path to Gaussian .log file')
    parser.add_argument('--frag_a', type=str, default=None,
                        help='Atom indices for fragment A (0-based)')
    parser.add_argument('--frag_b', type=str, default=None,
                        help='Atom indices for fragment B (0-based)')
    parser.add_argument('--sigma_a', type=int, default=1,
                        help='Rotational symmetry number for fragment A')
    parser.add_argument('--sigma_b', type=int, default=1,
                        help='Rotational symmetry number for fragment B')
    parser.add_argument('--n_remove', type=int, default=6,
                        help='Number of lowest frequencies to remove (default: 6)')
    parser.add_argument('--validate', action='store_true',
                        help='Run validation test against Gaussian output and exit with status')
    args = parser.parse_args()

    # ---------------------------------------------------------------
    # Parse log file
    # ---------------------------------------------------------------
    print(f"Parsing: {args.logfile}")
    data = parse_log(args.logfile)
    print(f"  Atoms: {len(data.atom_numbers)}")
    print(f"  Temperature: {data.temperature:.3f} K")
    print(f"  Pressure: {data.pressure:.5f} atm")
    print(f"  Molecular mass: {data.molecular_mass:.5f} amu")
    print(f"  SCF energy: {data.scf_energy:.8f} Hartree")
    print(f"  Frequencies: {len(data.frequencies)} total, {data.n_imaginary} imaginary")
    print(f"  Rotational temperatures: {data.rot_temperatures}")
    print(f"  Rotational symmetry number: {data.rot_symmetry_number}")
    print()

    # ---------------------------------------------------------------
    # Standard thermochemistry (Gaussian-equivalent)
    # ---------------------------------------------------------------
    # Get real frequencies only (exclude imaginary)
    real_freqs = [f for f in data.frequencies if f > 0]
    print(f"Real vibrational modes: {len(real_freqs)}")
    print()

    result = compute_thermochem(
        freqs=real_freqs,
        mass_amu=data.molecular_mass,
        rot_temps=data.rot_temperatures,
        sigma=data.rot_symmetry_number,
        T=data.temperature,
        P_atm=data.pressure,
    )

    # ---------------------------------------------------------------
    # Print E, Cv, S comparison
    # ---------------------------------------------------------------
    print("=" * 100)
    print("VALIDATION: Comparison with Gaussian output")
    print("=" * 100)
    print()

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

    # ---------------------------------------------------------------
    # Thermal corrections comparison
    # ---------------------------------------------------------------
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

    # ---------------------------------------------------------------
    # Final energies comparison
    # ---------------------------------------------------------------
    print("  Final energies (Hartree):")
    print(f"  {'Quantity':<45s} {'Ours':>14s}  {'Gaussian':>14s}  {'Difference':>12s}")
    print(f"  {'-'*45} {'-'*14}  {'-'*14}  {'-'*12}")
    E_elec = data.scf_energy
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

    # ZPE in kcal/mol
    print(f"  ZPE: {result['ZPE_kcal']:.5f} kcal/mol  (Gaussian: {data.zpe_kcal_per_mol:.5f})")
    print()

    # ---------------------------------------------------------------
    # --validate: check results against Gaussian and exit
    # ---------------------------------------------------------------
    if args.validate:
        tol = 5e-4  # tolerance in Hartree (~0.3 kcal/mol)
        checks = [
            ('ZPE correction', result['ZPE_hartree'], data.zpe_correction),
            ('Thermal correction to Energy', result['thermal_correction_energy'], data.thermal_correction_energy),
            ('Thermal correction to Enthalpy', result['thermal_correction_enthalpy'], data.thermal_correction_enthalpy),
            ('Thermal correction to Gibbs', result['thermal_correction_gibbs'], data.thermal_correction_gibbs),
            ('E_elec + ZPE', E_elec + result['ZPE_hartree'], data.sum_elec_zpe),
            ('E_elec + G', E_elec + result['thermal_correction_gibbs'], data.sum_elec_thermal_g),
        ]
        all_pass = True
        print("=" * 70)
        print("VALIDATION TEST")
        print("=" * 70)
        for label, ours, ref in checks:
            diff = abs(ours - ref)
            status = "PASS" if diff < tol else "FAIL"
            if status == "FAIL":
                all_pass = False
            print(f"  {status}  {label:<40s}  diff = {diff:.2e} Hartree")
        print()
        if all_pass:
            print("  All validation checks PASSED.")
            sys.exit(0)
        else:
            print("  Some validation checks FAILED.")
            sys.exit(1)

    # ---------------------------------------------------------------
    # Two-fragment correction
    # ---------------------------------------------------------------
    if args.frag_a and args.frag_b:
        print()
        print("=" * 100)
        print("TWO-FRAGMENT CORRECTED THERMOCHEMISTRY")
        print("=" * 100)
        print()

        idx_a = parse_index_spec(args.frag_a)
        idx_b = parse_index_spec(args.frag_b)

        print(f"  Fragment A: atoms {idx_a[0]}-{idx_a[-1]} ({len(idx_a)} atoms)")
        print(f"  Fragment B: atoms {idx_b[0]}-{idx_b[-1]} ({len(idx_b)} atoms)")
        print()

        # Compute fragment properties
        rot_a, moments_a, mass_a = compute_fragment_rot_temps(
            data.coordinates, data.atom_masses, idx_a)
        rot_b, moments_b, mass_b = compute_fragment_rot_temps(
            data.coordinates, data.atom_masses, idx_b)

        print(f"  Fragment A: mass = {mass_a:.5f} amu")
        print(f"    Rotational temperatures: {rot_a[0]:.5f}  {rot_a[1]:.5f}  {rot_a[2]:.5f} K")
        print(f"    Principal moments (amu·A²): {moments_a[0]:.3f}  {moments_a[1]:.3f}  {moments_a[2]:.3f}")
        print()
        print(f"  Fragment B: mass = {mass_b:.5f} amu")
        print(f"    Rotational temperatures: {rot_b[0]:.5f}  {rot_b[1]:.5f}  {rot_b[2]:.5f} K")
        print(f"    Principal moments (amu·A²): {moments_b[0]:.3f}  {moments_b[1]:.3f}  {moments_b[2]:.3f}")
        print()

        # Run two-fragment thermochemistry
        result2 = compute_two_fragment_thermochem(
            freqs_all_real=real_freqs,
            mass_A=mass_a,
            mass_B=mass_b,
            rot_temps_A=rot_a,
            rot_temps_B=rot_b,
            sigma_A=args.sigma_a,
            sigma_B=args.sigma_b,
            T=data.temperature,
            P_atm=data.pressure,
            n_remove=args.n_remove,
        )

        print(f"  Removed {args.n_remove} lowest frequencies (cm⁻¹):")
        for i, f in enumerate(result2['removed_freqs']):
            print(f"    {i+1}. {f:.4f}")
        print(f"  Remaining vibrational modes: {len(result2['kept_freqs'])}")
        print()

        # Print two-fragment results
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

        # Corrections
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

        # Final energies
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

        # Key result
        G_standard = E_elec + result['thermal_correction_gibbs']
        G_corrected = E_elec + result2['thermal_correction_gibbs']
        print(f"  *** E_elec + G (standard):      {G_standard:.6f} Hartree")
        print(f"  *** E_elec + G (two-fragment):   {G_corrected:.6f} Hartree")
        print(f"  *** Correction:                  {(G_corrected - G_standard):.6f} Hartree")
        print(f"  ***                              {(G_corrected - G_standard) / kcal_to_hartree:.3f} kcal/mol")
        print()


if __name__ == '__main__':
    main()
