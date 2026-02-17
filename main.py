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
    python main.py <logfile> --frag_a 0-31 --frag_b 31-84

Single fragment (other inferred as complement):
    python main.py <logfile> --frag_a 0-31

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
    freq_to_vib_temp, compute_fragment_rot_temps,
)
from output import (
    print_parse_summary,
    print_gaussian_comparison,
    run_validation_checks,
    print_validation_test,
    print_two_fragment_results,
)


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


def validate_fragments(idx_a, idx_b, n_atoms):
    """Validate fragment indices. Returns a list of error messages (empty if valid)."""
    errors = []
    set_a = set(idx_a)
    set_b = set(idx_b)
    all_indices = set(range(n_atoms))

    oob_a = set_a - all_indices
    oob_b = set_b - all_indices
    if oob_a or oob_b:
        errors.append("ERROR: Fragment indices are out of range.")
        errors.append(f"       This molecule has {n_atoms} atoms (valid indices: 0-{n_atoms - 1}).")
        if oob_a:
            errors.append(f"       Fragment A contains invalid indices: {sorted(oob_a)}")
        if oob_b:
            errors.append(f"       Fragment B contains invalid indices: {sorted(oob_b)}")
        return errors

    overlap = set_a & set_b
    if overlap:
        errors.append("ERROR: Fragments A and B share atoms.")
        errors.append(f"       Overlapping indices: {sorted(overlap)}")
        errors.append("       Each atom must belong to exactly one fragment.")
        return errors

    missing = all_indices - set_a - set_b
    if missing:
        errors.append("ERROR: Not all atoms are assigned to a fragment.")
        errors.append(f"       Missing indices: {sorted(missing)}")
        errors.append(f"       Fragments cover {len(set_a) + len(set_b)} of {n_atoms} atoms.")
        errors.append("       Every atom must belong to either fragment A or fragment B.")

    return errors


def parse_args():
    """Parse command-line arguments."""
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
    parser.add_argument('--override', action='store_true',
                        help='Continue even if validation against Gaussian fails')
    return parser.parse_args()


def run_standard_thermochem(data):
    """Compute standard thermochemistry from parsed data.

    Returns (result, real_freqs, E_elec, our_vib_temps).
    """
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

    E_elec = data.scf_energy
    our_vib_temps = [freq_to_vib_temp(f) for f in real_freqs]
    print_gaussian_comparison(result, data, E_elec, our_vib_temps)

    return result, real_freqs, E_elec, our_vib_temps


def run_validation(result, data, E_elec, our_vib_temps, args):
    """Run validation checks, exit if --validate or if failures without --override."""
    tol = 5e-4
    checks, failures = run_validation_checks(result, data, E_elec, our_vib_temps, tol)

    if args.validate:
        print_validation_test(checks, tol)
        if not failures:
            print("  All validation checks PASSED.")
            sys.exit(0)
        else:
            print("  Some validation checks FAILED.")
            sys.exit(1)

    if failures and not args.override:
        print("ERROR: Standard thermochemistry does not match Gaussian output.")
        print("       The following checks exceeded the tolerance of {:.0e} Hartree:".format(tol))
        for label, diff, unit in failures:
            print(f"         - {label}: off by {diff:.2e} {unit}")
        print()
        print("       This likely means the log file was not parsed correctly.")
        print("       Use --override to continue anyway, or --validate for details.")
        sys.exit(1)


def _format_index_range(indices):
    """Format a sorted list of indices as a compact range string (e.g. '31-83')."""
    if not indices:
        return "(none)"
    if indices[-1] - indices[0] + 1 == len(indices):
        return f"{indices[0]}-{indices[-1]}"
    return ",".join(str(i) for i in indices)


def resolve_fragments(args, n_atoms):
    """Parse fragment args, infer missing fragment, validate.

    Returns (idx_a, idx_b) or exits on error.
    """
    all_indices = set(range(n_atoms))

    if args.frag_a and args.frag_b:
        idx_a = parse_index_spec(args.frag_a)
        idx_b = parse_index_spec(args.frag_b)
    elif args.frag_a:
        idx_a = parse_index_spec(args.frag_a)
        idx_b = sorted(all_indices - set(idx_a))
        print(f"  Fragment B inferred: atoms {_format_index_range(idx_b)} ({len(idx_b)} atoms)")
    else:
        idx_b = parse_index_spec(args.frag_b)
        idx_a = sorted(all_indices - set(idx_b))
        print(f"  Fragment A inferred: atoms {_format_index_range(idx_a)} ({len(idx_a)} atoms)")

    frag_errors = validate_fragments(idx_a, idx_b, n_atoms)
    if frag_errors:
        for msg in frag_errors:
            print(msg)
        sys.exit(1)

    return idx_a, idx_b


def main():
    args = parse_args()

    # --- Parse log file ---
    print(f"Parsing: {args.logfile}")
    data = parse_log(args.logfile)
    print_parse_summary(data)

    # --- Standard thermochemistry ---
    result, real_freqs, E_elec, our_vib_temps = run_standard_thermochem(data)

    # --- Validation ---
    run_validation(result, data, E_elec, our_vib_temps, args)

    # --- Two-fragment correction ---
    if args.frag_a or args.frag_b:
        idx_a, idx_b = resolve_fragments(args, len(data.atom_numbers))

        rot_a, moments_a, mass_a = compute_fragment_rot_temps(
            data.coordinates, data.atom_masses, idx_a)
        rot_b, moments_b, mass_b = compute_fragment_rot_temps(
            data.coordinates, data.atom_masses, idx_b)

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

        print_two_fragment_results(result, result2, E_elec, {
            'idx_a': idx_a, 'idx_b': idx_b,
            'mass_a': mass_a, 'mass_b': mass_b,
            'rot_a': rot_a, 'rot_b': rot_b,
            'moments_a': moments_a, 'moments_b': moments_b,
            'n_remove': args.n_remove,
        })


if __name__ == '__main__':
    main()
