# Free Energy Calculator

A Python tool for recomputing Gibbs free energies from Gaussian frequency
calculations, with an optional **two-fragment correction** for transition states
where two molecular fragments are loosely associated.

## Why does this exist?

Gaussian's built-in thermochemistry treats the entire system as a single rigid
body.  That's fine for a normal molecule, but for a transition state where two
fragments are barely held together, the six lowest vibrational modes are really
just the fragments wobbling relative to each other.  Those low-frequency modes
dominate the entropy and can badly distort the free energy.

The two-fragment correction fixes this by:

1. Removing the six lowest vibrational frequencies (the interfragment modes).
2. Treating each fragment as its own independent translator and rotor.

The result is a free energy that better reflects the physics of a dissociative
transition state.

## What's in here

| File | Purpose |
|---|---|
| `main.py` | Command-line driver; parses arguments and prints results |
| `parse_gaussian.py` | Extracts all relevant data from a Gaussian `.log` file |
| `thermochem.py` | IGRRHO thermochemistry (translational, rotational, vibrational, electronic contributions) |
| `inertia.py` | Computes moments of inertia and rotational temperatures for molecular fragments |

## Requirements

- Python 3.7+
- NumPy

No other dependencies.  Install NumPy if you don't have it:

```bash
pip install numpy
```

## Quick start

### 1. Standard thermochemistry (validation mode)

This reproduces Gaussian's own thermochemistry from the `.log` file and compares
every quantity side-by-side.  Use this to confirm the code is working correctly:

```bash
python main.py your_file.log
```

You'll see tables comparing E, Cv, S, thermal corrections, and final energies
against Gaussian's printed values.  Everything should match to within rounding
(~10^-6 Hartree).

### 2. Validate against Gaussian (automated test)

To run a quick pass/fail check that the code reproduces Gaussian's output:

```bash
python main.py your_file.log --validate
```

This checks all thermal corrections and final energies against the values
Gaussian printed in the log file.  It exits with code 0 if everything matches
and code 1 if anything is off.

### 3. Two-fragment corrected free energy

Specify which atoms belong to each fragment using 0-based indices:

```bash
python main.py your_file.log --frag_a 0-30 --frag_b 31-83
```

The index syntax is flexible:
- Ranges: `0-30` (inclusive on both ends)
- Lists: `0,1,2,3,4`
- Mixed: `0-4,10,12-15`

This will print the standard thermochemistry first, then a full comparison
showing how each component (translation, rotation, vibration) changes under
the two-fragment treatment, and finally the corrected E_elec + G.

### Additional options

| Flag | Default | Description |
|---|---|---|
| `--sigma_a` | 1 | Rotational symmetry number for fragment A |
| `--sigma_b` | 1 | Rotational symmetry number for fragment B |
| `--n_remove` | 6 | Number of lowest frequencies to remove |

For most organic molecules with no special symmetry, the defaults are fine.

## How to figure out your fragment indices

Open your `.log` file and look for the "Input orientation" block near the end.
The atoms are listed in order starting from 1 (Gaussian numbering).  Decide
which atoms belong to each fragment, then subtract 1 to convert to 0-based
indices.

For example, if fragment A is Gaussian atoms 1-31, use `--frag_a 0-30`.

The two fragments should account for every atom in the molecule.  The code will
report each fragment's mass, and you can check that they sum to the total
molecular mass as a sanity check.

## Example output (key lines)

```
  *** E_elec + G (standard):      -1626.374561 Hartree
  *** E_elec + G (two-fragment):   -1626.392935 Hartree
  *** Correction:                  -0.018374 Hartree
  ***                              -11.530 kcal/mol
```

The correction is negative (stabilizing) because the two-fragment treatment
captures the additional translational and rotational entropy that the fragments
gain when treated as independent bodies.

## Theory

The standard thermochemistry follows the ideal-gas/rigid-rotor/harmonic-oscillator
(IGRRHO) formalism as described in:

- Ochterski, J. W. "Thermochemistry in Gaussian." (2000).
  Available at https://gaussian.com/thermo/
- McQuarrie, D. A. *Statistical Mechanics*. University Science Books, 2000.

The two-fragment correction replaces:

- 1 translational partition function (total mass) with 2 (one per fragment)
- 1 rotational partition function (whole-molecule inertia tensor) with 2 (one per fragment)
- Removes the *n* lowest vibrational modes that correspond to interfragment motion

All physical constants match the CODATA values used internally by Gaussian, so
the standard mode reproduces Gaussian's output to numerical precision.
