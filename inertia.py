"""
Moment of inertia calculator for molecular fragments.

Computes principal moments of inertia and converts to rotational
temperatures for use in partition function calculations.
"""

import numpy as np

# Physical constants
kB = 1.3806503e-23          # J/K
hbar = 1.054571628e-34      # J·s (h / 2π)
amu_to_kg = 1.66053878e-27  # kg/amu
angstrom_to_m = 1e-10       # m/Å


def compute_fragment_rot_temps(coordinates, masses, atom_indices):
    """Compute rotational temperatures for a molecular fragment.

    Parameters
    ----------
    coordinates : list of [x, y, z]
        All atomic coordinates in Angstroms.
    masses : list of float
        All atomic masses in amu.
    atom_indices : list of int
        0-based indices of atoms belonging to this fragment.

    Returns
    -------
    rot_temps : list of 3 floats
        Rotational temperatures [Θ_A, Θ_B, Θ_C] in Kelvin, sorted descending.
    principal_moments : list of 3 floats
        Principal moments of inertia in amu·Å², sorted ascending.
    fragment_mass : float
        Total mass of the fragment in amu.
    """
    coords = np.array([coordinates[i] for i in atom_indices])  # Å
    m = np.array([masses[i] for i in atom_indices])             # amu

    fragment_mass = np.sum(m)

    # Center of mass
    com = np.sum(m[:, None] * coords, axis=0) / fragment_mass
    coords_com = coords - com  # shift to COM frame

    # Build inertia tensor (in amu·Å²)
    I = np.zeros((3, 3))
    for i in range(len(m)):
        x, y, z = coords_com[i]
        r2 = x*x + y*y + z*z
        I[0, 0] += m[i] * (r2 - x*x)
        I[1, 1] += m[i] * (r2 - y*y)
        I[2, 2] += m[i] * (r2 - z*z)
        I[0, 1] -= m[i] * x * y
        I[0, 2] -= m[i] * x * z
        I[1, 2] -= m[i] * y * z
    I[1, 0] = I[0, 1]
    I[2, 0] = I[0, 2]
    I[2, 1] = I[1, 2]

    # Diagonalize
    eigenvalues = np.linalg.eigvalsh(I)  # sorted ascending
    principal_moments = sorted(eigenvalues)  # amu·Å²

    # Convert to rotational temperatures: Θ_rot = ℏ² / (2·I·kB)
    rot_temps = []
    for I_val in principal_moments:
        I_si = I_val * amu_to_kg * (angstrom_to_m ** 2)  # kg·m²
        if I_si > 0:
            theta = hbar**2 / (2.0 * I_si * kB)
        else:
            theta = 0.0
        rot_temps.append(theta)

    # Sort descending (Gaussian convention: largest Θ first)
    rot_temps = sorted(rot_temps, reverse=True)

    return rot_temps, principal_moments, fragment_mass
