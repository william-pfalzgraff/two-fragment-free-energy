"""
IGRRHO thermochemistry calculations (Ochterski / McQuarrie formulas).

All functions return values in the units noted. Final assembly converts
to Hartree for comparison with Gaussian output.
"""

import math
import numpy as np

# Physical constants (CODATA values used by Gaussian)
R_JmolK = 8.314472          # J/(mol·K)
R_calmolK = 1.987216        # cal/(mol·K)  (thermochemical calorie)
kB = 1.3806503e-23          # J/K
h = 6.62606896e-34          # J·s
NA = 6.02214199e23          # 1/mol
c = 2.99792458e10           # cm/s
amu_to_kg = 1.66053878e-27  # kg/amu
hbar = 1.054571628e-34      # J·s (h / 2π)
angstrom_to_m = 1e-10       # m/Å
hartree_to_joule = 4.35974394e-18  # J/Hartree
atm_to_Pa = 101325.0        # Pa/atm
cal_to_J = 4.184            # J/cal
kcal_to_hartree = 1.0 / 627.5095   # Hartree/kcal


# ---------------------------------------------------------------------------
# Translational contributions (ideal gas, 3D)
# ---------------------------------------------------------------------------

def e_trans(T):
    """Translational thermal energy in kcal/mol."""
    return 1.5 * R_JmolK * T / (cal_to_J * 1000.0)


def cv_trans():
    """Translational Cv in cal/(mol·K)."""
    return 1.5 * R_calmolK


def s_trans(mass_amu, T, P_atm):
    """Sackur-Tetrode translational entropy in cal/(mol·K).

    Parameters
    ----------
    mass_amu : float
        Molecular mass in amu.
    T : float
        Temperature in Kelvin.
    P_atm : float
        Pressure in atm.
    """
    m = mass_amu * amu_to_kg  # kg
    P = P_atm * atm_to_Pa     # Pa

    # q_trans/V = (2*pi*m*kB*T/h^2)^(3/2)
    # q_trans = (q_trans/V) * (kB*T/P)
    # S = R * [ln(q_trans) + 5/2]
    lambda_thermal_sq = (2.0 * math.pi * m * kB * T) / (h * h)
    q_trans = lambda_thermal_sq ** 1.5 * (kB * T / P)
    return R_calmolK * (math.log(q_trans) + 2.5)


# ---------------------------------------------------------------------------
# Rotational contributions (rigid rotor, nonlinear molecule)
# ---------------------------------------------------------------------------

def e_rot(T):
    """Rotational thermal energy in kcal/mol (nonlinear)."""
    return 1.5 * R_JmolK * T / (cal_to_J * 1000.0)


def cv_rot():
    """Rotational Cv in cal/(mol·K) (nonlinear)."""
    return 1.5 * R_calmolK


def s_rot(rot_temps, sigma, T):
    """Rotational entropy in cal/(mol·K) for a nonlinear molecule.

    Parameters
    ----------
    rot_temps : list of 3 floats
        Rotational temperatures [Θ_A, Θ_B, Θ_C] in Kelvin.
    sigma : int
        Rotational symmetry number.
    T : float
        Temperature in Kelvin.
    """
    theta_prod = rot_temps[0] * rot_temps[1] * rot_temps[2]
    q_rot = (math.sqrt(math.pi) / sigma) * (T ** 1.5) / math.sqrt(theta_prod)
    return R_calmolK * (math.log(q_rot) + 1.5)


# ---------------------------------------------------------------------------
# Vibrational contributions (harmonic oscillator)
# ---------------------------------------------------------------------------

def freq_to_vib_temp(freq_cm):
    """Convert frequency in cm^-1 to vibrational temperature Θ_v in K."""
    return h * c * freq_cm / kB


def zpe_mode(freq_cm):
    """Zero-point energy for one mode in Joules/mol."""
    return 0.5 * NA * h * c * freq_cm


def e_vib_mode(freq_cm, T):
    """Vibrational thermal energy for one mode in kcal/mol.

    Includes ZPE: E = R*Θ_v*(1/2 + 1/(exp(Θ_v/T)-1))
    """
    theta = freq_to_vib_temp(freq_cm)
    u = theta / T
    # E in J/mol = R_JmolK * theta * (0.5 + 1/(exp(u)-1))
    e_jmol = R_JmolK * theta * (0.5 + 1.0 / (math.exp(u) - 1.0))
    return e_jmol / (cal_to_J * 1000.0)  # kcal/mol


def cv_vib_mode(freq_cm, T):
    """Vibrational Cv for one mode in cal/(mol·K)."""
    theta = freq_to_vib_temp(freq_cm)
    u = theta / T
    eu = math.exp(u)
    return R_calmolK * u * u * eu / ((eu - 1.0) ** 2)


def s_vib_mode(freq_cm, T):
    """Vibrational entropy for one mode in cal/(mol·K)."""
    theta = freq_to_vib_temp(freq_cm)
    u = theta / T
    eu = math.exp(u)
    return R_calmolK * (u / (eu - 1.0) - math.log(1.0 - math.exp(-u)))


# ---------------------------------------------------------------------------
# Sum over vibrational modes
# ---------------------------------------------------------------------------

def e_vib_total(freqs, T):
    """Total vibrational thermal energy in kcal/mol (includes ZPE)."""
    return sum(e_vib_mode(f, T) for f in freqs)


def cv_vib_total(freqs, T):
    """Total vibrational Cv in cal/(mol·K)."""
    return sum(cv_vib_mode(f, T) for f in freqs)


def s_vib_total(freqs, T):
    """Total vibrational entropy in cal/(mol·K)."""
    return sum(s_vib_mode(f, T) for f in freqs)


def zpe_total(freqs):
    """Total ZPE in kcal/mol."""
    return sum(zpe_mode(f) for f in freqs) / (cal_to_J * 1000.0)


def zpe_total_hartree(freqs):
    """Total ZPE in Hartree."""
    return zpe_total(freqs) * kcal_to_hartree


# ---------------------------------------------------------------------------
# Electronic contributions
# ---------------------------------------------------------------------------

def e_elec():
    """Electronic thermal energy in kcal/mol (singlet ground state)."""
    return 0.0


def cv_elec():
    """Electronic Cv in cal/(mol·K)."""
    return 0.0


def s_elec(multiplicity=1):
    """Electronic entropy in cal/(mol·K)."""
    return R_calmolK * math.log(multiplicity)


# ---------------------------------------------------------------------------
# Full thermochemistry assembly
# ---------------------------------------------------------------------------

def compute_thermochem(freqs, mass_amu, rot_temps, sigma, T, P_atm,
                       multiplicity=1):
    """Compute standard IGRRHO thermochemistry matching Gaussian output.

    Parameters
    ----------
    freqs : list of float
        Real vibrational frequencies in cm^-1 (imaginary already excluded).
    mass_amu : float
        Total molecular mass in amu.
    rot_temps : list of 3 floats
        Rotational temperatures in Kelvin.
    sigma : int
        Rotational symmetry number.
    T : float
        Temperature in Kelvin.
    P_atm : float
        Pressure in atm.

    Returns
    -------
    dict with all thermochemical quantities.
    """
    # Individual components
    E_t = e_trans(T)
    E_r = e_rot(T)
    E_v = e_vib_total(freqs, T)
    E_e = e_elec()

    Cv_t = cv_trans()
    Cv_r = cv_rot()
    Cv_v = cv_vib_total(freqs, T)
    Cv_e = cv_elec()

    S_t = s_trans(mass_amu, T, P_atm)
    S_r = s_rot(rot_temps, sigma, T)
    S_v = s_vib_total(freqs, T)
    S_e = s_elec(multiplicity)

    E_total = E_t + E_r + E_v + E_e
    Cv_total = Cv_t + Cv_r + Cv_v + Cv_e
    S_total = S_t + S_r + S_v + S_e

    # ZPE
    zpe_kcal = zpe_total(freqs)
    zpe_hartree = zpe_kcal * kcal_to_hartree

    # Thermal corrections in Hartree
    E_thermal_hartree = E_total * kcal_to_hartree
    kT_hartree = R_JmolK * T / (NA * hartree_to_joule)
    H_correction = E_thermal_hartree + kT_hartree
    # S is in cal/(mol·K), so T*S is in cal/mol; divide by 1000 for kcal/mol
    TS_hartree = S_total * T / 1000.0 * kcal_to_hartree
    G_correction = H_correction - TS_hartree

    return {
        # Components in kcal/mol, cal/mol-K
        'E_trans': E_t, 'E_rot': E_r, 'E_vib': E_v, 'E_elec': E_e,
        'E_total': E_total,
        'Cv_trans': Cv_t, 'Cv_rot': Cv_r, 'Cv_vib': Cv_v, 'Cv_elec': Cv_e,
        'Cv_total': Cv_total,
        'S_trans': S_t, 'S_rot': S_r, 'S_vib': S_v, 'S_elec': S_e,
        'S_total': S_total,
        # ZPE
        'ZPE_kcal': zpe_kcal,
        'ZPE_hartree': zpe_hartree,
        # Corrections in Hartree
        'thermal_correction_energy': E_thermal_hartree,
        'thermal_correction_enthalpy': H_correction,
        'thermal_correction_gibbs': G_correction,
        'kT_hartree': kT_hartree,
    }


def compute_two_fragment_thermochem(
    freqs_all_real, mass_A, mass_B,
    rot_temps_A, rot_temps_B,
    sigma_A, sigma_B,
    T, P_atm, n_remove=6, multiplicity=1
):
    """Compute two-fragment corrected thermochemistry.

    Parameters
    ----------
    freqs_all_real : list of float
        All real vibrational frequencies in cm^-1 (sorted ascending).
    mass_A, mass_B : float
        Fragment masses in amu.
    rot_temps_A, rot_temps_B : list of 3 floats
        Rotational temperatures for each fragment in Kelvin.
    sigma_A, sigma_B : int
        Rotational symmetry numbers for each fragment.
    T : float
        Temperature in Kelvin.
    P_atm : float
        Pressure in atm.
    n_remove : int
        Number of lowest real frequencies to remove (default 6).

    Returns
    -------
    dict with all thermochemical quantities.
    """
    # Sort and split frequencies
    sorted_freqs = sorted(freqs_all_real)
    removed_freqs = sorted_freqs[:n_remove]
    kept_freqs = sorted_freqs[n_remove:]

    # Translation: two independent fragments
    E_t = 2 * e_trans(T)
    Cv_t = 2 * cv_trans()
    S_t = s_trans(mass_A, T, P_atm) + s_trans(mass_B, T, P_atm)

    # Rotation: two independent fragments
    E_r = 2 * e_rot(T)
    Cv_r = 2 * cv_rot()
    S_r = s_rot(rot_temps_A, sigma_A, T) + s_rot(rot_temps_B, sigma_B, T)

    # Vibration: only kept modes
    E_v = e_vib_total(kept_freqs, T)
    Cv_v = cv_vib_total(kept_freqs, T)
    S_v = s_vib_total(kept_freqs, T)

    # Electronic
    E_e = e_elec()
    Cv_e = cv_elec()
    S_e = s_elec(multiplicity)

    E_total = E_t + E_r + E_v + E_e
    Cv_total = Cv_t + Cv_r + Cv_v + Cv_e
    S_total = S_t + S_r + S_v + S_e

    # ZPE (only kept modes)
    zpe_kcal = zpe_total(kept_freqs)
    zpe_hartree = zpe_kcal * kcal_to_hartree

    # Also compute ZPE of removed modes for reporting
    zpe_removed_kcal = zpe_total(removed_freqs)

    # Thermal corrections in Hartree
    E_thermal_hartree = E_total * kcal_to_hartree
    kT_hartree = R_JmolK * T / (NA * hartree_to_joule)
    H_correction = E_thermal_hartree + kT_hartree
    # S is in cal/(mol·K), so T*S is in cal/mol; divide by 1000 for kcal/mol
    TS_hartree = S_total * T / 1000.0 * kcal_to_hartree
    G_correction = H_correction - TS_hartree

    return {
        # Components
        'E_trans': E_t, 'E_rot': E_r, 'E_vib': E_v, 'E_elec': E_e,
        'E_total': E_total,
        'Cv_trans': Cv_t, 'Cv_rot': Cv_r, 'Cv_vib': Cv_v, 'Cv_elec': Cv_e,
        'Cv_total': Cv_total,
        'S_trans': S_t, 'S_rot': S_r, 'S_vib': S_v, 'S_elec': S_e,
        'S_total': S_total,
        # ZPE
        'ZPE_kcal': zpe_kcal,
        'ZPE_hartree': zpe_hartree,
        'ZPE_removed_kcal': zpe_removed_kcal,
        # Corrections in Hartree
        'thermal_correction_energy': E_thermal_hartree,
        'thermal_correction_enthalpy': H_correction,
        'thermal_correction_gibbs': G_correction,
        'kT_hartree': kT_hartree,
        # Frequency info
        'removed_freqs': removed_freqs,
        'kept_freqs': kept_freqs,
    }


# ---------------------------------------------------------------------------
# Fragment moment of inertia / rotational temperatures
# ---------------------------------------------------------------------------

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
