"""
Parse Gaussian frequency log files to extract thermochemistry data.
"""

import re
from dataclasses import dataclass, field


@dataclass
class GaussianData:
    """All data extracted from a Gaussian frequency log file."""
    # Basic info
    temperature: float = 0.0          # Kelvin
    pressure: float = 0.0             # Atm
    molecular_mass: float = 0.0       # amu

    # Atoms
    atom_numbers: list = field(default_factory=list)     # atomic numbers
    atom_masses: list = field(default_factory=list)      # amu
    coordinates: list = field(default_factory=list)      # [[x,y,z], ...] in Angstroms

    # Frequencies
    frequencies: list = field(default_factory=list)      # all frequencies in cm^-1
    n_imaginary: int = 0

    # Rotational
    rot_temperatures: list = field(default_factory=list)  # Kelvin [Θ1, Θ2, Θ3]
    rot_constants_ghz: list = field(default_factory=list) # GHz
    rot_symmetry_number: int = 1

    # SCF energy
    scf_energy: float = 0.0           # Hartree

    # Gaussian's printed thermochemistry (for validation)
    zpe_joules_per_mol: float = 0.0
    zpe_kcal_per_mol: float = 0.0
    zpe_correction: float = 0.0       # Hartree
    thermal_correction_energy: float = 0.0     # Hartree
    thermal_correction_enthalpy: float = 0.0   # Hartree
    thermal_correction_gibbs: float = 0.0      # Hartree

    # Sum lines
    sum_elec_zpe: float = 0.0
    sum_elec_thermal_e: float = 0.0
    sum_elec_thermal_h: float = 0.0
    sum_elec_thermal_g: float = 0.0

    # E, Cv, S breakdown
    e_total: float = 0.0       # kcal/mol
    e_electronic: float = 0.0
    e_translational: float = 0.0
    e_rotational: float = 0.0
    e_vibrational: float = 0.0

    cv_total: float = 0.0      # cal/mol-K
    cv_electronic: float = 0.0
    cv_translational: float = 0.0
    cv_rotational: float = 0.0
    cv_vibrational: float = 0.0

    s_total: float = 0.0       # cal/mol-K
    s_electronic: float = 0.0
    s_translational: float = 0.0
    s_rotational: float = 0.0
    s_vibrational: float = 0.0


def parse_log(filepath: str) -> GaussianData:
    """Parse a Gaussian frequency log file."""
    with open(filepath, 'r') as f:
        text = f.read()

    data = GaussianData()

    # Temperature and Pressure
    m = re.search(r'Temperature\s+([0-9.]+)\s+Kelvin\.\s+Pressure\s+([0-9.]+)\s+Atm\.', text)
    if m:
        data.temperature = float(m.group(1))
        data.pressure = float(m.group(2))

    # Molecular mass
    m = re.search(r'Molecular mass:\s+([0-9.]+)\s+amu\.', text)
    if m:
        data.molecular_mass = float(m.group(1))

    # Atomic numbers and masses
    for m in re.finditer(r'Atom\s+\d+\s+has atomic number\s+(\d+)\s+and mass\s+([0-9.]+)', text):
        data.atom_numbers.append(int(m.group(1)))
        data.atom_masses.append(float(m.group(2)))

    # Last Input orientation or Standard orientation block
    orientation_pattern = re.compile(
        r'(?:Input orientation|Standard orientation):\s*\n'
        r'\s*-+\s*\n'
        r'\s*Center\s+Atomic\s+Atomic\s+Coordinates \(Angstroms\)\s*\n'
        r'\s*Number\s+Number\s+Type\s+X\s+Y\s+Z\s*\n'
        r'\s*-+\s*\n'
        r'(.*?)\n\s*-+',
        re.DOTALL
    )
    matches = list(orientation_pattern.finditer(text))
    if matches:
        block = matches[-1].group(1)
        data.coordinates = []
        for line in block.strip().split('\n'):
            parts = line.split()
            if len(parts) >= 6:
                data.coordinates.append([float(parts[3]), float(parts[4]), float(parts[5])])

    # Frequencies
    for m in re.finditer(r'Frequencies --\s+(.*)', text):
        freqs = m.group(1).split()
        for f in freqs:
            data.frequencies.append(float(f))

    # Count imaginary
    data.n_imaginary = sum(1 for f in data.frequencies if f < 0)

    # Rotational temperatures
    m = re.search(r'Rotational temperatures \(Kelvin\)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)', text)
    if m:
        data.rot_temperatures = [float(m.group(1)), float(m.group(2)), float(m.group(3))]

    # Rotational constants
    m = re.search(r'Rotational constants \(GHZ\):\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)', text)
    if m:
        data.rot_constants_ghz = [float(m.group(1)), float(m.group(2)), float(m.group(3))]

    # Rotational symmetry number
    m = re.search(r'Rotational symmetry number\s+(\d+)\.', text)
    if m:
        data.rot_symmetry_number = int(m.group(1))

    # SCF energy (last occurrence)
    scf_matches = list(re.finditer(r'SCF Done:\s+E\([A-Za-z0-9]+\)\s+=\s+(-?[0-9.]+)', text))
    if scf_matches:
        data.scf_energy = float(scf_matches[-1].group(1))

    # ZPE
    m = re.search(r'Zero-point vibrational energy\s+([0-9.]+)\s+\(Joules/Mol\)', text)
    if m:
        data.zpe_joules_per_mol = float(m.group(1))
    m = re.search(r'Zero-point vibrational energy.*?\n\s+([0-9.]+)\s+\(Kcal/Mol\)', text)
    if m:
        data.zpe_kcal_per_mol = float(m.group(1))

    # Thermal corrections
    m = re.search(r'Zero-point correction=\s+([0-9.]+)', text)
    if m:
        data.zpe_correction = float(m.group(1))
    m = re.search(r'Thermal correction to Energy=\s+([0-9.]+)', text)
    if m:
        data.thermal_correction_energy = float(m.group(1))
    m = re.search(r'Thermal correction to Enthalpy=\s+([0-9.]+)', text)
    if m:
        data.thermal_correction_enthalpy = float(m.group(1))
    m = re.search(r'Thermal correction to Gibbs Free Energy=\s+([0-9.]+)', text)
    if m:
        data.thermal_correction_gibbs = float(m.group(1))

    # Sum lines
    m = re.search(r'Sum of electronic and zero-point Energies=\s+(-?[0-9.]+)', text)
    if m:
        data.sum_elec_zpe = float(m.group(1))
    m = re.search(r'Sum of electronic and thermal Energies=\s+(-?[0-9.]+)', text)
    if m:
        data.sum_elec_thermal_e = float(m.group(1))
    m = re.search(r'Sum of electronic and thermal Enthalpies=\s+(-?[0-9.]+)', text)
    if m:
        data.sum_elec_thermal_h = float(m.group(1))
    m = re.search(r'Sum of electronic and thermal Free Energies=\s+(-?[0-9.]+)', text)
    if m:
        data.sum_elec_thermal_g = float(m.group(1))

    # E, Cv, S breakdown table
    thermo_table_pattern = re.compile(
        r'E \(Thermal\)\s+CV\s+S\s*\n'
        r'\s*KCal/Mol\s+Cal/Mol-Kelvin\s+Cal/Mol-Kelvin\s*\n'
        r'(.*?)(?:\n\s*\n|\n\s*Q\s)',
        re.DOTALL
    )
    m = thermo_table_pattern.search(text)
    if m:
        block = m.group(1)
        for line in block.strip().split('\n'):
            line = line.strip()
            if not line:
                continue
            # Match lines like: "Total   488.699   150.289   250.996"
            # or "Vibration  1   0.593   1.986   7.296"
            parts = line.split()
            if len(parts) < 4:
                continue
            label = parts[0]
            # Get the last 3 numbers
            try:
                e_val = float(parts[-3])
                cv_val = float(parts[-2])
                s_val = float(parts[-1])
            except (ValueError, IndexError):
                continue

            if label == 'Total':
                data.e_total = e_val
                data.cv_total = cv_val
                data.s_total = s_val
            elif label == 'Electronic':
                data.e_electronic = e_val
                data.cv_electronic = cv_val
                data.s_electronic = s_val
            elif label == 'Translational':
                data.e_translational = e_val
                data.cv_translational = cv_val
                data.s_translational = s_val
            elif label == 'Rotational':
                data.e_rotational = e_val
                data.cv_rotational = cv_val
                data.s_rotational = s_val
            elif label == 'Vibrational':
                data.e_vibrational = e_val
                data.cv_vibrational = cv_val
                data.s_vibrational = s_val

    return data
