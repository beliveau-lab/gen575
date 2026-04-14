"""Isotope distribution calculation library.

Implements the Kubinyi algorithm (Analytica Chimica Acta, 247, 1991, 107-119)
for predicting isotope distributions from elemental compositions. Natural
isotope abundances are from D.E. Matthews as implemented in IDCalc.
"""

import re

import numpy as np


# ---------------------------------------------------------------------------
# Monoisotopic masses (Da)
# ---------------------------------------------------------------------------
MONO_MASSES = {
    "C": 12.0,
    "H": 1.007825,
    "N": 14.003074,
    "O": 15.9949146,
    "S": 31.9720718,
    "P": 30.9737634,
}

# ---------------------------------------------------------------------------
# Natural isotope abundances (relative to most abundant = 100)
# Values from D.E. Matthews, as used in IDCalc
# ---------------------------------------------------------------------------
ISOTOPE_ABUNDANCES = {
    "C": [100.0, 1.0958793],
    "H": [100.0, 0.0142],
    "O": [100.0, 0.03799194, 0.20499609],
    "N": [100.0, 0.368351851],
    "S": [100.0, 0.789308, 4.430646, 0.0, 0.021048],
    "P": [100.0],
}

# ---------------------------------------------------------------------------
# Amino acid residue compositions: (C, H, N, O, S)
# These are residue masses (loss of water already accounted for in peptide
# assembly; water is added back for termini in sequence_to_composition).
# ---------------------------------------------------------------------------
AMINO_ACIDS = {
    "A": (3, 5, 1, 1, 0),
    "R": (6, 12, 4, 1, 0),
    "N": (4, 6, 2, 2, 0),
    "D": (4, 5, 1, 3, 0),
    "C": (3, 5, 1, 1, 1),
    "Q": (5, 8, 2, 2, 0),
    "E": (5, 7, 1, 3, 0),
    "G": (2, 3, 1, 1, 0),
    "H": (6, 7, 3, 1, 0),
    "I": (6, 11, 1, 1, 0),
    "L": (6, 11, 1, 1, 0),
    "K": (6, 12, 2, 1, 0),
    "M": (5, 9, 1, 1, 1),
    "F": (9, 9, 1, 1, 0),
    "P": (5, 7, 1, 1, 0),
    "S": (3, 5, 1, 2, 0),
    "T": (4, 7, 1, 2, 0),
    "W": (11, 10, 2, 1, 0),
    "Y": (9, 9, 1, 2, 0),
    "V": (5, 9, 1, 1, 0),
}

PROTON_MASS = 1.00727647  # Da


# ---------------------------------------------------------------------------
# Formula and sequence parsing
# ---------------------------------------------------------------------------

def parse_formula(formula: str) -> dict:
    """Parse an elemental formula string into a composition dictionary.

    Args:
        formula: Elemental formula such as "C60H86N16O15".

    Returns:
        Dictionary mapping element symbols to atom counts.

    Raises:
        ValueError: If the formula contains unsupported elements.
    """
    tokens = re.findall(r"([A-Z][a-z]?)(\d*)", formula)
    composition = {}
    for element, count_str in tokens:
        if not element:
            continue
        count = int(count_str) if count_str else 1
        if element not in MONO_MASSES:
            raise ValueError(
                f"Unsupported element '{element}'. "
                f"Supported: {', '.join(sorted(MONO_MASSES))}"
            )
        composition[element] = composition.get(element, 0) + count
    return composition


def sequence_to_composition(sequence: str) -> dict:
    """Convert a peptide/protein sequence to an elemental composition.

    Sums amino acid residue compositions and adds water (H2O) for the
    N- and C-termini.

    Args:
        sequence: Amino acid sequence in single-letter code (uppercase).

    Returns:
        Dictionary mapping element symbols to atom counts.

    Raises:
        ValueError: If the sequence contains invalid amino acid codes.
    """
    elements = ["C", "H", "N", "O", "S"]
    composition = {e: 0 for e in elements}

    for aa in sequence.upper():
        if aa not in AMINO_ACIDS:
            raise ValueError(f"Invalid amino acid '{aa}'.")
        residue = AMINO_ACIDS[aa]
        for i, element in enumerate(elements):
            composition[element] += residue[i]

    # Add water for termini: +2H, +1O
    composition["H"] += 2
    composition["O"] += 1

    # Remove elements with zero count
    return {k: v for k, v in composition.items() if v > 0}


def parse_input(text: str) -> dict:
    """Auto-detect input type and return elemental composition.

    If the string contains digits, it is treated as an elemental formula.
    Otherwise it is treated as a peptide/protein sequence.

    Args:
        text: Formula or amino acid sequence string.

    Returns:
        Dictionary mapping element symbols to atom counts.
    """
    text = text.strip()
    if any(ch.isdigit() for ch in text):
        return parse_formula(text)
    return sequence_to_composition(text)


# ---------------------------------------------------------------------------
# Kubinyi isotope distribution algorithm
# ---------------------------------------------------------------------------

def _convolve_element(pattern: np.ndarray, abundances: np.ndarray,
                      count: int, precision: float) -> np.ndarray:
    """Convolve an element's isotope pattern using binary expansion.

    Decomposes *count* into powers of two so that a molecule with hundreds
    of atoms of a given element requires only ~log2(count) convolutions
    instead of *count* convolutions.

    Args:
        pattern: Current cumulative isotope pattern.
        abundances: Isotope abundance array for this element.
        count: Number of atoms of this element.
        precision: Minimum relative abundance to retain.

    Returns:
        Updated isotope pattern after convolving *count* atoms.
    """
    if count == 0 or len(abundances) <= 1:
        return pattern

    # Build the single-atom pattern (normalized)
    single = abundances / abundances.max()

    power_pattern = single.copy()
    remaining = count

    while remaining > 0:
        if remaining & 1:
            pattern = np.convolve(pattern, power_pattern)
            max_val = pattern.max()
            if max_val > 0:
                pattern /= max_val

            # Trim trailing values below precision
            above = np.where(pattern > precision)[0]
            if len(above) > 0:
                pattern = pattern[above[0]:above[-1] + 1]

        remaining >>= 1
        if remaining > 0:
            power_pattern = np.convolve(power_pattern, power_pattern)
            max_val = power_pattern.max()
            if max_val > 0:
                power_pattern /= max_val
            above = np.where(power_pattern > precision)[0]
            if len(above) > 0:
                power_pattern = power_pattern[above[0]:above[-1] + 1]

    return pattern


def calculate_isotope_distribution(composition: dict,
                                   precision: float = 1e-7) -> np.ndarray:
    """Calculate the isotope distribution for a given elemental composition.

    Implements the Kubinyi algorithm with binary expansion for performance.

    Args:
        composition: Element symbol to atom count mapping.
        precision: Minimum relative abundance to retain (default 1e-7).

    Returns:
        1-D array of relative abundances normalized so the maximum is 1.0.
        Index 0 corresponds to the monoisotopic peak, index i corresponds
        to monoisotopic mass + i Da.
    """
    pattern = np.array([1.0])

    for element, count in composition.items():
        if element not in ISOTOPE_ABUNDANCES:
            continue
        abundances = np.array(ISOTOPE_ABUNDANCES[element], dtype=np.float64)
        pattern = _convolve_element(pattern, abundances, count, precision)

    # Final normalization
    max_val = pattern.max()
    if max_val > 0:
        pattern /= max_val

    return pattern


# ---------------------------------------------------------------------------
# Mass calculations
# ---------------------------------------------------------------------------

def monoisotopic_mass(composition: dict) -> float:
    """Calculate the neutral monoisotopic mass of a composition.

    Args:
        composition: Element symbol to atom count mapping.

    Returns:
        Monoisotopic mass in Daltons.
    """
    mass = 0.0
    for element, count in composition.items():
        if element in MONO_MASSES:
            mass += count * MONO_MASSES[element]
    return mass


def average_mass(composition: dict, distribution: np.ndarray) -> float:
    """Calculate the abundance-weighted average mass.

    Args:
        composition: Element symbol to atom count mapping.
        distribution: Isotope distribution array from
            calculate_isotope_distribution.

    Returns:
        Average mass in Daltons.
    """
    mono = monoisotopic_mass(composition)
    indices = np.arange(len(distribution), dtype=np.float64)
    total = distribution.sum()
    if total == 0:
        return mono
    return np.sum(distribution * (mono + indices)) / total


def mz_values(mono_mass: float, num_peaks: int, charge: int) -> np.ndarray:
    """Calculate m/z values for each isotopic peak at a given charge state.

    Args:
        mono_mass: Neutral monoisotopic mass in Daltons.
        num_peaks: Number of isotopic peaks.
        charge: Charge state (must be >= 1).

    Returns:
        Array of m/z values for each isotopic peak.
    """
    indices = np.arange(num_peaks, dtype=np.float64)
    return (mono_mass + indices + charge * PROTON_MASS) / charge


# ---------------------------------------------------------------------------
# Instrument-dependent resolving power
# ---------------------------------------------------------------------------

# Instrument types and their resolution scaling behavior:
#   TOF:        constant resolving power across m/z
#   Quad/Trap:  constant peak width (FWHM), so R scales linearly with m/z
#   Orbitrap:   R proportional to 1/sqrt(m/z)
#   FT-ICR:     R proportional to 1/m/z

INSTRUMENT_TYPES = {
    "TOF": "TOF",
    "Quad/Ion Trap": "quad",
    "Orbitrap": "orbitrap",
    "FT-ICR": "fticr",
}


def effective_resolving_power(mz: float, instrument: str,
                              rp_ref: float, mz_ref: float = 200.0) -> float:
    """Calculate the effective resolving power at a given m/z.

    Args:
        mz: The m/z value at which to calculate resolving power.
        instrument: Instrument type key from INSTRUMENT_TYPES values
            ("TOF", "quad", "orbitrap", "fticr").
        rp_ref: Resolving power specified at the reference m/z.
        mz_ref: Reference m/z at which rp_ref is defined (default 200).

    Returns:
        Effective resolving power at the given m/z.
    """
    if instrument == "TOF":
        return rp_ref
    elif instrument == "quad":
        # Constant FWHM: FWHM = mz_ref / rp_ref; R(mz) = mz / FWHM
        return rp_ref * mz / mz_ref
    elif instrument == "orbitrap":
        # R proportional to 1/sqrt(m/z)
        return rp_ref * np.sqrt(mz_ref / mz)
    elif instrument == "fticr":
        # R proportional to 1/m/z
        return rp_ref * mz_ref / mz
    else:
        return rp_ref


# ---------------------------------------------------------------------------
# Profile generation
# ---------------------------------------------------------------------------

def generate_profile(mz_centers: np.ndarray, abundances: np.ndarray,
                     resolving_power: float,
                     num_points: int = 5000,
                     instrument: str = "TOF",
                     mz_ref: float = 200.0) -> tuple:
    """Generate a continuous peak profile from discrete isotopic peaks.

    Each peak is modeled as a Gaussian with FWHM determined by the
    instrument-dependent resolving power at that peak's m/z.

    Args:
        mz_centers: m/z values of the isotopic peaks.
        abundances: Relative abundances of the peaks.
        resolving_power: Resolving power at the reference m/z.
        num_points: Number of points in the output profile.
        instrument: Instrument type ("TOF", "quad", "orbitrap", "fticr").
        mz_ref: Reference m/z at which resolving_power is specified.

    Returns:
        Tuple of (x_array, y_array) for plotting.
    """
    # FWHM to sigma conversion: sigma = FWHM / (2 * sqrt(2 * ln(2)))
    fwhm_to_sigma = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))

    # Determine x range with padding based on the widest peak
    # (lowest effective R gives the widest peak)
    rp_values = np.array([
        effective_resolving_power(mz, instrument, resolving_power, mz_ref)
        for mz in mz_centers
    ])
    fwhm_values = mz_centers / rp_values
    max_fwhm = fwhm_values.max()
    padding = 3.0 * max_fwhm
    x = np.linspace(mz_centers.min() - padding,
                     mz_centers.max() + padding,
                     num_points)

    y = np.zeros_like(x)
    for mz, abundance, rp in zip(mz_centers, abundances, rp_values):
        fwhm = mz / rp
        sigma = fwhm * fwhm_to_sigma
        y += abundance * np.exp(-0.5 * ((x - mz) / sigma) ** 2)

    # Normalize to max = 100
    max_y = y.max()
    if max_y > 0:
        y = y / max_y * 100.0

    return x, y


# ---------------------------------------------------------------------------
# Convenience: format composition as formula string
# ---------------------------------------------------------------------------

def composition_to_formula(composition: dict) -> str:
    """Format a composition dictionary as a formula string.

    Args:
        composition: Element symbol to atom count mapping.

    Returns:
        Formula string like "C60H86N16O15".
    """
    # Standard order: C, H, N, O, S, P, then alphabetical
    order = ["C", "H", "N", "O", "S", "P"]
    parts = []
    for element in order:
        if element in composition and composition[element] > 0:
            count = composition[element]
            parts.append(f"{element}{count}" if count > 1 else element)
    for element in sorted(composition):
        if element not in order and composition[element] > 0:
            count = composition[element]
            parts.append(f"{element}{count}" if count > 1 else element)
    return "".join(parts)


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # Water: H2O
    comp = parse_formula("H2O")
    mono = monoisotopic_mass(comp)
    print(f"Water monoisotopic mass: {mono:.4f} Da (expected ~18.0106)")

    # C60H86N16O15 from lecture slide
    comp = parse_formula("C60H86N16O15")
    mono = monoisotopic_mass(comp)
    dist = calculate_isotope_distribution(comp)
    avg = average_mass(comp, dist)
    print(f"C60H86N16O15 monoisotopic mass: {mono:.4f} Da")
    print(f"C60H86N16O15 average mass: {avg:.4f} Da")
    print(f"  Distribution ({len(dist)} peaks): {dist[:6].round(4)}")

    # MYFAVRITEPEPTIDE
    comp = sequence_to_composition("MYFAVRITEPEPTIDE")
    formula = composition_to_formula(comp)
    mono = monoisotopic_mass(comp)
    dist = calculate_isotope_distribution(comp)
    avg = average_mass(comp, dist)
    print(f"\nMYFAVRITEPEPTIDE: {formula}")
    print(f"  Monoisotopic mass: {mono:.4f} Da")
    print(f"  Average mass: {avg:.4f} Da")

    # Myoglobin formula
    comp = parse_formula("C774H1224N210O222S5")
    mono = monoisotopic_mass(comp)
    dist = calculate_isotope_distribution(comp)
    avg = average_mass(comp, dist)
    print(f"\nMyoglobin (C774H1224N210O222S5):")
    print(f"  Monoisotopic mass: {mono:.4f} Da")
    print(f"  Average mass: {avg:.4f} Da")
    print(f"  Distribution has {len(dist)} peaks")
    print(f"  Most abundant peak: M+{np.argmax(dist)}")
