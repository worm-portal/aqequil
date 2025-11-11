"""
Shared utility functions for parsing chemical formulas and oxidation states.

This module consolidates formula parsing logic used across multiple modules
to avoid code duplication and ensure consistency.
"""

import re
import pandas as pd


# Common regex patterns for formula parsing
# Updated to handle elements of any length (e.g., "Edta", "Na", "Mg")
# [A-Z][a-z]* matches: capital letter followed by zero or more lowercase letters
FORMULA_OX_PATTERN = r'([0-9]*\.?[0-9]+)?([A-Z][a-z]*)([\+\-][0-9]*)?'
FORMULA_OX_WITH_PSEUDOELEM_PATTERN = r'([0-9]*\.?[0-9]+)?([A-Z][a-z]*j?[ivxlcdm]*[pnz]?)([\+\-][0-9]*)?'


def normalize_charge_string(charge_str):
    """
    Normalize charge string format.

    Converts '+' to '+1' and '-' to '-1' for consistency.

    Parameters
    ----------
    charge_str : str or None
        Charge string (e.g., '+', '-', '+2', '-3')

    Returns
    -------
    str or None
        Normalized charge string, or None if input was None

    Examples
    --------
    >>> normalize_charge_string('+')
    '+1'
    >>> normalize_charge_string('-')
    '-1'
    >>> normalize_charge_string('+3')
    '+3'
    >>> normalize_charge_string(None)
    None
    """
    if charge_str is None:
        return None
    if charge_str == '+':
        return '+1'
    elif charge_str == '-':
        return '-1'
    else:
        return charge_str


def parse_formula_ox_part(part, include_pseudoelements=False):
    """
    Parse a single part of a formula_ox string.

    Parameters
    ----------
    part : str
        Single component from formula_ox (e.g., "2Fe+3", "S-2", "Fejiiip+3")
    include_pseudoelements : bool
        Whether to recognize pseudoelement names (e.g., "Fejiiip")

    Returns
    -------
    tuple or None
        (coefficient, element, charge_str) or None if parsing fails

    Examples
    --------
    >>> parse_formula_ox_part("2Fe+3")
    (2.0, 'Fe', '+3')
    >>> parse_formula_ox_part("S-2")
    (1.0, 'S', '-2')
    >>> parse_formula_ox_part("Fejiiip+3", include_pseudoelements=True)
    (1.0, 'Fejiiip', '+3')
    """
    if include_pseudoelements:
        pattern = FORMULA_OX_WITH_PSEUDOELEM_PATTERN
    else:
        pattern = FORMULA_OX_PATTERN

    match = re.match(pattern, part)
    if not match:
        return None

    coeff_str, element, charge_str = match.groups()
    coeff = float(coeff_str) if coeff_str else 1.0
    charge_str = normalize_charge_string(charge_str)

    return (coeff, element, charge_str)


def parse_formula_ox_string(formula_ox_str, include_pseudoelements=False):
    """
    Parse a complete formula_ox string into components.

    Parameters
    ----------
    formula_ox_str : str
        Formula with oxidation states (e.g., "Fe+3 2S-2", "Co+2 2Co+3 4S-2")
    include_pseudoelements : bool
        Whether to recognize pseudoelement names

    Returns
    -------
    dict
        Dictionary mapping "Element+charge" to coefficient
        Returns empty dict if parsing fails

    Examples
    --------
    >>> parse_formula_ox_string("Fe+3 2S-2")
    {'Fe+3': 1.0, 'S-2': 2.0}
    >>> parse_formula_ox_string("Co+2 2Co+3")
    {'Co+2': 1.0, 'Co+3': 2.0}
    """
    if not formula_ox_str or pd.isna(formula_ox_str) or formula_ox_str == 'nan' or formula_ox_str == '':
        return {}

    result = {}
    parts = str(formula_ox_str).strip().split()

    for part in parts:
        parsed = parse_formula_ox_part(part, include_pseudoelements)
        if parsed:
            coeff, element, charge_str = parsed
            elem_with_charge = element + (charge_str if charge_str else '')
            result[elem_with_charge] = coeff

    return result


def get_element_oxidation_states_from_formula_ox(formula_ox_str, element):
    """
    Extract oxidation state(s) of a specific element from formula_ox string.

    Parameters
    ----------
    formula_ox_str : str
        Formula with oxidation states (e.g., "Co+2 2Co+3 4S-2")
    element : str
        Element to extract (e.g., "Co")

    Returns
    -------
    list of tuples
        List of (oxidation_state, count) tuples
        e.g., for Co in "Co+2 2Co+3 4S-2": [(2, 1.0), (3, 2.0)]

    Examples
    --------
    >>> get_element_oxidation_states_from_formula_ox("Co+2 2Co+3 4S-2", "Co")
    [(2, 1.0), (3, 2.0)]
    >>> get_element_oxidation_states_from_formula_ox("Fe+3 2S-2", "Fe")
    [(3, 1.0)]
    """
    if not formula_ox_str or pd.isna(formula_ox_str) or formula_ox_str == 'nan' or formula_ox_str == '':
        return []

    result = []
    parts = str(formula_ox_str).strip().split()

    for part in parts:
        # Check if this part contains the element (but not as part of another element)
        if not re.search(f'{element}([^a-z]|$)', part):
            continue

        parsed = parse_formula_ox_part(part, include_pseudoelements=False)
        if not parsed:
            continue

        coeff, elem, charge_str = parsed
        if elem != element:
            continue

        # Parse charge to get oxidation state
        if charge_str:
            ox_state = int(charge_str)
        else:
            ox_state = 0

        result.append((ox_state, coeff))

    return result


def calculate_average_oxidation_state(formula_ox_str, element):
    """
    Calculate average oxidation state of an element in a formula_ox string.

    Parameters
    ----------
    formula_ox_str : str
        Formula with oxidation states
    element : str
        Element to calculate average for

    Returns
    -------
    float
        Average oxidation state

    Examples
    --------
    >>> calculate_average_oxidation_state("Co+2 2Co+3", "Co")
    2.666...
    >>> calculate_average_oxidation_state("Fe+3", "Fe")
    3.0
    """
    ox_states = get_element_oxidation_states_from_formula_ox(formula_ox_str, element)
    if not ox_states:
        return 0.0

    total_charge = sum(ox * count for ox, count in ox_states)
    total_atoms = sum(count for ox, count in ox_states)

    return total_charge / total_atoms if total_atoms > 0 else 0.0
