"""
Module for handling redox suppression in thermodynamic databases.

This module implements the logic to prevent equilibration between different
oxidation states of specified elements by creating "pseudoelements" that
trick EQ3/6 into treating different oxidation states as different elements.

For example, if Fe redox is suppressed:
- Fe+2 becomes the element "Fejiip"
- Fe+3 becomes the element "Fejiiip"
- Fe0 becomes the element "Fejz"

Pseudoelement naming convention:
- Element name (e.g., "Fe") is separated from charge by 'j'
- The charge magnitude is in lowercase roman numerals
- The letter at the end is either 'p' (positive), 'n' (negative), or 'z' (zero)
- Examples: Fe+3 -> Fejiiip, C-4 -> Cjivn, S0 -> Sjz
"""

import pandas as pd
import numpy as np
import pychnosz
import roman
import re
from .formula_parser import parse_formula_ox_part


def create_pseudoelement_name(element, oxidation_state):
    """
    Create a pseudoelement name from an element and its oxidation state.

    Parameters
    ----------
    element : str
        Element symbol (e.g., "Fe", "S", "C")
    oxidation_state : int or str
        Oxidation state as integer or string with sign (e.g., +3, "-2", 0)

    Returns
    -------
    str
        Pseudoelement name following the convention (e.g., "Fejiiip" for Fe+3)

    Examples
    --------
    >>> create_pseudoelement_name("Fe", 3)
    'Fejiiip'
    >>> create_pseudoelement_name("Fe", "+3")
    'Fejiiip'
    >>> create_pseudoelement_name("S", -2)
    'Sjiin'
    >>> create_pseudoelement_name("S", 0)
    'Sjz'
    """
    # Convert oxidation state to integer if it's a string
    if isinstance(oxidation_state, str):
        oxidation_state = int(oxidation_state.replace('+', '').replace('-', ''))
        if '-' in str(oxidation_state):
            oxidation_state = -abs(oxidation_state)

    # Determine charge sign suffix
    if oxidation_state > 0:
        suffix = 'p'  # positive
        magnitude = abs(oxidation_state)
    elif oxidation_state < 0:
        suffix = 'n'  # negative
        magnitude = abs(oxidation_state)
    else:
        suffix = 'z'  # zero
        magnitude = 0

    # Convert magnitude to lowercase roman numerals
    if magnitude == 0:
        roman_mag = ''
    else:
        roman_mag = roman.toRoman(magnitude).lower()

    # Construct pseudoelement name
    pseudoelement = f"{element}j{roman_mag}{suffix}"

    return pseudoelement


def parse_pseudoelement_name(pseudoelement):
    """
    Parse a pseudoelement name to extract the element and oxidation state.

    Parameters
    ----------
    pseudoelement : str
        Pseudoelement name (e.g., "Fejiiip", "Sjiin", "Sjz")

    Returns
    -------
    tuple
        (element, oxidation_state) where element is str and oxidation_state is int

    Examples
    --------
    >>> parse_pseudoelement_name("Fejiiip")
    ('Fe', 3)
    >>> parse_pseudoelement_name("Sjiin")
    ('S', -2)
    >>> parse_pseudoelement_name("Sjz")
    ('S', 0)
    """
    if 'j' not in pseudoelement or len(pseudoelement) <= 2:
        # Not a pseudoelement, return as-is
        return (pseudoelement, None)

    # Split on 'j'
    parts = pseudoelement.split('j')
    element = parts[0]
    charge_part = parts[1]

    # Determine sign from last character
    if charge_part[-1] == 'p':
        sign = 1
        roman_part = charge_part[:-1]
    elif charge_part[-1] == 'n':
        sign = -1
        roman_part = charge_part[:-1]
    elif charge_part[-1] == 'z':
        return (element, 0)
    else:
        raise ValueError(f"Invalid pseudoelement name: {pseudoelement}")

    # Convert roman numerals to integer
    if roman_part:
        magnitude = roman.fromRoman(roman_part.upper())
        oxidation_state = sign * magnitude
    else:
        oxidation_state = 0

    return (element, oxidation_state)


def extract_oxidation_states_from_thermo_df(thermo_df, element):
    """
    Extract all oxidation states for a given element from the thermodynamic database.

    Parameters
    ----------
    thermo_df : pd.DataFrame
        Thermodynamic database with 'formula_ox' column
    element : str
        Element symbol to search for

    Returns
    -------
    list
        List of unique oxidation states found for this element (as strings like "Fe+3", "Fe+2")
    """
    oxidation_states = set()

    # Iterate through all formula_ox entries
    for idx, row in thermo_df.iterrows():
        formula_ox = row.get('formula_ox', 'nan')

        if pd.isna(formula_ox) or formula_ox == 'nan' or formula_ox == '':
            continue

        # Parse the formula_ox string using shared function
        parts = formula_ox.strip().split()

        for part in parts:
            parsed = parse_formula_ox_part(part, include_pseudoelements=False)
            if parsed:
                _, elem, charge = parsed
                if elem == element:
                    if charge:
                        oxidation_states.add(f"{elem}{charge}")
                    else:
                        # No charge means oxidation state 0
                        oxidation_states.add(f"{elem}")

    return sorted(list(oxidation_states))


def create_redox_element_states_dict(thermo_df, suppress_redox):
    """
    Create a dictionary mapping elements to their pseudoelement names for each oxidation state.

    Only includes oxidation states that have aqueous basis species available.
    Species with oxidation states that lack aqueous basis species will be rejected.

    Parameters
    ----------
    thermo_df : pd.DataFrame
        Thermodynamic database with 'formula_ox' column
    suppress_redox : list
        List of element symbols to suppress redox for (e.g., ["Fe", "S"])

    Returns
    -------
    tuple
        (redox_elem_states, rejected_species) where:
        - redox_elem_states: dict {element: {oxidation_state_string: pseudoelement_name}}
        - rejected_species: list of species names that were rejected
    """
    redox_elem_states = {}
    rejected_species = []

    for element in suppress_redox:
        # Extract all oxidation states for this element
        ox_states = extract_oxidation_states_from_thermo_df(thermo_df, element)

        if len(ox_states) == 0:
            continue

        # Create pseudoelement names for each oxidation state that has an aqueous basis species
        redox_entry = {}
        ox_states_without_basis = []

        for ox_state_str in ox_states:
            # Parse oxidation state from string like "Fe+3" or just "Fe" (for zero)
            match = re.match(r'([A-Z][a-z]?)([\+\-])?([0-9]+)?', ox_state_str.strip())
            if match:
                elem, sign, magnitude = match.groups()
                if sign and magnitude:
                    ox_state = int(sign + magnitude)
                else:
                    # No sign/magnitude means oxidation state 0
                    ox_state = 0
            else:
                continue

            # Check if there's an aqueous or gaseous basis species for this oxidation state
            # EQ3/6 basis species must be aqueous or gaseous (e.g., O2(g), S2(g))
            # Look for species containing this oxidation state with tag='basis' or 'aux' and state in ['aq', 'gas']
            has_valid_basis = False
            for idx, row in thermo_df.iterrows():
                formula_ox_val = row.get('formula_ox', '')
                # Handle NaN and non-string values
                if pd.notna(formula_ox_val) and isinstance(formula_ox_val, str):
                    # Check if this species contains the element in this oxidation state
                    # Parse formula_ox to find all element-oxidation state pairs
                    contains_ox_state = False
                    parts = formula_ox_val.strip().split()
                    for part in parts:
                        # Use shared parsing function
                        parsed = parse_formula_ox_part(part, include_pseudoelements=False)
                        if parsed:
                            _, elem, charge = parsed
                            # Reconstruct the element+charge string
                            part_normalized = elem + (charge if charge else '')

                            if part_normalized == ox_state_str:
                                contains_ox_state = True
                                break

                    if (contains_ox_state and
                        (row.get('tag') == 'basis' or row.get('tag') == 'aux') and
                        row.get('state') in ['aq', 'gas']):
                        has_valid_basis = True
                        break

            if has_valid_basis:
                # Create pseudoelement name
                pseudoelement = create_pseudoelement_name(element, ox_state)
                redox_entry[ox_state_str] = pseudoelement
            else:
                ox_states_without_basis.append(ox_state_str)

        # Reject species containing oxidation states without aqueous basis species
        if ox_states_without_basis:
            for ox_state_str in ox_states_without_basis:
                # Find all species with this oxidation state
                for idx, row in thermo_df.iterrows():
                    formula_ox = row.get('formula_ox', '')
                    if pd.notna(formula_ox) and formula_ox != 'nan' and formula_ox != '':
                        # Check if this oxidation state appears in formula_ox
                        parts = formula_ox.strip().split()
                        for part in parts:
                            # Extract element from part
                            elem_match = re.match(r'[0-9]*\.?[0-9]*([A-Z][a-z]?)', part)
                            if elem_match and elem_match.group(1) == element:
                                # This part contains our element - check oxidation state
                                # For exact matching: "Fe+2" matches "Fe+2", "Fe" matches "Fe", but "Fe" does NOT match "Fe+2"
                                if ox_state_str.strip() == part.strip():
                                    if row['name'] not in rejected_species:
                                        rejected_species.append(row['name'])

        if redox_entry:
            redox_elem_states[element] = redox_entry

    return redox_elem_states, rejected_species


def add_pseudoelements_to_chnosz(redox_elem_states, element_df=None):
    """
    Add pseudoelements to CHNOSZ's element database.

    Parameters
    ----------
    redox_elem_states : dict
        Dictionary mapping elements to pseudoelement names
        Format: {element: {ox_state_string: pseudoelement_name}}
    element_df : pd.DataFrame, optional
        Element database to use as template (if None, uses CHNOSZ default)

    Returns
    -------
    None
        Modifies CHNOSZ's global element database
    """
    # Get current element database from CHNOSZ
    current_elements = pychnosz.thermo().element

    # Add each pseudoelement
    for element in redox_elem_states:
        # Get properties of the original element
        elem_row = current_elements[current_elements['element'] == element]

        if elem_row.empty:
            print(f"Warning: Element {element} not found in CHNOSZ element database")
            continue

        # Use the first matching row
        elem_row = elem_row.iloc[0]

        # Create pseudoelement entries
        for ox_state_str, pseudoelement in redox_elem_states[element].items():
            # Create new element entry with pseudoelement name
            new_elem = pd.DataFrame({
                'element': [pseudoelement],
                'state': [elem_row['state']],
                'source': [f"redox suppression for {ox_state_str}"],
                'mass': [elem_row['mass']],
                's': [elem_row['s']],
                'n': [elem_row['n']]
            })

            # Add to CHNOSZ element database
            updated_elements = pd.concat([current_elements, new_elem], ignore_index=True)
            current_elements = updated_elements

    # Update CHNOSZ with the new element database
    pychnosz.thermo(element=current_elements)


def modify_formulas_for_redox_suppression(thermo_df, redox_elem_states, verbose=1):
    """
    Modify chemical formulas in thermo_df to use pseudoelements for redox-suppressed elements.

    This creates 'formula_modded' and 'formula_ox_modded' columns where redox-suppressed
    elements are replaced with their corresponding pseudoelements.

    Parameters
    ----------
    thermo_df : pd.DataFrame
        Thermodynamic database with 'formula' and 'formula_ox' columns
    redox_elem_states : dict
        Dictionary mapping elements to pseudoelement names
    verbose : int
        Verbosity level

    Returns
    -------
    pd.DataFrame
        Modified thermodynamic database with 'formula_modded' and 'formula_ox_modded' columns
    """
    thermo_df = thermo_df.copy()

    # Initialize modded columns if they don't exist
    if 'formula_modded' not in thermo_df.columns:
        thermo_df['formula_modded'] = thermo_df['formula']

    if 'formula_ox_modded' not in thermo_df.columns:
        thermo_df['formula_ox_modded'] = thermo_df.get('formula_ox', 'nan')

    # Process each species
    for idx, row in thermo_df.iterrows():
        species_name = row['name']
        formula = row['formula']
        formula_ox = row.get('formula_ox', 'nan')

        # Handle simple ion formulas like "Fe+2" which should become "Fejiip+2"
        # Check if the formula matches the pattern: Element+charge or Element-charge
        import re
        simple_ion_match = re.match(r'^([A-Z][a-z]?)([+\-]\d+)$', formula)
        if simple_ion_match:
            element = simple_ion_match.group(1)
            charge_str = simple_ion_match.group(2)

            # Check if this element is redox-suppressed
            if element in redox_elem_states:
                # Parse formula_ox to get the oxidation state
                if pd.notna(formula_ox) and formula_ox != 'nan' and formula_ox != '':
                    # formula_ox for simple ions is typically same as formula, e.g., "Fe+2"
                    ox_state_str = formula_ox.strip()
                    if ox_state_str in redox_elem_states[element]:
                        pseudoelement = redox_elem_states[element][ox_state_str]
                        # Replace formula with pseudoelement version
                        thermo_df.at[idx, 'formula_modded'] = f"{pseudoelement}{charge_str}"
                        thermo_df.at[idx, 'formula_ox_modded'] = f"{pseudoelement}{charge_str}"
                        if verbose >= 2:
                            print(f"  Modified simple ion {species_name}: {formula} -> {pseudoelement}{charge_str}")
                        continue

        # Skip if no formula_ox
        if pd.isna(formula_ox) or formula_ox == 'nan' or formula_ox == '':
            thermo_df.at[idx, 'formula_modded'] = formula
            thermo_df.at[idx, 'formula_ox_modded'] = formula_ox
            continue

        # Parse the formula to get charge
        try:
            formula_makeup = pychnosz.makeup(formula)
            charge = formula_makeup.get('Z', 0)
        except:
            charge = 0

        # Build modified formula and formula_ox
        modified_formula_parts = {}
        modified_formula_ox_parts = []

        # Parse formula_ox using shared function
        parts = formula_ox.strip().split()

        for part in parts:
            # Use shared parsing function
            parsed = parse_formula_ox_part(part, include_pseudoelements=False)
            if not parsed:
                continue

            coeff, element, ox_charge = parsed

            # Check if this element is redox-suppressed
            if element in redox_elem_states:
                # Build oxidation state string
                if ox_charge:
                    ox_state_str = f"{element}{ox_charge}"
                else:
                    # No charge means oxidation state 0
                    ox_state_str = element

                if ox_state_str in redox_elem_states[element]:
                    pseudoelement = redox_elem_states[element][ox_state_str]

                    # Add to modified formula
                    if pseudoelement in modified_formula_parts:
                        modified_formula_parts[pseudoelement] += coeff
                    else:
                        modified_formula_parts[pseudoelement] = coeff

                    # Add to modified formula_ox
                    if ox_charge:
                        # Has oxidation state charge
                        if coeff == 1.0:
                            modified_formula_ox_parts.append(f"{pseudoelement}{ox_charge}")
                        else:
                            modified_formula_ox_parts.append(f"{coeff:.4g}{pseudoelement}{ox_charge}")
                    else:
                        # No oxidation state charge (oxidation state 0)
                        if coeff == 1.0:
                            modified_formula_ox_parts.append(pseudoelement)
                        else:
                            modified_formula_ox_parts.append(f"{coeff:.4g}{pseudoelement}")
                else:
                    # Oxidation state not in our mapping, keep original
                    if element in modified_formula_parts:
                        modified_formula_parts[element] += coeff
                    else:
                        modified_formula_parts[element] = coeff

                    if coeff == 1.0:
                        modified_formula_ox_parts.append(part)
                    else:
                        modified_formula_ox_parts.append(f"{coeff:.4g}{element}{ox_charge}")
            else:
                # Not redox-suppressed, keep original element
                if element in modified_formula_parts:
                    modified_formula_parts[element] += coeff
                else:
                    modified_formula_parts[element] = coeff

                modified_formula_ox_parts.append(part)

        # Build final formula string
        formula_parts = []
        for elem, count in modified_formula_parts.items():
            if count == 1.0:
                formula_parts.append(elem)
            else:
                formula_parts.append(f"{elem}{count:.4g}")

        # Add charge
        if charge > 0:
            if charge == 1:
                formula_parts.append('+')
            else:
                formula_parts.append(f'+{int(charge)}')
        elif charge < 0:
            if charge == -1:
                formula_parts.append('-')
            else:
                formula_parts.append(f'{int(charge)}')

        modified_formula = ''.join(formula_parts)
        modified_formula_ox = ' '.join(modified_formula_ox_parts)

        # Update dataframe
        if modified_formula:
            thermo_df.at[idx, 'formula_modded'] = modified_formula
        else:
            thermo_df.at[idx, 'formula_modded'] = formula

        if modified_formula_ox_parts:
            thermo_df.at[idx, 'formula_ox_modded'] = modified_formula_ox
        else:
            thermo_df.at[idx, 'formula_ox_modded'] = formula_ox

    return thermo_df
