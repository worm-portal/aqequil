"""
Module for checking and balancing dissociation reactions.

This module provides functionality to check if dissociation reactions are balanced
and to balance them if needed using CHNOSZ basis species.
"""

import pychnosz
import pandas as pd
import numpy as np


def check_dissrxn_balanced(dissrxn_str, species_name=None, verbose=False):
    """
    Check if a dissociation reaction is balanced.

    This uses CHNOSZ's info() and makeup() to check if the reaction balances,
    similar to the R version's subcrt_bal() function.

    Parameters
    ----------
    dissrxn_str : str
        Dissociation reaction string in format: "-1.0000 species1 1.0000 species2 ..."
    species_name : str, optional
        Name of the species being checked (for error messages)
    verbose : bool
        Whether to print detailed information

    Returns
    -------
    bool
        True if reaction is balanced, False otherwise
    """
    if not dissrxn_str or dissrxn_str.strip() == '' or dissrxn_str == 'nan':
        return False

    try:
        # Parse the dissociation reaction
        parts = dissrxn_str.strip().split()
        if len(parts) < 2:  # Need at least coeff1 species1
            return False

        coeffs = []
        species_names = []

        for i in range(0, len(parts), 2):
            try:
                coeffs.append(float(parts[i]))
                species_names.append(parts[i+1])
            except (ValueError, IndexError):
                return False

        # Get species indices from CHNOSZ
        try:
            ispecies = pychnosz.info(species_names, messages=False)
        except Exception as e:
            if verbose:
                print(f"Error getting species info: {e}")
            return False

        # Calculate mass balance using CHNOSZ's makeup
        # This is equivalent to R's makeup(ispecies, coeff, sum=TRUE)
        # Use species indices instead of names to properly look up from database
        total_makeup = {}
        for i, (idx, coeff) in enumerate(zip(ispecies, coeffs)):
            try:
                # Use index to get makeup from database
                sp_makeup = pychnosz.makeup(idx)
                for elem, count in sp_makeup.items():
                    if elem not in total_makeup:
                        total_makeup[elem] = 0.0
                    total_makeup[elem] += coeff * count
            except Exception as e:
                if verbose:
                    print(f"Error getting makeup for {species_names[i]} (index {idx}): {e}")
                return False

        # Check if all elements balance (should be ~0)
        # Take out very small numbers (matching R's abs(mss) < 1e-7)
        tolerance = 1e-7
        for elem in total_makeup:
            if abs(total_makeup[elem]) < tolerance:
                total_makeup[elem] = 0.0

        # Check if any element is non-zero (unbalanced)
        if any(abs(v) > tolerance for v in total_makeup.values()):
            if verbose:
                print(f"Unbalanced reaction for {species_name}: {total_makeup}")
            return False

        return True

    except Exception as e:
        if verbose:
            print(f"Error checking balance for {species_name}: {e}")
        return False


def balance_dissrxn(species_name, species_list, basis_species, fixed_species=None):
    """
    Balance a dissociation reaction by adding basis species.

    Parameters
    ----------
    species_name : str
        Name of the species to balance
    species_list : list
        List of species names already in the reaction
    basis_species : list
        List of basis species names to use for balancing
    fixed_species : list, optional
        List of fixed species (H2O, H+, O2, etc.)

    Returns
    -------
    dict
        Dictionary with 'species' and 'coeffs' keys, or None if balancing fails
    """
    if fixed_species is None:
        fixed_species = ["H2O", "H+", "O2(g)", "water", "e-"]

    try:
        # Set up basis
        pychnosz.basis(basis_species, messages=False)

        # Use subcrt to get balanced reaction
        result = pychnosz.subcrt([species_name], [-1], messages=False)

        if result is None or not hasattr(result, 'reaction'):
            return None

        # Extract coefficients and species names from the reaction DataFrame
        coeffs = result.reaction['coeff']
        names = result.reaction['name']

        # Replace 'water' with 'H2O' for EQ3 compatibility
        names = ['H2O' if n == 'water' else n for n in names]

        return {
            'species': names,
            'coeffs': [float(c) for c in coeffs]
        }

    except Exception as e:
        return None


def format_dissrxn_string(species_names, coeffs):
    """
    Format a dissociation reaction as a string.

    Parameters
    ----------
    species_names : list
        List of species names
    coeffs : list
        List of coefficients

    Returns
    -------
    str
        Formatted dissociation reaction string
    """
    parts = []
    for coeff, name in zip(coeffs, species_names):
        parts.append(f"{coeff:.4f}")
        parts.append(name)

    return " ".join(parts)
