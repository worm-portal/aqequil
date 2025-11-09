"""
check_TP_grid.py - Python conversion of check_TP_grid.r

Checks that the temperature(s) and pressure(s) of an aqueous speciation calculation
allow for a liquid water phase, and calculates polynomial fits.
"""

import numpy as np
from numpy.polynomial import polynomial as P
import copy


def vmessage(m, vlevel, verbose):
    """Print messages if 'verbose' setting >= vlevel of message."""
    if verbose >= vlevel:
        print(m)


def roundup(x, n):
    """
    Round up the last digit of a number.
    e.g., 52.6112 becomes 52.6113 with n=4.
    This function is useful for reporting rounded PSAT pressures that keeps water in the liquid phase.
    """
    return round(x, n) + 10**-n


def check_TP_grid(grid_temps, grid_press, P1, water_model="SUPCRT92", check_for_errors=True, verbose=1):
    """
    Check temperature-pressure grid for aqueous speciation calculations.

    Parameters
    ----------
    grid_temps : array-like
        Temperature grid in degrees C
    grid_press : array-like or str
        Pressure grid in bar, or "psat" for saturation pressure
    P1 : float
        Reference pressure for PSAT calculations
    water_model : str, optional
        Water model to use (default: "SUPCRT92")
    check_for_errors : bool, optional
        Whether to check for errors in T-P grid (default: True)
    verbose : int, optional
        Verbosity level (default: 1)

    Returns
    -------
    dict
        Dictionary with keys: "grid_temps", "grid_press", "poly_coeffs_1", "poly_coeffs_2"
    """
    import pychnosz

    # Convert to numpy array and numeric
    grid_temps = np.array(grid_temps, dtype=float)

    # Set water model to SUPCRT92 to calculate Psat because other water models like DEW won't
    pychnosz.water("SUPCRT92", messages=False)

    # Round grid temperatures to four decimal places
    grid_temps = np.round(grid_temps, 4)

    # Calculate PSAT pressure if specified by user
    if isinstance(grid_press, str) and grid_press.lower() == "psat":
        vmessage("Calculating pressure grid along liquid-vapor saturation curve...", 2, verbose)
        grid_press = []
        for T in grid_temps:
            psat_result = pychnosz.water("Psat", T=T+273.15, Psat_floor=P1, messages=False)
            grid_press.append(psat_result[0])
        grid_press = np.array(grid_press)
    else:
        grid_press = np.array(grid_press, dtype=float)

    # Check TP polynomial
    if len(grid_temps) >= 8:

        if len(grid_temps) % 2 == 0:
            # even-length TP grid
            n_mid1 = len(grid_temps) // 2 + 1
            n_mid2 = n_mid1 - 1
        else:
            # odd-length TP grid
            n_mid1 = len(grid_temps) // 2 + 1
            n_mid2 = n_mid1 - 1

        # Polynomial for the first T-P range
        temps_1 = grid_temps[0:n_mid1]
        press_1 = grid_press[0:n_mid1]
        degree_1 = len(temps_1) - 1
        poly_coeffs_1 = np.polyfit(temps_1, press_1, degree_1)
        # Reverse to match R convention (intercept first)
        poly_coeffs_1 = poly_coeffs_1[::-1]

        # Polynomial for the second T-P range
        temps_2 = grid_temps[n_mid1-1:]
        press_2 = grid_press[n_mid1-1:]
        degree_2 = len(temps_2) - 1
        poly_coeffs_2 = np.polyfit(temps_2, press_2, degree_2)
        # Reverse to match R convention (intercept first)
        poly_coeffs_2 = poly_coeffs_2[::-1]

        # Test polynomial coefficients for first half of temperature grid
        for T in grid_temps[0:n_mid2]:
            test_sum_1 = 0
            for i in range(len(poly_coeffs_1)):
                test_sum_1 = test_sum_1 + poly_coeffs_1[i] * T**(i)
            if np.isnan(test_sum_1):
                raise ValueError(f"Error: Could not compute the coefficients of an interpolating polynomial "
                               f"for the first half of the values of the temperature grid: "
                               f"[{', '.join(map(str, grid_temps[0:n_mid2]))}].")

            test_sum_2 = 0
            for i in range(len(poly_coeffs_2)):
                test_sum_2 = test_sum_2 + poly_coeffs_2[i] * T**(i)
            if np.isnan(test_sum_2):
                raise ValueError(f"Error: Could not compute the coefficients of an interpolating polynomial "
                               f"for the last half of the values of the temperature grid: "
                               f"[{', '.join(map(str, grid_temps[n_mid2:]))}].")

        # Convert to list for consistency with R output
        poly_coeffs_1 = poly_coeffs_1.tolist()
        poly_coeffs_2 = poly_coeffs_2.tolist()

    else:
        poly_coeffs_1 = 'None'
        poly_coeffs_2 = 'None'

    # reset water model
    pychnosz.water(water_model, messages=False)
    
    if check_for_errors:
        # Check that water is a liquid at each T-P point
        TP_grid_errors = []
        for i in range(len(grid_temps)):
            if len(grid_press) == 1 and i == 1:
                # break on 2nd T-P point in dynamic db
                break

            psat_press = None
            try:
                psat_press = pychnosz.water("Psat", T=grid_temps[i]+273.15, Psat_floor=P1, messages=False)
            except:
                psat_press = None

            if psat_press is None or np.isnan(psat_press):
                rho_val = None
                try:
                    rho_val = pychnosz.water("rho", T=grid_temps[i]+273.15, P=grid_press[i], Psat_floor=P1, messages=False)
                except:
                    rho_val = None

                if rho_val is None or np.isnan(rho_val):
                    TP_grid_errors.append(
                        f"\nWater density could not be calculated at "
                        f"{grid_temps[i]} degrees C and "
                        f"{grid_press[i]} bar with the water model {water_model}"
                    )
                elif rho_val <= 310:
                    TP_grid_errors.append(
                        f"\nWater density is {roundup(rho_val, 3)} kg m^3 at "
                        f"{grid_temps[i]} degrees C and "
                        f"{grid_press[i]} bar. This is too low (< 310 kg m^3) "
                        f"for this calculation. Increase pressure or decrease temperature."
                    )
            elif grid_press[i] < psat_press:
                TP_grid_errors.append(
                    f"\n{grid_press[i]} bar is below liquid-vapor "
                    f"saturation pressure {roundup(psat_press, 4)} "
                    f"bar at {grid_temps[i]} degrees C."
                )

        if len(TP_grid_errors) > 0:
            error_msg = "\n".join(TP_grid_errors)
            error_msg += "\n\nIncrease the pressure at these temperature points to keep water in a liquid state."
            raise ValueError(error_msg)
    
    return {
        "grid_temps": grid_temps.tolist(),
        "grid_press": grid_press.tolist(),
        "poly_coeffs_1": poly_coeffs_1,
        "poly_coeffs_2": poly_coeffs_2
    }
