"""
HKF and CGL equation of state wrapper for AqEquil.

This module imports equation of state functions from the chnosz package
rather than redefining them. This eliminates code duplication and ensures
we're using the well-tested, vectorized implementations from chnosz.

All functions previously defined in this file are now imported from:
- chnosz.models.hkf (hkf, gfun, convert_cm3bar)
- chnosz.models.cgl (cgl, quartz_coesite)
- chnosz.models.berman (Berman)
- chnosz.models.hkf_helpers (calc_logK, calc_G_TP, G2logK, dissrxn2logK, OBIGT2eos)

Note: The chnosz implementations are vectorized and handle edge cases better
than the original code. They are also actively maintained.
"""

# Import all required functions from chnosz
from pychnosz import (
    # HKF equation of state functions
    hkf,
    gfun,

    # CGL equation of state functions
    cgl,
    quartz_coesite,

    # Berman equation of state functions
    Berman,

    # Helper functions for logK calculations
    calc_logK,
    calc_G_TP,
    G2logK,
    dissrxn2logK,

    # OBIGT database processing
    OBIGT2eos,

    # Water properties
    water,

    # Formula utilities
    entropy
)

# For backwards compatibility, ensure convert_cm3bar is available
def convert_cm3bar(value):
    """Convert cm3*bar to calories (SUPCRT92 convention)."""
    return value * 4.184 * 10

# Re-export all functions so existing imports continue to work
__all__ = [
    'hkf',
    'gfun',
    'cgl',
    'quartz_coesite',
    'calc_logK',
    'calc_G_TP',
    'G2logK',
    'dissrxn2logK',
    'OBIGT2eos',
    'convert_cm3bar',
    'water',
    'entropy',
    'Berman'
]
