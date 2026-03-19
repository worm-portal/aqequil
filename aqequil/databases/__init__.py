"""
Bundled thermodynamic databases for aqequil.

This directory contains thermodynamic database files downloaded from WORM-db
(https://github.com/worm-portal/WORM-db) during the wheel build process:

- data0.wrm: Text-format thermodynamic database (source)
- data1.wrm: Binary-format database compiled by EQPT (used by EQ3NR)
- wrm_data_latest.csv: Main thermodynamic data in CSV format
- elements.csv: Element properties database
- solid_solutions.csv: Solid solution parameters
- wrm_data_logK.csv: Equilibrium constants (logK) database
- wrm_data_logK_S.csv: Equilibrium constants (logK) and entropy (deltaS) for aqueous species
- speciation_groups_WORM.txt: Speciation group definitions for plotting

The data1.wrm file is generated from data0.wrm by running EQPT.
Most files are downloaded directly from WORM-db during wheel builds.
"""

import os


def get_database_path(filename=None):
    """
    Get the absolute path to the databases directory or a specific database file.

    Parameters
    ----------
    filename : str, optional
        Name of the database file (e.g., "data1.wrm").
        If None, returns the path to the databases directory.

    Returns
    -------
    str
        Absolute path to the databases directory or file.
    """
    if filename is None:
        return os.path.dirname(__file__)
    return os.path.join(os.path.dirname(__file__), filename)
