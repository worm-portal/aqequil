"""
Bundled thermodynamic databases for aqequil.

This directory contains pre-compiled data1 files for use with EQ3/6.
The data1.wrm file is generated from data0.wrm during the wheel build process.
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
