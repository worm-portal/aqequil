"""
Test data for aqequil integration tests.

This directory contains sample input files used for testing the package.
"""

import os
from pathlib import Path


def get_test_data_path(filename):
    """
    Get the absolute path to a test data file.

    Parameters
    ----------
    filename : str
        Name of the test data file.

    Returns
    -------
    str
        Absolute path to the test data file.
    """
    return os.path.join(os.path.dirname(__file__), filename)
