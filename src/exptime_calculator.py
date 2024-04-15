#!/bin/ebv python
"""
This module contains functions to calculate exposure time for a given magnitude
for S-PLUS
Developed by: Andre Santos, 2024
"""

import numpy as np


def get_mag(texp: float, filt: str) -> float:
    """
    Given a exposure time and filter, returns a value for magnitude.

    Parameters:
    -----------

    """

    # Limiting Magnitudes for SNR > 3
    # Time per exposure. maglim is achieved after coadding 3 exposures
    maglim = {
        'u': (21.0, 3 * 227.0),
        'g': (21.3, 3 * 33.0),
        'r': (21.3, 3 * 40.0),
        'i': (20.9, 3 * 46.0),
        'z': (20.1, 3 * 56.0),
    }

    m0 = maglim[filt][0]
    t0 = maglim[filt][1]

    maglim = m0 + 2.5 * np.log10(texp/t0)

    return maglim


def get_exp_time(mag: float, filt: str) -> float:
    """
    Given a magnitude and filter, returns an exposure time.

    Parameters:
    -----------
        mag: float
            magtinude one want to achieve.
        filt: str
            telescope filter to take in account.

    Returns:
    -------
        Exposure time.
    """

    maglim = {
        'u': (21.0, 227.0),
        'g': (21.3, 33.0),
        'r': (21.3, 40.0),
        'i': (20.9, 46.0),
        'z': (20.1, 56.0),
    }

    m0 = maglim[filt][0]
    t0 = maglim[filt][1]

    t1 = 3 * t0 * 10 ** ((mag - m0) / 2.5)

    return t1


if __name__ == "__main__":
    print(get_mag(138, 'i'))
    print(get_exp_time(21.5, 'i')/3)
