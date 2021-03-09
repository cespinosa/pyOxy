#!/usr/bin/env python3

import numpy as np

def U_Mor16_O23_ts_cal(O23):
    '''
    Ionization Parameter calibrator based on OH abundance

    ref: 2016Morisset_aa594

    Parameters
    ----------
    O23 : float
        OIIl3727/OIIIl5007 flux

    Returns
    -------
    U : float
        Ionization Parameter
    '''
    if isinstance(O23, (int, float)):
        if O23 > 0:
            U = -2.74 - 1.00 * np.log10(O23)
        else:
            U = np.nan
        return U
    mask = O23 > 0
    U = np.empty(O23.shape)
    U[:] = np.nan
    U[mask] = -2.74 - 1.00 * np.log10(O23[mask])
    return U

def U_Mor16_O23_fs_cal(O23):
    '''
    Ionization Parameter calibrator based on OH abundance

    ref: 2016Morisset_aa594

    Parameters
    ----------
    O23 : float
        OIIl3727/OIIIl5007 flux

    Returns
    -------
    U : float
        Ionization Parameter
    '''
    if isinstance(O23, (int, float)):
        if O23 > 0:
            U = -2.38 - 2.36 * np.log10(O23)
        else:
            U = np.nan
        return U
    mask = O23 > 0
    U = np.empty(O23.shape)
    U[:] = np.nan
    U[mask] = -2.38 - 2.36 * np.log10(O23[mask])
    return U
