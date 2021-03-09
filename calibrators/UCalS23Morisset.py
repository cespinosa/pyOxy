#!/usr/bin/env python3

import numpy as np

def U_Mor16_S23_cal(S23):
    '''
    Ionization Parameter calibrator based on OH abundance

    ref: 2016Morisset_aa594

    Parameters
    ----------
    S23 : float
        SIIl6716+31/SIIIl6312 flux

    Returns
    -------
    U : float
        Ionization Parameter
    '''
    if isinstance(S23, (int, float)):
        if S23 > 0:
            U = -2.62 - 1.22 * np.log10(S23)
        else:
            U = np.nan
        return U
    mask = S23 > 0
    U = np.empty(S23.shape)
    U[:] = np.nan
    U[mask] = -2.62 - 1.22 * np.log10(S23[mask])
    return U
