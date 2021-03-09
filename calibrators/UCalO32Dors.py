#!/usr/bin/env python3

import numpy as np

def U_Dors_O32_cal(O32):
    '''
    Ionization Parameter calibrator based on OIII/OII ratio

    ref: 2011Dors_mnras415

    Parameters
    ----------
    O32 : float
        O32 ratio

    Returns
    -------
    U : float
        Ionization Parameter
    '''
    U_f = lambda l_O32: 1.22 * l_O32 - 2.25
    if isinstance(O32, (int, float)):
        if O32 > 0:
            l_O32 = np.log10(O32)
            U = U_f(l_O32)
        else:
            U = np.nan
        return U
    mask_O32 = O32 > 0
    U = np.empty(O32.shape)
    U[:] = np.nan
    U[mask_O32] = U_f(np.log10(O32[mask_O32]))
    return U
