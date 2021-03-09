#!/usr/bin/env python3

import numpy as np

def U_Dors_S_cal(SIIHa):
    '''
    Ionization Parameter calibrator based on OIII/OII ratio

    ref: 2011Dors_mnras415

    Parameters
    ----------
    SIIHa : float
        (SIIl6716 + SIIl6730) / Ha ratio

    Returns
    -------
    U : float
        Ionization Parameter
    '''
    if isinstance(SIIHa, (int, float)):
        if SIIHa > 0:
            U = -1.66 * np.log10(SIIHa) - 4.13
        else:
            U = np.nan
        return U

    mask_SIIHa = SIIHa > 0
    U = np.empty(SIIHa.shape)
    U[:] = np.nan
    U[mask_SIIHa] = -1.66 * np.log10(SIIHa[mask_SIIHa]) - 4.13
    return U
