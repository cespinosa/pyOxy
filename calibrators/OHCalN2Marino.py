#!/usr/bin/env python3

import numpy as np

def OH_M13_N2_cal(N2):
    '''
    Oxygen Abundance calibrator based on N2 ratio
    N2 = NIIl6583/Ha

    ref: 2013Marino_AA559

    Parameters
    ----------
    N2 : float
        NIIl6583/Ha ratio

    Returns
    -------
    OH : float
        Oxygen Abundace: 12 + log(O/H)
    '''
    OH_f = lambda l_N2: 8.743 + 0.462 * l_N2
    if isinstance(N2, (int, float)):
        if (np.log10(N2) > -1.6) & (np.log10(N2) < -0.2):
            l_N2 = np.log10(N2)
            OH = OH_f(l_N2)
        else:
            OH = np.nan
        return OH
    mask_N2_g = np.log10(N2) > -1.6
    mask_N2_l = np.log10(N2) < -0.2
    mask = mask_N2_g & mask_N2_l
    N2 = np.array(N2)
    OH = np.empty(N2.shape)
    OH[:] = np.nan
    OH[mask] = OH_f(np.log10(N2[mask]))
    return OH
