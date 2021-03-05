#!/usr/bin/env python3

import numpy as np

def OH_M13_O3N2_cal(O3N2):
    '''
    Oxygen Abundance calibrator based on O3N2 ratio
    O3N2 = OIIIl5007/Hb * Ha/NIIl6583

    ref: 2013Marino_AA559

    Parameters
    ----------
    O3N2 : float
        OIII/Hb * Ha/NIIl6583

    Returns
    -------
    OH : float
        Oxygen Abundace: 12 + log(O/H)
    '''
    OH_f = lambda l_O3N2: 8.533 - 0.214 * l_O3N2
    if isinstance(O3N2, (int, float)):
        if (np.log10(O3N2) > -1.1) & (np.log10(O3N2) < 1.7):
            l_O3N2 = np.log10(O3N2)
            OH = OH_f(l_O3N2)
        else:
            OH = np.nan
        return OH
    mask_O3N2_g = np.log10(O3N2) > -1.1
    mask_O3N2_l = np.log10(O3N2) < 1.7
    mask = mask_O3N2_g & mask_O3N2_l
    O3N2 = np.array(O3N2)
    OH = np.empty(O3N2.shape)
    OH[:] = np.nan
    OH[mask] = OH_f(np.log10(O3N2[mask]))
    return OH
