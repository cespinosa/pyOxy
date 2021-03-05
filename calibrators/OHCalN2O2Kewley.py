#!/usr/bin/env python3

import numpy as np

def OH_Kew02_N2O2_cal(R):
    '''
    Oxygen Abundance calibrator base on the Kewley & Dopita's work (2002)

    ref: 2002Kewley_ApJSS142

    Parameters
    ----------
    R : float
       NIIl6564/OIIl3727+29

    Returns
    -------
    OH : float
        Oxygen Abundance: 12 + log(O/H)
    '''
    if isinstance(R, (int, float)):
        if R > 0:
            l_R = np.log10(R)
            tmp = 1.54020 + 1.26602 * l_R + 0.167977 * l_R**2
            if tmp > 0:
                OH = 8.93 + np.log10(tmp)
            else:
                OH = np.nan
            return OH
    R = np.array(R)
    l_R = np.log10(R)
    tmp = 1.54020 + 1.26602 * l_R + 0.167977 * l_R**2
    tmp = np.array(tmp)
    mask = tmp > 0
    OH = np.empty(R.shape)
    OH[:] = np.nan
    OH[mask] = 8.93 + np.log10(tmp[mask])
    return OH
