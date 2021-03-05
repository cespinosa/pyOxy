#!/usr/bin/env python3

import numpy as np
from scipy.interpolate import interp1d

def OH_C20_R3_cal(R3):
    '''
    Oxygen Abundnace calibrator based on R3 ratio

    ref = 2020Curti_mnras491

    Parameters
    ----------
    R3 : float
        OIII5007/Hb ratio

    Returns
    -------
    OH : float
        Oxygen Abundance: 12 + log(O/H)
    '''
    R = lambda x: -0.277 - 3.549*x - 3.593*x**2 - 0.981*x**3
    npts = 1000
    z = np.linspace(7.65, 8.85, npts)
    x = z-8.69
    y = R(x)
    f = interp1d(y, x)
    if isinstance(R3, (int, float)):
        if R3 > 0:
            l_R3 = np.log10(RS32)
            OH = f(l_R3) + 8.69
        else:
            OH = np.nan
        return OH
    mask_R3 = R3 > 0
    mask_max = np.log10(R3) <= y.max()
    mask_min = np.log10(R3) >= y.min()
    mask = mask_max & mask_min
    OH = np.empty(R3.shape)
    OH[:] = np.nan
    OH[mask_R3 & mask] = f(np.log10(R3[mask_R3 & mask])) + 8.69
    return OH
