#!/usr/bin/env python3

import numpy as np
from scipy.interpolate import interp1d

def OH_C20_O3S2_cal(O3S2):
    '''
    Oxygen Abundnace calibrator based on O3S2 ratio

    ref = 2020Curti_mnras491

    Parameters
    ----------
    O3S2 : float
        (OIII5007/Hb)/(SII6717,31/Ha) ratio

    Returns
    -------
    OH : float
        Oxygen Abundance: 12 + log(O/H)
    '''
    R = lambda x:  0.191 - 4.292*x - 2.538*x**2 + 0.053*x**3 + 0.332*x**4
    npts = 1000
    z = np.linspace(7.65, 8.85, npts)
    x = z-8.69
    y = R(x)
    f = interp1d(y, x)
    if isinstance(O3S2, (int, float)):
        if O3S2 > 0:
            l_O3S2 = np.log10(O3S2)
            OH = f(l_O3S2) + 8.69
        else:
            OH = np.nan
        return OH
    mask_O3S2 = O3S2 > 0
    mask_max = np.log10(O3S2) <= y.max()
    mask_min = np.log10(O3S2) >= y.min()
    mask = mask_max & mask_min
    OH = np.empty(O3S2.shape)
    OH[:] = np.nan
    OH[mask_O3S2 & mask] = f(np.log10(O3S2[mask_O3S2 & mask])) + 8.69
    return OH
