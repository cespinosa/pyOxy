#!/usr/bin/env python3

import numpy as np
from scipy.interpolate import interp1d

def OH_C20_R2_cal(R2):
    '''
    Oxygen Abundnace calibrator based on R2 ratio

    ref = 2020Curti_mnras491

    Parameters
    ----------
    R2 : float
        OII3727,29/Hb ratio

    Returns
    -------
    OH : float
        Oxygen Abundance: 12 + log(O/H)
    '''
    R = lambda x: 0.435 - 1.362*x - 5.655*x**2 - 4.851*x**3 - 0.478*x**4 + 0.736*x**5
    npts = 1000
    z = np.linspace(7.65, 8.85, npts)
    x = z-8.69
    y = R(x)
    f = interp1d(y, x)
    if isinstance(R2, (int, float)):
        if R2 > 0:
            l_R2 = np.log10(R2)
            OH = f(l_R2) + 8.69
        else:
            OH = np.nan
        return OH
    mask_R2 = R2 > 0
    mask_max = np.log10(R2) <= y.max()
    mask_min = np.log10(R2) >= y.min()
    mask = mask_max & mask_min
    OH = np.empty(R2.shape)
    OH[:] = np.nan
    OH[mask_R2 & mask] = f(np.log10(R2[mask_R2 & mask])) + 8.69
    return OH
