#!/usr/bin/env python3

import numpy as np
from scipy.interpolate import interp1d

def OH_C20_S2_cal(S2):
    '''
    Oxygen Abundnace calibrator based on S2 ratio

    ref = 2020Curti_mnras491

    Parameters
    ----------
    S2 : float
        SII6717,31/Ha ratio

    Returns
    -------
    OH : float
        Oxygen Abundance: 12 + log(O/H)
    '''
    R = lambda x: -0.442 - 0.360*x - 6.271*x**2 - 8.339*x**3 - 3.559*x**4
    npts = 1000
    z = np.linspace(7.65, 8.85, npts)
    x = z-8.69
    y = R(x)
    f = interp1d(y, x)
    if isinstance(S2, (int, float)):
        if S2 > 0:
            l_S2 = np.log10(S2)
            OH = f(l_S2) + 8.69
        else:
            OH = np.nan
        return OH
    mask_S2 = S2 > 0
    mask_max = np.log10(S2) <= y.max()
    mask_min = np.log10(S2) >= y.min()
    mask = mask_max & mask_min
    OH = np.empty(S2.shape)
    OH[:] = np.nan
    OH[mask_S2 & mask] = f(np.log10(S2[mask_S2 & mask])) + 8.69
    return OH
