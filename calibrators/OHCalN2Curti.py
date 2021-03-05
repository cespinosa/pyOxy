#!/usr/bin/env python3

import numpy as np
from scipy.interpolate import interp1d

def OH_C20_N2_cal(N2):
    '''
    Oxygen Abundnace calibrator based on N2 ratio

    ref = 2020Curti_mnras491

    Parameters
    ----------
    N2 : float
        NII6584/Ha ratio

    Returns
    -------
    OH : float
        Oxygen Abundance: 12 + log(O/H)
    '''
    R = lambda x: -0.489 + 1.513*x - 2.554*x**2 - 5.293*x**3 - 2.867*x**4
    npts = 1000
    z = np.linspace(7.65, 8.85, npts)
    x = z-8.69
    y = R(x)
    f = interp1d(y, x)
    if isinstance(N2, (int, float)):
        if N2 > 0:
            l_N2 = np.log10(N2)
            OH = f(l_N2) + 8.69
        else:
            OH = np.nan
        return OH
    mask_N2 = N2 > 0
    mask_max = np.log10(N2) <= y.max()
    mask_min = np.log10(N2) >= y.min()
    mask = mask_max & mask_min
    OH = np.empty(N2.shape)
    OH[:] = np.nan
    OH[mask_N2 & mask] = f(np.log10(N2[mask_N2 & mask])) + 8.69
    return OH
