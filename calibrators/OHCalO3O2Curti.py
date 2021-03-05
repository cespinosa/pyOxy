#!/usr/bin/env python3

import numpy as np
from scipy.interpolate import interp1d

def OH_C20_O3O2_cal(O3O2):
    '''
    Oxygen Abundnace calibrator based on O3O2 ratio

    ref = 2020Curti_mnras491

    Parameters
    ----------
    O3O2 : float
        (OIII5007/OII3727,29) ratio

    Returns
    -------
    OH : float
        Oxygen Abundance: 12 + log(O/H)
    '''
    R = lambda x: -0.691 - 2.944*x - 1.308*x**2
    npts = 1000
    z = np.linspace(7.65, 8.85, npts)
    x = z-8.69
    y = R(x)
    f = interp1d(y, x)
    if isinstance(O3O2, (int, float)):
        if O3O2 > 0:
            l_O3O2 = np.log10(O3O2)
            OH = f(l_O3O2) + 8.69
        else:
            OH = np.nan
        return OH
    mask_O3O2 = O3O2 > 0
    mask_max = np.log10(O3O2) <= y.max()
    mask_min = np.log10(O3O2) >= y.min()
    mask = mask_max & mask_min
    OH = np.empty(O3O2.shape)
    OH[:] = np.nan
    OH[mask_O3O2 & mask] = f(np.log10(O3O2[mask_O3O2 & mask])) + 8.69
    return OH
