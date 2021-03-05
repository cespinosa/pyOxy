#!/usr/bin/env python3

import numpy as np
from scipy.interpolate import interp1d

def OH_C20_RS32_cal(RS32):
    '''
    Oxygen Abundnace calibrator based on RS32 ratio
    This is a bivaluated function with a maximun at z=7.999337,
    x=z-8.69=-0.690663 (RS32=0.788207)

    ref = 2020Curti_mnras491

    Parameters
    ----------
    RS32 : float
        (SII6717,31/Ha)+(OIII5007/Hb) ratio

    Returns
    -------
    OH : float
        Oxygen Abundance: 12 + log(O/H)
    '''
    R = lambda x: -0.054 - 2.546*x - 1.970*x**2 + 0.082*x**3 + 0.222*x**4
    npts = 1000
    z = np.linspace(7.65, 8.85, npts)
    x = z-8.69
    y = R(x)
    f = interp1d(y, x)
    if isinstance(RS32, (int, float)):
        if RS32 > 0:
            l_RS32 = np.log10(RS32)
            OH = f(l_RS32) + 8.69
        else:
            OH = np.nan
        return OH
    mask_RS32 = RS32 > 0
    mask_max = np.log10(RS32) <= y.max()
    mask_min = np.log10(RS32) >= y.min()
    mask = mask_max & mask_min
    OH = np.empty(RS32.shape)
    OH[:] = np.nan
    OH[mask_RS32 & mask] = f(np.log10(RS32[mask_RS32 & mask])) + 8.69
    return OH
