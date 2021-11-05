#!/usr/bin/env python3

import numpy as np
from scipy.interpolate import interp1d


def OH_C20_RS32_cal(RS32, N2O2=None):
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
    def R(x): return -0.054 - 2.546*x - 1.970*x**2 + 0.082*x**3 + 0.222*x**4
    npts = 1000
    if isinstance(N2O2, (int, float)):
        if N2O2 is None:
            lower_lim = 7.65
            upper_lim = 8.85
        else:
            log_n2o2 = np.log10(N2O2)
            if log_n2o2 < -1.2:
                lower_lim = 7.65
                upper_lim = 7.999337
            else:
                lower_lim = 7.999337
                upper_lim = 8.85
        z = np.linspace(lower_lim, upper_lim, npts)
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
    shape_map = RS32.shape
    OH = np.array([])
    for rs23_value, n2o2_value in zip(RS32.ravel(), N2O2.ravel()):
        if n2o2_value is None:
            lower_lim = 7.65
            upper_lim = 8.85
        else:
            log_n2o2 = np.log10(n2o2_value)
            if log_n2o2 < -1.2:
                lower_lim = 7.65
                upper_lim = 7.999337
            else:
                lower_lim = 7.999337
                upper_lim = 8.85
        z = np.linspace(lower_lim, upper_lim, npts)
        x = z-8.69
        y = R(x)
        line_upper_lim = y.max()
        line_lower_lim = y.min()
        f = interp1d(y, x)
        if rs23_value > 0:
            l_rs23 = np.log10(rs23_value)
            if ((l_rs23 >= line_lower_lim) and (l_rs23 <= line_upper_lim)):
                OH_value = f(l_rs23) + 8.69
            else:
                OH_value = np.nan
        else:
            OH_value = np.nan
        OH = np.append(OH, OH_value)
    # mask_RS32 = RS32 > 0
    # mask_max = np.log10(RS32) <= y.max()
    # mask_min = np.log10(RS32) >= y.min()
    # mask = mask_max & mask_min
    # OH = np.empty(RS32.shape)
    # OH[:] = np.nan
    # OH[mask_RS32 & mask] = f(np.log10(RS32[mask_RS32 & mask])) + 8.69
    OH_ima = OH.reshape(shape_map)
    return OH_ima
