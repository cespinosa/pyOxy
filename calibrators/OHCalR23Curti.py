#!/usr/bin/env python3

import numpy as np
from scipy.interpolate import interp1d


def OH_C20_R23_cal(R23, N2O2=None):
    '''
    Oxygen Abundnace calibrator based on R23 ratio

    ref = 2020Curti_mnras491

    Parameters
    ----------
    R23 : float
        (OII3727,29 + OIII4959,5007)/Hb ratio

    Returns
    -------
    OH : float
        Oxygen Abundance: 12 + log(O/H)
    '''
    def R(x): return 0.527 - 1.569*x - 1.652*x**2 - 0.421*x**3
    npts = 1000
    if isinstance(N2O2, (int, float)):
        if N2O2 is None:
            lower_lim = 7.65
            upper_lim = 8.85
        else:
            log_n2o2 = np.log10(N2O2)
            if log_n2o2 < -1.2:
                lower_lim = 7.65
                upper_lim = 8.066528864124953
            else:
                lower_lim = 8.066528864124953
                upper_lim = 8.85
        z = np.linspace(lower_lim, upper_lim, npts)
        x = z-8.69
        y = R(x)
        f = interp1d(y, x)
    if isinstance(R23, (int, float)):
        if R23 > 0:
            l_R23 = np.log10(R23)
            OH = f(l_R23) + 8.69
        else:
            OH = np.nan
        return OH
    shape_map = R23.shape
    OH = np.array([])
    for r23_value, n2o2_value in zip(R23.ravel(), N2O2.ravel()):
        if n2o2_value is None:
            lower_lim = 7.65
            upper_lim = 8.85
        else:
            log_n2o2 = np.log10(n2o2_value)
            if log_n2o2 < -1.2:
                lower_lim = 7.65
                upper_lim = 8.066528864124953
            else:
                lower_lim = 8.066528864124953
                upper_lim = 8.85
        z = np.linspace(lower_lim, upper_lim, npts)
        x = z-8.69
        y = R(x)
        line_upper_lim = y.max()
        line_lower_lim = y.min()
        f = interp1d(y, x)
        if r23_value > 0:
            l_r23 = np.log10(r23_value)
            if ((l_r23 >= line_lower_lim) and (l_r23 <= line_upper_lim)):
                OH_value = f(l_r23) + 8.69
            else:
                OH_value = np.nan
        else:
            OH_value = np.nan
        OH = np.append(OH, OH_value)
    # mask_R23 = R23 > 0
    # mask_max = np.log10(R23) <= y.max()
    # mask_min = np.log10(R23) >= y.min()
    # mask = mask_max & mask_min
    # OH = np.empty(R23.shape)
    # OH[:] = np.nan
    # OH[mask_R23 & mask] = f(np.log10(R23[mask_R23 & mask])) + 8.69
    OH_ima = OH.reshape(shape_map)
    return OH_ima
