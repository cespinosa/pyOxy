#!/usr/bin/env python3

import numpy as np
from scipy.interpolate import interp1d


def OH_C20_S2_cal(S2, N2O2=None):
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
    def R(x): return -0.442 - 0.360*x - 6.271*x**2 - 8.339*x**3 - 3.559*x**4
    npts = 1000
    if isinstance(N2O2, (int, float)):
        if N2O2 is None:
            lower_lim = 7.65
            upper_lim = 8.85
        else:
            log_n2o2 = np.log10(N2O2)
            if log_n2o2 < -1.2:
                lower_lim = 7.65
                upper_lim = 8.659469502260592
            else:
                lower_lim = 8.659469502260592
                upper_lim = 8.85
        z = np.linspace(lower_lim, upper_lim, npts)
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
    shape_map = S2.shape
    OH = np.array([])
    for s2_value, n2o2_value in zip(S2.ravel(), N2O2.ravel()):
        if n2o2_value is None:
            lower_lim = 7.65
            upper_lim = 8.85
        else:
            log_n2o2 = np.log10(n2o2_value)
            if log_n2o2 < -1.2:
                lower_lim = 7.65
                upper_lim = 8.659469502260592
            else:
                lower_lim = 8.659469502260592
                upper_lim = 8.85
        z = np.linspace(lower_lim, upper_lim, npts)
        x = z-8.69
        y = R(x)
        line_upper_lim = y.max()
        line_lower_lim = y.min()
        f = interp1d(y, x)
        if s2_value > 0:
            l_s2 = np.log10(s2_value)
            if ((l_s2 >= line_lower_lim) and (l_s2 <= line_upper_lim)):
                OH_value = f(l_s2) + 8.69
            else:
                OH_value = np.nan
        else:
            OH_value = np.nan
        OH = np.append(OH, OH_value)
    # mask_S2 = S2 > 0
    # mask_max = np.log10(S2) <= y.max()
    # mask_min = np.log10(S2) >= y.min()
    # mask = mask_max & mask_min
    # OH = np.empty(S2.shape)
    # OH[:] = np.nan
    # OH[mask_S2 & mask] = f(np.log10(S2[mask_S2 & mask])) + 8.69
    OH_ima = OH.reshape(shape_map)
    return OH_ima
