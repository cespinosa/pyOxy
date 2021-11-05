#!/usr/bin/env python3

import numpy as np
from scipy.interpolate import interp1d


def OH_C20_R3_cal(R3, N2O2=None):
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
    def R(x): return -0.277 - 3.549*x - 3.593*x**2 - 0.981*x**3
    npts = 1000
    if isinstance(N2O2, (int, float)):
        if N2O2 is None:
            lower_lim = 7.65
            upper_lim = 8.85
        else:
            log_n2o2 = np.log10(N2O2)
            if log_n2o2 < -1.2:
                lower_lim = 7.65
                upper_lim = 8.002610729099287
            else:
                lower_lim = 8.002610729099287
                upper_lim = 8.85
        z = np.linspace(lower_lim, upper_lim, npts)
        x = z-8.69
        y = R(x)
        f = interp1d(y, x)
    if isinstance(R3, (int, float)):
        if R3 > 0:
            l_R3 = np.log10(R3)
            OH = f(l_R3) + 8.69
        else:
            OH = np.nan
        return OH
    shape_map = R3.shape
    OH = np.array([])
    for r3_value, n2o2_value in zip(R3.ravel(), N2O2.ravel()):
        if n2o2_value is None:
            lower_lim = 7.65
            upper_lim = 8.85
        else:
            log_n2o2 = np.log10(n2o2_value)
            if log_n2o2 < -1.2:
                lower_lim = 7.65
                upper_lim = 8.002610729099287
            else:
                lower_lim = 8.002610729099287
                upper_lim = 8.85
        z = np.linspace(lower_lim, upper_lim, npts)
        x = z-8.69
        y = R(x)
        line_upper_lim = y.max()
        line_lower_lim = y.min()
        f = interp1d(y, x)
        if r3_value > 0:
            l_r3 = np.log10(r3_value)
            if ((l_r3 >= line_lower_lim) and (l_r3 <= line_upper_lim)):
                OH_value = f(l_r3) + 8.69
            else:
                OH_value = np.nan
        else:
            OH_value = np.nan
        OH = np.append(OH, OH_value)
    # mask_R3 = R3 > 0
    # mask_max = np.log10(R3) <= y.max()
    # mask_min = np.log10(R3) >= y.min()
    # mask = mask_max & mask_min
    # OH = np.empty(R3.shape)
    # OH[:] = np.nan
    # OH[mask_R3 & mask] = f(np.log10(R3[mask_R3 & mask])) + 8.69
    OH_ima = OH.reshape(shape_map)
    return OH_ima
