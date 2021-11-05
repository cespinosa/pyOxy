#!/usr/bin/env python3

import numpy as np
from scipy.interpolate import interp1d


def OH_C20_R2_cal(R2, N2O2=None):
    '''
    Oxygen Abundnace calibrator based on R2 ratio

    ref = 2020Curti_mnras491

    Parameters
    ----------
    R2 : float
        OII3727,29/Hb ratio

    N2O2 : float
        NII6584/OII3727 ratio

    Returns
    -------
    OH : float
        Oxygen Abundance: 12 + log(O/H)
    '''
    def R(x): return 0.435 - 1.362*x - 5.655 * \
        x**2 - 4.851*x**3 - 0.478*x**4 + 0.736*x**5
    npts = 1000
    if isinstance(N2O2, (int, float)):
        if N2O2 is None:
            lower_lim = 7.65
            upper_lim = 8.85
        else:
            log_n2o2 = np.log10(N2O2)
            if log_n2o2 < -1.2:
                lower_lim = 7.65
                upper_lim = 8.542149770955938
            else:
                lower_lim = 8.542149770955938
                upper_lim = 8.85
        z = np.linspace(lower_lim, upper_lim, npts)
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
    shape_map = R2.shape
    OH = np.array([])
    for r2_value, n2o2_value in zip(R2.ravel(), N2O2.ravel()):
        if n2o2_value is None:
            lower_lim = 7.65
            upper_lim = 8.85
        else:
            log_n2o2 = np.log10(n2o2_value)
            if log_n2o2 < -1.2:
                lower_lim = 7.65
                upper_lim = 8.542149770955938
            else:
                lower_lim = 8.542149770955938
                upper_lim = 8.85
        z = np.linspace(lower_lim, upper_lim, npts)
        x = z-8.69
        y = R(x)
        line_upper_lim = y.max()
        line_lower_lim = y.min()
        f = interp1d(y, x)
        if r2_value > 0:
            l_r2 = np.log10(r2_value)
            if ((l_r2 >= line_lower_lim) and (l_r2 <= line_upper_lim)):
                OH_value = f(l_r2) + 8.69
            else:
                OH_value = np.nan
        else:
            OH_value = np.nan
        OH = np.append(OH, OH_value)
    # mask_R2 = R2 > 0
    # mask_max = np.log10(R2) <= y.max()
    # mask_min = np.log10(R2) >= y.min()
    # mask = mask_max & mask_min
    # OH = np.empty(R2.shape)
    # OH[:] = np.nan
    # OH[mask_R2 & mask] = f(np.log10(R2[mask_R2 & mask])) + 8.69
    OH_ima = OH.reshape(shape_map)
    return OH_ima
