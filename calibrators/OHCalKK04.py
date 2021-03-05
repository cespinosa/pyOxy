#!/usr/bin/env python3

import numpy as np
from scipy.constants import c # m/s

def KK_cal_l_q(O32, OH):
    y = np.log10(O32)
    z = OH
    l_q = (32.81 - 1.153 * y**2 \
           + z * (-3.396 - 0.025 * y + 0.1444 * y**2)) \
           * (4.603 - 0.3119 * y - 0.163 * y**2 \
              + z * (-0.48 + 0.0271 * y + 0.02037 * y**2))**(-1)
    return l_q

def KK_cal_OH_lower(R23, l_q):
    x = np.log10(R23)
    OH = 9.40 + 4.65 * x - 3.17 * x**2 \
        - l_q * (0.272 + 0.547 * x - 0.513 * x**2)
    return OH

def KK_cal_OH_upper(R23, l_q):
    x = np.log10(R23)
    OH = 9.72 - 0.777 * x - 0.951 * x**2 \
        - 0.072 * x**3 - 0.811 * x**4 \
        - l_q * (0.0737 - 0.0713 * x \
                 - 0.141 * x**2 + 0.0373 * x**3 - 0.058 * x**4)
    return OH

def Get_OH_lq(N2O2, R23, O32):
    delta_iter = 1e12
    converg = 0.001
    n_iter = 0
    if N2O2 > 0:
        l_N2O2 = np.log10(N2O2)
        if l_N2O2 < -1.2:
            OH_0 = 8.2
        elif l_N2O2 >= -1.2:
            OH_0 = 8.7
        if (R23 > 0) & (O32 > 0):
            while (delta_iter>converg):
                lq = KK_cal_l_q(O32, OH_0)
                if OH_0 > 8.4:
                    OH_1 = KK_cal_OH_upper(R23, lq)
                else:
                    OH_1 = KK_cal_OH_lower(R23, lq)
                delta_iter = np.abs(OH_1 - OH_0)
                n_iter += 1
                OH_0 = OH_1
                if n_iter > 100:
                    break
        else:
            OH_1 = np.nan
    else:
        OH_1 = np.nan
    return OH_1

def OH_KK04_cal(N2O2, R23, O32):
    '''
    Oxygen Abundances and Ionization Parameter calibrator

    ref: 2004Kobulnicky_ApJ617

    Parameters
    ----------
    N2O2 : float
        NIIl6564/OII3727

    R23 : float
        R23 ratio where R32=R3+R2, R2 = OII3227/Hb and
        R3 = OIII5007/Hb

    O23 : float
        OIIIl5007/OIIl3727 ratio

    Returns
    -------
    OH : float
        Oxygen Abundace: 12 + log(O/H)

    Log q : float
        Ionization Parameter
    '''
    if isinstance(N2O2, (int, float)) and isinstance(R23, (int, float)) \
            and isinstance(O32, (int, float)):
        print('Enter if 1')
        OH = Get_OH_lq(N2O2, R23, O32)
        return OH
    N2O2 = np.array(N2O2)
    R23 = np.array(R23)
    O32 = np.array(O32)
    OH = []
    # print(len(N2O2.ravel()))
    for i in range(len(N2O2.ravel())):
        OH_kk = Get_OH_lq(N2O2.ravel()[i], R23.ravel()[i], O32.ravel()[i])
        OH.append(OH_kk)
    OH = np.array(OH)
    OH = OH.reshape(N2O2.shape)
    return OH
