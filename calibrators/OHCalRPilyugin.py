#!/usr/bin/env python3

import numpy as np


def OH_Pil16_R_cal(R3, N2, R2):
    '''
    Oxygen Abundance calibrator based on Pilyugin's work (2016)
    N2 = NII(l6548 + l6584)/Hb
    R3 = OIII(l4959 + l5007)/Hb
    R2 = OII(l3727 + l3729)/Hb

    ref: 2016Pilyugin_MNRAS457

    Parameters
    ----------
    R3 : float
        R3 ratio
    N2 : float
        N2 ratio
    R2 : float
        R2 ratio
    Returns
    -------
    OH : float
        Oxygen Abundace: 12 + log(O/H)
    '''
    if isinstance(N2, (int, float)) and isinstance(R3, (int,float)) \
            and isinstance(R2, (int, float)):
        if np.log10(N2) >= -0.6:
            OH = 8.589  + 0.022 * np.log10(R3/R2) + 0.399 * np.log10(N2) + \
                (-0.137 + 0.164 * np.log10(R3/R2) + \
                0.589 * np.log10(N2)) * np.log10(R2)
        elif np.log10(N2) < -0.6:
            OH = 7.932 + 0.944 * np.log10(R3 / R2) + 0.695 * np.log10(N2) + \
                (0.970 - 0.291 * np.log10(R3 / R2) - \
                0.019 * np.log10(N2)) * np.log10(R2)
        else:
            OH = np.nan
        return OH
    if isinstance(N2, (int, float)): N2 = np.full(R3.shape, N2, dtype=R3.dtype)
    if isinstance(R2, (int, float)): R2 = np.full(R3.shape, R2, dtype=R3.dtype)
    if isinstance(R3, (int, float)): R3 = np.full(R2.shape, R3, dtype=R2.dtype)

    mask_N2_g = np.log10(N2) >= -0.6
    mask_N2_l = np.log10(N2) < -0.6

    OH = np.empty(N2.shape)
    OH[:] = np.nan
    OH[mask_N2_g] = 8.589  + 0.022 * np.log10(R3[mask_N2_g]/R2[mask_N2_g]) + \
                    0.399 * np.log10(N2[mask_N2_g]) + \
                    (-0.137 + 0.164 * np.log10(R3[mask_N2_g]/R2[mask_N2_g]) + \
                    0.589 * np.log10(N2[mask_N2_g])) * np.log10(R2[mask_N2_g])
    OH[mask_N2_l] = 7.932 + 0.944 * np.log10(R3[mask_N2_l]/R2[mask_N2_l]) + \
                    0.695 * np.log10(N2[mask_N2_l]) + \
                    (0.970 - 0.291 * np.log10(R3[mask_N2_l]/R2[mask_N2_l]) - \
                    0.019 * np.log10(N2[mask_N2_l])) * np.log10(R2[mask_N2_l])
    return OH
