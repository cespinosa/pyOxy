#!/usr/bin/env python3

import numpy as np

def OH_Pil16_S_cal(R3, N2, S2):
    '''
    Oxygen Abundance calibrator based on Pilyugin's work (2016)
    N2 = NII(l6548 + l6584)/Hb
    R3 = OIII(l4959 + l5007)/Hb
    S2 = SII(l6717 + l6731)/Hb

    ref: 2016Pilyugin_MNRAS457

    Parameters
    ----------
    R3 : float
        R3 ratio
    N2 : float
        N2 ratio
    S2 : float
        S2 ratio
    Returns
    -------
    OH : float
        Oxygen Abundace: 12 + log(O/H)
    '''
    if isinstance(R3, (int, float)) and isinstance(N2, (int, float)) \
            and isinstance(S2, (int, float)):
        if np.log10(N2) >= -0.6:
            OH = 8.424  + 0.030 * np.log10(R3/S2) + 0.751 * np.log10(N2) + \
                (-0.349 + 0.182 * np.log10(R3/S2) + \
                 0.508 * np.log10(N2)) * np.log10(S2)
        elif np.log10(N2) < -0.6:
            OH = 8.072 + 0.789 * np.log10(R3 / S2) + 0.726 * np.log10(N2) + \
                (1.069 - 0.170 * np.log10(R3 / S2) + \
                 0.022 * np.log10(N2)) * np.log10(S2)
        else:
            OH = np.nan
        return OH
    if isinstance(N2, (int, float)): N2 = np.full(R3.shape, N2, dtype=R3.dtype)
    if isinstance(S2, (int, float)): S2 = np.full(R3.shape, S2, dtype=R3.dtype)
    if isinstance(R3, (int, float)): R3 = np.full(S2.shape, R3, dtype=S2.dtype)

    mask_N2_g = np.log10(N2) >= -0.6
    mask_N2_l = np.log10(N2) < -0.6

    OH = np.empty(N2.shape)
    OH[:] = np.nan
    OH[mask_N2_g] = 8.424  + 0.030 * np.log10(R3[mask_N2_g]/S2[mask_N2_g]) + \
                    0.751 * np.log10(N2[mask_N2_g]) + \
                    (-0.349 + 0.182 * np.log10(R3[mask_N2_g]/S2[mask_N2_g]) + \
                    0.508 * np.log10(N2[mask_N2_g])) * np.log10(S2[mask_N2_g])
    OH[mask_N2_l] = 8.072 + 0.789 * np.log10(R3[mask_N2_l]/S2[mask_N2_l]) + \
                    0.726 * np.log10(N2[mask_N2_l]) + \
                    (1.069 - 0.170 * np.log10(R3[mask_N2_l]/S2[mask_N2_l]) + \
                    0.022 * np.log10(N2[mask_N2_l])) * np.log10(S2[mask_N2_l])
    return OH
