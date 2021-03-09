#!/usr/bin/env python3

import numpy as np

def NH_Pil16_R_cal(R3, N2, R2):
    '''
    Nitrogen Abundance calibrator based on Pilyugin's work (2016)
    N2 = NII(l6548 + l6584)/Hb
    R2 = OII(l3727 + l3729)/Hb
    R3 = OIII(l4959 + l5007)/Hb

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
    NH : float
       Nitrogen Abundace: 12 + log(N/H)
    '''
    if isinstance(N2, (int, float)) and isinstance(R3, (int,float)) \
            and isinstance(R2, (int, float)):
        if np.log10(N2) >= -0.6:
            NH = 7.939  + 0.135 * np.log10(R3/R2) + 1.217 * np.log10(N2) \
                + (-0.765 + 0.166 * np.log10(R3/R2) + 0.449 * np.log10(N2)) \
                * np.log10(R2)
        elif np.log10(N2) < -0.6:
            NH = 7.476  + 0.879 * np.log10(R3/R2) + 1.451 * np.log10(N2) \
                + (-0.011 - 0.327 * np.log10(R3/R2) - 0.064 * np.log10(N2)) \
                * np.log10(R2)
        else:
            NH = np.nan
        return NH

    if isinstance(N2, (int, float)): N2 = np.full(R3.shape, N2, dtype=R3.dtype)
    if isinstance(R2, (int, float)): R2 = np.full(R3.shape, R2, dtype=R3.dtype)
    if isinstance(R3, (int, float)): R3 = np.full(R2.shape, R3, dtype=R2.dtype)

    mask_N2_g = np.log10(N2) >= -0.6
    mask_N2_l = np.log10(N2) < -0.6

    NH = np.empty(N2.shape)
    NH[:] = np.nan
    NH[mask_N2_g] = 7.939  + 0.135 * np.log10(R3[mask_N2_g]/R2[mask_N2_g]) + \
                    1.217 * np.log10(N2[mask_N2_g]) + \
                    (-0.765 + 0.166 * np.log10(R3[mask_N2_g]/R2[mask_N2_g]) + \
                    0.449 * np.log10(N2[mask_N2_g])) * np.log10(R2[mask_N2_g])
    NH[mask_N2_l] = 7.476  + 0.879 * np.log10(R3[mask_N2_l]/R2[mask_N2_l]) + \
                    1.451 * np.log10(N2[mask_N2_l]) + \
                    (-0.011 - 0.327 * np.log10(R3[mask_N2_l]/R2[mask_N2_l]) - \
                    0.064 * np.log10(N2[mask_N2_l])) * np.log10(R2[mask_N2_l])
    return NH
