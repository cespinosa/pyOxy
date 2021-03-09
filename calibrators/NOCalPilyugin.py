#!/usr/bin/env python3

import numpy as np

def NO_Pil16_cal(N2, R2):
    '''
    Nitrogen Oxygen ratio calibrator based on Pilyugin's work (2016)
    N2 = NII(l6548 + l6584)/Hb
    R2 = OII(l3727 + l3729)/Hb

    ref: 2016Pilyugin_MNRAS457

    Parameters
    ----------
    N2 : float
        N2 ratio
    R2 : float
        R2 ratio

    Returns
    -------
    NO : float
       Nitrogen Abundace: log(N/O)
    '''
    if isinstance(N2, (int, float)) and isinstance(R2, (int, float)):
        if -1.6 < np.log10(N2/R2) < 0.6:
            if N2 > 0:
                NO = -0.657 - 0.201 * np.log10(N2) + \
                    (0.742-0.075*np.log10(N2))*np.log10(N2/R2)
            else:
                NO = np.nan
        else:
            NO = np.nan
        return NO
    if isinstance(N2, (int, float)): N2 = np.full(R2.shape, N2, dtype=R2.dtype)
    if isinstance(R2, (int, float)): R2 = np.full(N2.shape, R2, dtype=N2.dtype)
    mask_g = np.log10(N2/R2) >= -1.6
    mask_l = np.log10(N2/R2) <= 0.6
    mask_N2 = N2 > 0
    NO = np.empty(N2.shape)
    NO[:] = np.nan
    NO[mask_g & mask_l & mask_N2] = -0.657 - \
        0.201 * np.log10(N2[mask_g & mask_l & mask_N2]) + \
        (0.742-0.075*np.log10(N2[mask_g & mask_l & mask_N2]))*\
        np.log10(N2[mask_g & mask_l & mask_N2]/R2[mask_g & mask_l & mask_N2])
    return NO
