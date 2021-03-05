#!/usr/bin/env python3
import numpy as np

def OH_T04_cal(R23, NIIHa):
    '''
    Oxygen Abundance calibrator for the upper branch based on R23 (R2 + R3) ratio
    R2 = OII3227+29/Hb
    R3 = OIII5007/Hb

    ref: 2004Tremonti_ApJ613

    Parameters
    ----------
    R23 : float
        R23 = (R2 + R3)

    NIIHa : float
        NII6584/Ha ratio

    Returns
    -------
    OH : float
        Oxygen Abundace: 12 + log(O/H)
    '''
    OH_f = lambda l_R23: 9.185 - 0.313 * l_R23 - 0.264 * (l_R23**2) - \
        0.321 * (l_R23**3)
    if isinstance(R23, (int, float)) and isinstance(NIIHa, (int, float)):

        if R23 > 0:
            l_R23 = np.log10(R23)
            if np.log10(NIIHa) >= -1.1:
                OH = OH_f(l_R23)
            else:
                OH = np.nan
        else:
            OH = np.nan
        return OH
    if isinstance(R23, (int, float)):
        R23 = np.full(R23.shape, R23, dtype=R23.dtype)
    if isinstance(NIIHa, (int, float)):
        NIIHa = np.full(NIIHa.shape, NIIHa, dtype=NIIHa.dtype)
    R23 = np.array(R23)
    mask_R23 = R23 > 0
    mask_NII = np.log10(NIIHa) >= -1.1
    OH = np.empty(R23.shape)
    OH[:] = np.nan
    OH[mask_R23 & mask_NII] = OH_f(np.log10(R23[mask_R23 & mask_NII]))
    return OH
