#!/usr/bin/env python3
import numpy as np

def OH_Pet04_N2_poly_cal(N2):
    '''
    Oxygen Abundance calibrator base on the Pettini & Pagel's work (2004)

    ref: 2004Pettini_MNRAS348

    Parameters
    ----------
    N2 : float
       NII6584)/Ha

    Returns
    -------
    OH : float
        Oxygen Abundance: 12 + log(O/H)
    '''
    OH_f = lambda l_N2: 9.37 + 2.03 * l_N2 + 1.26 * l_N2 * l_N2 + \
                    0.32 * l_N2 * l_N2 * l_N2
    if isinstance(N2, (int, float)):
        if N2 > 0:
            l_N2 = np.log10(N2)
            if (l_N2 > -2.5) & (l_N2 < -0.3):
                OH = OH_f(l_N2)
            else:
                OH = np.nan
        else:
            OH = np.nan
        return OH
    N2 = np.array(N2)
    l_N2 = np.log10(N2)
    mask1 = l_N2 > -2.5
    mask2 = l_N2 < -0.3
    mask = mask1 & mask2
    OH = np.empty(N2.shape)
    OH[:] = np.nan
    OH[mask] = OH_f(l_N2[mask])
    return OH
