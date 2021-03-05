#!/usr/bin/env python3
import numpy as np

def OH_Pet04_N2_lin_cal(N2):
    '''
    Oxygen Abundance calibrator base on the Pettini & Pagel's work (2004)

    ref: 2004Pettini_MNRAS348

    Parameters
    ----------
    N2 : float
       NII6584/Ha

    Returns
    -------
    OH : float
        Oxygen Abundance: 12 + log(O/H)
    '''
    OH_f = lambda l_N2: 8.9 + 0.57 * l_N2
    if isinstance(N2, (int, float)):
        if N2 > 0:
            l_N2 = np.log10(N2)
            if -2.5 < l_N2 < -0.3:
                OH = OH_f(l_N2)
            else:
                OH = np.nan
        else:
            OH = np.nan
        return OH
    N2 = np.array(N2)
    l_N2 = np.log10(N2)
    mask_N2_l = l_N2 > -2.5
    mask_N2_u = l_N2 < -0.3
    mask = mask_N2_u & mask_N2_l
    OH = np.empty(N2.shape)
    OH[:] = np.nan
    OH[mask] = OH_f(l_N2[mask])
    return OH
