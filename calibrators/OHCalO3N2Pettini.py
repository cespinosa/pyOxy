#!/usr/bin/env python3
import numpy as np

def OH_Pet04_O3N2_cal(O3N2):
    '''
    Oxygen Abundance calibrator base on the Pettini & Pagel's work (2004)

    ref: 2004Pettini_MNRAS348

    Parameters
    ----------
    O3N2 : float
        OIIIl5007/Hb * Ha/NIIl6583

    Returns
    -------
    OH : float
        Oxygen Abundance: 12 + log(O/H)
    '''
    OH_f = lambda l_O3N2: 8.73 - 0.32 * l_O3N2
    if isinstance(O3N2, (int, float)):
        if O3N2 > 0:
            l_O3N2 = np.log10(O3N2)
            if -1 < l_O3N2 < 1.9:
                OH = OH_f(l_O3N2)
            else:
                OH = np.nan
        else:
            OH = np.nan
        return OH
    O3N2 = np.array(O3N2)
    l_O3N2 = np.log10(O3N2)
    mask_O3N2_l = l_O3N2 > -1
    mask_O3N2_u = l_O3N2 < 1.9
    mask = mask_O3N2_l & mask_O3N2_u
    OH = np.empty(O3N2.shape)
    OH[:] = np.nan
    OH[mask] = OH_f(l_O3N2[mask])
    return OH
