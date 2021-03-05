#!/usr/bin/env python3

import numpy as np
from OHCalO3N2Marino import OH_M13_O3N2_cal

def OH_M08_cal(R23, O3N2, N2O2):
    '''
    Oxygen Abundance calibrator based on Maiolino's work (2008)
    If log(N2O2) > 1.2, then OH is calculated by the Marino's calibrator
    The OH default value is nan

    ref: 2008Maiolino_AA488

    Parameters
    ----------
    R23 : float
        R23 ratio where R32=R3+R2, R2 = OII3227/Hb and
        R3 = OIIIl5007/Hb

    O3N2 : float
        OIIIl5007/Hb * Ha/NIIl6583

    N2O2 : float
        NIIl6564/OII3727

    Returns
    -------
    OH : float
        Oxygen Abundace: 12 + log(O/H)
    '''
    if isinstance(R23, (int, float)) and isinstance(O3N2, (int, float)) \
       and isinstance(N2O2, (int, float)):
        l_R23 = np.log10(R23)
        l_O3N2 = np.log10(O3N2)
