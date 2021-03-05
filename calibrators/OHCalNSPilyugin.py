#!/usr/bin/env python3

import numpy as np

def OH_Pil11_NS_cal(R3, N2, S2):
    '''
    Oxygen Abundance calibrator based on Pilyugin's work (2011)
    N2 = NII(l6548 + l6584)/Hb
    S2 = SII(l6717 + l6731)/Hb
    R3 = OIII(l4959 + l5007)/Hb

    ref: 2011Pilyugin_MNRAS412

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
            if np.log10(N2) > -0.1:
                OH = 8.454  - 0.216 * np.log10(R3) - 0.362 * np.log10(S2) + \
                    (-0.101 * np.log10(N2/S2))
            elif np.log10(N2) <= -0.1:
                if (np.log10(N2/S2) > -0.25):
                    OH = 8.456 + 0.082 * np.log10(R3) + 0.391 * np.log10(N2) + \
                        0.290 * np.log10(N2/S2)
                elif (np.log10(N2/S2) <= -0.25):
                    OH = 7.881 + 0.929 * np.log10(R3) + 0.650 * np.log10(N2) + \
                        0.025 * np.log10(N2/S2)
                else:
                    OH = np.nan
            else:
                OH = np.nan
            return OH
    if isinstance(R3, (int, float)): R3 = np.full(R2.shape, R3, dtype=R3.dtype)
    if isinstance(N2, (int, float)): N2 = np.full(R3.shape, N2, dtype=R3.dtype)
    if isinstance(S2, (int, float)): S2 = np.full(S2.shape, S2, dtype=R3.dtype)
    mask1 = np.log10(N2) > -0.1
    mask2 = np.log10(N2) <= -0.1
    mask3 = np.log10(N2/S2) > -0.25
    mask4 = np.log10(N2/S2) <= -0.25
    OH = np.empty(R3.shape)
    OH[:] = np.nan
    OH[mask1] = 8.454  - 0.216 * np.log10(R3[mask1]) - \
        0.362 * np.log10(S2[mask1]) + \
        (-0.101 * np.log10(N2[mask1]/S2[mask1]))
    OH[mask2 & mask3] = 8.456 + 0.082 * np.log10(R3[mask2 & mask3]) + \
        0.391 * np.log10(N2[mask2 & mask3]) + \
        0.290 * np.log10(N2[mask2 & mask3]/S2[mask2 & mask3])
    OH[mask2 & mask4] = 7.881 + 0.929 * np.log10(R3[mask2 & mask4]) + \
        0.650 * np.log10(N2[mask2 & mask4]) + \
        0.025 * np.log10(N2[mask2 & mask4]/S2[mask2 & mask4])
    return OH
