#!/usr/bin/env python3
import numpy as np

def OH_Pil10_ONS_cal(R3, N2, R2, S2, P):
    '''
    Oxygen Abundance calibrator based on Pilyugin's work (2010)
    N2 = NII(l6548 + l6584)/Hb
    R2 = OII(l3727 + l3729)/Hb
    S2 = SII(l6717 + l6731)/Hb
    R3 = OIII(l4959 + l5007)/Hb
    P  = R3/(R3 + R2)

    ref: 2010Pilyugin_ApJ720

    Parameters
    ----------
    R3 : float
        R3 ratio
    N2 : float
        N2 ratio
    R2 : float
        R2 ratio
    S2 : float
        S2 ratio
    P  : float
        P ratio

    Returns
    -------
    OH : float
        Oxygen Abundace: 12 + log(O/H)
    '''
    if isinstance(R3, (int, float)) and isinstance(N2, (int, float)) \
       and isinstance(R2, (int, float)) and isinstance(S2, (int, float)) \
       and isinstance(P, (int, float)):
            if np.log10(N2) > -0.1:
                OH = 8.277 + 0.657 * P - 0.399 * np.log10(R3) + \
                    (-1)*0.061 * np.log10(N2 / R2) + 0.005 * np.log10(S2 / R2)
            elif np.log10(N2) <= -0.1:
                if (np.log10(N2/S2) > -0.25):
                    OH = 8.816 - 0.733 * P + 0.454 * np.log10(R3) + \
                        0.710 * np.log10(N2 / R2) - 0.337 * np.log10(S2 / R2)
                elif (np.log10(N2/S2) < -0.25):
                    OH = 8.774 - 1.855 * P + 1.517 * np.log10(R3) + \
                        0.304 * np.log10(N2 / R2) + 0.328 * np.log10(S2 / R2)
                else:
                    OH = np.nan
            else:
                OH = np.nan
            return OH
    if isinstance(R3, (int, float)): R3 = np.full(R2.shape, R3, dtype=R3.dtype)
    if isinstance(N2, (int, float)): N2 = np.full(R3.shape, N2, dtype=R3.dtype)
    if isinstance(R2, (int, float)): R2 = np.full(R3.shape, R2, dtype=R3.dtype)
    if isinstance(S2, (int, float)): S2 = np.full(S2.shape, S2, dtype=R3.dtype)
    if isinstance(P,  (int, float)): P  = np.full(P.shape,  P,  dtype=R3.dtype)
    mask1 = np.log10(N2) > -0.1
    mask2 = np.log10(N2) <= -0.1
    mask3 = np.log10(N2/S2) > -0.25
    mask4 = np.log10(N2/S2) <= -0.25
    OH = np.empty(R3.shape)
    OH[:] = np.nan
    OH[mask1] = 8.277 + 0.657 * P[mask1] - 0.399 * np.log10(R3[mask1]) + \
                    (-1)*0.061 * np.log10(N2[mask1]/R2[mask1]) + \
                    0.005 * np.log10(S2[mask1]/R2[mask1])
    OH[mask2 & mask3] = 8.816 - 0.733 * P[mask2 & mask3] + \
        0.454 * np.log10(R3[mask2 & mask3]) + \
        0.710 * np.log10(N2[mask2 & mask3]/R2[mask2 & mask3]) - \
        0.337 * np.log10(S2[mask2 & mask3]/R2[mask2 & mask3])
    OH[mask2 & mask4] = 8.774 - 1.855 * P[mask2 & mask4] + \
        1.517 * np.log10(R3[mask2 & mask4]) + \
        0.304 * np.log10(N2[mask2 & mask4]/R2[mask2 & mask4]) + \
        0.328 * np.log10(S2[mask2 & mask4]/R2[mask2 & mask4])
    return OH
