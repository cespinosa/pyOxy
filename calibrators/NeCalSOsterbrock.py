#!/usr/bin/env python3

import numpy as np

def Ne_Oster_S_cal(SII15, SII30):
    '''
    Electronic density calculation based on SII16SII30 ratio

    SII15SII30 = SIIl6716/SIIl6730 ratio

    ref: 1989Osterbrock_book

    Parameters
    ----------
    SII15 : float
       SIIl6716 flux
    SII30 : float
       SIIl6730 flux

    Returns
    -------
    Ne : float
        Electronic density
    '''
    b = 1.49
    c = 3.77
    d = 12.8
    T=1
    if isinstance(SII15, (int, float)) and isinstance(SII30, (int, float)):
        SII16SII30 = SII15 / SII30
        x = (SII16SII30/b-1)/(c-SII16SII30*d/b)
        Ne = np.nan
        if x > 0:
            Ne = np.log10(x / (10 ** (-4) * T ** (-0.5)))
        return Ne
    if isinstance(SII15, (int, float)): SII15 = np.full(SII30.shape, SII15,
            dtype=SII30.dtype)
    if isinstance(SII30, (int, float)): SII30 = np.full(SII15.shape, SII30,
            dtype=SII15.dtype)
    SII16SII30 = SII15 / SII30
    x = (SII16SII30/b-1)/(c-SII16SII30*d/b)
    mask_x = x > 0
    Ne = np.empty(SII30.shape)
    Ne[:] = np.nan
    Ne[mask_x] = np.log10(x[mask_x] / (10 ** (-4) * T ** (-0.5)))
    return Ne
