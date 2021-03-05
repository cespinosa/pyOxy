#!/usr/bin/env python3

import OxygenMLP
import numpy as np

def OH_Ho_cal(OIIHb, OIIIHb, NIIHb, SIIHb):
    oxygenClass = OxygenMLP.OxygenMLP()
    try:
        oxygenClass.ingestLines(OIIHb.ravel(), OIIIHb.ravel(),
                                NIIHb.ravel(), SIIHb.ravel())
        OH, eOH = oxygenClass.predictZ()
        OH = OH.reshape(OIIHb.shape)
    except:
        OH = np.empty(OIIHb.shape)
        OH[:] = np.nan
    return OH
