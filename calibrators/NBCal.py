#!/usr/bin/env python3

import numpy as np
import pandas as pd
from NebulaBayes import NB_Model

def NB_cal(OII, Hb, OIII5007, OI, Ha, NII6583, SII16, SII31,
           eOII, eHb, eOIII5007, eOI, eHa, eNII6583, eSII16, eSII31):
    varNames = {"12 + log O/H":"OH_NB", "log P/k":"l_Pk_NB",
                "log U":"U_NB"}
    varTypes = ["Estimate", "CI68_low", "CI68_high"]
    linelist_NB = ["OII3726_29", "Hbeta", "OIII5007", "OI6300", "Halpha",
                   "NII6583", "SII6716", "SII6731"]
    obs_flux = [OII, Hb, OIII5007, OI, Ha, NII6583, SII16, SII31]
    obs_errs = [eOII, eHb, eOIII5007, eOI, eHa, eNII6583, eSII16, eSII31]
    dicOut = {}
    if np.isfinite(obs_flux).all() and np.isfinite(obs_errs).all():
        NB_Model_HII = NB_Model("HII", line_list=linelist_NB,
                                interpd_grid_shape=(100, 100, 100))
        try:
            Result_HII = NB_Model_HII(obs_flux, obs_errs, linelist_NB)
            Estimate_table = Result_HII.Posterior.DF_estimates
            for varNBName, varOutName in varNames.items():
                dicOut[varOutName] = Estimate_table.loc[varNBName, varTypes[0]]
                dicOut['e1'+varOutName] = dicOut[varOutName] - \
                    Estimate_table.loc[varNBName,varTypes[1]]
                dicOut['e2'+varOutName] = Estimate_table.loc[varNBName,
                                                             varTypes[2]] - \
                                                             dicOut[varOutName]
        except:
            for varNBName, varOutName in varNames.items():
                dicOut[varOutName] = np.nan
                dicOut['e1'+varOutName] = np.nan
                dicOut['e2'+varOutName] = np.nan
    else:
        for varNBName, varOutName in varNames.items():
            dicOut[varOutName] = np.nan
            dicOut['e1'+varOutName] = np.nan
            dicOut['e2'+varOutName] = np.nan
    return pd.Series(dicOut)
