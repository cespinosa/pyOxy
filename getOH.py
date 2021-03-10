#!/usr/bin/env python3
# ./getOH.py /data/CALIFA_DATA/dataproducts/fe_files/flux_elines.NGC5947.cube.fits.gz /home/espinosa/tmp/testOH
# --indexLines 45 28 27 26 47 0 46 49 50 42 198
# --lineIDs Ha Hb OIII4959 OIII5007 NII6548 OII3727 NII6583 SII6717 SII6731 SIII6312 EWHa
# --indexeLines 249 232 231 230 251 204 250 253 246 254
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from extinction import Av_calculation, make_correction
from calibrators.OHCalN2Marino    import OH_M13_N2_cal
from calibrators.OHCalO3N2Marino  import OH_M13_O3N2_cal
from calibrators.OHCalR23Tremonti import OH_T04_cal
from calibrators.OHCalN2Pettini   import OH_Pet04_N2_lin_cal
from calibrators.OHCalN2bPettini  import OH_Pet04_N2_poly_cal
from calibrators.OHCalO3N2Pettini import OH_Pet04_O3N2_cal
from calibrators.OHCalN2O2Kewley  import OH_Kew02_N2O2_cal
from calibrators.OHCalONSPilyugin import OH_Pil10_ONS_cal
from calibrators.OHCalONPilyugin  import OH_Pil10_ON_cal
from calibrators.OHCalNSPilyugin  import OH_Pil11_NS_cal
from calibrators.OHCalRPilyugin   import OH_Pil16_R_cal
from calibrators.OHCalSPilyugin   import OH_Pil16_S_cal
from calibrators.OHCalRS32Curti   import OH_C20_RS32_cal
from calibrators.OHCalR3Curti     import OH_C20_R3_cal
from calibrators.OHCalO3O2Curti   import OH_C20_O3O2_cal
from calibrators.OHCalS2Curti     import OH_C20_S2_cal
from calibrators.OHCalR2Curti     import OH_C20_R2_cal
from calibrators.OHCalN2Curti     import OH_C20_N2_cal
from calibrators.OHCalR23Curti    import OH_C20_R23_cal
from calibrators.OHCalO3N2Curti   import OH_C20_O3N2_cal
from calibrators.OHCalO3S2Curti   import OH_C20_O3S2_cal
from calibrators.OHCalKK04        import OH_KK04_cal
from calibrators.OHCalHo          import OH_Ho_cal
from calibrators.UCalO32Dors      import U_Dors_O32_cal
from calibrators.UCalS2Dors       import U_Dors_S_cal
from calibrators.UCalO23Morisset  import U_Mor16_O23_fs_cal
from calibrators.UCalO23Morisset  import U_Mor16_O23_ts_cal
from calibrators.UCalS23Morisset  import U_Mor16_S23_cal
from calibrators.NHCalRPilyugin   import NH_Pil16_R_cal
from calibrators.NOCalPilyugin    import NO_Pil16_cal
from calibrators.NeCalSOsterbrock import Ne_Oster_S_cal

wl_lines = {'Ha':6562.68, 'Hb':4861.32, 'OIII4959':4958.91, 'OIII5007':5006.84,
            'NII6548':6548.08, 'OII3727':3727.4, 'NII6583':6583.41,
            'SII6717':6716.39, 'SII6731':6730.74, 'SIII6312':6312.40}

def getAllOH(dicIn):
    dicOut = {}
    dicOut['OH_Mar13_N2'] = OH_M13_N2_cal(dicIn['N2Ha'])
    dicOut['OH_Mar13_O3N2'] = OH_M13_O3N2_cal(dicIn['O3N2'])
    dicOut['OH_T04'] = OH_T04_cal(dicIn['R23'], dicIn['N2Ha'])
    dicOut['OH_Pet04_N2_lin'] = OH_Pet04_N2_lin_cal(dicIn['N2Ha'])
    dicOut['OH_Pet04_N2_poly'] = OH_Pet04_N2_poly_cal(dicIn['N2Ha'])
    dicOut['OH_Pet04_O3N2'] = OH_Pet04_O3N2_cal(dicIn['O3N2'])
    dicOut['OH_Kew02_N2O2'] = OH_Kew02_N2O2_cal(dicIn['N2O2Kew'])
    dicOut['OH_Pil10_ONS']  = OH_Pil10_ONS_cal(dicIn['R3'], dicIn['N2_Hb'],
                                                dicIn['R2'], dicIn['S2'],
                                                dicIn['P'])
    dicOut['OH_Pil10_ON']  = OH_Pil10_ON_cal(dicIn['R3'], dicIn['N2_Hb'],
                                              dicIn['R2'], dicIn['S2'])
    dicOut['OH_Pil11_NS']  = OH_Pil11_NS_cal(dicIn['R3'], dicIn['N2_Hb'],
                                              dicIn['S2'])
    dicOut['OH_Cur20_RS32'] = OH_C20_RS32_cal(dicIn['RS32'])
    dicOut['OH_Cur20_R3'] = OH_C20_R3_cal(dicIn['O3Hb'])
    dicOut['OH_Cur20_O3O2'] = OH_C20_O3O2_cal(dicIn['O3O2'])
    dicOut['OH_Cur20_S2'] = OH_C20_S2_cal(dicIn['S2Ha'])
    dicOut['OH_Cur20_R2'] = OH_C20_R2_cal(dicIn['R2'])
    dicOut['OH_Cur20_N2'] = OH_C20_N2_cal(dicIn['N2Ha'])
    dicOut['OH_Cur20_R23'] = OH_C20_R23_cal(dicIn['R23'])
    dicOut['OH_Cur20_O3N2'] = OH_C20_O3N2_cal(dicIn['O3N2'])
    dicOut['OH_Cur20_O3S2'] = OH_C20_O3S2_cal(dicIn['O3S2'])
    dicOut['OH_KK04'] = OH_KK04_cal(dicIn['N2O2'], dicIn['R23'],
                                     dicIn['O23'])
    dicOut['OH_Pil16_R'] = OH_Pil16_R_cal(dicIn['R3'], dicIn['N2'],
                                          dicIn['R2'])
    dicOut['OH_Pil16_S'] = OH_Pil16_S_cal(dicIn['R3'], dicIn['N2'],
                                          dicIn['S2'])
    dicOut['OH_Ho'] = OH_Ho_cal(dicIn['R2'], dicIn['O3Hb'],
                                dicIn['N2Hb'], dicIn['S2'])
    dicOut['U_Dors_O32']     = U_Dors_O32_cal(dicIn['O32'])
    dicOut['U_Dors_S']     = U_Dors_S_cal(dicIn['S2Ha'])
    dicOut['U_Mor16_O23_fs'] = U_Mor16_O23_fs_cal(dicIn['O23'])
    dicOut['U_Mor16_O23_ts'] = U_Mor16_O23_ts_cal(dicIn['O23'])
    dicOut['U_Mor16_S23'] = U_Mor16_S23_cal(dicIn['S23'])
    dicOut['NH_Pil16_R'] = NH_Pil16_R_cal(dicIn['R3'], dicIn['N2_Hb'],
                                           dicIn['R2'])
    dicOut['NO_Pil16_R'] = dicOut['NH_Pil16_R'] - dicOut['OH_Pil16_R']
    dicOut['NO_Pil16_Ho_R'] = dicOut['NH_Pil16_R'] - dicOut['OH_Ho']
    dicOut['NO_Pil16_N2_R2'] = NO_Pil16_cal(dicIn['N2_Hb'], dicIn['R2'])
    dicOut['Ne_Oster_S'] = Ne_Oster_S_cal(dicIn['SII6717_cor'],
                                           dicIn['SII6731_cor'])
    return dicOut

def main(InputPath, OutputPath, indices, eindices, lineIDs, nMC=1,
         EWHaCut=6, SFHCube=None, AgeCut=0.3548, fyCut=0.04,
         MaskMaps=0):
    CreateDic  = True
    data = fits.getdata(InputPath)
    if 'EWHa' in lineIDs:
        EWHaIDX = indices[lineIDs.index('EWHa')]
        IDX2Remove = lineIDs.index('EWHa')
        lineIDs.pop(IDX2Remove)
        indices.pop(IDX2Remove)
    else:
        EWHaIDX = None
    dicInOri = {label: [data[index], data[eindex]]
                for index, eindex, label in zip(indices, eindices,lineIDs)}
    for i in np.arange(nMC):
        dicIn = {key: np.random.normal(dicInOri[key][0],
                                       dicInOri[key][1])
                 for key in dicInOri.keys()}
        # Make exctintion correction
        dicIn['Av'] = Av_calculation(dicIn['Ha']/dicIn['Hb'])
        for label in lineIDs:
            dicIn['{}_cor'.format(label)] = make_correction(dicIn[label],
                                                            dicIn['Av'],
                                                            wl_lines[label])
        dicIn['O32'] = 1.333 * dicIn['OIII5007_cor'] / dicIn['OII3727_cor']
        dicIn['R3'] = (dicIn['OIII4959_cor'] + dicIn['OIII5007_cor'])/dicIn['Hb_cor']
        dicIn['N2'] = (dicIn['NII6548_cor'] + dicIn['NII6583_cor'])/dicIn['Hb_cor']
        dicIn['N2Ha'] = (dicIn['NII6583_cor']/dicIn['Ha_cor'])
        dicIn['N2Hb'] = (dicIn['NII6583_cor']/dicIn['Hb_cor'])
        dicIn['O3N2'] = (dicIn['OIII5007_cor']/dicIn['Hb_cor'] * dicIn['Ha_cor']/dicIn['NII6583_cor'])
        dicIn['R2'] = dicIn['OII3727_cor']/dicIn['Hb_cor']
        dicIn['S2'] = (dicIn['SII6717_cor'] + dicIn['SII6731_cor'])/dicIn['Hb_cor']
        dicIn['S23'] = (dicIn['SII6717_cor'] + dicIn['SII6731_cor']) / dicIn['SIII6312_cor']
        dicIn['R23']   = dicIn['R2'] + dicIn['R3']
        dicIn['N2O2Kew']  = dicIn['NII6583_cor'] / dicIn['OII3727_cor']
        dicIn['N2_Hb'] = 1.333 * dicIn['NII6583_cor'] / dicIn['Hb_cor']
        dicIn['N2_Ha'] = 1.333 * dicIn['NII6583_cor'] / dicIn['Ha_cor']
        dicIn['P']     = dicIn['R3']/(dicIn['R3']+dicIn['R2'])
        dicIn['RS32']  = ((dicIn['SII6717_cor']+\
                            dicIn['SII6731_cor'])/dicIn['Ha_cor']) + \
                            (dicIn['OIII5007_cor']/dicIn['Hb_cor'])
        dicIn['O3Hb']  = dicIn['OIII5007_cor'] / dicIn['Hb_cor']
        dicIn['O3O2']   = dicIn['OIII5007_cor'] / dicIn['OII3727_cor']
        dicIn['S2Ha']  = (dicIn['SII6717_cor']+\
                            dicIn['SII6731_cor'])/dicIn['Ha_cor']
        dicIn['O3S2']   = dicIn['OIII5007_cor'] / dicIn['S2Ha']
        dicIn['N2O2']  = 1.333 * dicIn['NII6583_cor'] / dicIn['OII3727_cor']
        dicIn['O23']   = dicIn['OII3727_cor'] / dicIn['OIII5007_cor']

        dicTMP = getAllOH(dicIn)
        if CreateDic:
            dicTMP2 = {key: dicTMP[key][np.newaxis,...]
                      for key in dicTMP.keys()}
            CreateDic=False
        else:
            for key in dicTMP2.keys():
                dicTMP2[key] = np.vstack([(dicTMP2[key]),
                                         (dicTMP[key])[np.newaxis,...]])
    # Mask EWHa:
    mask = np.ones_like(dicIn['Ha_cor'], dtype=bool)
    if not EWHaIDX is None:
        # print('EWHa mask')
        EWHaMap = data[EWHaIDX]
        EWHaMap = -EWHaMap
        maskEWHa =  EWHaMap > EWHaCut
        mask = mask * maskEWHa
    if not SFHCube is None:
        # print('SFH mask')
        SFHData = fits.getdata(SFHCube)
        SFHHeader = fits.getheader(SFHCube)
        idxs = []
        ages = []
        for key in SFHHeader.keys():
            if key.startswith('DESC'):
                if SFHHeader[key].startswith('Luminosity'):
                    if ('age' in SFHHeader[key]) and (not 'met' in SFHHeader[key]):
                        idx = key.split('_')[1]
                        age = (SFHHeader[key].split())[-2]
                        idxs.append(int(idx))
                        ages.append(float(age))
        fy = np.zeros_like(dicIn['Ha_cor'])
        for idx, age in zip(idxs, ages):
            if age > AgeCut:
                break
            fy += SFHData[idx]
        maskSFH = fy > fyCut
        mask = mask * maskSFH
    dicOut = {}
    for key in dicTMP2.keys():
        dicOut[key] = np.nanmean(dicTMP2[key], axis=0)
        # mask = (dicOut[key] < 7) | (dicOut[key] > 9.5)
        # dicOut[key][mask] = np.nan
        dicOut['e_'+key] = np.nanstd(dicTMP2[key], axis=0)
        # dicOut['e_'+key][mask] = np.nan
    for key in dicOut.keys():
        if not ('Ne' in key):
            dicOut[key][~mask] = np.nan
    if MaskMaps == 1:
        dicOut['maskEWHa'] = maskEWHa
        dicOut['maskSFH'] = maskSFH
        dicOut['mask'] = mask
    llist = []
    for i, key in enumerate(dicOut.keys()):
        llist.append(dicOut[key])
    CubeOut = np.array(llist)
    # print(dicIn['Ha_cor'])
    hdu = fits.PrimaryHDU(CubeOut)
    hdu.header['COMMENTARY1'] = 'FORMAT DESCXX HEADERS: OH_AUTHOR_INDEX_DESC'
    for i, key in enumerate(dicOut.keys()):
        hdu.header['DESC{}'.format(i+1)] = key
    hdul = fits.HDUList([hdu])
    hdul.writeto('{}.fits'.format(OutputPath), overwrite=True)
    # OH = getAllOH(OII3727, OIII4959, OIII5007, NII6548, NII6583,
    #             SII6717, SII6731, Ha, Hb)

if __name__ == "__main__":
    # Initiate the parser
    desc = 'Get abundance maps from abundance emission maps'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-V", "--version", help="Show program version",
                        action="store_true")
    parser.add_argument("inPath", help="Input file path", type=str)
    parser.add_argument("outPath", help="Output Path", type=str)
    parser.add_argument("--indexLines", help="Emission Line Indices in file",
                        nargs='+', type=int, default=None)
    parser.add_argument("--indexeLines", help="Emission Line Error Indices in file",
                        nargs='+', type=int, default=None)
    parser.add_argument("--lineIDs", help="Emission Line Name",
                        nargs='+', type=str, default=None)
    parser.add_argument("--nMC", help="MC iterations",
                        type=int, default=10)
    parser.add_argument("--EWHaCut", help="EWHa cut value for filter the map if the EWHa index is given",
                        type=float, default=6.0)
    parser.add_argument("--SFHCube", help="Path of SFHCube for stellar populations filter",
                        type=str,default=None)
    parser.add_argument("--AgeCut", help="EWHa cut value for filter the map if the EWHa index is given",
                        type=float, default=0.3548)
    parser.add_argument("--fyCut", help="EWHa cut value for filter the map if the EWHa index is given",
                        type=float, default=0.04)
    parser.add_argument("--MaskMaps", help="flag to get the mask maps used,  1-True, 0-False",
                        type=int, default=0)

    # Read arguments from the command line
    args = parser.parse_args()
    # Print version
    if args.version:
        print("pyOxy version 0.1")
    inPath = args.inPath
    outPath = args.outPath
    indices = args.indexLines
    eindices = args.indexeLines
    lineIDs = args.lineIDs
    nMC = args.nMC
    EWHaCut = args.EWHaCut
    SFHCube = args.SFHCube
    AgeCut = args.AgeCut
    fyCut = args.fyCut
    MaskMaps = args.MaskMaps
    if indices is None:
        print("Indices are required")
        sys.exit()
    if eindices is None:
        print("Indices are required")
        sys.exit()
    if lineIDs is None:
        print("Line IDs are required")
        sys.exit()
    main(InputPath=inPath, OutputPath=outPath, indices=indices, lineIDs=lineIDs,
         eindices=eindices, nMC=nMC, EWHaCut=EWHaCut, SFHCube=SFHCube,
         AgeCut=AgeCut, fyCut=fyCut, MaskMaps=MaskMaps)
