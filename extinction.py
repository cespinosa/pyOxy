import numpy as np

def A_l(R_v, lw):
    #  From Cardelli,1989
    #  F_cor = F * 10 ***(0.4*Av*A_l(R_v,l))
    lw = lw / 10000
    x = 1 / lw
    if x > 1.1:
        y = x - 1.82
        a_x = 1.0 + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 \
            + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7
        b_x = 1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 \
            - 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7
    else:
        a_x = 0.574 * x ** 1.61
        b_x = -0.527 * x ** 1.61
    A_l_ = a_x + b_x/R_v
    return A_l_

def make_correction(flux, Av, wl, Rv=3.1):
    Al = A_l(Rv, wl)
    flux_cor = flux * 10 ** (0.4 * Av * Al)
    return flux_cor
    
def Av_calculation(HaHb, Rv=3.1):
    a1 = A_l(Rv, 4861)
    a2 = A_l(Rv, 6562)
    if isinstance(HaHb, (int, float)):
        if HaHb >= 2.86:
            Av = (2.5 * np.log10(HaHb/2.86))/(a1 - a2)
        elif HaHb < 2.86:
            Av = 0 
        return Av
    Av = np.empty(HaHb.shape)
    Av[:] = np.nan
    mask_HaHb_g = HaHb >= 2.86
    mask_HaHb_l = HaHb < 2.86
    Av[mask_HaHb_g] = (2.5 * np.log10(HaHb[mask_HaHb_g]/2.86))/(a1 - a2)
    Av[mask_HaHb_l] = 0
    return Av
