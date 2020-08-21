
# file for non-dimension calculations

from CoolProp.CoolProp import PropsSI as SI
from math import log, sqrt

G = 9.81  # gravity

    
def film_coeff(dh, p_secc_i, T_ref_i, Nu):
    """
    Heat transfer coefficient

    :param float dh: Hydraulic diameter
    :param float p_secc_i: pressure along pipe length
    :param float T_ref_i: reference temperature
    :param float Nu: nusselt number
    :rtype float:
    """
    # added total contact area to find W/(m2K)
    film = (Nu*SI('L', 'T', T_ref_i, 'P', p_secc_i, 'helium') / dh)

    return film
    

def nusselt(MassFlow, T_ref_i, A1, p_secc_i, dh, epsi, deltaz):

    v_s = MassFlow/(SI('D', 'T', T_ref_i, 'P', p_secc_i, 'helium')*A1)
    re = SI('D', 'T', T_ref_i, 'P', p_secc_i, 'helium') \
            * v_s * dh / \
                (SI('V', 'T', T_ref_i, 'P', p_secc_i, 'helium') /\
                 (SI('D', 'T', T_ref_i, 'P', p_secc_i, 'helium')))

    a = 2.0/(log(10))
    d = log(10)*re/5.2
    b = epsi/(dh*3.7)
    s = b*d + log(d)
    r = s**(s/(s+1))
    m = b*d + log(d/r)
    p = log(r/m)
    DLA = p*(m/(m+1))
    DCFA = DLA*(1+((p/2)/(((m+1)**2)+((p/3)*(2*m-1)))))
    RHS = a*((log(d/r)) + DCFA)
    fricc = (1/RHS)**2

    h_1 = (fricc * deltaz * (v_s)**2) / (2 * float(G) * dh) #in [m]

    # smooth tube:
    # friction=(1.82*log10(re)-1.64)**(-2.)

    pr = (SI('V', 'T', T_ref_i, 'P', p_secc_i, 'helium')/\
          SI('D', 'T', T_ref_i, 'P', p_secc_i, 'helium')) \
        * SI('C', 'T', T_ref_i, 'P', p_secc_i, 'helium') \
        / SI('L', 'T', T_ref_i, 'P', p_secc_i, 'helium')

    # Gnielinski (Nu correlation, no Petukhov correction):
    return ( (fricc/8.)*re*pr/(1.07+12.7*sqrt(fricc/8.) \
            * (pr**(2./3.)-1)) ), re, pr, h_1, v_s