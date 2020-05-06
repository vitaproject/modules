from CoolProp.CoolProp import PropsSI as SI
from math import*

g = 9.81 # gravity

class non_dimensional:
    
    def __init__(self):
        
        return
    
    def film_coeff(dh, input_pressure, T_ref_i, Nu):
        # added total contact area to find W/(m2K)
        film = (Nu*SI('L', 'T', T_ref_i, 'P', input_pressure, 'helium') / dh) 
        
        return film
    
    def nusselt(MassFlow, T_ref_i, Ag, input_pressure, dh, epsi, deltaz):
            # Nu=f(Re,Pr) ,, Depends on correlation
            # Reynold's number, Re = rho*v_s*Dh/mu
            # v_s = q_r/(rho*A) -> Bulk velocity (away from boundary layer)
        
        v_s = MassFlow/(SI('D', 'T', T_ref_i, 'P', input_pressure, 'helium')*Ag)
        re = SI('D', 'T', T_ref_i, 'P', input_pressure, 'helium') \
                * v_s * dh / \
                    (SI('V', 'T', T_ref_i, 'P', input_pressure, 'helium') /\
                     (SI('D', 'T', T_ref_i, 'P', input_pressure, 'helium')))
          
          # fricc = 0.25 / (log10 ((epsi /(3.71*dh())) + (5.74 / (re)**0.9)))**2
        a = 2/(log(10))
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
        
        h_1 = (fricc * deltaz * (v_s)**2) / (2 * float(g) * dh) #in [m] 
        
        # smooth tube:
        # friction=(1.82*log10(re)-1.64)**(-2.)
        pr = (SI('V', 'T', T_ref_i, 'P', input_pressure, 'helium')/\
              SI('D', 'T', T_ref_i, 'P', input_pressure, 'helium')) \
            * SI('C', 'T', T_ref_i, 'P', input_pressure, 'helium') \
            / SI('L', 'T', T_ref_i, 'P', input_pressure, 'helium')
        
        # Gnielinski (Nu correlation, no Petukhov correction):
        return ( (fricc/8.)*re*pr/(1.07+12.7*sqrt(fricc/8.) \
                * (pr**(2./3.)-1)) ), re, pr, h_1, v_s
        
        # def iterable_film_coeff(dh, P_secc_i, T_ref_i, Nu, iternal_radius):
        # # added total contact area to find W/K
        #     film = Nu*SI('L', 'T', T_ref_i, 'P', P_secc_i, 'helium') \
        #             / dh \
        #     # * (2 * pi * internal_radius * 0.0591)
            
        # return film
    
        # def nusselt_iterable(MassFlow, T_ref_i, Ag, P_secc_i, dh, epsi, deltaz):
        #     # Nu=f(Re,Pr) ,, Depends on correlation
        #     # Reynold's number, Re = rho*v_s*Dh/mu
        #     # v_s = q_r/(rho*A) -> Bulk velocity (away from boundary layer)
        
        #     v_s = MassFlow/(SI('D', 'T', T_ref_i, 'P', P_secc_i, 'helium')*Ag)
        #     re = SI('D', 'T', T_ref_i, 'P', P_secc_i, 'helium') \
        #             * v_s * dh / SI('V', 'T', T_ref_i, 'P', P_secc_i, 'helium')
              
        #       # fricc = 0.25 / (log10 ((epsi /(3.71*dh())) + (5.74 / (re)**0.9)))**2
        #     a = 2/(log(10))
        #     d = log(10)*re/5.2
        #     b = epsi/(dh*3.7)
        #     s = b*d + log(d)
        #     r = s**(s/(s+1))
        #     m = b*d + log(d/r)
        #     p = log(r/m)
        #     DLA = p*(m/(m+1))
        #     DCFA = DLA*(1+((p/2)/(((m+1)**2)+((p/3)*(2*m-1)))))
        #     RHS = a*((log(d/r)) + DCFA)
        #     fricc = (1/RHS)**2
            
        #     h_1 = (fricc * deltaz * (v_s)**2) / (2 * float(g) * dh)
            
        #     moodyf.append(fricc)
            
        #     # smooth tube:
        #     # friction=(1.82*log10(re)-1.64)**(-2.)
        #     pr = SI('V', 'T', T_ref_i, 'P', P_secc_i, 'helium') \
        #         * SI('C', 'T', T_ref_i, 'P', P_secc_i, 'helium') \
        #         / SI('L', 'T', T_ref_i, 'P', P_secc_i, 'helium')
            
        # # Gnielinski (Nu correlation, no Petukhov correction):
        # return ( (fricc/8.)*re*pr/(1.07+12.7*sqrt(fricc/8.) \
        #         * (pr**(2./3.)-1)) ), re, pr, h_1, v_s