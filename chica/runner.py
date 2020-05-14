from CoolProp.CoolProp import PropsSI as SI
from non_dimensional import non_dimensional
from math import *

g = 9.81 # gravity
R = 2077.1
gamma = 1.66

class solver:
    
    def __init__(self):
        
        return
    
    def initial(T_ref, section_0, section_1, MassFlow, input_power, input_pressure, A1, dh, \
                epsi, deltaz, T_metal, Re, Pr, h_f, hf_tot, P_secc, v_secc, Nusselt, htc_0, A2, Ma):
        
        ## 1st loop
        for i in range( len(section_1)-1 ) :
            
            T_ref[i+1] = T_ref[i] + (input_power[i] / (MassFlow*SI('C', 'T', T_ref[i], 'P', \
            input_pressure, 'helium'))) 
            
            # Computing the temperature of interface:
            
            Nu, re, pr, h_1, v_s = non_dimensional.nusselt(MassFlow, T_ref[i+1], A1, \
            input_pressure, dh, epsi, deltaz)
            
            T_metal[i+1] = T_ref[i+1] + input_power[i]/(non_dimensional.film_coeff(dh,\
            input_pressure, T_ref[i+1], Nu)*A2)
            
            Ma.append(v_s / sqrt(gamma * R * T_ref[i+1]))
            Nusselt.append(Nu)
            Re.append(re)
            Pr.append(pr)
            h_f.append(h_1)
            v_secc.append(v_s)
            htc_0.append(non_dimensional.film_coeff(dh, input_pressure, T_ref[i+1], Nu) \
            * A2)
            
            P_secc[i+1] = P_secc[i] + (0.5*SI('D', 'T', T_ref[i+1], 'P', P_secc[i], 'helium') \
            * (v_secc[i]**2 - v_secc[i+1]**2))\
            + (SI('D', 'T', T_ref[i+1], 'P', P_secc[i], 'helium')*float(g))\
            * float(section_1[i+1]-section_1[i])\
            - (SI('D', 'T', T_ref[i+1], 'P', P_secc[i], 'helium')*float(g)*float(h_f[i+1]))
            
            # print([(0.5*SI('D', 'T', T_ref[i+1], 'P', P_secc[i], 'helium') * (v_secc[i]**2 - v_secc[i+1]**2)), \
            # (SI('D', 'T', T_ref[i+1], 'P', P_secc[i], 'helium')*float(g)) * float(section_1[i+1]-section_1[i]), \
            # - (SI('D', 'T', T_ref[i+1], 'P', P_secc[i], 'helium')*float(g)*float(h_f[i+1]))])
            
        hf_tot[0] = h_f[0]
        
        ## 2nd loop
        for x in range(len(h_f)-1) :
            hf_tot[x+1] = float(h_f[x]) + hf_tot[x]
        
        return P_secc, v_secc, hf_tot, T_metal, T_ref
        
    def looper(T_ref, section_0, section_1, MassFlow, input_power, Ag, dh, epsi, \
               deltaz, T_metal, Re, Pr, h_f, hf_tot, P_secc, v_secc, Nusselt, htc_0, Alist, Ma):
        
        P_secc_old = [10E9] # just assigns a value to P_secc_old to initialise
        count = 0
        
        while sqrt((1-(P_secc[len(P_secc)-1]/P_secc_old[len(P_secc_old)-1]))**2) > 10E-8:
        
            P_secc_old = P_secc
            
            ## 1st loop
            for i in range( len(section_1)-1 ):
                
                T_ref[i+1] = T_ref[i] + input_power[i] / (MassFlow*SI('C', 'T', T_ref[i], 'P', P_secc[i], 'helium'))
                
                # Computing the temperature of interface:
                
                Nu, re, pr, h_1, v_s = non_dimensional.nusselt(MassFlow, T_ref[i+1], Ag[i+1], \
                P_secc[i+1], dh, epsi, deltaz[i])
                
                T_metal[i+1] = T_ref[i+1] + input_power[i]/(non_dimensional.film_coeff(dh,\
                P_secc[i+1], T_ref[i+1], Nu)*Alist[i])
            
                Ma.append(v_s / sqrt(gamma * R * T_ref[i+1]))
                Nusselt.append(Nu)
                Re.append(re)
                Pr.append(pr)
                h_f.append(h_1)
                v_secc.append(v_s)
                htc_0.append(non_dimensional.film_coeff(dh, P_secc[i], T_ref[i], Nu) \
                * Alist[i]) # 9 values in htc_0 output list
                
                P_secc[i+1] = P_secc[i] + (0.5*SI('D', 'T', T_ref[i+1], 'P', P_secc[i], 'helium') \
                * (v_secc[i]**2 - v_secc[i+1]**2))\
                + (SI('D', 'T', T_ref[i+1], 'P', P_secc[i], 'helium')*float(g))\
                * float(section_1[i]-section_1[i+1])\
                - (SI('D', 'T', T_ref[i+1], 'P', P_secc[i], 'helium')*float(g)*float(h_f[i+1]))
                    
            hf_tot[0] = h_f[0]
            
            
            ## 2nd loop
            for x in range(len(h_f)-1) :
                hf_tot[x+1] = float(h_f[x]) + hf_tot[x]
            
            count += 1
            print(P_secc[len(P_secc)-1]/P_secc_old[len(P_secc_old)-1])
        
        return P_secc, v_secc, hf_tot, T_metal, T_ref
    


