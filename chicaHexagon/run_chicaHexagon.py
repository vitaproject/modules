
from chicaHexagon.geometry_setup import initial_domain_space, \
                                               refined_domain_space, \
                                               first_row, \
                                               next_rows
                                   
from chicaHexagon.flow_properties import mach_number_sympy, \
                                                mass_flow_sympy, \
                                                check_cases, \
                                                Reynolds_sympy, \
                                                Euler_sympy, \
                                                pdrop_sympy, \
                                                coolant_temperature_sympy, \
                                                metal_temperature_sympy, \
                                                mstar_sympy, \
                                                taus_sympy, \
                                                massflow_nextrows

from chicaHexagon.misc_tools import write_CAD_input_files, \
                                    plot_temperature, \
                                    input_data

from numpy import linspace, sqrt
from CoolProp.CoolProp import PropsSI as SI
import time

# ------------------------------- centre point lists ------------------------------- #

start_time = time.time()

n_plates, width, inner_radius, outer_radius, x_displacement_from_origin, \
           y_displacement_from_origin, gap, Rstrike, channel_width, massFlow, \
               Ti, Pi, D, t, HFf, sbar, heatflux = input_data("inputs.asc", "q.asc")

height, theta, point_x1, point_y1, line_angle, no_points_x,  \
    cp, A = initial_domain_space(width, n_plates, gap, inner_radius, outer_radius)

hcp_final, tcp_final = refined_domain_space(line_angle, cp, point_x1, \
point_y1, width, inner_radius, outer_radius, x_displacement_from_origin, \
y_displacement_from_origin, Rstrike)
    
#--------------------------------------- flow properties -------------------------------------#

system_level_results = []
runs = []
massflow_split = 0.9 # massflow split, weighting to outboard
Mach_number = linspace(0.15, 0.3, 4)
No_hexagons_1st_row = linspace(200, 700, 6)

for Ma_in in Mach_number:
    for n_in in No_hexagons_1st_row:
        
        Ma, u_in = mach_number_sympy(Ti, Ma = Ma_in)
        massflow_in, u, Ajet, density_in = mass_flow_sympy(Ti, Pi, u = u_in)
        mstar_in, massflow_in, cp_in = mstar_sympy(Ti, Pi, A, massflow = massflow_in)
        
        massflow_total = n_in * massflow_in # number of hexagons
        
        channels, hcp_inboard, hcp_outboard = first_row(hcp_final, n_in)
        
        hcp = [hcp_inboard, hcp_outboard]
        
        for j, direction in enumerate(channels):
            
            Tinput = Ti
            Pinput = Pi
            u = u_in
            density = density_in
            massflow = massflow_in
            mstar = mstar_in
            n = n_in
            specific_heat_capacity = cp_in
            mu = SI("V", "T", Tinput, "P", Pinput, "helium")
            
            for i, channel in enumerate(direction):
                
                for coordinates in channel:
                    Ri = sqrt(coordinates[0]**2 + coordinates[1]**2)
                    HFi = HFf(Ri) * A
                    Re, u = Reynolds_sympy(density, mu, D, u = u)
                    Eu, Re = Euler_sympy(Re = Re)
                    dp, Eu, u = pdrop_sympy(density, Eu = Eu, u = u)
                    Tnew, massflow = coolant_temperature_sympy(HFi, specific_heat_capacity, Tinput, massflow = massflow)
                    mstar, massflow, specific_heat_capacity = mstar_sympy(Tinput, Pinput, A, massflow = massflow)
                    taus, mstar = taus_sympy(mstar = mstar)
                    Tmetal, Tnew = metal_temperature_sympy(Tinput, taus, Tnew = Tnew)
                    coordinates.append([Tnew, dp, Tmetal, massflow, mstar, \
                                        Re, density, u, Ma])

                Tinput = sum([temp[-1][0] for temp in direction[i]]) / len(direction[i])
                Pinput += -dp
                
                mu = SI("V", "T", Tinput, "P", Pinput, "helium")
                density = SI("D", "T", Tinput, "P", Pinput, "helium")
                specific_heat_capacity = SI("C", "T", Tinput, "P", Pinput, "helium")
                
                n = int(len(direction[i])/2.0)
                Re, u = Reynolds_sympy(density, mu, D, u = u)
                check_cases(Re, mstar, Ma)
            
            massflow_actual = [massflow_total * (1-massflow_split), massflow_total * massflow_split]
            ni = 0
            Tpeak = max([i[-1][2] for i in channels[0][0]]) # Kelvin
            
            while ni <= len(hcp[j]):
                
                direction.append([])
                massflow, ni = massflow_nextrows(hcp[j][n+1], massflow, specific_heat_capacity, A, Tpeak, Tinput, HFf, massflow_actual[j])
                
                if massflow == 0.0:
                    flag = 1.0
                    ni = len(hcp[j])
                    massflow = (massflow_actual[j]) / (ni - n)
                elif ni <= len(hcp[j]):
                    flag = 0.0
                    ni = int((massflow_actual[j]) / massflow)
                    ni = int(n + ni)
                else:
                    ni = len(hcp[j])
                    massflow = massflow_actual[j] / (ni - n)
                
                direction[len(direction)-1] = next_rows(hcp[j], n, ni)
                
                for coordinates in direction[len(direction)-1]:
                    Ri = sqrt(coordinates[0]**2 + coordinates[1]**2)
                    HFi = HFf(Ri) * A
                    Tnew, massflow = coolant_temperature_sympy(HFi, specific_heat_capacity, Tinput, massflow = massflow)
                    mstar, massflow, specific_heat_capacity = mstar_sympy(Tinput, Pinput, A, massflow = massflow)
                    taus, mstar = taus_sympy(mstar = mstar)
                    Tmetal, Tnew = metal_temperature_sympy(Tinput, taus, Tnew = Tnew)
                    massflow, u, Ajet, density = mass_flow_sympy(Tinput, Pinput, massflow = massflow)
                    Re, u = Reynolds_sympy(density, mu, D, u = u)
                    Ma, u = mach_number_sympy(Tinput, u = u)
                    Eu, Re = Euler_sympy(Re = Re)
                    dp, Eu, u = pdrop_sympy(density, Eu = Eu, u = u)
                    coordinates.append([Tnew, dp, Tmetal, massflow, mstar, \
                                        Re, density, u, Ma])
                Tinput = sum([temp[-1][0] for temp in direction[i]]) / len(direction[i])
                Pinput += -dp
                Re, u = Reynolds_sympy(density, mu, D, u = u)
                check_cases(Re, mstar, Ma)
                
                mu = SI("V", "T", Tinput, "P", Pinput, "helium")
                density = SI("D", "T", Tinput, "P", Pinput, "helium")
                specific_heat_capacity = SI("C", "T", Tinput, "P", Pinput, "helium")
                n = ni
                
                if flag == 1.0:
                    ni += 1
                elif ni == len(hcp[j]):
                    ni += 1
                else:
                    pass
        
        runs.append(channels)
        T_max = max([i[-1][2] for i in channels[1][0]])
        
        system_level_results.append([massflow_total, massflow_in, n_in, Pinput, \
                                     T_max, len(channels[-1]), Ma, u, Eu, Re, mstar])
        plot_temperature(channels, sbar, heatflux)
        print("--- %s seconds ---" % (time.time() - start_time))

# write_CAD_input_files(height, width, hcp_final, tcp_final, channels, cp)

print("--- %s seconds ---" % (time.time() - start_time))
  