
from chicaHexagon.geometry_setup_update import initial_domain_space, \
                                               refined_domain_space, \
                                               toroidal_flow_channel_assignment
                                   
from chicaHexagon.flow_properties_update import coolant_temperature

from chicaHexagon.misc_tools import write_CAD_input_files, \
                                    plot_temperature, \
                                    input_data

import time

# ------------------------------- centre point lists ------------------------------- #

start_time = time.time()

n_plates, width, inner_radius, outer_radius, x_displacement_from_origin, \
           y_displacement_from_origin, gap, Rstrike, \
           channel_width, massFlow, Ti, Pi, D, t, HFf, sbar, heatflux = input_data("inputs.asc", "q.asc")

height, theta, point_x1, point_y1, line_angle, no_points_x,  \
    cp, A = initial_domain_space(width, n_plates, gap, inner_radius, outer_radius)
    
hcp_final, tcp_final = refined_domain_space(line_angle, cp, point_x1, \
    point_y1, width, inner_radius, outer_radius, x_displacement_from_origin, \
    y_displacement_from_origin)

#--------------------------------------- flow properties -------------------------------------#
channels, RN = toroidal_flow_channel_assignment(theta, hcp_final, Rstrike, inner_radius, \
                                            outer_radius, channel_width)

channels, sbar, heatflux, values_per_channel = coolant_temperature(massFlow, channels, Ti, Pi, \
                              inner_radius, outer_radius, cp, hcp_final, A, D, \
                                  heatflux, sbar, t)

plot_temperature(cp, channels, hcp_final, sbar, heatflux, RN, theta)

write_CAD_input_files(height, width, hcp_final, tcp_final, channels, cp, \
                        heatflux, values_per_channel)

plot_temperature(channels, sbar, heatflux)
write_CAD_input_files(height, width, hcp_final, tcp_final, channels, cp)

print("--- %s seconds ---" % (time.time() - start_time))



    