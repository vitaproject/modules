
from chicaHexagon.geometry_setup import initial_domain_space, \
                                           refined_domain_space, \
                                           plot_centre_points, \
                                           write_CAD_input_files
from chicaHexagon.utility import get_example_data_path
import time

# ------------------------------- centre point lists ------------------------------- #

values = []
cp = [] ## a list for assigning centre points of hexagon array, (x, y)
cp_remove = [] ## list of values to remove from cp as are out of bounds
cp_mirror = []
hcp_final = [] ## list of remaining central points
hcp_x = [] ## a list of x values for plotting
hcp_y = [] ## a list of y values for plotting
hcp_ID = [] ## a list of the matrix ID's left after removing out of bounds points
hexagon = [] ## a list of coordinates describing the points of the regular hexagon
tcp_final = []
tcp_x = [] ## list of centre points of triangles that make up the repeating hexagon
tcp_y = []
triangle = [] ## a list of coordinates describing the points of the equilateral triangle

start_time = time.time()

input_data = get_example_data_path("inputs.asc")

with open(input_data) as inputs:
    for line in inputs:
        words = line.split()
        values.append(float(words[len(words)-1:][0]))

n, width, inner_radius, outer_radius, x_displacement_from_origin, \
    y_displacement_from_origin, gap = values

height, theta, point_x1, point_y1, line_angle, no_points, cp = \
    initial_domain_space(width, n, gap, inner_radius, outer_radius, cp)
    
hcp_final, tcp_final = refined_domain_space(line_angle, cp, point_x1, \
    point_y1, cp_remove, width, inner_radius, outer_radius, tcp_final, \
    hcp_final, x_displacement_from_origin, y_displacement_from_origin)

plot_centre_points(hcp_final, tcp_final)
    
write_CAD_input_files(height, width, hcp_ID, hexagon, hcp_final, tcp_final, \
                          cp_mirror)

print("--- %s seconds ---" % (time.time() - start_time))