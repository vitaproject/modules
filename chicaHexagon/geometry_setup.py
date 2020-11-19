## hexagon centre point map, uses a symmetry line at the origin

from math import sin, cos, tan, pi
from matplotlib.pyplot import scatter
from sympy import symbols

# ------------------------------- initial geometry calcs ------------------------------- #

def initial_domain_space(width, n, gap, inner_radius, outer_radius, cp):
    height = width/2 * tan(60*pi/180);                                         ## height of a hexagon (also the displacement between centre points) from Qdot, m
    theta = (360/n/2 - gap) * pi/180                                           ## half swept angle of the plate, rads
    point_x1 = inner_radius*cos(theta)
    point_y1 = inner_radius*sin(theta)
    line_angle = tan(theta)
    no_points = int(round((outer_radius/min(width, height)), 0))

# ------------------------------- point map ------------------------------- #

    for i in range(no_points):
        for j in range(no_points):
            cp.append([i*1.5*width, j*height])
            cp.append([i*1.5*width + 0.75*width, j*height + 0.5*height])
    
    return height, theta, point_x1, point_y1, line_angle, no_points, cp

def refined_domain_space(line_angle, cp, point_x1, point_y1, cp_remove, \
                         width, inner_radius, outer_radius, tcp_final, \
                         hcp_final, x_displacement_from_origin, \
                         y_displacement_from_origin):

# ------------------------------- expressions defining the bounding area ------------------------------- #

    x, y, r, k, h, x1, y1, angle = symbols("x y r k h x1 y1 angle")
    topy = (angle * (x-x1)) + y1
    arc = ((r**2 - ((y-k)**2))**0.5) + h
    
# ------------------------------- remove out of bounds values ------------------------------- #
    
    for i, z in enumerate(cp):
        
        if topy.subs([(x, z[0]), (x1, point_x1), (y1, point_y1), \
                      (angle, line_angle)]) <= z[1]: 
            cp_remove.append(z) 
    
        ## inner arc
        try:
            float(arc.subs([(r, inner_radius),(h, x_displacement_from_origin),\
                            (k, y_displacement_from_origin),(y, z[1])]))
        except:
            None
        else:
            if arc.subs([(r, inner_radius),(h, x_displacement_from_origin),\
                         (k, y_displacement_from_origin),(y, z[1])]) >= z[0]:
                cp_remove.append(z)
        
        ## outer arc
        try:
            float(arc.subs([(r, outer_radius),(h, x_displacement_from_origin),\
                            (k, y_displacement_from_origin),(y, z[1])]))
        except:
            None
        else:
            if arc.subs([(r, outer_radius),(h, x_displacement_from_origin),\
                         (k, y_displacement_from_origin),(y, z[1])]) <= z[0]:
                cp_remove.append(z)
    
    ## remove values from the total list and creates a new "final" list
    for i in cp:
        if i not in cp_remove: 
            hcp_final.append(i)

    # ------------------------------- create centre points for triangles ------------------------------- #
            
    R = (width / 4) / cos(pi/6)
    for i in hcp_final:
        for j in range(7):
            tcp_final.append([(i[0] + R*cos(j * pi/3 + pi/6))*1000, \
                              0, \
                              (i[1] + R*sin(j * pi/3 + pi/6))*1000])
    
    return hcp_final, tcp_final

# ------------------------------- plot triangle and hexagon centre points ------------------------------- #

def plot_centre_points(hcp_final, tcp_final):
    
    hcp_x = hcp_y = tcp_x = tcp_y = []
    
    for i in hcp_final:
        hcp_x.append(i[0])
        hcp_y.append(i[1])
    
    for i in tcp_final:
        tcp_x.append(i[0])
        tcp_y.append(i[1])
    
    #scatter(hcp_x, hcp_y, s = 1, c = "magenta")
    scatter(tcp_x, tcp_y, s = 1)

# ------------------------------- create vertex coordinates for first hexagon ------------------------------- #

def write_CAD_input_files(height, width, hcp_ID, hexagon, hcp_final, tcp_final, \
                          cp_mirror):
    
    R = (height/2) / sin(pi/3)
    for j in range(7):
            hexagon.append([(hcp_final[0][0] + R*cos(j * pi/3))*1000, \
                            0, \
                            (hcp_final[0][1] + R*sin(j * pi/3))*1000])
    
    # ------------------------------- create ID list for centre points ------------------------------- #
    
    ## assign ID's to centre points
    for i in hcp_final:
        hcp_ID.append([i[0]/width - hcp_final[0][0]/width,0, i[1]/height])
    
    for i in hcp_ID:
        if i[2] == 0:
            None
        else:
            cp_mirror.append([i[0], 0, -i[2]])
    
    for i in cp_mirror:
        hcp_ID.append(i)
    
    # ------------------------------- write files for FreeCAD ------------------------------- #
       
    f = open("results/hexagon.asc", "w")
    
    for line in hexagon:
        for i in line:
            f.write(str(i) + "\t")
        f.write("\n")
    
    f.close()
    
    f = open("results/triangle.asc", "w")
    
    for line in tcp_final:
        for i in line:
            f.write(str(i) + "\t")
        f.write("\n")
    
    f.close()
    
    f = open("results/IDs.asc", "w")
    
    for line in hcp_ID:
        for i in line:
            f.write(str(i) + "\t")
        f.write("\n")
    
    f.close()
    
    f = open("results/singleHexagon.asc", "w")
    
    geometry = str(width) + " " + str(height)
    f.write(geometry)
    
    f.close()