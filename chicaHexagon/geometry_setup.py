## hexagon centre point map, uses a symmetry line at the origin

from math import sin, cos, tan, atan, pi, ceil, sqrt
from numpy import linspace
from matplotlib.pyplot import scatter
from sympy import symbols
# from string import ascii_lowercase

def initial_domain_space(width, n, gap, inner_radius, outer_radius, cp):
    
    """
    Defining the initial point map of tesselated hexagon centre points 

    :param float width: Width "vertex to vertex" of a single hexagon
    :param float height: Height "flat to flat" of a single hexagon
    :param float theta: Swept angle of a single panel
    :param float gap: Swept angle defining the gap between plates
    :param float point_x1: x coordinate defining most extreme point on the arc \
        defined by inner_radius and theta
    :param float point_y1: y coordinate defining most extreme point on the arc \
        defined by inner_radius and theta
    :param float line_angle: gradient of the line joining the most extreme \
        points on the arcs defined by inner_radius, outer_radius and theta
        
    :param int no_points: Number of points in one direction to define the \
        point map of hexagon centre points
        
    :param list cp: Hexagon centre points map
    :rtype float:
    """
    
    height = width/2 * tan(60*pi/180)
    theta = (360/n/2 - gap) * pi/180
    point_x1 = inner_radius*cos(theta)
    point_y1 = inner_radius*sin(theta)
    line_angle = tan(theta)
    no_points = int(round((outer_radius/min(width, height)), 0))

    for i in range(no_points):
        for j in range(no_points):
            cp.append([i*1.5*width, j*height])
            cp.append([i*1.5*width + 0.75*width, j*height + 0.5*height])
    
    return height, theta, point_x1, point_y1, line_angle, no_points, cp

def refined_domain_space(angle, cp, x1, y1, cp_remove, width, inner_radius, \
                         outer_radius, tcp_final, hcp_final, \
                         x_displacement_from_origin, y_displacement_from_origin):

    """
    Refining the initial point map of tesselated hexagon centre points by \
        introducing a bounding box, where the internal points are kept

    :param sympy.core.symbol.Symbol x: x coorindate of hexagon centre \
        point being assessed 
    :param sympy.core.symbol.Symbol y: y coorindate of hexagon centre \
        point being assessed
    :param sympy.core.symbol.Symbol r: Radius defining a given boundaing arc
    :param sympy.core.add.Add topy: Line equation defining plate bounding edge
    :param sympy.core.add.Add arc: Arc equation defining plate inner and outer \
        bounding edge
        
    :param float R: Radius defining the coordinates of the hexagon vertices
    
    :param list hcp_final: Refined list of coordinates for every hexagon \
        centre point
    :param list tcp_final: Coordinates of the regular triangles that \
        constitute the tesselated hexagon array
    :rtype float:
    """

    x, y, r = symbols("x y r")
    topy = (angle * (x-x1)) + y1
    arc = ((r**2 - ((y-y_displacement_from_origin)**2))**0.5) + \
        x_displacement_from_origin

    for i, coordinates in enumerate(cp):
        
        if topy.subs([(x, coordinates[0])]) <= coordinates[1]: 
            cp_remove.append(coordinates) 
    
        ## inner arc
        try:
            float(arc.subs([(r, inner_radius), (y, coordinates[1])]))
        except:
            None
        else:
            if arc.subs([(r, inner_radius), (y, coordinates[1])]) >= coordinates[0]:
                cp_remove.append(coordinates)
        
        ## outer arc
        try:
            float(arc.subs([(r, outer_radius), (y, coordinates[1])]))
        except:
            None
        else:
            if arc.subs([(r, outer_radius) ,(y, coordinates[1])]) <= coordinates[0]:
                cp_remove.append(coordinates)
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
    
    """
    Plot of hexagon and triangle centre point arrays for error checking

    :param list hcp_x: x coordinates of hexagon array
    :param list hcp_y: y coordinates of hexagon array
    :param list tcp_x: x coordinates of triangle array
    :param list tcp_y: y coordinates of triangle array
    :rtype float:
    """
    
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
                          hcp_mirror):

    """
    Writes the relevant output files for importing in to FreeCAD

    :param float width: Width "vertex to vertex" of a single hexagon
    :param float height: Height "flat to flat" of a single hexagon
    
    :param list hexagon: Vertex coorindates of the first hexagon
    :param list hcp_final: Refined list of coordinates for every hexagon \
        centre point
    :param list hcp_mirror: Mirrors coordinates of hcp_final about the line \
        x = 0, ensuring not to repeat coordinates on the line x = 0
    :param list tcp_final: Coordinates of the regular triangles that \
        constitute the tesselated hexagon array
        
    :param file hexagon.asc: Height and width of a single hexagon
    :param file triangle.asc: Coordinates if tesselated triangle centre points 
    :param file IDs.asc: Height and width multipliers for every hexagon centre point \
        defining the position of each hexagon reletive to the centre point of \
            the hexagon defined in singleHexagon.asc
    :param file singleHexagon.asc: Vertex coordinates of the first hexagon
    :rtype float:
    """
    
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
            hcp_mirror.append([i[0], 0, -i[2]])
    
    for i in hcp_mirror:
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
    
def flow_channel_assignment(number_of_channels, theta, hcp_final, R, inner_radius, outer_radius):             ## flow channel and radial position assignment for every hexagon centre point
    
    """
    Splits the hexagon centre point list into "channels"", defined as the flow \
        path of the coolant from a user defined injection point R. Will be \
            used to characterise the coolant flow properties.

    :param float R: Radius of the injection point of the coolant into the panel
    :param float inner_radius: Radius defining the inner bounding arc of the \
        divertor target plate
    :param float outer_radius: Radius defining the outer bounding arc of the \
        divertor target plate
    :param float number_of_channels: Number of channels to be considered, \
        specified by the user
    :param float theta: Swept angle of a given target plate
    :param float hcp_final: Refined list of coordinates for every hexagon \
        centre point
    
    :param float n_rounded: Number of channels within the half swept angle of \
        a given plate
    :param float coordinate_angle: Angle of the centre point defined by the \
        origin (0,0)
    :param list radial_location: Radial location of the centre point defined \
        by the origin (0,0)
    :param list RN: Radial bounds to determine hexagon centre point channel \
        alocation
    :param list thetaN: Angular bounds of each channel to determine hexagon \
        centre point channel alocation
    :param list channels: List of lists containing hexagon centre point \ 
        coordinates determined to be within each channel
    
    :param file channels.asc: Coordinates
        
    :rtype float:
    """
    
    if R <= inner_radius or R >= outer_radius:
        n_rounded = ceil(number_of_channels/2)
        channels = [[] for i in range(n_rounded)]
        thetaN = [i for i in linspace(-theta, theta, int(number_of_channels)+1) if i >= 0]
        
        if thetaN[0] != 0:
            thetaN.insert(0, 0)

        coordinate_angle = [atan(i[1]/i[0]) for i in hcp_final]
    
        for j, theta in enumerate(coordinate_angle):
            for i, x in enumerate(thetaN):
                if theta >= x and theta <= thetaN[i+1]:
                    channels[i].append(hcp_final[j])

    else:
        n_rounded = ceil(number_of_channels/2)*2
        channels = [[] for i in range(n_rounded)]
        thetaN = [i for i in linspace(-theta, theta, int(number_of_channels)+1) if i >= 0]

        if thetaN[0] != 0:
            thetaN.insert(0, 0)
            thetaN.insert(int(number_of_channels)-1, 0) 

        RN = [inner_radius, R, outer_radius]

        coordinate_angle = [atan(i[1]/i[0]) for i in hcp_final]
        radial_location = [sqrt(i[1]**2+i[0]**2) for i in hcp_final]
        
        for j, theta in enumerate(coordinate_angle):
            for i, x in enumerate(thetaN):
                for k, Radius in enumerate(RN):
                    if theta >= x and theta <= thetaN[i+1]:
                        if radial_location[j] >= Radius and radial_location[j] <= RN[k+1]:
                            channels[i + int((n_rounded/2)) * k].append(hcp_final[j])
  
    for i, x in enumerate(channels):
        for j in x:
            j.append(sqrt((R-sqrt(j[0]**2+j[1]**2))**2))
        channels[i] = sorted(x, key = lambda x:x[2], reverse=False)
          
    f = open("results/channels.asc", "w")
    
    for line in channels:
        for i in line:
            for j in i:
                f.write(str(j) + "\t")
            f.write("\n")
        f.write("\n")
    
    f.write(str(thetaN) + "\n")
    
    f.close()
    
