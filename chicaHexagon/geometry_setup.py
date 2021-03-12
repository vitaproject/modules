
from math import sin, cos, tan, atan, pi, ceil, sqrt
from numpy import linspace
from sympy import symbols
from tqdm import tqdm
from operator import itemgetter
from copy import deepcopy

def initial_domain_space(width, n, gap, inner_radius, outer_radius):
    
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
    """

    cp = [] ## a list for assigning centre points of hexagon array, (x, y)
    height = width/2 * tan(60*pi/180)
    theta = (360/n/2 - gap) * pi/180
    point_x1 = inner_radius*cos(theta)
    point_y1 = inner_radius*sin(theta)
    line_angle = tan(theta)
    
    topy = outer_radius * sin(theta)
    no_points_y = int(round((topy/height), 0))
    no_points_x = int(round(((outer_radius-point_x1)/(1.5*width)), 0))
    
    a = width/2
    A = (3 * sqrt(3) / 2) * a ** 2

    for i in range(no_points_x):
        for j in range(no_points_y):
            cp.append([i*1.5*width + point_x1, j*height])
            cp.append([i*1.5*width + 0.75*width + point_x1, j*height + 0.5*height])
    
    return height, theta, point_x1, point_y1, line_angle, no_points_x, \
        cp, A

def refined_domain_space(angle, cp, x1, y1, width, inner_radius, \
                         outer_radius, x_displacement_from_origin, \
                         y_displacement_from_origin, Rstrike):

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
    """

    cp_remove = [] ## list of values to remove from cp as are out of bounds
    hcp_final = [] ## list of remaining central points
    tcp_final = []
    
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
    
    for i, value in enumerate(hcp_final):
        hcp_final[i].append(sqrt(value[0]**2 + value[1]**2)-Rstrike)
    
    return hcp_final, tcp_final
    
def radial_flow_channel_assignment(number_of_channels, theta, hcp_final, R, inner_radius, outer_radius):             ## flow channel and radial position assignment for every hexagon centre point
    
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
    
    return channels
    
def reorder(channels, R, RN):
    
    """
    Need to reorder the channels list to account for a strike point between \
        Rmin and Rmax
    
    :param float R: Radius of the injection point of the coolant into the panel
    
    :param list RN: Radial bounds to determine hexagon centre point channel \
        alocation
    :param list channels: List of lists containing hexagon centre point \ 
        coordinates determined to be within each channel
    :param channels_with_direction: Adds the channels to "direction" \
        designations. i.e. for a given R, the flow will iterate in 2 \
            directions, hence this distinction is required for calculating \
                flow properties in each direction
    """
    
    for i, Ri in enumerate(RN):
        if R >= Ri and R <= RN[i+1]:
            first_channel = i
    
    channels_with_direction = []
    
    channels_with_direction.insert(0, list(reversed(channels[:first_channel+1])))
    channels_with_direction.insert(1, channels[first_channel-1:])
    
    return channels_with_direction

def first_row(hcp_final, number_of_hexagons):
    
    """
    Divides the hexagon list about the strike point radius into inboard and \
        outboard categories. Then adds hexagons to the first row. The \
            inboard/outboard assigned hexagons will be used to populate the \
                remaining rows as they are created.

    :param float number_of_hexagons: Number of hexagons to be added to the \
        first row
    
    :param list hcp_inboard: Hexagon centre points inboard of strike point
    :param list hcp_outboard: Hexagon centre points outboard of strike point
    :param list hcp_first_row: Hexagon centre points in first row
    :param list channels: list containing lists of hexagon centre points \
        separated by row designation
    """

    
    hcp_inboard = []
    hcp_outboard = []
    
    for i in hcp_final:
        if i[2] >= 0:
            hcp_outboard.append(i)
        else:
            hcp_inboard.append(i)
    
    hcp_inboard = sorted(hcp_inboard, key=itemgetter(2), reverse = True)
    hcp_outboard = sorted(hcp_outboard, key=itemgetter(2))
    hcp_first_row = []
    
    for coordinates in hcp_inboard[:int(number_of_hexagons/2)]:
        hcp_first_row.append(coordinates)
    for coordinates in hcp_outboard[:int(number_of_hexagons/2)]:
        hcp_first_row.append(coordinates)
    
    channels = [[],[]] # two sets of lists of lists, indicating channels going inboard/outboard
    
    for direction in channels:
        direction.insert(0, deepcopy(hcp_first_row))
    
    return channels, hcp_inboard, hcp_outboard

def next_rows(hcp, n, ni):
    
    """
    Adds hexagons to the next row.
    
    :param float n: Number of hexagons to be added to the already assigned to \
        a row
    :param float ni: n + the number of hexagons to be added to the current row
    
    :param list hcp_next_row: Hexagon centre points outboard of strike point
    :param list hcp: Hexagon centre points separated by inboard/outboard \
        assignment
    """

    hcp_next_row = []
    
    for coordinates in hcp[n:ni]:
        hcp_next_row.append(deepcopy(coordinates))
    
    return hcp_next_row




