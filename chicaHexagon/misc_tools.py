             
from scipy.interpolate import interp1d                   
from chicaHexagon.utility import get_example_data_path
import matplotlib.pyplot as plt
from matplotlib.pyplot import scatter
from math import sin, cos, pi

# ------------------------------- plot triangle and hexagon centre points ------------------------------- #
def input_data(input_variable, input_heatflux):
    
    """
    Takes the input files containing user run options and heatflux data from \
        VITA and returns this data in a format useable by CHICA
    
    :param float n_plates: Number of divertor plates in the Spherical Tokamak
    :param float width: Individual hexagon width
    :param float inner_radius: Radius defining the inboard divertor plate edge
    :param float outer_radius: Radius defining the outboard divertor plate \
        edge
    :param float x_displacement_from_origin: Displacement along the x axis of \
        the centre point defining the plate sweept angle
    :param float y_displacement_from_origin: Displacement along the y axis of \
        the centre point defining the plate sweept angle
    :param float gap: Swept angle in degrees of the gap between plates
    :param float number_of_channels: Number of channels defining the active \
        cooling layup
    :param float Rstrike: Strike point radius
    :param float channel_width: Width of a given channel
    :param float massFlow: Total massflow to a single plate
    :param float Ti: Coolant input temperature
    :param float Pi: Coolant input pressure
    :param float D: Jet diameter for Double Layer Hybrid design
    :param float t: Plate armour thickness
        
    :param file input_variable: File containing user input data
    :param file input_heatflux: File containing the heat profile output from \
        VITA
        
    :param list values: User input data from input_variable 
    :param list q: Heat flux data from input_variable
    :param list s: Distance data from input_variable, relative to q
    :param list heatflux: Adjusted q to account for adjusted s
    :param list sbar: Adjusted from s data, s starts at 0 m relative to the \
        strike point, for the purpose of considering heat fluxes inboard of \
            the strike point, s data adjusted to be relative to the Tokamak \
                centre line.
    """
    
    values = []
    q = [] # power input from VITA
    s = [] # distance from strike point
    
    user_data = get_example_data_path(input_variable)
    heatflux_data = get_example_data_path(input_heatflux)
    
    with open(user_data) as inputs:
        for line in inputs:
            words = line.split()
            values.append(float(words[len(words)-1:][0]))
    
    with open(heatflux_data) as geometry:
        for line in geometry:
            heatflux = line.split()
            s.append(float(heatflux[0])/100.0)
            q.append(float(heatflux[1]))
    
    n_plates, width, inner_radius, outer_radius, x_displacement_from_origin, \
        y_displacement_from_origin, gap, Rstrike, channel_width, massFlow, Ti, \
            Pi, D, t = values
    
    sbar = [i + Rstrike for i in s]
    heatflux = [i * 1e6 for i in q]
    
    sbar.insert(0, Rstrike - (sbar[1] - sbar[0]))
    sbar.insert(0, 0)
    
    heatflux.insert(0,0)
    heatflux.insert(0,0)
    
    HFf = interp1d(sbar, heatflux)
    
    return n_plates, width, inner_radius, outer_radius, \
           x_displacement_from_origin, y_displacement_from_origin, gap, \
           Rstrike, channel_width, massFlow, Ti, Pi, D, t, HFf, sbar, heatflux


def plot_centre_points(hcp_final, tcp_final):
    
    """
    Plot of hexagon and triangle centre point arrays for error checking

    :param list hcp_x: x coordinates of hexagon array
    :param list hcp_y: y coordinates of hexagon array
    :param list tcp_x: x coordinates of triangle array
    :param list tcp_y: y coordinates of triangle array
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

def plot_temperature(channels, sbar, heatflux):
    
    """

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
    
    xyT = []
    
    for direction in channels:
        for channel in direction:
            for coordinates in channel:
                xyT.append(coordinates)
    
    stored_level = -1
        
    x = [item[0] for item in xyT]
    y = [item[1] for item in xyT]
    Tcoolant = [item[stored_level][0] for item in xyT]
    p = [item[stored_level][1] for item in xyT]
    Tmetal = [item[stored_level][2] for item in xyT]
    massflow = [item[stored_level][3] for item in xyT]
    mstar = [item[stored_level][4] for item in xyT]
    Re = [item[stored_level][5] for item in xyT]
    density = [item[stored_level][6] for item in xyT]
    velocity = [item[stored_level][7] for item in xyT]
    Mach = [item[stored_level][8] for item in xyT]
    
    fig = plt.figure()
    fig1 = plt.figure()
    fig2 = plt.figure()
    fig3 = plt.figure()
    fig4 = plt.figure()
    fig5 = plt.figure()
    fig6 = plt.figure()
    fig7 = plt.figure()
    fig8 = plt.figure()
    
    ax = fig.add_subplot(111, projection='3d')
    ax1 = fig1.add_subplot(111, projection='3d')
    ax2 = fig2.add_subplot(111, projection='3d')
    ax3 = fig3.add_subplot(111, projection='3d')
    ax4 = fig4.add_subplot(111, projection='3d')
    ax5 = fig5.add_subplot(111, projection='3d')
    ax6 = fig6.add_subplot(111, projection='3d')
    ax7 = fig7.add_subplot(111, projection='3d')
    ax8 = fig8.add_subplot(111, projection='3d')
    
    ax.set_title("Coolant Temperature")
    ax1.set_title("Pressure Drop")
    ax2.set_title("Metal Temperature")
    ax3.set_title("Mass Flow per Hexagon")
    ax4.set_title("Non Dimensional Mass Flow")
    ax5.set_title("Reynolds Number")
    ax6.set_title("Density")
    ax7.set_title("Velocity")
    ax8.set_title("Mach")
    
    ax.set_zlabel("Temperature [K]")
    ax1.set_zlabel("Pressure [Pa]")
    ax2.set_zlabel("Temperature [K]")
    ax3.set_zlabel("Mass flow [kg/s]")
    ax4.set_zlabel("[]")
    ax5.set_zlabel("[]")
    ax6.set_zlabel("Density, [kg/m3]")
    ax7.set_zlabel("Velocity, [m/s]")
    ax8.set_zlabel("[]")
    
    ax.scatter(x, y, Tcoolant, s=0.5, c = Tcoolant, depthshade = True)
    ax1.scatter(x, y, p, s=0.5, c = p, depthshade = True)
    ax2.scatter(x, y, Tmetal, s=0.5, c = Tmetal, depthshade = True)
    ax3.scatter(x, y, massflow, s=0.5, c = massflow, depthshade = True)
    ax4.scatter(x, y, mstar, s=0.5, c = mstar, depthshade = True)
    ax5.scatter(x, y, Re, s=0.5, c = Re, depthshade = True)
    ax6.scatter(x, y, density, s=0.5, c = density, depthshade = True)
    ax7.scatter(x, y, velocity, s=0.5, c = velocity, depthshade = True)
    ax8.scatter(x, y, Mach, s=0.5, c = Mach, depthshade = True)
    
    plt.plot(sbar, heatflux)
    plt.show()
    
def write_CAD_input_files(height, width, hcp_final, tcp_final, channels, cp):

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
    """
    
    hcp_ID = [] ## a list of the matrix ID's left after removing out of bounds points
    hcp_mirror = []
    hexagon = [] ## a list of coordinates describing the points of the regular hexagon
    hcp_final_y = []
    hcp_final_y_mirror = []
    tcp_final_y = []
    tcp_final_y_mirror = []
    
    R = (height/2) / sin(pi/3)
    for j in range(7):
            hexagon.append([(hcp_final[0][0] + R*cos(j * pi/3))*1000, \
                            0, \
                            (hcp_final[0][1] + R*sin(j * pi/3))*1000])
    
    # ------------------------------- create ID list for centre points ------------------------------- #
    
    ## assign ID's to centre points
    for i in hcp_final:
        hcp_ID.append([i[0]/width - hcp_final[0][0]/width,0, i[1]/height])
        hcp_final_y.append([i[0]*1000, 0, i[1]*1000])
    
    for i in hcp_final_y:
        if i[2] == 0:
            None
        else:
            hcp_final_y_mirror.append([i[0], 0, -i[2]])
    
    for i in hcp_final_y_mirror:
        hcp_final_y.append(i)
    
    for i in tcp_final:
        tcp_final_y.append([i[0], 0, i[2]])
    
    for i in tcp_final_y:
        if i[2] == 0:
            None
        else:
            tcp_final_y_mirror.append([i[0], 0, -i[2]])
    
    for i in tcp_final_y_mirror:
        tcp_final_y.append(i)
    print(tcp_final_y)
    for i in hcp_ID:
        if i[2] == 0:
            None
        else:
            hcp_mirror.append([i[0], 0, -i[2]])
    
    for i in hcp_mirror:
        hcp_ID.append(i)
        
    # ------------------------------- write files for FreeCAD ------------------------------- #
    
    file_path = "results/"
    data_path = ""
    
    f = open(data_path + file_path + "hexagon.asc", "w")
    
    for line in hexagon:
        for i in line:
            f.write(str(i) + "\t")
        f.write("\n")
    
    f.close()
    
    f = open(data_path + file_path + "IDs.asc", "w")
    
    for line in hcp_ID:
        for i in line:
            f.write(str(i) + "\t")
        f.write("\n")
    
    f.close()
    
    f = open(data_path + file_path + "singleHexagon.asc", "w")
    
    geometry = str(width) + " " + str(height)
    f.write(geometry)
    
    f.close()
    
    f = open(data_path + file_path + "channels.asc", "w")
    
    for direction in channels:
        for line in direction:
            for i in line:
                for j in i:
                    f.write(str(j) + "\t")
                f.write("\n")
            f.write("\n")
    
    f.close()
    
    f = open(data_path + file_path + "hexagon_centre_points.asc", "w")
    
    for line in hcp_final_y:
        for i in line:
            f.write(str(i) + "\t")
        f.write("\n")
    
    f.close()
    
    f = open(data_path + file_path + "triangle_centre_points.asc", "w")
    
    for line in tcp_final_y:
        for i in line:
            f.write(str(i) + "\t")
        f.write("\n")
    
    f.close()