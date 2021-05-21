from scipy.interpolate import interp1d
from math import sqrt, sin, cos, tan, pi
from CHICA.utility import get_example_data_path
import xml.etree.ElementTree as ET
from CoolProp.CoolProp import PropsSI as SI
from sympy import symbols
from numpy import linspace
import time
from CHICA.flow_properties import mach_number_sympy, \
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

class Material:

    def helium(self, input_pressure, input_temperature):
        
        self.viscosity = SI("V", "T", input_temperature, "P", input_pressure, "helium")
        self.density = SI("D", "T", input_temperature, "P", input_pressure, "helium")
        self.specific_heat_capacity = SI("C", "T", input_temperature, "P", input_pressure, "helium")

class Hexagon:
    
    ID = 0
    
    def __init__(self, hexagon_data, HFf, cross_sectional_area, specific_heat_capacity, \
                 input_temperature, mass_flow, non_dimensional_temperature):

        self.x = hexagon_data[0]
        self.y = hexagon_data[2]
        self.R = sqrt((self.x**2) + (self.y**2))
        heatflux = HFf(self.R) * cross_sectional_area
        self.developed_coolant_temperature, self.mass_flow = coolant_temperature_sympy(heatflux, \
            specific_heat_capacity, input_temperature, massflow = mass_flow) 
        self.metal_temperature, self.developed_coolant_temperature = \
            metal_temperature_sympy(input_temperature, non_dimensional_temperature, \
            Tnew = self.developed_coolant_temperature) 
        
        self.ID = self.identification()
        
    @classmethod
    def identification(cls):
        cls.ID += 1
        return cls.ID

class Channel:
    
    def __init__(self, input_pressure, input_temperature, panel_hexagons, \
                 HFf, cross_sectional_area, Ma = 0.0):
        
        ## panel level properties i.e. only get assigned/changed at this level
        self.panel_hexagons = panel_hexagons # needs to be accessed at hexagon level
        self.HFf = HFf # needs to be accessed at hexagon level
        self.cross_sectional_area = cross_sectional_area # needs to be accessed at hexagon level
        self.mass_flow_total: float = None # calculated in the first channel
        self.jet_cross_sectional_area: float = None # needs to be accessed at channel level, calculated in the first channel
        self.jet_diameter: float = None # needs to be accessed at channel level, calculated in the first channel
        
        ## channel level properties i.e. only get assigned/changed at this level
        self.input_pressure = input_pressure # needs to be accessed at hexagon level
        self.input_temperature = input_temperature # needs to be accessed at hexagon level
        self.mach_number = Ma # needs to be accessed at channel level
        self.reynolds_number: float = None # needs to be accessed at channel level
        self.jet_velocity: float = None # needs to be accessed at channel level
        self.euler_number: float = None # needs to be accessed at channel level
        self.pressure_drop: float = None # needs to be accessed at channel level
        self.non_dimensional_temperature: float = None # needs to be accessed at channel level
        self.non_dimensional_mass_flow: float = None # needs to be accessed at channel level
        self.mass_flow: float = None # needs to be accessed at hexagon level
        self.fluid_properties = Material()
        self.fluid_properties.helium(input_pressure, input_temperature) # needs to be accessed at hexagon level
        self.hexagons: list = [] # needs to be accessed at hexagon level
        self.output_temperature: float = None # needs to be accessed at channel level
        self.output_pressure = input_pressure # needs to be accessed at channel level
        self.max_temperature: float = None # needs to be accessed at channel level
        self.channel_ID: int = None
        
    def repeated_operations(self):
        
        # potential properties: 
        # fluid input properties calculated:
        # viscosity, Reynolds number, Euler number, 
        # fluid developed properties calculated:
        # pressure drop
        
        self.non_dimensional_mass_flow, self.mass_flow, self.fluid_properties.specific_heat_capacity = \
            mstar_sympy(self.input_temperature, self.input_pressure, self.cross_sectional_area, \
            massflow = self.mass_flow)
        self.reynolds_number, self.jet_velocity = Reynolds_sympy(self.fluid_properties.density, \
            self.fluid_properties.viscosity, self.jet_diameter, u = self.jet_velocity)
        self.euler_number, self.reynolds_number = Euler_sympy(Re = self.reynolds_number)
        self.pressure_drop, self.euler_number, self.jet_velocity = pdrop_sympy(\
            self.fluid_properties.density, Eu = self.euler_number, u = self.jet_velocity)
        self.non_dimensional_temperature, self.non_dimensional_mass_flow = \
            taus_sympy(mstar = self.non_dimensional_mass_flow)
            
    def get_hexagons(self, hexagon_selection):
        
        for hexagon in hexagon_selection:
            
            self.hexagons.append(Hexagon(hexagon, self.HFf, self.cross_sectional_area, \
                self.fluid_properties.specific_heat_capacity, self.input_temperature, self.mass_flow, \
                self.non_dimensional_temperature))
                
            self.panel_hexagons.remove(hexagon)
    
    @classmethod
    def create_first_channel(cls, input_pressure, input_temperature, panel_hexagons, \
                 HFf, cross_sectional_area, number_of_hexagons, Ma=0.0):


        obj = cls(input_pressure, input_temperature, panel_hexagons, \
                         HFf, cross_sectional_area, Ma)

        obj.mach_number, obj.jet_velocity = mach_number_sympy(obj.input_temperature, \
                                                                Ma=obj.mach_number)
        obj.mass_flow, obj.jet_velocity, obj.jet_cross_sectional_area, \
        obj.density, obj.jet_diameter = mass_flow_sympy(obj.input_temperature, \
                                                          obj.input_pressure, u=obj.jet_velocity)

        obj.mass_flow_total = number_of_hexagons * obj.mass_flow

        obj.repeated_operations()

        reordered_hexagons = sorted(obj.panel_hexagons, key=lambda x:x[3], reverse=False)
        hexagon_selection = reordered_hexagons[:number_of_hexagons]

        obj.get_hexagons(hexagon_selection)

        obj.output_temperature = sum([hexagon.developed_coolant_temperature for hexagon in obj.hexagons]) \
                                  / number_of_hexagons
        obj.output_pressure -= obj.pressure_drop
        obj.max_temperature = max([hexagon.metal_temperature for hexagon in obj.hexagons])
        
        return obj
        
    @classmethod
    def create_channel(cls, panel_hexagons, HFf, cross_sectional_area, \
                     remaining_hexagons, previous_channel, mass_flow_actual, Ma=0.0):

        obj = cls(previous_channel.output_pressure, previous_channel.output_temperature, panel_hexagons, \
                         HFf, cross_sectional_area)

        obj.mass_flow = previous_channel.mass_flow
        obj.max_temperature_from_previous_channel = previous_channel.max_temperature
        
        obj.mass_flow, number_of_hexagons_required = massflow_nextrows( \
            remaining_hexagons[0], obj.mass_flow, obj.fluid_properties.specific_heat_capacity, \
            obj.cross_sectional_area, obj.max_temperature_from_previous_channel, \
            obj.input_temperature, obj.HFf, mass_flow_actual)

        if obj.mass_flow == 0.0:
            number_of_hexagons_required = len(remaining_hexagons)
            obj.mass_flow = mass_flow_actual / number_of_hexagons_required

        elif number_of_hexagons_required <= len(remaining_hexagons):
            number_of_hexagons_required = int(number_of_hexagons_required)
            obj.mass_flow = mass_flow_actual / number_of_hexagons_required

        else:
            number_of_hexagons_required = len(remaining_hexagons)
            obj.mass_flow = mass_flow_actual / number_of_hexagons_required

        obj.mass_flow, obj.jet_velocity, obj.jet_cross_sectional_area, \
        obj.fluid_properties.density, obj.jet_diameter = mass_flow_sympy(obj.input_temperature, \
                                                                           obj.input_pressure, massflow=obj.mass_flow)

        obj.mach_number, obj.jet_velocity = mach_number_sympy(obj.input_temperature, \
                                                                u=obj.jet_velocity)

        obj.repeated_operations()

        obj.get_hexagons(remaining_hexagons[:number_of_hexagons_required])
        del remaining_hexagons[:number_of_hexagons_required]

        obj.output_temperature = sum([hexagon.developed_coolant_temperature for hexagon in obj.hexagons]) \
                                  / number_of_hexagons_required
        obj.output_pressure = obj.input_pressure - obj.pressure_drop
        obj.max_temperature = previous_channel.max_temperature
        
        return obj

class Panel:
    
    def __init__(self, mach_number, number_of_hexagons_in_first_channel, mass_split, \
                 input_file, heat_flux_file):
        
        self.channels = []
        self.panel_hexagon = []
        self.mach_number = mach_number
        self.number_of_hexagons_in_first_channel = number_of_hexagons_in_first_channel
        self.mass_flow_split = mass_split
        self.mass_flow_total:float = None
        self.read_inputs(input_file, heat_flux_file)

    def configure(self):
        self.populate_hexagons()
        self.add_first_channel()
        self.add_subsequent_channels(self.channels[0])

    def add_first_channel(self):
        
        self.channels.append(Channel.create_first_channel(self.input_data["input_pressure"], 
            self.input_data["input_temperature"], self.panel_hexagons, self.HFf, self.cross_sectional_area, \
            self.number_of_hexagons_in_first_channel, Ma = self.mach_number))
        
    def add_subsequent_channels(self, first_channel):

        new_channels = [[], []]
        
        reordered_hexagons = sorted(self.panel_hexagons, key = lambda x:x[3], reverse = False)
        inboard_hexagons, outboard_hexagons = [hexagons for hexagons in reordered_hexagons if hexagons[-1] <= 0.0], \
            [hexagons for hexagons in reordered_hexagons if hexagons[-1] >= 0.0]
        all_remaining_hexagons = [inboard_hexagons, outboard_hexagons]
        
        self.mass_flow_total = first_channel.mass_flow_total
        
        self.mass_flow_actual = [self.mass_flow_total * (1-self.mass_flow_split), self.mass_flow_total * \
            self.mass_flow_split]
        
        for index, direction in enumerate(all_remaining_hexagons):
            
            mass_flow_actual = self.mass_flow_actual[index]
            
            new_channels[index].append(Channel.create_channel(self.panel_hexagons, \
                    self.HFf, self.cross_sectional_area, direction, first_channel, \
                    mass_flow_actual))
            
            while len(direction) >= 1.0:
                new_channels[index].append(Channel.create_channel(self.panel_hexagons, \
                    self.HFf, self.cross_sectional_area, direction, new_channels[index][-1], \
                    mass_flow_actual))
        
        new_channels[0].reverse()
        for channel in new_channels[0]:
            self.channels.insert(0, channel)
        
        for channel in new_channels[1]:
            self.channels.append(channel)
        
    def read_inputs(self, inputs_filename, q_filename):
        # read asc files and store results to member variables (self...)
        
        user_data = get_example_data_path(inputs_filename)
        heatflux_data = get_example_data_path(q_filename)
        
        input_data = dict()
        tree = ET.parse(user_data)
        trunk = tree.getroot()
        
        for root in trunk:
            input_data[root.tag] = float(root.text)
        
        self.input_data = input_data
        
        q = [] # power input from VITA
        s = [] # distance from strike point
        
        with open(heatflux_data) as geometry:
            for line in geometry:
                heatflux = line.split()
                s.append(float(heatflux[0])/100.0)
                q.append(float(heatflux[1])) # /6.87 * 30.0)
        
        sbar = [i + self.input_data["strike_radius"] for i in s]
        heatflux = [i * 1e6 for i in q]
        
        sbar.insert(0, self.input_data["strike_radius"] - (sbar[1] - sbar[0]))
        sbar.insert(0, 0)
        
        heatflux.insert(0,0)
        heatflux.insert(0,0)
        
        self.HFf = interp1d(sbar, heatflux)
        
    def populate_hexagons(self):
        
        cp = []  ## a list for assigning centre points of hexagon array, (x, y)
        cp_remove = [] ## list of values to remove from cp as are out of bounds
        hcp_final = [] ## list of remaining central points
        hcp_ID = [] ## a list of the matrix ID's left after removing out of bounds points
        hcp_final_y = []
        hcp_final_y_mirror = []
        
        hexagon_height = self.input_data["hexagon_width"] / 2 * tan(60 * pi / 180)
        theta = (360.0 / self.input_data["number_of_plates"] / 2 - self.input_data["swept_angle_gap_between_plates"]) * pi / 180
        
        point_x1 = self.input_data["inner_radius"]*cos(theta)
        point_y1 = self.input_data["inner_radius"]*sin(theta)
        line_angle = tan(theta)
        
        topy = self.input_data["outer_radius"] * sin(theta)
        no_points_y = int(round((topy/hexagon_height), 0))
        no_points_x = int(round(((self.input_data["outer_radius"]-point_x1)/(1.5*self.input_data["hexagon_width"])), 0))
        
        a = self.input_data["hexagon_width"]/2
        self.cross_sectional_area = (3 * sqrt(3) / 2) * a ** 2
    
        for i in range(no_points_x):
            for j in range(no_points_y):
                cp.append([i*1.5*self.input_data["hexagon_width"] + point_x1, j*hexagon_height])
                cp.append([i*1.5*self.input_data["hexagon_width"] + 0.75*self.input_data["hexagon_width"] + \
                            point_x1, j*hexagon_height + 0.5*hexagon_height])
    
        x, y, r = symbols("x y r")
        topy = (line_angle * (x-point_x1)) + point_y1
        arc = ((r**2 - ((y-self.input_data["y_displacement_from_origin"])**2))**0.5) + \
            self.input_data["x_displacement_from_origin"]
        
        for i, coordinates in enumerate(cp):
            
            if topy.subs([(x, coordinates[0])]) <= coordinates[1]: 
                cp_remove.append(coordinates) 
        
            ## inner arc
            try:
                float(arc.subs([(r, self.input_data["inner_radius"]), (y, coordinates[1])]))
            except:
                None
            else:
                if arc.subs([(r, self.input_data["inner_radius"]), (y, coordinates[1])]) >= coordinates[0]:
                    cp_remove.append(coordinates)
            
            ## outer arc
            try:
                float(arc.subs([(r, self.input_data["outer_radius"]), (y, coordinates[1])]))
            except:
                None
            else:
                if arc.subs([(r, self.input_data["outer_radius"]) ,(y, coordinates[1])]) <= coordinates[0]:
                    cp_remove.append(coordinates)
                                                                                    ## remove values from the total list and creates a new "final" list
        for i in cp:
            if i not in cp_remove: 
                hcp_final.append(i)
        
        # ------------------------------- create centre points for triangles ------------------------------- #
        
        for i in hcp_final:
            hcp_ID.append([i[0]/self.input_data["hexagon_width"] - hcp_final[0][0]/ \
                           self.input_data["hexagon_width"],0, i[1]/hexagon_height])
            hcp_final_y.append([i[0], 0, i[1]])
    
        for i in hcp_final_y:
            if i[2] == 0:
                None
            else:
                hcp_final_y_mirror.append([i[0], 0, -i[2]])
        
        for i in hcp_final_y_mirror:
            hcp_final_y.append(i)
        
        for i, value in enumerate(hcp_final_y):
            hcp_final_y[i].append(sqrt((sqrt(value[0]**2 + value[1]**2)-self.input_data["strike_radius"])**2))
            hcp_final_y[i].append(sqrt(value[0]**2 + value[1]**2)-self.input_data["strike_radius"])
        
        self.panel_hexagons = hcp_final_y
        
if __name__ == "__main__":
    # start_time = time.time()
    
    ## BUSINESS CODE
    
    results = dict()
    number_of_hexagons_range = linspace(500, 2000, 1) #4
    mach_number_range = linspace(0.15, 0.3, 1) #4
    mass_split_range = linspace(0.5, 0.9, 1) #3
    
    for number_of_hexagons in number_of_hexagons_range:
        for mach_number in mach_number_range:
            for mass_split in mass_split_range:
        
                print([number_of_hexagons, mach_number, mass_split])
                
                key = "panel" + "_" + \
                    str(number_of_hexagons) + "_" + \
                    str(mach_number) + "_" + \
                    str(mass_split)
                
                panel = Panel(float(mach_number), int(number_of_hexagons), float(mass_split), "inputs.xml", "q_adjusted.asc")
                panel.configure()
        
                results[key] = panel
    
    keys = results.keys()
            
    f = open("hexagon_centre_points.asc", "w")
    for key in keys:
        for channel in results[key].channels:
            for hexagon in channel.hexagons:
                f.write(str(hexagon.x) + "\t" + \
                    "0" + "\t" + \
                        str(hexagon.y))
                f.write("\n")
    f.close()
    
    f = open("hexagon_data.asc", "w")
    for key in keys:
        for channel in results[key].channels:
            for hexagon in channel.hexagons:
                f.write(str(hexagon.x) + "\t" + \
                    str(hexagon.y) + "\t" + \
                        str(hexagon.metal_temperature) + "\t" + \
                            str(hexagon.developed_coolant_temperature))
                f.write("\n")
        f.write(key)
        f.write("\n")
    f.close()
    
    f = open("channel_data.asc", "w")
    for key in keys:
        for index, channel in enumerate(results[key].channels):
            f.write(str(key) + "\t" + \
                    str(index) + "\t" + \
                    str(len(channel.hexagons)) + "\t" + \
                    str(channel.output_pressure) + "\t" + \
                    str(channel.mach_number) + "\t" + \
                    str(channel.mass_flow) + "\t" + \
                    str(channel.mass_flow * len(channel.hexagons)) + "\t" + \
                    str(channel.reynolds_number) + "\t" + \
                    str(channel.euler_number) + "\t" + \
                    str(channel.jet_velocity) + "\t" + \
                    str(channel.non_dimensional_temperature) + "\t" + \
                    str(channel.non_dimensional_mass_flow) + "\t" + \
                    str(channel.max_temperature) + "\t" + \
                    str(channel.fluid_properties.density) + "\t" + \
                    str(channel.fluid_properties.viscosity) + "\t" + \
                    str(channel.fluid_properties.specific_heat_capacity))
            f.write("\n")
    f.close()