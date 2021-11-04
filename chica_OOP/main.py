from scipy.interpolate import interp1d
from math import sqrt, sin, cos, tan, pi, atan
from CHICA.utility import get_example_data_path
import xml.etree.ElementTree as ET
from CoolProp.CoolProp import PropsSI as SI
from sympy import symbols
from numpy import linspace, array
from copy import deepcopy
import time
import matplotlib.pyplot as plt
from CHICA.flow_properties import mach_number_sympy, \
                                        mass_flow_sympy, \
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
        self.viscosity = SI("V", "T", input_temperature, "P", input_pressure,
                            "helium")
        self.density = SI("D", "T", input_temperature, "P", input_pressure,
                          "helium")
        self.specific_heat_capacity = SI("C", "T", input_temperature, "P",
                                         input_pressure, "helium")

class VJ_cell:

    VJ_cell_ID = 0

    def __init__(self, x, y, HFf, cross_sectional_area, specific_heat_capacity, \
                 input_temperature, mass_flow, non_dimensional_temperature):

        self.x = x
        self.y = y
        self.R = sqrt(x**2 * y**2)

        heatflux = HFf(self.R) * cross_sectional_area
        self.developed_coolant_temperature, self.mass_flow = \
            coolant_temperature_sympy(heatflux, specific_heat_capacity,
                                      input_temperature, massflow=mass_flow)
        self.metal_temperature, self.developed_coolant_temperature = \
            metal_temperature_sympy(input_temperature,
                                    non_dimensional_temperature,
                                    Tnew=self.developed_coolant_temperature)

        self.VJ_cell_ID = self.identification()

    @classmethod
    def identification(cls):
        cls.VJ_cell_ID += 1
        return cls.VJ_cell_ID

class Hexagon:

    hexagon_ID = 0

    def __init__(self, hexagon_data, HFf, cross_sectional_area,
                 specific_heat_capacity, input_temperature, mass_flow,
                 non_dimensional_temperature):

        self.x = hexagon_data[0]
        self.y = hexagon_data[2]
        self.R = sqrt((self.x**2) + (self.y**2))
        heatflux = HFf(self.R) * cross_sectional_area
        self.developed_coolant_temperature, self.mass_flow = \
            coolant_temperature_sympy(heatflux, specific_heat_capacity,
                                      input_temperature, massflow=mass_flow)
        self.metal_temperature, self.developed_coolant_temperature = \
            metal_temperature_sympy(input_temperature,
                                    non_dimensional_temperature,
                                    Tnew=self.developed_coolant_temperature)
        self.hexagon_ID = self.identification()

    @classmethod
    def identification(cls):
        cls.hexagon_ID += 1
        return cls.hexagon_ID

class Channel:

    def __init__(self, input_pressure, input_temperature, panel_hexagons,
                 HFf, cross_sectional_area, Ma=0.0):

        # panel level properties i.e. only get assigned/changed at this level
        self.panel_hexagons = panel_hexagons  # needs to be accessed at hexagon
        # level
        self.HFf = HFf  # needs to be accessed at hexagon level
        self.cross_sectional_area = cross_sectional_area  # needs to be
        # accessed at hexagon level
        self.mass_flow_total: float = None  # calculated in the first channel
        self.jet_cross_sectional_area: float = None  # needs to be accessed at
        # channel level, calculated in the first channel
        self.jet_diameter: float = None  # needs to be accessed at channel
        # level, calculated in the first channel

        # channel level properties i.e. only get assigned/changed at this level
        self.input_pressure = input_pressure  # needs to be accessed at hexagon
        # level
        self.input_temperature = input_temperature  # needs to be accessed at
        # hexagon level
        self.mach_number = Ma  # needs to be accessed at channel level
        self.reynolds_number: float = None  # needs to be accessed at channel
        #  level
        self.jet_velocity: float = None  # needs to be accessed at channel
        # level
        self.euler_number: float = None  # needs to be accessed at channel
        # level
        self.pressure_drop: float = None  # needs to be accessed at channel
        # level
        self.non_dimensional_temperature: float = None  # needs to be accessed
        # at channel level
        self.non_dimensional_mass_flow: float = None  # needs to be accessed at
        # channel level
        self.mass_flow: float = None  # needs to be accessed at hexagon level
        self.fluid_properties = Material()
        self.fluid_properties.helium(input_pressure, input_temperature)
        # needs to be accessed at hexagon level
        self.hexagons: list = []  # needs to be accessed at hexagon level
        self.output_temperature: float = None  # needs to be accessed at
        # channel level
        self.output_pressure = input_pressure  # needs to be accessed at
        # channel level
        self.max_temperature: float = None  # needs to be accessed at channel
        # level
        self.break_flag: float = 0

    def get_hexagons(self, hexagon_selection):

        for hexagon in hexagon_selection:

            self.hexagons.append(
                Hexagon(hexagon, self.HFf, self.cross_sectional_area,
                        self.fluid_properties.specific_heat_capacity,
                        self.input_temperature, self.mass_flow,
                        self.non_dimensional_temperature))

            self.panel_hexagons.remove(hexagon)

    @classmethod
    def create_first_channel(cls, input_pressure, input_temperature,
                             panel_hexagons, HFf, cross_sectional_area,
                             number_of_hexagons, distribution_type, layup_type,
                             Ma=0.0):

        obj = cls(input_pressure, input_temperature, panel_hexagons,
                  HFf, cross_sectional_area, Ma)

        obj.mach_number, obj.jet_velocity = \
            mach_number_sympy(obj.input_temperature, Ma=obj.mach_number)
        obj.mass_flow, obj.jet_velocity, obj.jet_cross_sectional_area, \
            obj.density, obj.jet_diameter = mass_flow_sympy(
                obj.input_temperature, obj.input_pressure, u=obj.jet_velocity)

        obj.mass_flow_total = number_of_hexagons * obj.mass_flow

        obj.repeated_operations()

        if distribution_type == "strikepoint":
            reordered_hexagons = sorted(
                obj.panel_hexagons, key=lambda x: x[3], reverse=False)
            hexagon_selection = reordered_hexagons[:number_of_hexagons]
        elif distribution_type == "inout" and layup_type == "HF_specific":
            reordered_hexagons = sorted(
                obj.panel_hexagons, key=lambda x: x[3], reverse=False)
            hexagon_selection = reordered_hexagons[:number_of_hexagons]
        elif distribution_type == "inout":
            reordered_hexagons = sorted(
                obj.panel_hexagons, key=lambda x: x[5], reverse=False)
            hexagon_selection = reordered_hexagons[:number_of_hexagons]

        obj.get_hexagons(hexagon_selection)

        obj.output_temperature = \
            sum([hexagon.developed_coolant_temperature for hexagon in
                 obj.hexagons]) / number_of_hexagons

        obj.output_pressure -= obj.pressure_drop
        obj.max_temperature = \
            max([hexagon.metal_temperature for hexagon in obj.hexagons])

        return obj

    @classmethod
    def create_HF_specific_channel(cls, panel_hexagons, HFf,
                                   cross_sectional_area, remaining_hexagons,
                                   previous_channel, mass_flow_actual, Ma=0.0):

        obj = cls(previous_channel.output_pressure,
                  previous_channel.output_temperature, panel_hexagons,
                  HFf, cross_sectional_area)

        obj.mass_flow = previous_channel.mass_flow
        obj.max_temperature_from_previous_channel = \
            previous_channel.max_temperature

        obj.mass_flow, obj.number_of_hexagons_required = massflow_nextrows(
            remaining_hexagons[0], obj.mass_flow,
            obj.fluid_properties.specific_heat_capacity,
            obj.cross_sectional_area,
            obj.max_temperature_from_previous_channel,
            obj.input_temperature, obj.HFf, mass_flow_actual)

        remainging_hexagons, obj.break_flag = \
            obj.repeated_operations(1, remaining_hexagons, mass_flow_actual)

        obj.max_temperature = previous_channel.max_temperature

        return obj

    @classmethod
    def create_structured_channel(cls, panel_hexagons, HFf,
                                  cross_sectional_area, remaining_hexagons,
                                  previous_channel, mass_flow_actual,
                                  specified_number_of_hexagons=None, Ma=0.0):

        obj = cls(previous_channel.output_pressure,
                  previous_channel.output_temperature, panel_hexagons, HFf,
                  cross_sectional_area)

        obj.mass_flow = previous_channel.mass_flow
        obj.max_temperature_from_previous_channel = \
            previous_channel.max_temperature

        if specified_number_of_hexagons is None:
            obj.number_of_hexagons_required = len(previous_channel.hexagons)
        else:
            obj.number_of_hexagons_required = specified_number_of_hexagons

        remainging_hexagons, obj.break_flag = \
            obj.repeated_operations(1, remaining_hexagons, mass_flow_actual)

        obj.max_temperature = \
            max([hexagon.metal_temperature for hexagon in obj.hexagons])

        return obj

    def repeated_operations(self, flag=None, remaining_hexagons=None,
                            mass_flow_actual=None):

        if flag == 1:
            if self.mass_flow == 0.0:
                self.number_of_hexagons_required = len(remaining_hexagons)
                self.mass_flow = \
                    mass_flow_actual / self.number_of_hexagons_required

            elif self.number_of_hexagons_required <= len(remaining_hexagons):
                self.number_of_hexagons_required = \
                    int(self.number_of_hexagons_required)
                self.mass_flow = \
                    mass_flow_actual / self.number_of_hexagons_required

            else:
                self.number_of_hexagons_required = len(remaining_hexagons)
                self.mass_flow = \
                    mass_flow_actual / self.number_of_hexagons_required

            self.mass_flow, self.jet_velocity, self.jet_cross_sectional_area, \
                self.fluid_properties.density, self.jet_diameter = \
                mass_flow_sympy(self.input_temperature,
                                self.input_pressure,
                                massflow=self.mass_flow)

            self.mach_number, self.jet_velocity = \
                mach_number_sympy(self.input_temperature, u=self.jet_velocity)

        self.non_dimensional_mass_flow, self.mass_flow, \
            self.fluid_properties.specific_heat_capacity = \
            mstar_sympy(self.input_temperature, self.input_pressure,
                        self.cross_sectional_area, massflow=self.mass_flow)
        self.reynolds_number, self.jet_velocity = \
            Reynolds_sympy(self.fluid_properties.density,
                           self.fluid_properties.viscosity, self.jet_diameter,
                           u=self.jet_velocity)
        self.euler_number, self.reynolds_number = \
            Euler_sympy("DLH", Re=self.reynolds_number)
        self.pressure_drop, self.euler_number, self.jet_velocity = \
            pdrop_sympy(self.fluid_properties.density, Eu=self.euler_number,
                        u=self.jet_velocity)
        self.non_dimensional_temperature, self.non_dimensional_mass_flow = \
            taus_sympy("DLH", mstar=self.non_dimensional_mass_flow)

        if flag == 1:
            self.get_hexagons(
                remaining_hexagons[:self.number_of_hexagons_required])
            del remaining_hexagons[:self.number_of_hexagons_required]

            self.output_temperature = \
                sum([hexagon.developed_coolant_temperature for hexagon in
                     self.hexagons]) / self.number_of_hexagons_required
            self.output_pressure = self.input_pressure - self.pressure_drop

            if self.output_pressure <= 0.0:
                self.break_flag = 1.0

            return remaining_hexagons, self.break_flag


class Panel:

    def __init__(self, input_file, heat_flux_file, mach_number,
                 number_of_hexagons_in_first_channel, mass_split, layup_type,
                 distribution_type):

        self.channels = []
        self.panel_hexagon = []
        self.mass_flow_total: float = None
        self.read_inputs(input_file, heat_flux_file)
        self.mach_number = mach_number
        self.number_of_hexagons_in_first_channel = \
            number_of_hexagons_in_first_channel
        self.mass_flow_split = mass_split
        self.layup_type = layup_type
        self.distribution_type = distribution_type

    def __call__(self, cell_type):

        self.cell_type = cell_type

        if self.cell_type == "DLH":       
            self.populate_cells()
            self.configure_hexagons()

        elif self.cell_type == "JIVC":       
            self.populate_cells()
            # self.configure_squares()

        else:
            raise ValueError("Error, no cell type assigned")

    def configure_squares(self):
        pass

    def configure_hexagons(self):

        # adds hexagons to middle channel and calcs flow properties
        self.channels.append(Channel.create_first_channel(
            self.input_data["input_pressure"],
            self.input_data["input_temperature"], self.panel_hexagons,
            self.HFf, self.cross_sectional_area,
            self.number_of_hexagons_in_first_channel, self.distribution_type,
            self.layup_type, Ma=self.mach_number))

        self.mass_flow_total = self.channels[0].mass_flow_total

        implementation_options = {"strikepoint": self.strikepoint_injection,
                                  "inout": self.inout_injection}
        implementation_choice = self.distribution_type

        run = implementation_options[implementation_choice]
        run(self.channels[0])

    def strikepoint_injection(self, first_channel):

        new_channels = [[], []]

        reordered_hexagons = sorted(
            self.panel_hexagons, key=lambda x: x[3], reverse=False)
        inboard_hexagons, outboard_hexagons = \
            [hexagons for hexagons in reordered_hexagons if hexagons[4]
             <= 0.0], [hexagons for hexagons in reordered_hexagons if
                       hexagons[4] >= 0.0]

        all_remaining_hexagons = [inboard_hexagons, outboard_hexagons]

        self.mass_flow_actual = [self.mass_flow_total *
                                 (1-self.mass_flow_split),
                                 self.mass_flow_total *
                                 self.mass_flow_split]

        self.add_channels(all_remaining_hexagons, new_channels, first_channel)

        for channel in new_channels[0]:
            self.channels.insert(0, channel)

        for channel in new_channels[1]:
            self.channels.append(channel)

    def inout_injection(self, first_channel):

        new_channels = [[]]

        if self.layup_type == "structured":

            reordered_hexagons = sorted(
                self.panel_hexagons, key=lambda x: x[5], reverse=False)

            all_remaining_hexagons = [reordered_hexagons]

            self.mass_flow_actual = [self.mass_flow_total]

            self.add_channels(all_remaining_hexagons, new_channels,
                              first_channel)

            for direction in new_channels:
                for channel in direction:
                    self.channels.append(channel)

        elif self.layup_type == "HF_specific":

            first_round = [[first_channel]]
            self.new_channels = [[]]

            reordered_hexagons = sorted(
                self.panel_hexagons, key=lambda x: x[3], reverse=False)
            inboard_hexagons, outboard_hexagons = \
                [hexagons for hexagons in reordered_hexagons if hexagons[4]
                 <= 0.0], [hexagons for hexagons in reordered_hexagons if
                           hexagons[4] >= 0.0]

            self.dummy_set_all_remaining_hexagons = \
                [deepcopy(inboard_hexagons)]
            all_remaining_hexagons = [inboard_hexagons, outboard_hexagons]

            self.mass_flow_actual = [self.mass_flow_total]

            self.add_channels(all_remaining_hexagons, first_round,
                              first_channel)

            self.channels = []

            for direction in self.new_channels:
                for channel in direction:
                    self.channels.append(channel)

    def add_channels(self, all_remaining_hexagons, new_channels,
                     first_channel):

        layup_options = {"structured": Channel.create_structured_channel,
                         "HF_specific": Channel.create_HF_specific_channel}

        layup_selection = layup_options[self.layup_type]

        if self.layup_type == "structured" and \
            self.distribution_type == "strikepoint" \
            or self.layup_type == "structured" and \
            self.distribution_type == "inout" \
                or self.layup_type == "HF_specific" and \
                self.distribution_type == "strikepoint":

            for index, direction in enumerate(all_remaining_hexagons):

                mass_flow_actual = self.mass_flow_actual[index]

                new_channels[index].append(
                    layup_selection(self.panel_hexagons, self.HFf,
                                    self.cross_sectional_area, direction,
                                    first_channel, mass_flow_actual))

                while len(direction) >= 1.0:
                    new_channels[index].append(
                        layup_selection(self.panel_hexagons, self.HFf,
                                        self.cross_sectional_area, direction,
                                        new_channels[index][-1],
                                        mass_flow_actual))
                    if new_channels[index][-1].break_flag == 1:
                        break

        else:

            # HF_specific from middle out to get channel sizes
            dummy_hexagons = deepcopy(self.panel_hexagons)

            for direction in self.dummy_set_all_remaining_hexagons:

                mass_flow_actual = self.mass_flow_actual[0]

                new_channels[0].append(
                    Channel.create_HF_specific_channel(
                        dummy_hexagons, self.HFf, self.cross_sectional_area,
                        direction, first_channel, mass_flow_actual))

                while len(direction) >= 1.0:

                    new_channels[0].append(
                        Channel.create_HF_specific_channel(
                            dummy_hexagons, self.HFf,
                            self.cross_sectional_area, direction,
                            new_channels[0][-1], mass_flow_actual))
                    if new_channels[0][-1].break_flag == 1:
                        break

            new_channels[0].reverse()
            self.populate_hexagons()

            self.layup_type = "structured"

            self.number_of_hexagons_in_first_channel = \
                len(new_channels[0][0].hexagons)
            self.mach_number = new_channels[0][0].mach_number

            self.new_channels[0].append(
                Channel.create_first_channel(
                    self.input_data["input_pressure"],
                    self.input_data["input_temperature"], self.panel_hexagons,
                    self.HFf, self.cross_sectional_area,
                    self.number_of_hexagons_in_first_channel,
                    self.distribution_type, self.layup_type,
                    Ma=self.mach_number))

            self.mass_flow_total = self.new_channels[0][0].mass_flow_total

            reordered_hexagons = sorted(
                self.panel_hexagons, key=lambda x: x[5], reverse=False)
            all_remaining_hexagons = [reordered_hexagons]

            del new_channels[0][0]

            for i, direction in enumerate(all_remaining_hexagons):

                mass_flow_actual = self.mass_flow_total

                for index, channel in enumerate(new_channels[0]):

                    self.new_channels[0].append(
                        Channel.create_structured_channel(
                            self.panel_hexagons, self.HFf,
                            self.cross_sectional_area, direction,
                            self.new_channels[0][-1], mass_flow_actual,
                            specified_number_of_hexagons=
                            len(channel.hexagons)))

                    if self.new_channels[0][-1].break_flag == 1:
                        break

                while len(direction) >= 1.0:

                    self.new_channels[0].append(
                        Channel.create_HF_specific_channel(
                            self.panel_hexagons, self.HFf,
                            self.cross_sectional_area, direction,
                            self.new_channels[0][-1], mass_flow_actual))

                    if self.new_channels[0][-1].break_flag == 1:
                        break

    def read_inputs(self, inputs_filename, q_filename):
        # read asc files and store results to member variables (self...)

        user_data = get_example_data_path(inputs_filename)
        heatflux_data = get_example_data_path(q_filename)

        input_data = dict()
        tree = ET.parse(user_data)
        trunk = tree.getroot()

        for root in trunk:
            try:
                input_data[root.tag] = float(root.text)
            except ValueError:
                input_data[root.tag] = list(root.text.split(" "))

        self.input_data = input_data

        q = []  # power input from VITA
        s = []  # distance from strike point

        with open(heatflux_data) as geometry:
            for line in geometry:
                heatflux = line.split()
                s.append(float(heatflux[0]) / 100.0)  # converting to m from mm, is this correct, might be cm
                q.append(float(heatflux[1])) # / 6.87) * 30.0)

        sbar = [i + self.input_data["strike_radius"] for i in s]
        heatflux = [i * 1e6 for i in q]
        for index, point in enumerate(heatflux):
            if point <= 30 * 1e6 * 0.2:
                heatflux[index] = 30 * 1e6 * 0.2

        sbar.insert(0, sbar[0])
        sbar.insert(0, 0)
        sbar.append(self.input_data["outer_radius"])
        self.sbar = sbar

        heatflux.insert(0, 30 * 1e6 * 0.2)
        heatflux.insert(0, 30 * 1e6 * 0.2)
        heatflux.append(30 * 1e6 * 0.2)

        plt.plot(sbar, heatflux)
        plt.xlabel("distance, cm")
        plt.ylabel("Heat Flux, MW/m2")

        self.HFf = interp1d(sbar, heatflux)
    
    def hexagon_starter(self, cp, point_x1, line_angle, topy):
        
        width = self.input_data["width"]
        height = width / 2 * tan(60 * pi / 180)
        
        no_points_y = int(round(topy/height, 0))
        no_points_x = int(round(((self.input_data["outer_radius"]-point_x1) /
                                 (1.5*width)), 0))

        a = width/2
        self.cross_sectional_area = (3 * sqrt(3) / 2) * a ** 2

        for i in range(no_points_x):
            for j in range(no_points_y):
                cp.append([i*1.5*width + point_x1, j*height])
                cp.append([i*1.5*width + 0.75*width + point_x1, 
                           j*height + 0.5*height])
        
        return cp, height, width

    def squares_starter(self, cp, point_x1, line_angle, topy):
        
        height = self.input_data["height"]
        width = self.input_data["width"]
        self.cross_sectional_area = height * width
        
        no_points_y = int(round(topy/height, 0))
        no_points_x = int(round((self.input_data["outer_radius"]-point_x1) /
                                 (width), 0))

        for i in range(no_points_x):
            for j in range(no_points_y):
                cp.append([i*width + point_x1, j*height])
                
        return cp, width, height
    
    def bounding_box(self, first_radius, second_radius, cp_final_y, rotation_theta_set, split_theta):

        cp = []  # a list for assigning centre points of hexagon array, (x, y)
        cp_remove = []  # list of values to remove from cp as are out of bounds
        cp_final = []  # list of remaining central points
        cp_ID = []  # a list of the matrix ID's left after removing out of
        # bounds points
        cp_final_y_mirror = []        
        cp_final_new = []

        point_x1 = first_radius*cos(split_theta)
        point_y1 = first_radius*sin(split_theta)
        line_angle = tan(split_theta)
        topy = second_radius * sin(split_theta)

        cell_type_options = {"DLH":self.hexagon_starter, \
                             "JIVC":self.squares_starter}

        cell_type_selection = cell_type_options[self.cell_type]
        cp, width, height = cell_type_selection(cp, point_x1, line_angle, topy)

        x, y, r = symbols("x y r")
        topy = (line_angle * (x-point_x1)) + point_y1
        arc = ((r**2 - ((y-self.input_data["y_displacement_from_origin"])**2))
                ** 0.5) + self.input_data["x_displacement_from_origin"]

        for i, coordinates in enumerate(cp):

            if topy.subs([(x, coordinates[0])]) <= coordinates[1]:
                cp_remove.append(coordinates)

            # inner arc
            try:
                float(arc.subs([(r, first_radius),
                                (y, coordinates[1])]))
            except:
                None
            else:
                if arc.subs([(r, first_radius),
                              (y, coordinates[1])]) >= coordinates[0]:
                    cp_remove.append(coordinates)

            # outer arc
            try:
                float(arc.subs([(r, second_radius),
                                (y, coordinates[1])]))
            except:
                None
            else:
                if arc.subs([(r, second_radius),
                              (y, coordinates[1])]) <= coordinates[0]:
                    cp_remove.append(coordinates)

        for i in cp:
            if i not in cp_remove:
                cp_final.append(i)

        for i in cp_final:
            cp_ID.append(
                [i[0] / width - cp_final[0][0] / width, 0, i[1] / height])
            cp_final_new.append([i[0], 0, i[1]])

        for i in cp_final_new:
            if i[2] == 0:
                pass
            else:
                cp_final_y_mirror.append([i[0], 0, -i[2]])

        # - need to copy and rotate the copies to 

        for i in cp_final_y_mirror:
            cp_final_new.append(i)

        for i, value in enumerate(cp_final_new):
            cp_final_new[i].append(
                sqrt(value[0]**2 + value[2]**2))
            # absolute radius
            cp_final_new[i].append(
                atan(value[2]/value[0]))
            # angle from the origin for each point

        cpsets = []

        for i in rotation_theta_set:
            cpsets.append(array(deepcopy(cp_final_new)))

        cp_final_new = []
        
        for i, cpset in enumerate(cpsets):
            for x, y, z, R, theta in cpset:
                cp_final_new.append([R * cos(rotation_theta_set[i] + theta), \
                                    y, \
                                    R * sin(rotation_theta_set[i] + theta)])

        for i, value in enumerate(cp_final_new):
            cp_final_new[i].append(
                sqrt((sqrt(value[0]**2 + value[2]**2) - self.sbar[2])**2))
            # distance from strikepoint, absolute
            cp_final_new[i].append(
                sqrt(value[0]**2 + value[2]**2) - self.sbar[2])
            # distance from strikepoint, direction dependant
            cp_final_new[i].append(
                sqrt(value[0]**2 + value[2]**2))
            # absolute radius
            cp_final_new[i].append(
                atan(value[2]/value[0]))
            # angle from the origin for each point

        for i in cp_final_new:
            cp_final_y.append(deepcopy(i))

        return cp_final_y

    def populate_cells(self):

        # - centre point definitions, same for both codes

        cp_final_y = []
        rotation_set = []
        theta_set = []

        n = [int(ni) for ni in self.input_data["n"]]
        m = int(self.input_data["m"])
        
        if len(n) != m:
            raise Exception(""" the number of radial slices is not equal to the number of provided lateral slices """)

        for ni in n:

            single_n_angle = self.input_data["number_of_plates"] * ni * \
                self.input_data["number_of_carriers"]

            theta_set.append((360.0 / single_n_angle / 2 - self.input_data["swept_angle_gap_between_plates"]) * pi / 180)
            rotation_set.append(theta_set[-1] * ni * 2 * linspace((ni-1)*(1/(2*ni)), - (ni-1)*(1/(2*ni)), ni))

        rset = linspace(self.input_data["inner_radius"], \
                        self.input_data["outer_radius"], \
                        m+1)

        for i, r in enumerate(rset):
            if len(rset) - 1 == i:
                break
            else:
                cp_final_y = self.bounding_box(r, rset[i+1], cp_final_y, rotation_set[i], theta_set[i])

        self.panel_hexagons = cp_final_y

if __name__ == "__main__":
    start_time = time.time()

    # BUSINESS CODE

    results = dict()
    number_of_hexagons_range = linspace(1000, 4000, 1)  # 4
    mach_number_range = linspace(0.2, 0.3, 1)  # 4
    mass_split_range = linspace(0.5, 0.9, 1)  # 2
    distribution_options = ["strikepoint"]  # strikepoint inout
    layup_options = ["HF_specific"]  # structured HF_specific

    for layup in layup_options:
        for distribution in distribution_options:
            for number_of_hexagons in number_of_hexagons_range:
                for mach_number in mach_number_range:
                    for mass_split in mass_split_range:

                        key = "panel" + "_" + \
                            str(number_of_hexagons) + "_" + \
                            str(mach_number) + "_" + \
                            str(mass_split) + "_" + \
                            str(layup) + "_" + \
                            str(distribution)

                        panel = Panel("inputs.xml", "q_adjusted.asc",
                                      float(mach_number),
                                      int(number_of_hexagons),
                                      float(mass_split), layup, distribution)
                        panel("JIVC")

                        results[key] = panel

    keys = results.keys()

    # f = open("hexagon_centre_points.asc", "w")
    # for key in keys:
    #     for channel in results[key].channels:
    #         for hexagon in channel.hexagons:
    #             f.write(str(hexagon.x) + "\t" + "0" + "\t" + str(hexagon.y))
    #             f.write("\n")
    # f.close()

    f = open("hexagon_centre_points.asc", "w")
    for key in keys:
        for hexagon in results[key].panel_hexagons:
            f.write(str(hexagon[0]) + "\t" + "0" + "\t" + str(hexagon[2]))
            f.write("\n")
    f.close()

    f = open("hexagon_data.asc", "w")
    for key in keys:
        for channel in results[key].channels:
            for hexagon in channel.hexagons:
                f.write(str(hexagon.x) + "\t" + str(hexagon.y) + "\t" +
                        str(hexagon.metal_temperature) + "\t" +
                        str(hexagon.developed_coolant_temperature))
                f.write("\n")
        f.write(key)
        f.write("\n")
    f.close()

    f = open("channel_data.asc", "w")
    for key in keys:
        for index, channel in enumerate(results[key].channels):
            f.write(str(key) + "\t" + \
                    # str(channel.channel_ID) + "\t" + \
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
                    str(channel.fluid_properties.specific_heat_capacity) + \
                    "\t" + str(channel.break_flag))
            f.write("\n")
    f.close()
