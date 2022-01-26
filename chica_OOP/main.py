from scipy.interpolate import interp1d
from math import sqrt, sin, cos, tan, pi, atan
from CHICA.utility import get_example_data_path
from CoolProp.CoolProp import PropsSI as SI
from numpy import linspace, array, average, flip, delete, append
from numpy import concatenate as concat
from pandas import DataFrame
from copy import deepcopy
from os import mkdir
import xml.etree.ElementTree as ET
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
                                        massflow_nextrows_DLH, \
                                        mass_flow_next_row_JIVC

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
        self.R = sqrt(x**2 + y**2)
        heatflux = HFf(self.R) * cross_sectional_area

        self.metal_temperature = metal_temperature_sympy(input_temperature,
                                    non_dimensional_temperature, mass_flow,
                                    heatflux, specific_heat_capacity)

        self.developed_coolant_temperature = \
            coolant_temperature_sympy(non_dimensional_temperature,
                                      input_temperature, self.metal_temperature)

        self.VJ_cell_ID = self.identification()

        # store variables, heavy handed as on channel created yet, will be deleted, just for checking
        self.taus = non_dimensional_temperature
        self.heatflux = heatflux
        self.input_temperature = input_temperature
        self.mass_flow = mass_flow
        self.Cp = specific_heat_capacity

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
        
        self.metal_temperature = metal_temperature_sympy(input_temperature,
                                    non_dimensional_temperature, mass_flow,
                                    heatflux, specific_heat_capacity)

        self.developed_coolant_temperature = \
            coolant_temperature_sympy(non_dimensional_temperature,
                                      input_temperature, self.metal_temperature)

        self.hexagon_ID = self.identification()

        self.taus = non_dimensional_temperature
        self.heatflux = heatflux
        self.input_temperature = input_temperature
        self.mass_flow = mass_flow
        self.Cp = specific_heat_capacity

    @classmethod
    def identification(cls):
        cls.hexagon_ID += 1
        return cls.hexagon_ID

class Channel:

    def __init__(self, input_pressure, input_temperature, HFf, \
                 cross_sectional_area, Ma=0.0, panel_hexagons = None):

        # # panel level properties i.e. only get assigned/changed at this level
        self.panel_hexagons = panel_hexagons  # needs to be accessed at hexagon
        # # level
        self.HFf = HFf  # needs to be accessed at hexagon level
        self.cross_sectional_area = cross_sectional_area  # needs to be
        # # accessed at hexagon level
        self.mass_flow_total: float = None  # calculated in the first channel
        self.jet_cross_sectional_area: float = None  # needs to be accessed at
        # # channel level, calculated in the first channel
        self.jet_diameter: float = None  # needs to be accessed at channel
        # # level, calculated in the first channel

        # # channel level properties i.e. only get assigned/changed at this level
        self.input_pressure = input_pressure  # needs to be accessed at hexagon
        # # level
        self.input_temperature = input_temperature  # needs to be accessed at
        # # hexagon level
        self.mach_number = Ma  # needs to be accessed at channel level
        self.reynolds_number: float = None  # needs to be accessed at channel
        # #  level
        self.jet_velocity: float = None  # needs to be accessed at channel
        # # level
        self.euler_number: float = None  # needs to be accessed at channel
        # # level
        self.pressure_drop: float = None  # needs to be accessed at channel
        # # level
        self.non_dimensional_temperature: float = None  # needs to be accessed
        # # at channel level
        self.non_dimensional_mass_flow: float = None  # needs to be accessed at
        # # channel level
        self.mass_flow: float = None  # needs to be accessed at hexagon level
        self.fluid_properties = Material()
        self.fluid_properties.helium(input_pressure, input_temperature)
        # # needs to be accessed at hexagon level
        self.hexagons: list = []  # needs to be accessed at hexagon level
        self.output_temperature: float = None  # needs to be accessed at
        # # channel level
        self.output_pressure = input_pressure  # needs to be accessed at
        # # channel level
        self.max_temperature: float = None  # needs to be accessed at channel
        # # level
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

        obj = cls(input_pressure, input_temperature, HFf, cross_sectional_area, \
                  Ma = Ma, panel_hexagons = panel_hexagons)
        
        obj.mach_number, obj.jet_velocity = \
            mach_number_sympy(obj.input_temperature, Ma=obj.mach_number)
        
        obj.mass_flow, obj.jet_velocity, obj.jet_cross_sectional_area, \
            obj.density, obj.jet_diameter = mass_flow_sympy(
                obj.input_temperature, obj.input_pressure, u=obj.jet_velocity)
        
        obj.mass_flow_total = number_of_hexagons * obj.mass_flow

        obj.repeated_operations_DLH()

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
                  previous_channel.output_temperature, HFf,
                  cross_sectional_area, panel_hexagons = panel_hexagons)
        
        obj.mass_flow = previous_channel.mass_flow
        obj.max_temperature_from_previous_channel = \
            previous_channel.max_temperature

        obj.mass_flow, obj.number_of_hexagons_required = massflow_nextrows_DLH(
            remaining_hexagons[0], obj.mass_flow,
            obj.fluid_properties.specific_heat_capacity,
            obj.cross_sectional_area,
            obj.max_temperature_from_previous_channel,
            obj.input_temperature, obj.HFf, mass_flow_actual)

        remainging_hexagons, obj.break_flag = \
            obj.repeated_operations_DLH(1, remaining_hexagons, mass_flow_actual)

        obj.max_temperature = previous_channel.max_temperature

        return obj

    @classmethod
    def create_structured_channel(cls, panel_hexagons, HFf,
                                  cross_sectional_area, remaining_hexagons,
                                  previous_channel, mass_flow_actual,
                                  specified_number_of_hexagons=None, Ma=0.0):
        
        obj = cls(previous_channel.output_pressure,
                  previous_channel.output_temperature, HFf,
                  cross_sectional_area, panel_hexagons = panel_hexagons)
        
        obj.mass_flow = previous_channel.mass_flow
        obj.max_temperature_from_previous_channel = \
            previous_channel.max_temperature

        if specified_number_of_hexagons is None:
            obj.number_of_hexagons_required = len(previous_channel.hexagons)
        else:
            obj.number_of_hexagons_required = specified_number_of_hexagons
        
        remainging_hexagons, obj.break_flag = \
            obj.repeated_operations_DLH(1, remaining_hexagons, mass_flow_actual)

        obj.max_temperature = \
            max([hexagon.metal_temperature for hexagon in obj.hexagons])

        return obj

    def repeated_operations_JIVC(self, assessment_set, n, peak_metal_temperature, \
        input_pressure, input_temperature, mass_flow_total, HFf, cross_sectional_area):
        
        cell_temperatures = []
        VJ_cells = []
        
        fluid_properties = Material()
        fluid_properties.helium(input_pressure, input_temperature)
        
        assessment_set_total_mass_flow = mass_flow_total / n
        jet_mass_flow = assessment_set_total_mass_flow / len(assessment_set)
        
        jet_mass_flow, jet_velocity, Ajet, density, jet_diameter \
            = mass_flow_sympy(input_temperature, input_pressure, massflow = jet_mass_flow)
        
        Reynolds_number, jet_velocity = Reynolds_sympy(density, fluid_properties.viscosity, jet_diameter, \
                       u = jet_velocity)
        
        Euler_number, Reynolds_number = Euler_sympy("JIVC", Re = Reynolds_number)
        
        pressure_drop, Euler_number, jet_velocity = \
            pdrop_sympy(density, Eu = Euler_number, u = jet_velocity)
        
        mstar, jet_mass_flow, specific_heat_capacity = mstar_sympy(\
            input_temperature, input_pressure, cross_sectional_area, massflow = jet_mass_flow)
        
        taus, mstar = taus_sympy("JIVC", mstar = mstar)
        
        mach_number, jet_velocity = mach_number_sympy(input_temperature, \
            u = jet_velocity)
        
        for cell in assessment_set:
            VJ_cell_i = VJ_cell(cell[0], cell[2], HFf, cross_sectional_area, 
                fluid_properties.specific_heat_capacity, \
                input_temperature, jet_mass_flow, taus)
            VJ_cell_i.jet_mass_flow = jet_mass_flow
            VJ_cells.append(VJ_cell_i)
            cell_temperatures.append(VJ_cell_i.developed_coolant_temperature)
            if peak_metal_temperature <= VJ_cell_i.metal_temperature:
                peak_metal_temperature = VJ_cell_i.metal_temperature
        
        new_input_temperature = average(cell_temperatures)
        
        return pressure_drop, new_input_temperature, mach_number, peak_metal_temperature, VJ_cells, fluid_properties, jet_velocity

    def repeated_operations_DLH(self, flag=None, remaining_hexagons=None,
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
                self.mass_flow = mass_flow_actual / self.number_of_hexagons_required

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

    def __init__(self, mach_number, input_file, heat_flux_file,
                 number_of_hexagons_in_first_channel = None, mass_split = None, 
                 layup_type = None, distribution_type = None, n = None, m1 = None, \
                 m2 = None, mass_flow_total = None, centre_channel_half_width = None):

        self.channels = []
        self.panel_hexagon = []
        self.read_inputs(input_file, heat_flux_file)
        self.mach_number = mach_number
        self.number_of_hexagons_in_first_channel = \
            number_of_hexagons_in_first_channel
        self.mass_flow_split = mass_split
        self.layup_type = layup_type
        self.distribution_type = distribution_type
        self.n = n
        self.m1 = m1
        self.m2 = m2
        self.mass_flow_total = mass_flow_total
        self.half_width = centre_channel_half_width

    def __call__(self, cell_type):

        self.cell_type = cell_type

        if self.cell_type == "DLH":       
            self.populate_cells()
            self.configure_DLH()

        elif self.cell_type == "JIVC":       
            self.populate_cells()
            self.configure_JIVC()

        else:
            raise ValueError("Error, no cell type assigned")

    def configure_DLH(self):
        
        # adds hexagons to middle channel and calcs flow properties
        self.channels.append(Channel.create_first_channel(
            self.input_data["input_pressure"],
            self.input_data["input_temperature"], self.panel_hexagons,
            self.HFf, self.cross_sectional_area,
            self.number_of_hexagons_in_first_channel, self.distribution_type,
            self.layup_type, Ma=self.mach_number))

        self.mass_flow_total = self.channels[0].mass_flow_total

        implementation_options = {"strikepoint": self.strikepoint_DLH,
                                  "inout": self.inout_DLH}
        implementation_choice = self.distribution_type

        run = implementation_options[implementation_choice]
        run(self.channels[0])

        output_pressures = []
        output_metal_temperatures = []
        output_coolant_temperatures = []
        number_of_channels = 0
        break_flag = []
        for channel in self.channels:
            output_pressures.append(channel.output_pressure)
            number_of_channels += 1
            output_metal_temperatures.append(channel.max_temperature)
            output_coolant_temperatures.append(channel.output_temperature)
            break_flag.append(channel.break_flag)
        
        self.output_pressure = min(output_pressures)
        self.number_of_channels = number_of_channels
        self.peak_metal_temperature = max(output_metal_temperatures)
        self.peak_coolant_temperature = max(output_coolant_temperatures)
        self.break_flag = max(break_flag)

    def configure_JIVC(self):
        
        implementation_options = {"strikepoint_structured": self.strikepoint_structured_JIVC,
                                  "strikepoint_HF_specific": self.strikepoint_HF_specific_JIVC,
                                  "inout_structured": self.inout_structured_JIVC, 
                                  "inout_HF_specific": self.inout_HF_specific_JIVC}
        
        implementation_choice = self.distribution_type + "_" + self.layup_type

        run = implementation_options[implementation_choice]
        run()

    def strikepoint_DLH(self, first_channel):

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

    def inout_DLH(self, first_channel):

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

    def inout_HF_specific_JIVC(self):
        pass

    def inout_structured_JIVC(self):
        
        mach_numbers = []
        mass_flow_totals = []
        mach_number_iter = 0
        number_of_cells = []
        error = lambda Ma_iter, Ma_target: sqrt((Ma_iter - Ma_target)**2)
        peak_metal_temperature = 0
        break_flag = 0
        output_temperatures = []
        output_pressures = []
        manifold_pipe_diameters = []
        delta_break_flag = 0

        while mach_number_iter != self.mach_number:
            
            for get_plotting_data in range(2):
                input_pressure = self.input_data["input_pressure"]
                input_temperature = self.input_data["input_temperature"]

                for i in range(self.m1+1):
                    pressure_drop, input_temperature, mach_number_iter, peak_metal_temperature, VJ_cells, fluid_properties, jet_velocity = \
                        Channel.repeated_operations_JIVC(self, self.assessment_set[i], self.n[i], \
                        peak_metal_temperature, input_pressure, input_temperature, self.mass_flow_total, self.HFf, self.cross_sectional_area)
                    input_pressure += - pressure_drop
                    if input_pressure <= 0:
                        break_flag = 1
                        break

                mach_numbers.append(mach_number_iter)
                mass_flow_totals.append(self.mass_flow_total)
                
                if mach_number_iter == self.mach_number:
                    break
                elif error(mach_number_iter, self.mach_number) <= 1E-3:
                    break
                elif break_flag == 1:
                    break
                elif mach_number_iter <= self.mach_number:
                    self.mass_flow_total += 2.0
                elif mach_number_iter >= self.mach_number:
                    self.mass_flow_total -= 0.5
            
            if error(mach_number_iter, self.mach_number) <= 1E-3:
                break
            elif break_flag == 1:
                break
            
            mach_mass = interp1d(mach_numbers, mass_flow_totals)
            try:
                self.mass_flow_total = float(mach_mass(self.mach_number))
            except ValueError:
                continue

        input_pressure = self.input_data["input_pressure"]
        input_temperature = self.input_data["input_temperature"]

        mach_numbers = []
        peak_metal_temperature = 0
        VJ_cells_panel = []
        
        for i in range(self.m1 + self.m2 + 1):
            
            pressure_drop, input_temperature, mach_number_iter, peak_metal_temperature, VJ_cells, fluid_properties, jet_velocity = \
                Channel.repeated_operations_JIVC(self, self.assessment_set[i], self.n[i], \
                peak_metal_temperature, input_pressure, input_temperature, self.mass_flow_total, self.HFf, self.cross_sectional_area)
            input_pressure += - pressure_drop
            
            manifold_pipe_diameters.append([sqrt(((mass_flow_total * self.mass_flow_split * 4) / \
                                             self.n[self.m1]) / \
                                            (fluid_properties.density * jet_velocity * pi)), i])
            
            output_pressures.append(input_pressure)
            output_temperatures.append(input_temperature)
            
            if input_pressure <= 0:
                break_flag = 1
                break
            
            mach_numbers.append(mach_number_iter)
            number_of_cells.append(len(self.assessment_set[i]) * self.n[i])
            VJ_cells_panel.append(VJ_cells)

        self.break_flag = break_flag
        
        manifold_pipe_diameters.sort(key = lambda x:x[1])
        manifold_pipe_diameters = array(next(zip(*manifold_pipe_diameters)))
        
        try:
            pipe_width_delta = self.microchannel_half_width - manifold_pipe_diameters
            for delta in pipe_width_delta:
                if delta <= 0.0:
                    delta_break_flag = 1.0
        except ValueError:
            pass
        
        if break_flag == 0:
            self.VJ_cells = VJ_cells_panel
            self.peak_metal_temperature = peak_metal_temperature
            self.peak_coolant_temperature = min(output_temperatures)
            self.number_of_cells = number_of_cells
            self.mach_numbers = mach_numbers
            self.number_of_channels = len(VJ_cells_panel)
            self.output_pressure = min(output_pressures)
            self.manifold_pipe_diameters = manifold_pipe_diameters
            self.delta_break_flag = delta_break_flag
        else:
            self.VJ_cells = VJ_cells_panel
            self.peak_metal_temperature = 0
            self.peak_coolant_temperature = 0
            self.number_of_cells = number_of_cells
            self.mach_numbers = mach_numbers
            self.number_of_channels = len(VJ_cells_panel)
            self.output_pressure = 0
            self.manifold_pipe_diameters = manifold_pipe_diameters
            self.delta_break_flag = delta_break_flag

    def strikepoint_HF_specific_JIVC(self):
        
        self.m1 = 1
        self.m2 = 1
        self.n = array([10 for i in range(m1 + m2 + 1)])
        
        number_of_cells = []
        mach_numbers = []        
        VJ_cells_panel = []
        output_pressures = []
        output_temperatures = []
        manifold_pipe_diameters = []
        delta_break_flag = 0
        n_jets_per_slab = array([])
        
        strikepoint_pressure = self.input_data["input_pressure"]
        input_temperature = self.input_data["input_temperature"]
        peak_metal_temperature = 0
        
        mach_number, jet_velocity = mach_number_sympy(input_temperature, Ma = self.mach_number)
        jet_mass_flow, jet_velocity, jet_cross_sectional_area, density, jet_diameter\
            = mass_flow_sympy(input_temperature, strikepoint_pressure, u = jet_velocity)
        mass_flow_total = len(self.assessment_set[self.m1]) * jet_mass_flow * self.n[self.m1]
        
        pressure_drop, strikepoint_temperature, mach_number_iter, peak_metal_temperature, VJ_cells, fluid_properties, jet_velocity =\
            Channel.repeated_operations_JIVC(self, self.assessment_set[self.m1], self.n[self.m1],\
            peak_metal_temperature, strikepoint_pressure, input_temperature, mass_flow_total,\
            self.HFf, self.cross_sectional_area)
        
        manifold_pipe_diameters.append([sqrt(((mass_flow_total* self.mass_flow_split * 4) / \
                                             self.n[self.m1]) / \
                                            (density * jet_velocity * pi)), self.m1])
        
        strikepoint_pressure += - pressure_drop
        mach_numbers.append(mach_number_iter)
        VJ_cells_panel.append(VJ_cells)
        assessment_set = [flip(self.assessment_set[0], 0), self.assessment_set[1], self.assessment_set[2]]
        
        remaining_set = deepcopy(assessment_set)
        loop_arrangement = [0, 2]

        for loop in loop_arrangement:
            
            input_pressure = deepcopy(strikepoint_pressure)
            input_temperature = deepcopy(strikepoint_temperature)
            
            unit_break_flag = 0
            break_flag = 0
            previous_jets = 0
            
            while len(remaining_set[loop]) > 0:
                
                fluid_properties = Material()
                fluid_properties.helium(input_pressure, input_temperature)
                
                jet_mass_flow, n_jets_total = mass_flow_next_row_JIVC(\
                remaining_set[loop][0], self.cross_sectional_area, \
                self.HFf, jet_mass_flow, fluid_properties.specific_heat_capacity, \
                input_temperature, peak_metal_temperature, mass_flow_total*self.mass_flow_split / self.n[loop], \
                self.cell_type)
                
                n_slabs_iter = round(n_jets_total / self.n_jets_per_slab[loop], 0)
                required_jets = int(n_slabs_iter * self.n_jets_per_slab[loop])

                if required_jets + previous_jets <= len(assessment_set[loop]):
                    next_group = assessment_set[loop][previous_jets:\
                        previous_jets + required_jets]
                else:
                    next_group = remaining_set[loop]
                    previous_jets = len(assessment_set)
                    unit_break_flag = 1
                
                pressure_drop, input_temperature, mach_number_iter, peak_metal_temperature, VJ_cells, fluid_properties, jet_velocity =\
                    Channel.repeated_operations_JIVC(self, next_group, self.n[loop], \
                    peak_metal_temperature, input_pressure, input_temperature, \
                    mass_flow_total * self.mass_flow_split, self.HFf, self.cross_sectional_area)

                manifold_pipe_diameters.append([sqrt(((mass_flow_total * self.mass_flow_split * 4) / \
                                             self.n[self.m1]) / \
                                            (fluid_properties.density * jet_velocity * pi)), loop])

                input_pressure += - pressure_drop
                mach_numbers.append(mach_number_iter)
                VJ_cells_panel.append(VJ_cells)
                
                if input_pressure <= 0:
                    break_flag = 1
                    break
                
                if unit_break_flag == 0:
                    n_jets_per_slab = append(n_jets_per_slab, n_slabs_iter)
                    previous_jets += required_jets
                    remaining_set[loop] = delete(remaining_set[loop], range(0,required_jets), 0)
                    
                else:
                    jet_mass_flow = 0.002
                    break

            self.mass_flow_split = 1 - self.mass_flow_split
            output_pressures.append(input_pressure)
            output_temperatures.append(input_temperature)
        
        self.break_flag = break_flag
        
        self.microchannel_half_width = n_jets_per_slab * self.input_data["height"]
        
        manifold_pipe_diameters.sort(key = lambda x:x[1])
        manifold_pipe_diameters = array(next(zip(*manifold_pipe_diameters)))
        
        try:
            pipe_width_delta = self.microchannel_half_width - manifold_pipe_diameters
            for delta in pipe_width_delta:
                if delta <= 0.0:
                    delta_break_flag = 1.0
        except ValueError:
            pass
        
        if break_flag == 0:
            self.mass_flow_total = mass_flow_total
            self.VJ_cells = VJ_cells_panel
            self.peak_metal_temperature = peak_metal_temperature
            self.peak_coolant_temperature = min(output_temperatures)
            self.number_of_cells = number_of_cells
            self.mach_numbers = mach_numbers
            self.number_of_channels = len(VJ_cells_panel)
            self.output_pressure = min(output_pressures)
            self.manifold_pipe_diameters = manifold_pipe_diameters
            self.delta_break_flag = delta_break_flag
        else:
            self.mass_flow_total = mass_flow_total
            self.VJ_cells = VJ_cells_panel
            self.peak_metal_temperature = 0
            self.peak_coolant_temperature = 0
            self.number_of_cells = number_of_cells
            self.mach_numbers = mach_numbers
            self.number_of_channels = len(VJ_cells_panel)
            self.output_pressure = 0
            self.manifold_pipe_diameters = manifold_pipe_diameters
            self.delta_break_flag = delta_break_flag

    def strikepoint_structured_JIVC(self):

        number_of_cells = []
        mach_numbers = []        
        VJ_cells_panel = []
        output_pressures = []
        output_temperatures = []
        break_flag = 0
        delta_break_flag = 0
        manifold_pipe_diameters = []
        
        strikepoint_pressure = self.input_data["input_pressure"]
        input_temperature = self.input_data["input_temperature"]
        peak_metal_temperature = 0
        
        loop_arrangement = [range(0, self.m1).__reversed__(), \
            range(self.m1 + 1, self.m1 + self.m2 + 1)]
        
        mach_number, jet_velocity = mach_number_sympy(input_temperature, Ma = self.mach_number)
        jet_mass_flow, jet_velocity, self.jet_cross_sectional_area, density, jet_diameter \
            = mass_flow_sympy(input_temperature, strikepoint_pressure, u = jet_velocity)
        mass_flow_total = len(self.assessment_set[self.m1]) * jet_mass_flow * self.n[self.m1]
        
        pressure_drop, strikepoint_temperature, mach_number_iter, peak_metal_temperature, VJ_cells, fluid_properties, jet_velocity =\
            Channel.repeated_operations_JIVC(self, self.assessment_set[self.m1], self.n[self.m1], \
            peak_metal_temperature, strikepoint_pressure, input_temperature, mass_flow_total, \
            self.HFf, self.cross_sectional_area)
        
        manifold_pipe_diameters.append([sqrt(((mass_flow_total* self.mass_flow_split * 4) / \
                                             self.n[self.m1]) / \
                                            (density * jet_velocity * pi)), self.m1])
        
        strikepoint_pressure += - pressure_drop
        input_pressure = deepcopy(strikepoint_pressure)
        input_temperature = deepcopy(strikepoint_temperature)
        mach_numbers.append(mach_number_iter)
        VJ_cells_panel.append(VJ_cells)
        
        for loop in loop_arrangement:
            
            for channel in loop:
                
                pressure_drop, input_temperature, mach_number_iter, peak_metal_temperature, VJ_cells, fluid_properties, jet_velocity =\
                    Channel.repeated_operations_JIVC(self, self.assessment_set[channel], self.n[channel], \
                    peak_metal_temperature, input_pressure, input_temperature, \
                    mass_flow_total * self.mass_flow_split, self.HFf, self.cross_sectional_area)

                manifold_pipe_diameters.append([sqrt(((mass_flow_total * self.mass_flow_split * 4) / \
                                             self.n[self.m1]) / \
                                            (fluid_properties.density * jet_velocity * pi)), channel])

                input_pressure += - pressure_drop
                mach_numbers.append(mach_number_iter)
                VJ_cells_panel.append(VJ_cells)

                if input_pressure <= 0:
                    break_flag = 1
                    break

            if input_pressure <= 0:
                break_flag = 1
                break

            self.mass_flow_split = 1 - self.mass_flow_split
            output_pressures.append(input_pressure)
            output_temperatures.append(input_temperature)
            input_pressure = deepcopy(strikepoint_pressure)
            input_temperature = deepcopy(strikepoint_temperature)
        
        self.break_flag = break_flag
        
        manifold_pipe_diameters.sort(key = lambda x:x[1])
        manifold_pipe_diameters = array(next(zip(*manifold_pipe_diameters)))
        
        try:
            pipe_width_delta = self.microchannel_half_width - manifold_pipe_diameters
            for delta in pipe_width_delta:
                if delta <= 0.0:
                    delta_break_flag = 1.0
        except ValueError:
            pass
        
        if break_flag == 0:
            self.mass_flow_total = mass_flow_total
            self.VJ_cells = VJ_cells_panel
            self.peak_metal_temperature = peak_metal_temperature
            self.peak_coolant_temperature = min(output_temperatures)
            self.number_of_cells = number_of_cells
            self.mach_numbers = mach_numbers
            self.number_of_channels = len(VJ_cells_panel)
            self.output_pressure = min(output_pressures)
            self.manifold_pipe_diameters = manifold_pipe_diameters
            self.delta_break_flag = delta_break_flag
        else:
            self.mass_flow_total = mass_flow_total
            self.VJ_cells = VJ_cells_panel
            self.peak_metal_temperature = 0
            self.peak_coolant_temperature = 0
            self.number_of_cells = number_of_cells
            self.mach_numbers = mach_numbers
            self.number_of_channels = len(VJ_cells_panel)
            self.output_pressure = 0
            self.manifold_pipe_diameters = manifold_pipe_diameters
            self.delta_break_flag = delta_break_flag
        
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
            self.populate_cells()

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
                s.append(float(heatflux[0]) / 100.0)  # converting to m from cm
                q.append(float(heatflux[1])) # / 6.87) * 30.0)

        sbar = [i + self.input_data["strike_radius"] for i in s]
        heatflux = [(i/q[0] * self.input_data["HF_peak"]) * 1e6 for i in q]
        for index, point in enumerate(heatflux):
            if point <= self.input_data["HF_peak"] * 1e6 * 0.1:
                heatflux[index] = self.input_data["HF_peak"] * 1e6 * 0.1

        sbar.insert(0, sbar[0])
        sbar.insert(0, 0)
        sbar.append(self.input_data["outer_radius"])
        self.sbar = sbar

        heatflux.insert(0, self.input_data["HF_peak"] * 1e6 * 0.1)
        heatflux.insert(0, self.input_data["HF_peak"] * 1e6 * 0.1)
        heatflux.append(self.input_data["HF_peak"] * 1e6 * 0.1)

        # plt.plot(sbar, heatflux)
        # plt.xlabel("distance, cm")
        # plt.ylabel("Heat Flux, W/m2"

        self.HFf = interp1d(sbar, heatflux)

    def DLH_starter(self, cp, point_x1, line_angle, topy):

        width = self.input_data["width"] + 2 * self.input_data["thickness"]
        height = width / 2 * tan(60 * pi / 180)

        no_points_y = int(round(topy/height, 0))
        no_points_x = int(round(((self.input_data["outer_radius"]-point_x1) /
                                 (1.5*width)), 0))

        self.cross_sectional_area = (3 * sqrt(3) / 2) * (width/2) ** 2

        for i in range(no_points_x):
            for j in range(no_points_y):
                cp.append([i*1.5*width + point_x1, j*height])
                cp.append([i*1.5*width + 0.75*width + point_x1, 
                           j*height + 0.5*height])

        return cp, height, width

    def JIVC_starter(self, cp, point_x1, line_angle, topy):

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
    
    def bounding_box(self, first_radius, second_radius, cp_final_y, \
                     rotation_theta_set, split_theta, assessment_set):

        cp = []  # a list for assigning centre points of hexagon array, (x, y)
        cp_remove = []  # list of values to remove from cp as are out of bounds
        cp_final = []  # list of remaining central points
        cpsets = []

        point_x1 = first_radius*cos(split_theta)
        point_y1 = first_radius*sin(split_theta)
        line_angle = tan(split_theta)
        topy = second_radius * sin(split_theta)

        cell_type_options = {"DLH":self.DLH_starter, \
                             "JIVC":self.JIVC_starter}

        cell_type_selection = cell_type_options[self.cell_type]
        cp, width, height = cell_type_selection(cp, point_x1, line_angle, topy)

        topy = lambda x, x1 = point_x1, y1 = point_y1, m = line_angle: m * (x - x1) + y1
        arc = lambda r, y, y_dis = self.input_data["y_displacement_from_origin"], \
                x_dis = self.input_data["x_displacement_from_origin"]: \
                    ((r**2 - ((y-y_dis)**2)) ** 0.5) + x_dis

        for coordinates in cp:

            if topy(coordinates[0]) <= coordinates[1]:
                cp_remove.append(coordinates)

            # - inner arc
            try:
                arc(first_radius, coordinates[1])
            except:
                pass
            else:
                if arc(first_radius, coordinates[1]) >= coordinates[0]:
                    cp_remove.append(coordinates)

            # - outer arc
            try:
                arc(second_radius, coordinates[1])
            except:
                pass
            else:
                if arc(second_radius, coordinates[1]) <= coordinates[0]:
                    cp_remove.append(coordinates)

        # - remove centre points from cp_final list and copy about 0

        for i in cp:
            if i not in cp_remove:
                cp_final.append([i[0], 0, i[1], sqrt(i[0]**2 + i[1]**2), atan(i[1]/i[0])])
                if cp_final[-1][2] == 0:
                    pass
                else:
                    cp_final.append([i[0], 0, -i[1], sqrt(i[0]**2 + i[1]**2), atan(-i[1]/i[0])])

        # - counting the number of slabs

        count = 0
        count1 = 0
        count2 = 0
        for i in cp_final:
            if i[2] == 0:
                count += 1

        for i in cp_final:
            if i[0] == cp_final[-1][0]:
                count1 += 1
            if i[0] == cp_final[0][0]:
                count2 += 1
        
        self.n_slabs = append(self.n_slabs, count)
        self.n_jets_per_slab = append(self.n_jets_per_slab, (count1 + count2) / 2.0)
        
        assessment_set.append(array(deepcopy(cp_final)))

        for i in rotation_theta_set:
            cpsets.append(array(deepcopy(cp_final)))

        cp_final = []
        
        for i, cpset in enumerate(cpsets):
            for x, y, z, R, theta in cpset:
                cp_final.append([R * cos(rotation_theta_set[i] + theta), \
                                 y, \
                                 R * sin(rotation_theta_set[i] + theta), \
                                 sqrt((sqrt(x**2 + z**2) - self.sbar[2])**2), \
                                 sqrt(x**2 + z**2) - self.sbar[2], \
                                 sqrt(x**2 + z**2), \
                                 atan(z/x)])

        for i in cp_final:
            cp_final_y.append(deepcopy(i))

        return cp_final_y, assessment_set

    def populate_cells(self):

        # - centre point definitions, same for both codes

        cp_final_y = []
        assessment_set = []
        rotation_set = []
        theta_set = []
        self.n_slabs = array([])
        self.n_jets_per_slab = array([])
        
        R1 = self.input_data["strike_radius"] - self.half_width
        R2 = self.input_data["strike_radius"] + self.half_width

        if len(self.n) != self.m1 + self.m2 + 1:

            raise Exception(""" the number of radial slices is not equal to the number of provided lateral slices """)

        if self.cell_type == "DLH":
            self.rset = linspace(self.input_data["inner_radius"], \
                        self.input_data["outer_radius"], \
                        self.m1+1)
        elif self.cell_type == "JIVC":
            self.rset = concat([linspace(self.input_data["inner_radius"], R1, self.m1 + 1), \
                        linspace(R2, self.input_data["outer_radius"], self.m2 + 1)])

        rset = self.rset
        
        theta_set = ((((360.0 / (self.input_data["number_of_plates"] * self.n * \
                    self.input_data["number_of_carriers"])) / 2) - \
                    (self.input_data["swept_angle_gap_between_plates"] / 2)) * \
                    pi / 180.0)

        for i, ni in enumerate(self.n):
            rotation_set.append(theta_set[i] * ni * 2 * linspace((ni-1)*(1/(2*ni)), \
                - (ni-1)*(1/(2*ni)), ni))

        for i, r in enumerate(rset):
            if len(rset) - 1 == i:
                break
            else:
                cp_final_y, assessment_set = self.bounding_box(r, rset[i+1], \
                cp_final_y, rotation_set[i], theta_set[i], assessment_set)

        self.n_jets_per_channel = self.n_slabs * self.n_jets_per_slab
        self.n_jets_total = sum(self.n_jets_per_channel)
        self.microchannel_half_width = self.n_jets_per_slab * self.input_data["height"]
        
        self.assessment_set = assessment_set
        self.panel_hexagons = cp_final_y
    
def write_out(results, cell_type, layup):

    keys = results.keys()
    
    panel_headers = ["output_pressure", "mass_flow_total", "number_of_channels",\
                          "peak_metal_temperature", "peak_coolant_temperature", "break_flag",\
                          "layup_type", "distribution_type", "mass_flow_split", "delta_break_flag", \
                          "manifold_pipe_diameters", "microchannel_half_width", "m1", "m2", \
                          "half_width", "mach_number"]

    unit_headers = ["x", "y", "R", "metal_temperature", "developed_coolant_temperature", \
                      "taus", "heatflux", "input_temperature", "mass_flow", "Cp"]

    panel_data = [[] for i in panel_headers]
        
    try:
        mkdir("panel_data")
    except FileExistsError:
        pass

    try:
        mkdir("panel_data/%s" % cell_type)
    except FileExistsError:
        pass
    
    try:
        mkdir("panel_data/%s/%s" % (cell_type, layup))
    except FileExistsError:
        pass
    try:
        mkdir("panel_data/%s/%s/cell_data" % (cell_type, layup))
    except FileExistsError:
        pass

    for key in keys:
        for index, header in enumerate(panel_headers):
            panel_data[index].append(getattr(results[key], header))
    
    keys_and_data = {k:v for k,v in zip(panel_headers, panel_data)}
    panel_data = DataFrame.from_dict(keys_and_data)
    panel_data.to_excel("panel_data/%s/%s/top_level_data_%s.xlsx" % (cell_type, layup, cell_type))

    for key in keys:
        
        unit_data = [[] for i in unit_headers]
        
        file = "panel_data/%s/%s/cell_data/%s.xlsx" % (cell_type, layup, key)
        if cell_type == "DLH":
            for hexagon_group in results[key].channels:
                for hexagon in hexagon_group.hexagons:
                    for index, header in enumerate(unit_headers):
                        unit_data[index].append(getattr(hexagon, header))
        else:
            for VJ_cell_group in results[key].VJ_cells:
                for VJ_cell in VJ_cell_group:
                    for index, header in enumerate(unit_headers):
                        unit_data[index].append(getattr(VJ_cell, header))

        keys_and_data = {k:v for k,v in zip(unit_headers, unit_data)}
        unit_data = DataFrame.from_dict(keys_and_data)
        unit_data.to_excel(file)

    f = open("hexagon_data.asc", "w")
    for key in keys:
        for group in results[key].assessment_set:
            for cell in group:
                f.write(str(cell[0]) + "\t" + str(cell[2]))
                f.write("\n")
            f.write("\n")
    f.close()

if __name__ == "__main__":
    start_time = time.time()

    # BUSINESS CODE

    results = dict()
    number_of_units_range = linspace(75, 230, 1, dtype = int)  # n / 13, n [1000, 3000]
    mach_number_range = linspace(0.15, 0.2, 1)
    mass_split_range = linspace(0.5, 0.75, 1)
    mass_flow_total = 4 # kg/s
    distribution_option = ["inout"]# "strikepoint" # strikepoint inout
    layup_option = "structured" # "structured" # structured HF_specific
    centre_channel_half_width = linspace(0.06, 0.2, 1) # 0.06 : 1000 , 0.2 : 3000
    m1_range = linspace(2, 2, 1, dtype = int)
    m2_range = linspace(1, 2, 1, dtype = int)
    count = 0
    cell_type = "JIVC"

    print("layup option: %s" % layup_option)
    # print("distribution option: %s" % distribution_option)
    
    for distribution in distribution_option:
        print("distribution option: %s" % distribution)
        for mach_number in mach_number_range:
            print("mach number: %1.2f" % mach_number)
            for R in centre_channel_half_width:
                print("half width: %1.2f" % R)
                for m1 in m1_range:
                    for m2 in m2_range:
                        print("m1: %i, m2: %i" % (m1, m2))
                        for mass_split in mass_split_range:
                            print("mass flow split: %1.2f" % (mass_split))
                            
                            key = "panel" + str(count)
                            count += 1
                            
                            if distribution_option == "inout" and layup_option == "HF_specific":
                                pass
                            else:
                                panel = Panel(float(mach_number), "inputs.xml", "q.asc",
                                              n = array([4 for i in range(m1 + m2 + 1)]), m1 = m1, m2 = m2,
                                              # number_of_hexagons_in_first_channel = number_of_units,
                                              mass_flow_total = mass_flow_total, centre_channel_half_width = R,
                                              distribution_type = distribution, layup_type = layup_option, 
                                              mass_split = mass_split) # n, m = 1 for DLH + JIVC HF specific
        
                                panel(cell_type)
            
                                results[key] = panel
    
    write_out(results, cell_type, layup_option)