from numpy import sqrt, pi
from scipy.interpolate import interp1d
from numpy import array, flip, delete, append
from copy import deepcopy
import xml.etree.ElementTree as ET
from CellGroups import Group, JIVCCells, DLHCells
from Material import Material
from CHICA.utility import get_example_data_path
from CHICA.flow_properties import mach_number_sympy, mass_flow_sympy, \
                                    mass_flow_next_row

class Cassette:

    def __init__(self, input_file, heat_flux_file):

        self.channels = []
        self.panel_hexagon = []
        self.read_inputs(input_file, heat_flux_file)
        n = array([4 for i in range(int(self.input_data["m1"]) +
                                    int(self.input_data["m2"]) + 1)])
        self.input_data["n"] = n

    def __call__(self):
        
        if self.input_data["cell_type"] == "DLH":
            self.cells = DLHCells(self.input_data)
            self.configure_DLH()

        elif self.input_data["cell_type"] == "JIVC":
            self.cells = JIVCCells(self.input_data)
            self.configure_JIVC()

        else:
            raise ValueError("Error, no cell type assigned")

    def configure_DLH(self):
        
        # adds hexagons to middle channel and calcs flow properties
        self.channels.append(Group.create_DLH_first_channel(self.input_data, self.cells))
        # self.mass_flow_total = self.channels[0].mass_flow_total

        implementation_options = {"strikepoint": self.strikepoint_DLH,
                                  "inout": self.inout_DLH}
        implementation_choice = self.input_data["distribution_type"]

        run = implementation_options[implementation_choice]
        run(self.channels[0])

        # output_pressures = []
        # output_metal_temperatures = []
        # output_coolant_temperatures = []
        # number_of_channels = 0
        # break_flag = []

        # for channel in self.channels:
        #     output_pressures.append(channel.output_pressure)
        #     number_of_channels += 1
        #     output_metal_temperatures.append(channel.max_temperature)
        #     output_coolant_temperatures.append(channel.output_temperature)
        #     break_flag.append(channel.break_flag)
        
        # self.output_pressure = min(output_pressures)
        # self.number_of_channels = number_of_channels
        # self.peak_metal_temperature = max(output_metal_temperatures)
        # self.peak_coolant_temperature = max(output_coolant_temperatures)
        # self.break_flag = max(break_flag)

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
        
        reordered_hexagons = sorted(first_channel.cells.panel_hexagons, key=lambda x: x[3], reverse=False)
        inboard_hexagons, outboard_hexagons = [hexagons for hexagons in reordered_hexagons if hexagons[4] <= 0.0],\
                                                [hexagons for hexagons in reordered_hexagons if hexagons[4] >= 0.0]

        all_remaining_hexagons = [inboard_hexagons, outboard_hexagons]

        self.mass_flow_actual = [first_channel["mass_flow_total"] * (1-first_channel["mass_flow_split"]),
                                 first_channel["mass_flow_total"] * first_channel["mass_flow_split"]]

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
                    pressure_drop, input_temperature, mach_number_iter, peak_metal_temperature, VJ_cells, working_fluid, jet_velocity = \
                        Group.repeated_operations_JIVC(self, self.assessment_set[i], self.n[i], \
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
            
            pressure_drop, input_temperature, mach_number_iter, peak_metal_temperature, VJ_cells, working_fluid, jet_velocity = \
                Group.repeated_operations_JIVC(self, self.assessment_set[i], self.n[i], \
                peak_metal_temperature, input_pressure, input_temperature, self.mass_flow_total, self.HFf, self.cross_sectional_area)
            input_pressure += - pressure_drop
            
            manifold_pipe_diameters.append([sqrt(((self.mass_flow_total * self.mass_flow_split * 4) / \
                                             self.n[self.m1]) / \
                                            (working_fluid.density * jet_velocity * pi)), i])
            
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
        self.n = array([10 for i in range(self.m1 + self.m2 + 1)])
        
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
        
        pressure_drop, strikepoint_temperature, mach_number_iter, peak_metal_temperature, VJ_cells, working_fluid, jet_velocity =\
            Group.repeated_operations_JIVC(self, self.assessment_set[self.m1], self.n[self.m1],\
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
                
                working_fluid = Material("helium")
                working_fluid.fluid_properties(input_pressure, input_temperature)
                
                jet_mass_flow, n_jets_total = mass_flow_next_row(\
                remaining_set[loop][0], self.cross_sectional_area, \
                self.HFf, jet_mass_flow, working_fluid.specific_heat_capacity, \
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
                
                pressure_drop, input_temperature, mach_number_iter, peak_metal_temperature, VJ_cells, working_fluid, jet_velocity =\
                    Group.repeated_operations_JIVC(self, next_group, self.n[loop], \
                    peak_metal_temperature, input_pressure, input_temperature, \
                    mass_flow_total * self.mass_flow_split, self.HFf, self.cross_sectional_area)

                manifold_pipe_diameters.append([sqrt(((mass_flow_total * self.mass_flow_split * 4) / \
                                             self.n[self.m1]) / \
                                            (working_fluid.density * jet_velocity * pi)), loop])

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
        
        pressure_drop, strikepoint_temperature, mach_number_iter, peak_metal_temperature, VJ_cells, working_fluid, jet_velocity =\
            Group.repeated_operations_JIVC(self, self.assessment_set[self.m1], self.n[self.m1], \
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
                
                pressure_drop, input_temperature, mach_number_iter, peak_metal_temperature, VJ_cells, working_fluid, jet_velocity =\
                    Group.repeated_operations_JIVC(self, self.assessment_set[channel], self.n[channel], \
                    peak_metal_temperature, input_pressure, input_temperature, \
                    mass_flow_total * self.mass_flow_split, self.HFf, self.cross_sectional_area)

                manifold_pipe_diameters.append([sqrt(((mass_flow_total * self.mass_flow_split * 4) / \
                                             self.n[self.m1]) / \
                                            (working_fluid.density * jet_velocity * pi)), channel])

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
        # doesn't work for inout + HF specific
        layup_options = {"structured": Group.create_DLH_structured_channel,
                         "HF_specific": Group.create_DLH_HF_specific_channel}

        layup_selection = layup_options[self.layup_type]

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

    def add_channels_inout_FIX(self):

        # HF_specific from middle out to get channel sizes
        dummy_hexagons = deepcopy(self.panel_hexagons)

        for direction in self.dummy_set_all_remaining_hexagons:

            mass_flow_actual = self.mass_flow_actual[0]

            new_channels[0].append(
                Group.create_HF_specific_channel(
                    dummy_hexagons, self.HFf, self.cross_sectional_area,
                    direction, first_channel, mass_flow_actual))

            while len(direction) >= 1.0:

                new_channels[0].append(
                    Group.create_HF_specific_channel(
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
            Group.create_first_channel(
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
                    Group.create_structured_channel(
                        self.panel_hexagons, self.HFf,
                        self.cross_sectional_area, direction,
                        self.new_channels[0][-1], mass_flow_actual,
                        specified_number_of_hexagons=
                        len(channel.hexagons)))

                if self.new_channels[0][-1].break_flag == 1:
                    break

            while len(direction) >= 1.0:

                self.new_channels[0].append(
                    Group.create_HF_specific_channel(
                        self.panel_hexagons, self.HFf,
                        self.cross_sectional_area, direction,
                        self.new_channels[0][-1], mass_flow_actual))

                if self.new_channels[0][-1].break_flag == 1:
                    break

    def read_inputs(self, inputs_filename, q_filename):
        # read asc files and store results to member variables (self...)

        user_data = open(inputs_filename)
        heatflux_data = get_example_data_path(q_filename)

        input_data = dict()
        tree = ET.parse(user_data)
        root = tree.getroot()

        keys = [child.keys() for child in root]

        for index, key in enumerate(keys):
            try:
                input_data[key[0]] = float(root[index].attrib[key[0]])
            except ValueError:
                input_data[key[0]] = root[index].attrib[key[0]]

        q = []  # power input from VITA
        s = []  # distance from strike point

        with open(heatflux_data) as geometry:
            for line in geometry:
                heatflux = line.split()
                s.append(float(heatflux[0]) / 100.0)  # converting to m from cm
                q.append(float(heatflux[1])) # / 6.87) * 30.0)

        sbar = [i + input_data["strike_radius"] for i in s]
        heatflux = [(i/q[0] * input_data["HF_peak"]) * 1e6 for i in q]
        for index, point in enumerate(heatflux):
            if point <= input_data["HF_peak"] * 1e6 * 0.1:
                heatflux[index] = input_data["HF_peak"] * 1e6 * 0.1

        sbar.insert(0, sbar[0])
        sbar.insert(0, 0)
        sbar.append(input_data["outer_radius"])
        # self.sbar = sbar

        heatflux.insert(0, input_data["HF_peak"] * 1e6 * 0.1)
        heatflux.insert(0, input_data["HF_peak"] * 1e6 * 0.1)
        heatflux.append(input_data["HF_peak"] * 1e6 * 0.1)

        # plt.plot(sbar, heatflux)
        # plt.xlabel("distance, cm")
        # plt.ylabel("Heat Flux, W/m2"

        # self.HFf = interp1d(sbar, heatflux)
        HFf = interp1d(sbar, heatflux)
        
        input_data["sbar"] = sbar
        input_data["HFf"] = HFf
        
        self.input_data = input_data