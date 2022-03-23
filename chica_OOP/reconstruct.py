# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 11:33:43 2022

@author: jack.taylor
"""

############################### DLH ###############################

def create_DLH_first_channel(cls, input_data, cells):

    obj = cls(input_data, cells)
    
    mach_number_sympy(obj.group_data, trigger = "jet_velocity")
    mass_flow_sympy(obj.group_data)
    
    obj.group_data["mass_flow_total"] = \
        obj.group_data["number_of_hexagons_in_first_channel"] *\
            obj.group_data["mass_flow"]

    obj.repeated_operations_DLH(obj.group_data)

    if obj.group_data["distribution_type"] == "strikepoint":
        reordered_hexagons = sorted(obj.cells.panel_hexagons, key=lambda x: x[3], reverse=False)
        hexagon_selection = reordered_hexagons[:int(obj.group_data["number_of_hexagons_in_first_channel"])]

    elif obj.group_data["distribution_type"] == "inout" and obj.group_data["layup_type"] == "HF_specific":
        reordered_hexagons = sorted(obj.cells.panel_hexagons, key=lambda x: x[3], reverse=False)
        hexagon_selection = reordered_hexagons[:int(obj.group_data["number_of_hexagons_in_first_channel"])]

    elif obj.group_data["distribution_type"] == "inout":
        reordered_hexagons = sorted(obj.cells.panel_hexagons, key=lambda x: x[5], reverse=False)
        hexagon_selection = reordered_hexagons[:int(obj.group_data["number_of_hexagons_in_first_channel"])]

    obj.get_hexagons(obj.group_data, obj.cells, hexagon_selection)

    obj.group_data["output_temperature"] = \
        sum([hexagon.developed_coolant_temperature for hexagon in
             obj.group_data["hexagons"]]) / obj.group_data["number_of_hexagons_in_first_channel"]

    obj.group_data["output_pressure"] = 100E5
    obj.group_data["output_pressure"] -= obj.group_data["pressure_drop"]
    obj.group_data["max_temperature"] = \
        max([hexagon.metal_temperature for hexagon in obj.group_data["hexagons"]])

    return obj

def create_DLH_HF_specific_channel(cls, panel_hexagons, HFf,
                               cross_sectional_area, remaining_hexagons,
                               previous_channel, mass_flow_actual, Ma=0.0):
    
    obj = cls(previous_channel.output_pressure,
              previous_channel.output_temperature, HFf,
              cross_sectional_area, panel_hexagons = panel_hexagons)
    
    obj.mass_flow = previous_channel.mass_flow
    obj.max_temperature_from_previous_channel = \
        previous_channel.max_temperature

    obj.mass_flow, obj.number_of_hexagons_required = mass_flow_next_row(
        remaining_hexagons[0], obj.mass_flow,
        obj.working_fluid.specific_heat_capacity,
        obj.cross_sectional_area,
        obj.max_temperature_from_previous_channel,
        obj.input_temperature, obj.HFf, mass_flow_actual)

    remainging_hexagons, obj.break_flag = \
        obj.repeated_operations_DLH(1, remaining_hexagons, mass_flow_actual)

    obj.max_temperature = previous_channel.max_temperature

    return obj

def create_DLH_structured_channel(cls, panel_hexagons, HFf,
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

def strikepoint_DLH(self, first_channel):

    # this essentially acts as a pre conditioner for the add_channels function
    # doesn't change dependant on HF_specific or structured
    # this could be the first @decorator for a strike point DLH ... setup

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

# inout wrapper, as need to do the organising of the cells first, which is 
# dependant on structured vs HF specific. This reordering is done first

def inout_wrapper(self, func):

    def wrapper(*args, **kwargs):

        all_remaining_hexagons, first_channel = func(*args, **kwargs)

        new_channels = [[]]

        self.mass_flow_actual = [self.mass_flow_total]

        self.add_channels(all_remaining_hexagons, new_channels, first_channel)

        for direction in new_channels:
            for channel in direction:
                self.channels.append(channel)

    return wrapper

def structured_DLH(self, first_channel):

    reordered_hexagons = sorted(self.cells.heat_cells, key=lambda x: x[5], reverse=False)

    all_remaining_hexagons = [reordered_hexagons]

    return all_remaining_hexagons, first_channel


def HF_specific_DLH(self, first_channel):

    first_round = [[first_channel]] # try to work this without this statement

    reordered_hexagons = sorted(self.cells.heat_cells, key=lambda x: x[3], reverse=False)

    inboard_hexagons, outboard_hexagons = [hexagons for hexagons in reordered_hexagons if hexagons[4]
         <= 0.0], [hexagons for hexagons in reordered_hexagons if hexagons[4] >= 0.0]

    # self.dummy_set_all_remaining_hexagons = [deepcopy(inboard_hexagons)] # currently used, not ideal
    all_remaining_hexagons = [inboard_hexagons, outboard_hexagons]

    return all_remaining_hexagons, first_round

############################### DLH ###############################

# proposed implementation
first_channel = []

test = inout_wrapper(HF_specific_DLH(first_channel))
test1 = inout_wrapper(structured_DLH(first_channel))

test2 = strikepoint_DLH(first_channel) # doesn't need to distinguish between HF_specific and structured






############################### JIVC ###############################

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
    
    # output setup
    
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
    
    # setup first group
    
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

    # repeated sections ^^

    assessment_set = [flip(self.assessment_set[0], 0), self.assessment_set[1], self.assessment_set[2]]
    remaining_set = deepcopy(assessment_set)

    # loop

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
    
    # setup first group
    
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
    mach_numbers.append(mach_number_iter)
    VJ_cells_panel.append(VJ_cells)
    
    # repeated sections ^^

    input_pressure = deepcopy(strikepoint_pressure)
    input_temperature = deepcopy(strikepoint_temperature)
    
    # loop

    loop_arrangement = [range(0, self.m1).__reversed__(), \
        range(self.m1 + 1, self.m1 + self.m2 + 1)]

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

############################### JIVC ###############################








