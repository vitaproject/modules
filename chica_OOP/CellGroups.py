from numpy import average, pi, linspace, array, sin, cos, tan, sqrt, append
from math import atan
from numpy import concatenate as concat
from copy import deepcopy
from collections import namedtuple
from Material import Material
from HeatCells import HeatCell
from CHICA.flow_properties import mach_number_sympy, \
                                        mass_flow_sympy, \
                                        reynolds_sympy, \
                                        euler_sympy, \
                                        pdrop_sympy, \
                                        m_star_sympy, \
                                        taus_sympy, \
                                        mass_flow_next_row
class Group:

    def __init__(self, input_data, cells):

        self.group_data = deepcopy(input_data)
        # self.cells = deepcopy(cells)
        working_fluid = Material("helium")
        working_fluid.fluid_properties(self.group_data["input_pressure"],\
                                            self.group_data["input_temperature"])
        self.group_data["fluid_properties"] = working_fluid

    def get_hexagons(self, group_data, cells, hexagon_selection):

        group_data["hexagons"] = []

        for hexagon in hexagon_selection:

            group_data["hexagons"].append(HeatCell(hexagon, group_data))
            print(hexagon)
            cells.heat_cells[int(group_data["direction"])].remove(hexagon)

    @classmethod
    def create_DLH_first_channel(cls, input_data, cells):

        obj = cls(input_data, cells)

        mach_number_sympy(obj.group_data, trigger = "jet_velocity")
        mass_flow_sympy(obj.group_data)

        obj.group_data["mass_flow_total"] = \
            obj.group_data["number_of_hexagons_in_first_channel"] *\
                obj.group_data["mass_flow"]

        obj.repeated_operations_DLH(obj.group_data)

        if obj.group_data["distribution_type"] == "strikepoint":
            reordered_hexagons = sorted(cells.heat_cells[0], key=lambda x: x[3], reverse=False)
            hexagon_selection = reordered_hexagons[:int(obj.group_data["number_of_hexagons_in_first_channel"])]

        elif obj.group_data["distribution_type"] == "inout" and obj.group_data["layup_type"] == "HF_specific":
            reordered_hexagons = sorted(cells.heat_cells[0], key=lambda x: x[3], reverse=False)
            hexagon_selection = reordered_hexagons[:int(obj.group_data["number_of_hexagons_in_first_channel"])]

        elif obj.group_data["distribution_type"] == "inout":
            reordered_hexagons = sorted(cells.heat_cells[0], key=lambda x: x[5], reverse=False)
            hexagon_selection = reordered_hexagons[:int(obj.group_data["number_of_hexagons_in_first_channel"])]

        obj.get_hexagons(obj.group_data, cells, hexagon_selection)

        obj.group_data["output_temperature"] = \
            sum([hexagon.developed_coolant_temperature for hexagon in
                 obj.group_data["hexagons"]]) / obj.group_data["number_of_hexagons_in_first_channel"]

        obj.group_data["output_pressure"] = obj.group_data["input_pressure"]
        obj.group_data["output_pressure"] -= obj.group_data["pressure_drop"]
        obj.group_data["max_temperature"] = \
            max([hexagon.metal_temperature for hexagon in obj.group_data["hexagons"]])

        return obj

    @classmethod
    def create_DLH_structured_channel(cls, input_data, cells):

        obj = cls(input_data, cells)

        max_temperature_from_previous_channel = obj.group_data["max_temperature"]

        if obj.group_data["specified_number_of_hexagons"] is None:
            obj.group_data["specified_number_of_hexagons"] = len(cells.heat_cells)

        if obj.group_data["specified_number_of_hexagons"] >= len(cells.heat_cells):
            obj.group_data["specified_number_of_hexagons"] = len(cells.heat_cells)

        obj.repeated_operations_DLH(obj.group_data)

        hexagon_selection = cells.heat_cells[obj.group_data["direction"]][:int(obj.group_data["specified_number_of_hexagons"])]
        
        obj.get_hexagons(obj.group_data, cells, hexagon_selection[0])

        obj.group_data["output_pressure"] -= obj.group_data["pressure_drop"]
        obj.group_data["max_temperature"] = max([hexagon.metal_temperature for hexagon in obj.group_data["hexagons"]])

        if obj.group_data["output_pressure"] <= 0:
            obj.group_data["break_flag"] = 1
        else:
            obj.group_data["break_flag"] = 0
        
        return obj

    @classmethod
    def create_DLH_HF_specific_channel(cls, heat_cells, HFf,
                                   cross_sectional_area, remaining_hexagons,
                                   previous_channel, mass_flow_actual, Ma=0.0):
        
        obj = cls(previous_channel.output_pressure,
                  previous_channel.output_temperature, HFf,
                  cross_sectional_area, heat_cells = heat_cells)
        
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

    def repeated_operations_JIVC(self, assessment_set, n, peak_metal_temperature, \
        input_pressure, input_temperature, mass_flow_total, HFf, cross_sectional_area):
        
        cell_temperatures = []
        VJ_cells = []
        
        working_fluid = Material("helium")
        working_fluid.fluid_properties(input_pressure, input_temperature)
        
        assessment_set_total_mass_flow = mass_flow_total / n
        jet_mass_flow = assessment_set_total_mass_flow / len(assessment_set)
        
        jet_mass_flow, jet_velocity, Ajet, density, jet_diameter \
            = mass_flow_sympy(input_temperature, input_pressure, massflow = jet_mass_flow)
        
        Reynolds_number, jet_velocity = Reynolds_sympy(density, working_fluid.viscosity, jet_diameter, \
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
            VJ_cell_i = HeatCell(cell, HFf, cross_sectional_area, 
                working_fluid.specific_heat_capacity, \
                input_temperature, jet_mass_flow, taus)
            VJ_cell_i.jet_mass_flow = jet_mass_flow
            VJ_cells.append(VJ_cell_i)
            cell_temperatures.append(VJ_cell_i.developed_coolant_temperature)
            if peak_metal_temperature <= VJ_cell_i.metal_temperature:
                peak_metal_temperature = VJ_cell_i.metal_temperature
        
        new_input_temperature = average(cell_temperatures)
        
        return pressure_drop, new_input_temperature, mach_number, peak_metal_temperature, VJ_cells, working_fluid, jet_velocity

    def repeated_operations_DLH(self, group_data):

        m_star_sympy(group_data, trigger = "non_dimensional_mass_flow")
        reynolds_sympy(group_data, trigger = "reynolds_number")
        euler_sympy(group_data, trigger = "euler_number")
        pdrop_sympy(group_data, trigger = "pressure_drop")
        taus_sympy(group_data, trigger = "non_dimensional_temperature")

class CellManipulation:

    # - centre point definitions, same for both codes

    def __init__(self):

        self.count = None
        self.count1 = None
        self.count2 = None
        self.rset = None
        self.assessment_set = None
        self.heat_cells = None

    def _case_setup(self):
        raise NotImplementedError()

    def _starter(self):
        raise NotImplementedError()

    def _populate_cells(self, input_data): # this is essentially a manager function

        cp_final_y = []
        assessment_set = []
        rotation_set = []
        theta_set = []

        if len(input_data["n"]) != input_data["m1"] + input_data["m2"] + 1:
            raise Exception(""" the number of radial slices is not equal to the number of provided lateral slices """)

        # self.case_setup(input_data)

        rset = self.rset

        theta_set = ((((360.0 / (input_data["number_of_plates"] * input_data["n"] * \
                    input_data["number_of_carriers"])) / 2) - \
                    (input_data["swept_angle_gap_between_plates"] / 2)) * \
                    pi / 180.0)

        for i, ni in enumerate(input_data["n"]):
            rotation_set.append(theta_set[i] * ni * 2 * linspace((ni-1)*(1/(2*ni)), \
                - (ni-1)*(1/(2*ni)), ni))

        for i, r in enumerate(rset):
            if len(rset) - 1 == i:
                break
            else:
                cp_final_y, assessment_set = self._bounding_box(input_data, r, rset[i+1], \
                cp_final_y, rotation_set[i], theta_set[i], assessment_set, input_data["sbar"])

        # self.jets(input_data, count, count1, count2)

        self.assessment_set = assessment_set
        self.heat_cells = [cp_final_y]

    def _bounding_box(self, input_data, first_radius, second_radius, cp_final_y, \
                     rotation_theta_set, split_theta, assessment_set, sbar):

        cp = []  # a list for assigning centre points of hexagon array, (x, y)
        cp_remove = []  # list of values to remove from cp as are out of bounds
        cp_final = []  # list of remaining central points
        cpsets = []

        point_x1 = first_radius*cos(split_theta)
        point_y1 = first_radius*sin(split_theta)
        line_angle = tan(split_theta)
        topy = second_radius * sin(split_theta)

        cp = self._starter(input_data, cp, point_x1, line_angle, topy)

        topy = lambda x, x1 = point_x1, y1 = point_y1, m = line_angle: m * (x - x1) + y1
        arc = lambda r, y, y_dis = input_data["y_displacement_from_origin"], \
                x_dis = input_data["x_displacement_from_origin"]: \
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
        
        assessment_set.append(array(deepcopy(cp_final)))

        for i in rotation_theta_set:
            cpsets.append(array(deepcopy(cp_final)))

        cp_final = []
        coord = namedtuple("coordinate", ["x", "y", "z", "R1", "R2", "R3"])
        
        for i, cpset in enumerate(cpsets):
            for x, y, z, R, theta in cpset:
                # optional new method, uses
                cp_final.append(coord(R * cos(rotation_theta_set[i] + theta), \
                                 y, \
                                 R * sin(rotation_theta_set[i] + theta), \
                                 sqrt((sqrt(x**2 + z**2) - sbar[2])**2), \
                                 sqrt(x**2 + z**2) - sbar[2], \
                                 sqrt(x**2 + z**2)))

                # cp_final.append([R * cos(rotation_theta_set[i] + theta), \
                #                  y, \
                #                  R * sin(rotation_theta_set[i] + theta), \
                #                  sqrt((sqrt(x**2 + z**2) - sbar[2])**2), \
                #                  sqrt(x**2 + z**2) - sbar[2], \
                #                  sqrt(x**2 + z**2)])
                                 # atan(z/x)]) # don't think this is necessary

        for i in cp_final:
            cp_final_y.append(deepcopy(i))

        self.count = count
        self.count1 = count1
        self.count2 = count2

        return cp_final_y, assessment_set

class JIVCCells(CellManipulation):

    def __init__(self, input_data):
        # all member variables of class should be anncounced here!!!
        self.n_slabs = array([])
        self.n_jets_per_slab = array([])
        self.n_jets_per_channel = None
        self.n_jets_total = None
        self.microchannel_half_width = None
        self.input_data = input_data

    def start_JIVC(self):

        super().__init__()
        self._case_setup(self.input_data)
        self._populate_cells(self.input_data)
        self._jets(self.input_data)

    def _case_setup(self, input_data):
        R1 = input_data["strike_radius"] - input_data["half_width"]
        R2 = input_data["strike_radius"] + input_data["half_width"]

        self.rset = concat([linspace(input_data["inner_radius"], R1, int(input_data["m1"]) + 1), \
                        linspace(R2, input_data["outer_radius"], int(input_data["m2"]) + 1)])

    def _starter(self, input_data, cp, point_x1, line_angle, topy):
        height = input_data["height"]
        width = input_data["width"]
        input_data["cross_sectional_area"] = height * width

        no_points_y = int(round(topy/height, 0))
        no_points_x = int(round((input_data["outer_radius"]-point_x1) /
                                 (width), 0))

        for i in range(no_points_x):
            for j in range(no_points_y):
                cp.append([i*width + point_x1, j*height])

        return cp

    def _jets(self, input_data):

        self.n_slabs = append(self.n_slabs, self.count)
        self.n_jets_per_slab = append(self.n_jets_per_slab, min(self.count1, self.count2))

        self.n_jets_per_channel = self.n_slabs * self.n_jets_per_slab
        self.n_jets_total = sum(self.n_jets_per_channel)
        self.microchannel_half_width = self.n_jets_per_slab * input_data["height"]

class DLHCells(CellManipulation):

    def __init__(self, input_data):
        # all member variables of class should be anncounced here!!!
        self.input_data = input_data

    def start_DLH(self):

        super().__init__()
        self._case_setup(self.input_data)
        self._populate_cells(self.input_data)

    def _case_setup(self, input_data):

        self.rset = linspace(input_data["inner_radius"], \
                            input_data["outer_radius"], \
                            int(input_data["m1"]) + 1)    

    def _starter(self, input_data, cp, point_x1, line_angle, topy):

        width = input_data["width"] + 2 * input_data["thickness"]
        height = width / 2 * tan(60 * pi / 180)

        no_points_y = int(round(topy/height, 0))
        no_points_x = int(round(((input_data["outer_radius"]-point_x1) /
                                 (1.5*width)), 0))

        input_data["cross_sectional_area"] = (3 * sqrt(3) / 2) * (width/2) ** 2

        for i in range(no_points_x):
            for j in range(no_points_y):
                cp.append([i*1.5*width + point_x1, j*height])
                cp.append([i*1.5*width + 0.75*width + point_x1, 
                           j*height + 0.5*height])

        return cp
