from numpy import sqrt
from CHICA.flow_properties import metal_temperature_sympy, \
                                    coolant_temperature_sympy

class HeatCell:

    cell_ID = 0

    def __init__(self, cell_data, group_data):
                 # HFf, cross_sectional_area, specific_heat_capacity, \
                 # input_temperature, mass_flow, non_dimensional_temperature):

        self.x = cell_data[0]
        self.y = cell_data[2]
        self.R = sqrt((self.x**2) + (self.y**2))
        heatflux = group_data["HFf"](self.R) * group_data["cross_sectional_area"]

        self.metal_temperature = metal_temperature_sympy(
            group_data["input_temperature"],
            group_data["non_dimensional_temperature"],
            group_data["mass_flow"], heatflux,
            group_data["fluid_properties"].specific_heat_capacity)

        self.developed_coolant_temperature = coolant_temperature_sympy(
            group_data["non_dimensional_temperature"],
            group_data["input_temperature"], self.metal_temperature)

        self.cell_ID = self.identification()

    @classmethod
    def identification(cls):
        cls.cell_ID += 1
        return cls.cell_ID