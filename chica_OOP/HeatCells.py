from numpy import sqrt
from dataclasses import dataclass
from CHICA.flow_properties import metal_temperature_sympy, \
                                    coolant_temperature_sympy

class HeatCell:

    cell_ID = 0

    def __init__(self, cell_data, group_data):
    
        self.x = cell_data.x
        self.z = cell_data.z
        self.R = sqrt((self.x**2) + (self.z**2))
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

@dataclass
class HeatCellData:
    
    ID: float
    developed_coolant_temperature: float
    metal_temperature: float