from CoolProp.CoolProp import PropsSI as SI

class Material:

    def __init__(self, material):
        
        self.material = material

    def fluid_properties(self, input_pressure, input_temperature):
        self.viscosity = SI("V", "T", input_temperature, "P", 
            input_pressure, self.material)
        self.density = SI("D", "T", input_temperature, "P", 
            input_pressure, self.material)
        self.specific_heat_capacity = SI("C", "T", input_temperature, "P",
            input_pressure, self.material)
    
    def solid_properties(self):
        # find a solid properties library, similar to CoolProp
        pass