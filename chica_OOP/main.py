from Cassettes import Cassette
from numpy import linspace
from CHICA.utility import get_example_data_path
from itertools import product
from pandas import DataFrame
from os import mkdir
from time import time
import xml.etree.ElementTree as ET

class XML:

    @property
    def xml_data(self):
        return self._xml_data

    @xml_data.setter
    def xml_data(self, new_xml_file):
        self._xml_data = new_xml_file

    @xml_data.deleter
    def xml_data(self):
        del self._xml_data

class SetupFile:

    def save(self, log):
        
        self.save_file = log

    @staticmethod
    def new_input_file(directory):

        data = ET.Element("input_data")
        ET.SubElement(data, "cell_type", attrib = {"type":"DLH"})
        ET.SubElement(data, "number_of_plates", attrib = {"min":"12", "max":"12", "num":"1"})
        ET.SubElement(data, 'width', attrib = {"min":"0.0064", "max":"0.0064", "num":"1"})
        ET.SubElement(data, 'height', attrib = {"min":"0.00308", "max":"0.00308", "num":"1"})
        ET.SubElement(data, 'inner_radius', attrib = {"min":"1.093", "max":"1.093", "num":"1"})
        ET.SubElement(data, 'outer_radius', attrib = {"min":"2.212", "max":"2.212", "num":"1"})
        ET.SubElement(data, 'x_displacement_from_origin', attrib = {"min":"0", "max":"0", "num":"1"})
        ET.SubElement(data, 'y_displacement_from_origin', attrib = {"min":"0", "max":"0", "num":"1"})
        ET.SubElement(data, 'swept_angle_gap_between_plates', attrib = {"min":"0", "max":"0", "num":"1"})
        ET.SubElement(data, 'strike_radius', attrib = {"min":"1.65", "max":"1.65", "num":"1"})
        ET.SubElement(data, 'input_temperature', attrib = {"min":"373.15", "max":"373.15", "num":"1"})
        ET.SubElement(data, 'input_pressure', attrib = {"min":"100E5", "max":"100E5", "num":"1"})
        ET.SubElement(data, 'jet_diameter', attrib = {"min":"1E-3", "max":"1E-3", "num":"1"})
        ET.SubElement(data, 'number_of_carriers', attrib = {"min":"3", "max":"3", "num":"1"})
        ET.SubElement(data, 'HF_peak', attrib = {"min":"10", "max":"10", "num":"1"})
        ET.SubElement(data, 'mach_number', attrib = {"min":"0.2", "max":"0.2", "num":"1"})
        ET.SubElement(data, 'number_of_hexagons_in_first_channel', attrib = {"min":"3000", "max":"3000", "num":"1"})
        ET.SubElement(data, 'mass_split', attrib = {"min":"0.5", "max":"0.5", "num":"1"})
        ET.SubElement(data, 'layup_type', attrib = {"type":"HF_specific"})
        ET.SubElement(data, 'distribution_type', attrib = {"type":"inout"})
        ET.SubElement(data, 'm1', attrib = {"min":"1", "max":"1", "num":"1"})
        ET.SubElement(data, 'm2', attrib = {"min":"1", "max":"1", "num":"1"})
        ET.SubElement(data, 'mass_flow_total', attrib = {"min":"4", "max":"4", "num":"1"})
        ET.SubElement(data, 'half_width', attrib = {"min":"0.06", "max":"0.06", "num":"1"})
        ET.SubElement(data, 'thickness', attrib = {"min":"0", "max":"0", "num":"1"})
        ET.SubElement(data, "direction", attrib = {"type":"0"})

        mydata = ET.tostring(data)
        myfile = open(directory + "/parameter_sweep.xml", "wb")
        myfile.write(mydata)

    @staticmethod
    def setup_runs(config_file):

        config_file = get_example_data_path(config_file)
        tree = ET.parse(config_file)
        root = tree.getroot()

        key_ring = [i.keys() for i in root]

        variations = SetupFile.number_of_variations(root)
        config_ID, number_of_runs = SetupFile.create_folders(root, variations)
        all_parameters = SetupFile.create_run_parameters(root, key_ring)
        SetupFile.create_run_configs(all_parameters, config_ID, number_of_runs)
        results = SetupFile.run_cases(config_ID, number_of_runs)
        
        return results
        
        ## write the config files into the correct folders using the correct WriteFile
    
    @staticmethod
    def number_of_variations(root):
        ## find out the number of runs required based on "num", does not work for entries without "num"
        variations = 1
        for child in root:
            try:
                if int(child.attrib["num"]) > 1:
                    variations *= int(child.attrib["num"])
            except KeyError:
                pass
        
        return variations
        
    @staticmethod
    def create_folders(root, variations):
        ## create config folder, this will contain runs, and creates run folders
        config_num = 0
        while config_num >= 0:
            try:
                mkdir("config_%i" % config_num)
                break
            except FileExistsError:
                config_num += 1
        config_num = config_num
        
        for index, run in enumerate(range(variations)):
            mkdir ("config_%i/run_%i" % (config_num, index))

        return config_num, index

    @staticmethod
    def create_run_parameters(root, key_ring):

        all_parameters = [[]]

        # repeat product function
        for index, keys in enumerate(key_ring):
            if len(keys) == 1:
                
                all_parameters[0].append(root[index].attrib[keys[0]])

            if len(keys) == 3:

                nested_list = list(product(all_parameters, 
                                           linspace(float(root[index].attrib[keys[0]]), 
                                           float(root[index].attrib[keys[1]]),
                                           int(root[index].attrib[keys[2]]))))

                all_parameters = SetupFile.remove_tuples(nested_list)

        return all_parameters

    @staticmethod
    def remove_tuples(packed_list):

        for i, x in enumerate(packed_list):
            try:
                packed_list[i] = list(x)
            except TypeError:
                pass
        
        for nested_list in packed_list:
            for values in nested_list:
                if type(values) == list:
                    for listed_value in values:
                        nested_list.insert(len(nested_list) - 1, listed_value)
                    del(nested_list[0])
        
        unpacked_list = packed_list
        return unpacked_list

    @staticmethod
    def create_run_configs(all_parameters, config_ID, number_of_runs):
        for index, run in enumerate(all_parameters):
            WriteFile(run, "config_%i/run_%i" % (config_ID, index), index)

    @staticmethod
    def run_cases(config_ID, number_of_runs):
        
        results = []
        
        for run in range(number_of_runs + 1):
            cassette = Cassette("config_%i/run_%i/run_%i.xml" % (config_ID, run, run),
                                "q.asc")
            cassette()
            results.append(cassette)

        return results

## not sure I need the separate WriteFile, could just put it all in SetupFile

class WriteFile:

    def __init__(self, run_parameters, directory, run_ID):

        WriteFile.input_file(directory, run_parameters, run_ID)

    @staticmethod
    def input_file(directory, parameters, run_ID):

        WriteFile.default(directory, run_ID)
        run_file = open(directory + "/run_%i.xml" % run_ID)

        tree = ET.parse(run_file)
        root = tree.getroot()

        for index, child in enumerate(root):
            child.set(child.tag, "%s" % parameters[index])
        
        tree.write(directory + "/run_%i.xml" % run_ID)

    @staticmethod
    def default(directory, run_ID):

        data = ET.Element("input_data")
        ET.SubElement(data, "cell_type")
        ET.SubElement(data, "number_of_plates")
        ET.SubElement(data, 'width')
        ET.SubElement(data, 'height')
        ET.SubElement(data, 'inner_radius')
        ET.SubElement(data, 'outer_radius')
        ET.SubElement(data, 'x_displacement_from_origin')
        ET.SubElement(data, 'y_displacement_from_origin')
        ET.SubElement(data, 'swept_angle_gap_between_plates')
        ET.SubElement(data, 'strike_radius')
        ET.SubElement(data, 'input_temperature')
        ET.SubElement(data, 'input_pressure')
        ET.SubElement(data, 'jet_diameter')
        ET.SubElement(data, 'number_of_carriers')
        ET.SubElement(data, 'HF_peak')
        ET.SubElement(data, 'mach_number')
        ET.SubElement(data, 'number_of_hexagons_in_first_channel')
        ET.SubElement(data, 'mass_split')
        ET.SubElement(data, 'layup_type')
        ET.SubElement(data, 'distribution_type')
        ET.SubElement(data, 'm1')
        ET.SubElement(data, 'm2')
        ET.SubElement(data, 'mass_flow_total')
        ET.SubElement(data, 'half_width')
        ET.SubElement(data, 'thickness')
        ET.SubElement(data, "direction")

        mydata = ET.tostring(data)
        myfile = open(directory + "/run_%i.xml" % run_ID, "wb")
        myfile.write(mydata)

    
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
    start_time = time()

    log = SetupFile()
    log.new_input_file("CHICA")
    results = log.setup_runs("parameter_sweep.xml")
    # log.save(results)
    
    # write_out(results, "DLH", "structured")
