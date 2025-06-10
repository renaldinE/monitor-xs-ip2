#%% Relevant packages

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from os.path import join
import xlsxwriter

#%% Custom packages

from .core.utils import load_config

config = load_config()

#%% Print results of test targets in Excel file

def print_test_target_results(targets, dir_to_save):
    
    # Get radionuclides produced according to target material
    df = pd.read_excel(config['MATERIALS_DATA'], sheet_name='Material_NuclideInventory')
    
    # Define common headers
    headers = [
        'Target No.',
        'Irradiation datetime',
        'Proton energy (MeV)',
        ]
    
    # Define headers specific to target material
    head_nucs_eval = {material: [nuc+'_eval' for nuc in df[material].dropna().to_list()] for material in df.columns}
    head_nucs = {material: df[material].dropna().to_list() for material in df.columns}
    head_err_nucs = {material: ['err_'+nuc for nuc in df[material].dropna().to_list()] for material in df.columns}
    
    # Create workbook and add worksheet
    workbook = xlsxwriter.Workbook(join(dir_to_save, 'results.xlsx'))
    worksheets = {material: workbook.add_worksheet(material+'_foils') for material in df.columns}
    
    # Define cell formats
    head_fmt = workbook.add_format({'bold': True})
    num_value_fmt = workbook.add_format({'num_format': '0.00E+00'})
    datetime_fmt = workbook.add_format({'num_format': 'dd/mm/yyyy hh:mm:ss'})
    
    # Initialize row and column counters
    row = {material: 0 for material in df.columns}
    
    for material in df.columns:
        # Print headers
        worksheets[material].write_row(row[material], 0, headers + head_nucs_eval[material] + head_nucs[material] + head_err_nucs[material], head_fmt)
        row[material] += 1
        
        for t in targets:
            if t.target_material == material:
                # Restart column position counting
                col = 0
                
                # Print tager number
                worksheets[material].write_string(row[material], col, t.target_id)
                col += 1
                
                # Print irradiation datetime
                worksheets[material].write(row[material], col, t.irr_end, datetime_fmt)
                col += 1
                
                # Print proton energy
                worksheets[material].write_number(row[material], col, t.energy)
                col += 1
                
                # Print EoB activities for each radionuclide evaluate using monitor XS and BDSIM
                for r in t.radionuclides:
                    worksheets[material].write_number(row[material], col, r.act_eob_eval, num_value_fmt)
                    col += 1
                
                # Print mean EoB activity for each radionuclide
                for r in t.radionuclides:
                    worksheets[material].write_number(row[material], col, r.mean_act_eob, num_value_fmt)
                    col += 1
                    
                # Print error of mean EoB activity for each radionuclide
                for r in t.radionuclides:
                    worksheets[material].write_number(row[material], col, r.err_mean_act_eob, num_value_fmt)
                    col += 1
                
                # Increment row position counting
                row[material] += 1
                
    # Close workbook to save file
    workbook.close()