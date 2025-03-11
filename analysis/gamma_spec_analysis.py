#%% Relevant packages

import numpy as np
import pandas as pd

#%% Custom packages

# from foil_analysis.radionuclide import Radionuclide, RadionuclideList
from monitor_xs_ip2.target import Target, TargetList
from monitor_xs_ip2.utils import load_config
from monitor_xs_ip2.data_plots import print_test_target_results

# Load configuration file
config = load_config()

#%% Load data

# Load target data
df = pd.read_excel(
    config['IRRADIATIONS'],
    sheet_name='Irradiations',
    usecols='A:F,I:K',
    names=['dt_start','dt_end','target_no','material','mass','thickness','integral','current','degrader'],
    dtype={'dt_start': str, 'dt_end': str}
    )

targets = TargetList(
    [Target(
        dt_start,
        dt_end,
        mass,
        tar_id,
        tar_mat,
        tar_thick,
        0, # enrichement, 0 = natural
        deg,
        current
        ) for dt_start, dt_end, mass, tar_id, tar_mat, tar_thick, deg, current in zip(
            df.dt_start,
            df.dt_end,
            df.mass,
            df.target_no,
            df.material,
            df.thickness,
            df.degrader,
            df.current)
        ]
    )

# Read reports and load gamma acquisition data into measurement objects
targets.get_acquisition_data(config['DIR_FOIL_REPORTS'], 'InterWinner')

# Calculate the activities of each radionuclide
targets.calculate_activities()

# Calculate mean activities with respect to the selected gamma energy
targets.calculate_mean_activities()

# Evaluate EoB activities for each target based on proton beam characteristics based on BDSIM simulations
targets.evaluate_eob_activites()

# Print results in an Excel file based on target material
print_test_target_results(targets, config['RESULTS'])