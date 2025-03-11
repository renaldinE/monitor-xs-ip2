#%% Relevant packages

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from os import listdir
from os.path import join
import re

#%% Custom packages

from .core.utils import load_config

#%% Monitor cross sections

pattern = r'[a-zA-Z]+\-\d+'

def load_monitor_cross_section(nuc: str):
    '''
    

    Parameters
    ----------
    nuc : str
        Radionuclide whose monitor cross section was selected to perform the measurement.

    Raises
    ------
    ValueError
        Error related to wrong input or empty CSV file.
    FileNotFoundError
        CSV file of monitor cross section not found.
    IOError
        Generic error occurring while reading CSV file.

    Returns
    -------
    pandas.DataFrame
        Data: energy, data, uncert.

    '''
    
    if not re.fullmatch(pattern, nuc):
        raise ValueError(f"Input nuclide name '{nuc}' does not match the required input pattern: (Chemical element)-(Mass No.).") 
    
    # Load configuration file
    config = load_config()
    
    # Get file paths of monitor cross sections
    filepaths_xs = {f.split('.')[0]: join(config['DIR_XS_DATA'],f) for f in listdir(config['DIR_XS_DATA'])}
    
    # Check if the nuclide exists in the directory
    if nuc not in filepaths_xs:
        raise FileNotFoundError(f"Monitor cross section data for nuclide '{nuc}' not found.")
    
    try:
        # Attempt to read the CSV file
        return pd.read_csv(filepaths_xs[nuc])
    except pd.errors.EmptyDataError:
        raise ValueError(f"CSV file for nuclide '{nuc}' is empty.")
    except pd.errors.ParserError:
        raise ValueError(f"Error parsing the CSV file for nuclide '{nuc}'.")
    except Exception as e:
        raise IOError(f"An error occurred while reading the file for nuclide '{nuc}'.") from e

#%% Function to linearly interpolate cross section function

def interpolate_cross_section(nuclide: str):
    ''' Interpolate tabulated data of Padé approximants of monitor cross sections '''
    
    # Load cross section -> pandas.DataFrame: energy, data, uncert
    cross_section = load_monitor_cross_section(nuclide)
    
    # Interpolate cross section 
    f_xs = interp1d(cross_section.energy, cross_section.data, kind='linear')
    f_unc_xs = interp1d(cross_section.energy, cross_section.uncert, kind='linear')
    
    return f_xs, f_unc_xs