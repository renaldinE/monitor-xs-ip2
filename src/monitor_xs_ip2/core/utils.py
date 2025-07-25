#%% Relevant packages

import numpy as np
import pandas as pd
import os
from math import log
from datetime import datetime as dt
import json
from typing import Tuple

#%% Constants

#%% Time difference function

def time_difference(date_hour_in: str, date_hour_end: str) -> float:
    '''
    Function that calculates the elapsed time between two dates and respective hours

    Parameters
    ----------
    date_hour_in : string
        format = 'dd-mm-yy hh:mm:ss'.
    date_hour_end : string
        format = 'dd-mm-yy hh:mm:ss'.

    Returns
    -------
    elapsed_time : float
        Elapsed time (s).

    '''
    
    try:
        tmp_date_hour_in = dt.strptime(date_hour_in,'%Y-%m-%d %H:%M:%S')
    except ValueError as msg:
        print(msg)
        tmp_date_hour_in = dt.strptime(date_hour_in,'%Y-%m-%d %H:%M:%S.%f')
        
    try:
        tmp_date_hour_end = dt.strptime(date_hour_end,'%Y-%m-%d %H:%M:%S')
    except ValueError as msg:
        print(msg)
        tmp_date_hour_end = dt.strptime(date_hour_end,'%Y-%m-%d %H:%M:%S.%f')
    
    elapsed_time = tmp_date_hour_end - tmp_date_hour_in
    
    return float(elapsed_time.total_seconds())

#%% Function to calculate decay constant and its error

def get_lambda_err(hf: float, err_hf: float, uom: str) -> Tuple[float, float]:
    """Converts the half-life to decay constant in s^-1.

    Args:
        hf (float): half-life.
        err_hf (float): error of half-life.
        uom (str): unit of measure of half-life.

    Raises:
        ValueError: wrong unit of measure

    Returns:
        tuple: (decay constant, error of decay constant)
    """
    
    if uom == "s":
        tmp_hf = hf
        tmp_err_hf = err_hf
    elif uom == "m":
        tmp_hf = hf * 60
        tmp_err_hf = err_hf * 60
    elif uom == "h":
        tmp_hf = hf * 60 * 60
        tmp_err_hf = err_hf * 60 * 60
    elif uom == "d":
        tmp_hf = hf * 60 * 60 * 24
        tmp_err_hf = err_hf * 60 * 60 * 24
    elif uom == "y":
        tmp_hf = hf * 60 * 60 * 24 * 365
        tmp_err_hf = err_hf * 60 * 60 * 24 * 365
    else:
        raise ValueError('No unit of measure for the half-life.')
        
    Lambda = log(2) / tmp_hf
    err_Lambda = Lambda**2 * tmp_err_hf / log(2)
    
    return Lambda, err_Lambda

#%% Load configuration JSON file

def load_config() -> dict:
    ''' Load configuration JSON file '''
    
    # Root directory
    root_directory = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
    
    # Define the path to the configuration file located in the .config folder
    config_file = os.path.join(root_directory, '.config', 'config.json')
     
    # Ensure the config file exists before attempting to open it
    if not os.path.isfile(config_file):
        raise FileNotFoundError(f"Configuration file not found: {config_file}")
    
    with open(config_file, 'r') as file:
        config = json.load(file)
    
    # Define the project root directory
    project_root = os.path.dirname(os.path.dirname(config_file))
    
    # Convert relative paths to absolute paths by combining with project root
    for key in config.keys():
        config[key] = os.path.join(project_root, config[key])
        
    return config