#%% Relevant packages

import numpy as np
import pandas as pd
from typing import List
from datetime import datetime as dt

#%% Custom packages

from .radionuclide import Radionuclide, RadionuclideList
from .core.utils import time_difference, load_config

config = load_config()

#%% Functions related to Measurement class

def get_datetime_meas(measurement):
    ''' Function to help sorting the measurements list '''
    if isinstance(measurement, Measurement):
        return dt.strptime(measurement.meas_date, "%Y-%m-%d %H:%M:%S")
    else:
        raise TypeError(f'{measurement} is not a Measurement class.')

#%% Measurement class

class Measurement():
    ''' Model a spectrometry measurement '''
    
    def __init__(
        self,
        meas_date: str,
        level: str,
        live_time: float,
        real_time: float,
        detector: str,
        tar_id: str,
        tar_mat: str):

        '''

        Parameters
        ----------
        meas_date: str
            Date and time of gamma ray spectrum acquisition.
        level: str
            Detector level.
        live_time: float
            Live time of the acquisition.
        real_time: float
            Real time of the acquisition
        detector: str
            Detector name used for the gamma ray acquisition.
        tar_id: str
            Target ID.
        tar_mat: str
            Target material.

        Returns
        -------
        None.

        '''

        self.meas_date: str = meas_date
        self.level: str = level
        self.live_time: float = live_time
        self.real_time: float = real_time
        self.detector: str = detector
        self.target_id: str = tar_id
        
        # Check file existence before proceeding
        try:
            df = pd.read_excel(config['GL_FILEPATH'], sheet_name="NuclideData", usecols="A:E")
        except FileNotFoundError:
            raise FileNotFoundError(f"File {config['GL_FILEPATH']} not found.")
        except Exception as e:
            raise Exception(f"Error reading the Excel file: {str(e)}")
            
        # Load radionuclides present in a specific target material
        df_nuc_inventory = pd.read_excel(config['MATERIALS_DATA'], sheet_name='Material_NuclideInventory')
        
        # Initialize radionuclide list
        self.radionuclides = RadionuclideList(
            [Radionuclide(
                nuc,
                hf,
                err_hf,
                uom,
                Ey_s) for nuc, hf, err_hf, uom, Ey_s in zip(
                        df.nuclide,
                        df.half_life,
                        df.err_half_life,
                        df.uom,
                        df.selected_g_line) if nuc in df_nuc_inventory[tar_mat].to_list()
                ]
            )
        
        # Load gamma lines for each radionuclide
        for radionuclide in self.radionuclides:
            try:
                gamma_lines = pd.read_excel(config['GL_FILEPATH'], sheet_name=radionuclide.name)
                radionuclide.g_energies = gamma_lines.energy.to_numpy()
                radionuclide.intensities = gamma_lines.intensity.to_numpy()
                radionuclide.err_intensities = gamma_lines.err_intensity.to_numpy()
            except FileNotFoundError:
                raise FileNotFoundError(f"Gamma line data for {radionuclide.name} not found in the Excel file.")
            except Exception as e:
                raise Exception(f"Error reading gamma line data for {radionuclide.name}: {str(e)}")
            
    
    def calculate_cooling_time(self,datetime_irr_end):
        '''
        Calculates the cooling time.

        Parameters
        ----------
        targets : list of Target variables
            List of targets.

        Returns
        -------
        None.

        '''
        self.t_cool = time_difference(datetime_irr_end, self.meas_date)
        self.err_t_cool = 2 # (s)
    
    def get_net_counts(self, report):
        '''
        Retrieve the relevant quantities from the report file

        Parameters
        ----------
        report : Report
            Report variable containining the net counts and associated errors.

        Returns
        -------
        None.

        '''
        
        for r in self.radionuclides: # go through the radionuclide list
            for Ey_s in r.g_energies: # go through the selected gamma lines for each radionuclide
                ck_Ey_s_found = False # to know if the gamma line was found or not
                
                # To store only the net counts of the detected peak closest to the gamma line
                net_counts_matches = np.array([])
                err_net_counts_matches = np.array([])
                deviations = np.array([])
                
                for Ey_m,net_counts,err_net_counts in zip(report.energy,
                                                          report.net_counts,
                                                          report.err_net_counts):
                    if (Ey_m < Ey_s + 1.) and ( Ey_m > Ey_s - 1.): # ck if the measured gamma line is among one of the selected
                        
                        net_counts_matches = np.append(net_counts_matches, net_counts)
                        err_net_counts_matches = np.append(err_net_counts_matches, err_net_counts)
                        deviations = np.append(deviations, abs(Ey_m - Ey_s) )
                        
                        ck_Ey_s_found = True
                
                if len(deviations) == 1:
                    self.radionuclides[r].net_counts.append( net_counts_matches[0] )
                    self.radionuclides[r].err_net_counts.append( err_net_counts_matches[0] )
                elif len(deviations) > 1:
                    self.radionuclides[r].net_counts.append( float(net_counts_matches[ deviations == deviations.min() ] ) )
                    self.radionuclides[r].err_net_counts.append( float(err_net_counts_matches[ deviations == deviations.min() ] ) )
                    
                if not ck_Ey_s_found: # in case not found, append a Null activity
                    self.radionuclides[r].net_counts.append( 0.0 )
                    self.radionuclides[r].err_net_counts.append( 0.0 )
        
            # Convert net counts into numpy array
            self.radionuclides[r].net_counts = np.array(self.radionuclides[r].net_counts)
            self.radionuclides[r].err_net_counts = np.array(self.radionuclides[r].err_net_counts)
    
    def __eq__(self,datetime_to_compare):
        # Check wether 'datetime_to_compare' is of type datetime.datetime
        if isinstance(datetime_to_compare, dt):
            datetime_meas = dt.strptime(self.datetime_meas,"%Y-%m-%d %H:%M:%S")
            if datetime_to_compare == datetime_meas:
                return True
            else:
                return False
        elif isinstance(datetime_to_compare, str):
            datetime_meas = dt.strptime(self.datetime_meas,"%Y-%m-%d %H:%M:%S")
            datetime_to_compare_dt = dt.strptime(datetime_to_compare,"%Y-%m-%d %H:%M:%S")

            if datetime_to_compare_dt == datetime_meas:
                return True
            else:
                return False
    
    def __str__(self):
        return f"Measurement of {self.meas_date} (Target: {self.target_id})."
    
    def __repr__(self):
        return f"Measurement('{self.meas_date}', '{self.target_id}')"

#%% MeasurementList class

class MeasurementList(list):
    def __init__(self, list_):
        # Check whether the input is a list
        if isinstance(list_, list):
            # Check whether all the list items are of Measurement type
            if not all(isinstance(i, Measurement) for i in list_):
                raise ValueError("All items in the list must be of Measurement type.")
            super().__init__(list_)
        else:
            raise ValueError(f"Input {list_} is not a list.")

    def append(self, item):
        if isinstance(item, Measurement):
            super().append(item)
        else:
            raise ValueError(f'{item} is not of Measurement type.')