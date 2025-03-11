#%% Relevant packages

import numpy as np
import pandas as pd
from os import listdir
from os.path import isfile,join
import subprocess

#%% Custom packages

from .radionuclide import Radionuclide, RadionuclideList
from .measurement import Measurement, MeasurementList, get_datetime_meas
from .spectrometry import Report
from .proton_beam import ProtonBeam
from .activity_evaluation import evaluate_eob_activity
from .core.utils import time_difference, load_config

config = load_config()

#%% Target class

class Target():
    ''' Model an irradiated target bombarded with either protons or neutrons '''
    
    def __init__(
        self,
        irr_start: str,
        irr_end: str,
        mass: float,
        tar_id: str,
        tar_mat: str,
        tar_thick: float,
        enr: float,
        degrader: str,
        current: float):
        '''

        Parameters
        ----------
        irr_start : str
            Date and time of irradiation start.
        irr_end : str
            Date and time of irradiation end.
        mass : float
            Target mass.
        tar_id : str
            Target identifier.
        tar_mat : str
            Target material.
        tar_thick : float
            Target thickness (um).
        enr : float
            Enrichment (%).
        degrader : str
            Degrader tag.
        current : float
            Proton current extracted from main beam with the beam splitter.

        Raises
        ------
        FileNotFoundError
            Error raised when Excel file of nuclear data not found.
        Exception
            Return error.

        Returns
        -------
        None.

        '''

        self.irr_start: str = irr_start
        self.irr_end: str = irr_end
        self.mass: float = mass
        self.target_id: str = tar_id
        self.target_material: str = tar_mat
        self.target_thick: float = tar_thick
        self.enr: float = enr
        self.degrader: str = degrader
        self.t_irr: float = time_difference(self.irr_start, self.irr_end)
        
        # Proton beam characteristics derived from BDSIM simulations
        proton_beam = ProtonBeam(self.degrader)
        self.energy = proton_beam.energy
        self.actual_current = current / 50 * proton_beam.current
        
        # Check file existence before proceeding
        try:
            df = pd.read_excel(config['GL_FILEPATH'], sheet_name="NuclideData", usecols="A:E")
        except FileNotFoundError:
            raise FileNotFoundError(f"File {config['GL_FILEPATH']} not found.")
        except Exception as e:
            raise Exception(f"Error reading the Excel file: {str(e)}")
            
        # Load radionuclides present in a specific target material
        df_nuc_inventory = pd.read_excel(config['MATERIALS_DATA'], sheet_name='Material_NuclideInventory')
        df_mol_weights = pd.read_excel(config['MATERIALS_DATA'], sheet_name='MolecularWeights')
        
        # Define molecular weight of target material
        self.mol_weight = df_mol_weights['mol_weight'][df_mol_weights['material'] == self.target_material].iloc[0]
        
        # Define target material density
        self.density = df_mol_weights['density'][df_mol_weights['material'] == self.target_material].iloc[0]

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
        
        # Load the gamma lines of each radionuclide
        for nuc in self.radionuclides:
            gamma_lines = pd.read_excel(config['GL_FILEPATH'],sheet_name=nuc.name)
            self.radionuclides[nuc].g_energies = gamma_lines.energy.to_numpy()
            self.radionuclides[nuc].intensities = gamma_lines.intensity.to_numpy()
            self.radionuclides[nuc].err_intensities = gamma_lines.err_intensity.to_numpy()
         
        # **** Initialize useful attributes ****
        # Cooling time
        self.t_cool = None
        self.act_eob_eval = 0
        
        # Initialize the MeasurementList
        self.measurements = MeasurementList([])
    
    def load_measurements(self, dir_reports, software):
        
        # Load file paths of gamma measurement reports
        try:
            report_list = [join(dir_reports, f) for f in listdir(dir_reports) if isfile(join(dir_reports, f))]
        except FileNotFoundError as error:
            print(f"Error: {error}")
            print("Attempting to map the Genie drive...")
        
            try:
                # Map the network drive
                subprocess.run(['net', 'use', 'T:', r'\\fs02\Genie2kLCH'], check=True, shell=True)
                print("Drive T: mapped successfully. Retrying file loading...")
                
                # Retry loading file paths after mapping the drive
                report_list = [join(dir_reports, f) for f in listdir(dir_reports) if isfile(join(dir_reports, f))]
            except subprocess.CalledProcessError as e:
                print(f"Failed to map drive: {e}")
                report_list = []  # Default to an empty list if drive mapping fails
        
        for filepath in report_list:
            
            # Read report file
            report = Report(filepath)
            if software == 'InterWinner':
                report.get_report_InterWinner()
            elif software == 'Genie2K':
                report.get_report_Genie2K()
            else:
                raise ValueError(f"{software} is not 'Genie2K' or 'InterWinner'.")
            
            if self.target_id == report.tar_id:
                self.measurements.append(
                    Measurement(
                        report.datetime_meas,
                        report.detector_level,
                        report.live_time,
                        report.real_time,
                        report.detector,
                        report.tar_id,
                        self.target_material)
                    )
                
                # Store the net counts in each Radionuclide instance
                self.measurements[-1].get_net_counts(report)
        
        # Sort the measurements based on the acquisition date and time
        self.measurements.sort(key=get_datetime_meas)
        
    def calculate_mean_act_eob(self):
        '''
        Calculates the mean EoB activity of each radionuclide

        Parameters
        ----------
        measurements : list of Measurement
            List of Measurement variables.

        Returns
        -------
        None.

        '''
        
        # Store only valid EoB activity results in the radionuclide list of Target class
        for m in self.measurements:
            for r in m.radionuclides:
                self.radionuclides[r].act_eob.append(float(r.act_eob[r.g_energies == r.g_line]))
                self.radionuclides[r].err_act_eob.append(float(r.err_act_eob[r.g_energies == r.g_line]))
        
        for r in self.radionuclides:
            self.radionuclides[r].act_eob = np.array(self.radionuclides[r].act_eob)
            self.radionuclides[r].err_act_eob = np.array(self.radionuclides[r].err_act_eob)
        
        # Now the target.radionuclides[i].act_eob is a list of non-zero act_eob
        for r in self.radionuclides:
            if r.act_eob.size == 0:
                # The mean_act_eob is already initialized to 0
                continue
            elif r.act_eob.size == 1:
                # Only 1 result, no reason to perform weighed mean
                self.radionuclides[r].mean_act_eob = r.act_eob[0]
                self.radionuclides[r].err_mean_act_eob = r.err_act_eob[0]
            else:
                # Calculate the weights = inverse of error squared
                weights = 1 / r.err_act_eob**2
                
                # Perform weighted average of the EoB activities
                self.radionuclides[r].mean_act_eob = np.average(
                    r.act_eob,
                    weights=weights
                    )
                
                # Calculate the error
                self.radionuclides[r].err_mean_act_eob = np.sqrt(1/weights.sum())
    
    def calculate_tty(self):
        '''
        Calculates the thin-tagret yield (MBq/uA).

        Returns
        -------
        None.

        '''
        
        proton_beam = ProtonBeam(self.degrader, self.mass)
        corr_factor = self.proton_current / 50
        
        for i,nuclide in enumerate(self.radionuclides):
            self.radionuclides[i].thin_target_yield = nuclide.mean_act_eob / ( corr_factor * proton_beam.current )
            self.radionuclides[i].err_thin_target_yield = nuclide.err_mean_act_eob / ( corr_factor * proton_beam.current )
                
    def __str__(self):
        return f"Target: {self.target_id} -> {self.target_material} ({self.degrader})."
        
    def __repr__(self):
        return f"Target: ('{self.target_material}', '{self.degrader}', '{self.target_id}')"
        
    def __eq__(self,tar_id_to_test):
        if tar_id_to_test == self.target_id:
            return True
        else:
            return False

class TargetList(list):
    def __init__(self, list_):
        # Check whether the input is a list
        if isinstance(list_,list):
            # Check whether all the list items are of GammaMeasurement type
            ck = True
            for i in list_:
                if not isinstance(i,Target):
                    ck = False
                    break
                
            if ck:
                super().__init__(list_)
            else:
                raise ValueError(f"Item {i} is not of Target type.")
        else:
            raise ValueError(f"Input {list_} is not a list.")
            
    def get_acquisition_data(self, dir_reports: str, software: str):
        '''
        Load all data of gamma spectrometry for each target

        Parameters
        ----------
        dir_reports : str
            Directory where to find the txt files of reports.
        software : str
            Software used to analyze gamma spectra: Genie2K or InterWinner.

        Returns
        -------
        None.

        '''
        
        for target in self:
            # Initialize gamma measurement objects and net counts
            target.load_measurements(dir_reports, software)
            
            for i, _ in enumerate(target.measurements):
                self[target].measurements[i].calculate_cooling_time(target.irr_end)
    
    def calculate_activities(self):
        ''' Calculate measured activity and EoB activity '''
        for t in self:
            for i, m in enumerate(t.measurements):
                for r in m.radionuclides:
                    # Calculate activity at the beginning of acquisition
                    self[t].measurements[i].radionuclides[r].calculate_activity(
                        m.level,
                        m.detector,
                        m.real_time,
                        m.live_time
                        )
                    
                    # Calculate EoB activity
                    self[t].measurements[i].radionuclides[r].calculate_eob_act(m.t_cool)
    
    def calculate_mean_activities(self):
        ''' Evaluate mean activities for each target and radionuclide '''
        
        for t in self:
            self[t].calculate_mean_act_eob()
            
    def evaluate_eob_activites(self):
        
        for t in self:
            for r in t.radionuclides:
                self[t].radionuclides[r].act_eob_eval = evaluate_eob_activity(
                    r.name,
                    t.energy,
                    t.t_irr,
                    t.actual_current,
                    t.density,
                    t.target_thick,
                    r.Lambda,
                    t.mol_weight
                    )
    
    def print_eob_activities(self):
        for t in self:
            print(t)
            for m in t.measurements:
                print(m)
                for r in m.radionuclides:
                    print(r)
                    print(r.act_eob)
    
    def append(self, item):
        if isinstance(item, Target):
            super().append(item)
        else:
            raise ValueError(f'{item} is not of Target type.')
    
    def __getitem__(self, key):
        # In case I have a string
        if isinstance(key,str) | isinstance(key,Target):
            for item in self:
                if item == key:
                    return item
                
            raise KeyError(f"No target found with ID {key}")
        # In case key not string but, e.g., a number -> classic indexing of list
        return super().__getitem__(key)