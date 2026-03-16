#%% Relevant packages

import numpy as np
from typing import Iterable, List
from periodictable import elements

#%% Custom packages

from .spectrometry import efficiency_fun
from .core.utils import get_lambda_err

#%% Radionuclide and RadionuclideList classes

class Radionuclide:
    def __init__(
        self,
        name: str,
        half_life: float,
        err_half_life: float,
        uom: str,
        g_line: float):
        '''
        Class representing a radionuclide

        Parameters
        ----------
        name: str
            Radionuclide name.
        half_life: flaot
            Radionuclide's half-life.
        err_half_life: float
            Error of radionuclide's half-life.
        g_line: flaot
            Main gamma line of radionuclide.

        Returns
        -------
        None.

        '''
        self.name: str = name
        self.half_life: float = half_life
        self.err_half_life: float = err_half_life
        self.uom: str = uom
        self.g_line: float = g_line
        
        # Convert the half-life into decay constant in seconds
        self.Lambda, self.err_Lambda = get_lambda_err(half_life, err_half_life, uom)
        
        # Define later on in script
        self.g_energies: List[float] = []  # selected gamma-line energies
        self.intensities: List[float] = []  # gamma-line intensities
        self.err_intensities: List[float] = []  # errors of gamma-line intensities
        
        # ***** Additional variables *****
        # Initialization - Used only by Measurements
        self.net_counts: List[float] = []
        self.err_net_counts: List[float] = []
        
        # Activity measured at the beginning of the gamma measurement
        self.act: List[float] = []
        self.err_act: List[float] = []
        
        # Also used by Target.calculate_mean_act_eob() to store the non-zero
        # EoB activities of the main gamma line
        self.act_eob: List[float] = []  # End-of-bombardment activity
        self.err_act_eob: List[float] = []
        
        # Initialization - Used only by Targets
        self.mean_act_eob = 0
        self.err_mean_act_eob = 0
        
        # ***** Variables only for proton irradiation *****
        self.thin_target_yield = 0
        self.err_thin_target_yield = 0
    
    def calculate_activity(self, level, detector, real_time, live_time):
        '''
        Calculates the measured activity and its associated error.

        Parameters
        ----------
        C_net : list
            List of net counts.
        err_C_net : list
            List of errors of net counts.
        t_RT : float
            Real time (s).
        t_LT : float
            Live time (s).
        level : int
            Detector level.

        Returns
        -------
        None.

        '''
        if real_time <= 0 or live_time <= 0:
            raise ValueError("real_time and live_time must be positive.")

        err_real_time = err_live_time = 0.0
        
        # Correction factor due to dead time
        C_dt = real_time / live_time
        err_C_dt = C_dt * np.sqrt((err_live_time / live_time) ** 2 + (err_real_time / real_time) ** 2)
        
        # Correction factor due to radioactive decay of the source
        one_minus_exp = 1 - np.exp(-self.Lambda * real_time)
        C_meas = self.Lambda / one_minus_exp
        temp = 1 - np.exp(-self.Lambda * real_time) * (1 + self.Lambda * real_time)
        err_C_meas = temp / (one_minus_exp ** 2) * self.err_Lambda
        
        # Efficiency calculation as a function of the gamma lines energy
        if len(self.g_energies) == 0:
            raise ValueError(f"No gamma-line energies are loaded for radionuclide '{self.name}'.")

        eff, err_eff = zip(*[efficiency_fun(Ey, level, detector) for Ey in self.g_energies])
        eff = np.asarray(eff, dtype=float)
        err_eff = np.asarray(err_eff, dtype=float)
        net_counts = np.asarray(self.net_counts, dtype=float)
        err_net_counts = np.asarray(self.err_net_counts, dtype=float)
        intensities = np.asarray(self.intensities, dtype=float) / 100
        err_intensities = np.asarray(self.err_intensities, dtype=float) / 100
        
        # Activity calculation
        with np.errstate(divide="ignore", invalid="ignore"):
            self.act = C_dt * C_meas * net_counts / (eff * intensities)
        self.act = np.nan_to_num(self.act, nan=0.0, posinf=0.0, neginf=0.0)
        
        # Calculation of the error of activity
        numerator = np.array([
            err_net_counts,
            err_eff,
            err_intensities,
            ])
        
        denominator = np.array([
            net_counts,
            eff,
            intensities
            ])
        
        # Error calculation for the 1D arrays
        with np.errstate(divide="ignore", invalid="ignore"):
            error_contrib = np.sum((numerator / denominator) ** 2, axis=0)
        error_contrib = np.nan_to_num(error_contrib, nan=0.0, posinf=0.0, neginf=0.0)
        
        # Add the contributions from the scalar terms
        error_contrib += (err_C_dt / C_dt) ** 2 + (err_C_meas / C_meas) ** 2
        
        # Calculate the final error
        self.err_act = np.sqrt(error_contrib) * self.act
        
        # Handle NaN values
        self.err_act = np.nan_to_num(self.err_act, nan=0.0, posinf=0.0, neginf=0.0)

    def calculate_eob_act(self, t_cool):
        '''
        Calculates the end-of-bombardment activity

        Parameters
        ----------
        t_cool : float
            Cooling time.

        Returns
        -------
        None.

        '''
        
        # Correction due to cooling time
        if t_cool is None:
            raise ValueError("t_cool cannot be None.")

        C_cool = np.exp(self.Lambda * t_cool)
        err_C_cool = t_cool * C_cool * self.err_Lambda

        # Calculate the end-of-bombardment activity
        act = np.asarray(self.act, dtype=float)
        err_act = np.asarray(self.err_act, dtype=float)
        self.act_eob = act * C_cool
        
        # Pre-allocate error array
        self.err_act_eob = np.zeros_like(self.act_eob)

        # Calculate the error for each element
        non_zero_mask = act != 0  # Only calculate error for non-zero activity

        # Compute error only where activity is non-zero
        self.err_act_eob[non_zero_mask] = np.sqrt(
            (err_act[non_zero_mask] / act[non_zero_mask]) ** 2 +
            (err_C_cool / C_cool) ** 2
        ) * self.act_eob[non_zero_mask]
        
        # For zero activity, the error remains zero
        self.err_act_eob[~non_zero_mask] = 0.0
    
    def __eq__(self, other):
        if isinstance(other, str):
            return other == self.name
        elif isinstance(other, Radionuclide):
            return other.name == self.name
        return NotImplemented
    
    def __str__(self):
        return f"Radionuclide: {self.name} ({self.half_life} {self.uom})."
    
    def __repr__(self):
        return f"Radionuclide('{self.name}', {self.half_life})"


class RadionuclideList(list):
    def __init__(self, list_):
        # Check whether the input is a list
        if not isinstance(list_, list):
            raise ValueError(f"Input {list_} is not a list.")

        invalid_item = next((i for i in list_ if not isinstance(i, Radionuclide)), None)
        if invalid_item is not None:
            raise ValueError(f"Item {invalid_item} is not of Radionuclide type.")

        super().__init__(list_)
    
    def append(self, item):
        if isinstance(item, Radionuclide):
            super().append(item)
        else:
            raise ValueError(f'{item} is not of Radionuclide type.')
            
    def insert(self, index, item):
        if isinstance(item, Radionuclide):
            super().insert(index, item)
        else:
            raise TypeError(f"The {item} to add is not of Radionuclide type.")
    
    def extend(self, other):
        if not isinstance(other, RadionuclideList):
            raise TypeError(f"The {other} to extend is not a RadionuclideList.")
        
        for item in other:
            if not isinstance(item, Radionuclide):
                raise TypeError(f"The item {item} to extend is not a Radionuclide.")
        
        super().extend(other)
                
    @staticmethod
    def get_element_No(isotope):
        return elements.symbol(isotope.name.split('-')[0]).number
    
    def sort_elements(self, reverse=False):
        self.sort(key=self.get_element_No, reverse=reverse)
    
    def __getitem__(self, key):
        # In case I have a string
        if isinstance(key,str):
            for item in self:
                if item == key:
                    return item
            raise KeyError(f"No radionuclide found with name {key}")
            
        elif isinstance(key,Radionuclide):
            for item in self:
                if item == key.name:
                    return item
            raise KeyError(f"No radionuclide found with name {key}")
        
        # In case key not string but, e.g., a number
        return super().__getitem__(key)
