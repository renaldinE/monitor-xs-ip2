#%% Relevant packages

import numpy as np
import pandas as pd
from typing import List
from periodictable import elements

#%% Custom packages

from .core.utils import load_config

#%% Proton beam characteristics

class ProtonBeam():
    """ Proton beam using new target station with thermocouple-based beam monitor profile """
    
    def __init__(self, degrader: str):
        
        # Load configuration file
        config = load_config()
        
        # Load proton beam characteristics simulated with BDSIM
        df = pd.read_excel(config['PROTON_BEAM'], sheet_name='BeamCharacteristics')
        
        # Initialize relevant attributes
        self.energy = df.energy[df.degrader == degrader].iloc[0]
        self.err_energy = df.err_energy[df.degrader == degrader].iloc[0]
        self.current = df.current[df.degrader == degrader].iloc[0]
    
    def __str__(self):
        return f"Proton beam: {self.energy:.1f} MeV, {self.current:.1f} uA."
    
    def __repr__(self):
        return f"Proton beam: ({self.energy}, {self.current})"