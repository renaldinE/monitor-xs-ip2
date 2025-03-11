#%% Relevant packages

from io import BytesIO
import scipy
import pandas as pd
import os
import re
from pathlib import Path
from srim import Ion, TRIM, Target

#%% Custom packages


#%% Super class based on initialized Transmit class of PySRIM, but never implemented

class Transmit():
    """ Read TRANSMIT.txt file generate by pysrim TRIM.run() """
    
    def __init__(self, directory, filename='TRANSMIT.txt'):
        '''reads the file named TRANSMIT.txt in SRIM Outputs folder'''
        
        try:
            with open(os.path.join(directory, filename), 'rb') as f:
                output = f.read()
        except FileNotFoundError as e:
            print(e)
            
            # User might forgot to add folder 'SRIM Output' to the file path
            directory_corrected = os.path.join(directory, 'SRIM Outputs')
            print('Attempting searching for TRANSMIT.txt file in ' + directory_corrected + ' ...')
            with open(os.path.join(directory_corrected, filename), 'rb') as f:
                output = f.read()
            
            print('\nTRANSMIT.txt file found.')
                
        self._input_data = self._read_trim_inputs(output)
        self._output_data = self._read_data(output)

    def _read_trim_inputs(self, output):
        ''' Example line to read from the file:
            ====== TRIM Calc.=  H(10 MeV) ==> Ti_1+Ni_1(  50 um) ========================= '''
        pattern_trim_calc = r'=+\s+TRIM\s+Calc\.=\s+([a-zA-Z]+)\((\d+)\s*([a-zA-Z]+)\)\s+==>\s+(.+?)\s*\(\s+(\d+)\s+([a-zA-Z]+)\)'
        match = re.search(pattern_trim_calc.encode('utf-8'), output)
        
        if match:
            out_dict = {
                'projectile': match.group(1),
                'beam_energy': (float(match.group(2)), match.group(3)), # (absolute value, unit of measure)
                'layer_name': match.group(4),
                'layer_thickness': (float(match.group(5)), match.group(6)) # (absolute value, unit of measure)
                }
            
        return out_dict
    
    def _read_data(self, output):
        ''' TRANSMIT file header:
        ============================== SRIM-2013.00 ==============================
        ==============================================================================
        =================  TRANSMIT.txt : File of Transmitted Ions  ==================
        =  This file tabulates the kinetics of ions or atoms leaving the target.     =
        =  Column #1: S= Sputtered Atom, B= Backscattered Ion, T= Transmitted Ion.   =
        =  Col.#2: Ion Number, Col.#3: Z of atom leaving, Col.#4: Atom energy (eV).  =
        =  Col.#5-7: Last location:  X= Depth into target, Y,Z= Transverse axes.     =
        =  Col.#8-10: Cosines of final trajectory.                                   =
        = *** This data file is in the same format as TRIM.DAT (see manual for uses).=
        ====== TRIM Calc.=  <Projectile atom>(<Abs. value energy> MeV) ==> <Layers' names>(  <total thickness of layers> <unit of measure>) =========================
         Ion  Atom   Energy        Depth       Lateral-Position        Atom Direction      
         Numb Numb    (eV)          X(A)        Y(A)       Z(A)      Cos(X)  Cos(Y) Cos(Z) '''
        
        # Last header pattern
        pattern_headers = r'\s*Numb\s+Numb\s+\(eV\)\s+X\(A\)\s+Y\(A\)\s+Z\(A\)\s+Cos\(X\)\s+Cos\(Y\)\s+Cos\(Z\)\s*'
        
        # Find location of first data line
        headers_match = re.search(pattern_headers.encode('utf-8'), output)
        start_idx = headers_match.end()
        
        # Extract raw data
        rawdata = BytesIO(output[start_idx:]).read().decode('utf-8')
        
        # Process the data
        data = [list(re.split(r'\s+', line.strip())) for line in rawdata.split('\r\n') if len(line.strip()) > 0]
        
        df_output = pd.DataFrame(
            data,
            columns=['particle_type', 'ion_number', 'z_leaving_atom', 'atom_energy', 'depth', 'Y_axis', 'Z_axis', 'cos_x', 'cos_y', 'cos_z']
            )
        
        return df_output

#%% Run srim

def run_pysrim(energy: float, layers: list, number_ions: int, srim_executable_directory_user=None):
    '''
    

    Parameters
    ----------
    energy : float
        Energy of the ion beam impinging on the surface of the first layer in MeV.
    layers : list
        List of Layer objects constituting the Target object of the SRIM simulation.
    number_ions : int
        Number of ions to simulate.
    srim_executable_directory_user : r-str, optional
        Absolute file path to the directory of SRIM executable. The default is None.

    Raises
    ------
    an
        DESCRIPTION.
    FileNotFoundError
        SRIM directory not found.

    Returns
    -------
    transmit_file : Transmit
        Data collected from TRANSMIT.txt file.

    '''
    
    # Define ion beam
    ion = Ion('H', energy=energy * 10**6)
    
    # Construct the target
    target = Target(layers)
    
    # Initialize a TRIM calculation
    trim = TRIM(target, ion, number_ions=number_ions, calculation=1,transmit=True) 
    
    if not srim_executable_directory_user:
        # Get absolute path to the "Documents" folder
        if os.name == 'nt':  # For Windows
            srim_directory = Path(os.environ['USERPROFILE']) / 'Documents' / 'SRIM-2013'
        else:  # For Linux/macOS
            srim_directory = Path.home() / 'Documents' / 'SRIM-2013'
    else:
        srim_directory = srim_executable_directory_user
    
    # Check if the SRIM directory exists and raise an error if not
    if not os.path.isdir(srim_directory):
        raise FileNotFoundError(f"The recommended location of SRIM executable was not found ({srim_directory}).\nPlease, specify SRIM location using 'srim_executable_directory_user' input variable.")
    
    # Run SRIM
    trim.run(srim_directory)
    
    # Parse and load TRANSMIT.txt file to get the energies of the transmitted protons
    transmit = Transmit(srim_directory)
    
    return transmit

#%% Calculate average energy

def get_energy_out(transmit):
    '''
    Fit the energy distribution of transmitted projectile particles

    Parameters
    ----------
    transmit : Transmit
        Transmit output file of SRIM simulation.

    Returns
    -------
    energy_out : float
        Mean energy fitted on Gaussian distribution of the transmitted projectiles energy in MeV.
    err_energy_out : float
        Error of the mean energy in MeV.

    '''
    
    if transmit._output_data.atom_energy.size:
        energy_out = None
        err_energy_out = None
        
        message = 'No transmitted ' + transmit._input_data['projectile'] + \
            ' with energy ' + transmit._input_data['beam_energy'][0] + ' ' + \
            transmit._input_data['beam_energy'][1] + ' through the ' + \
            transmit._input_data['layer_name'] + ' layer.'
            
        print(message)
    else:
        energy_out, err_energy_out = scipy.stats.norm.fit(
            transmit._output_data.atom_energy * 10**(-6)
            )
    
    return float(energy_out), float(err_energy_out)