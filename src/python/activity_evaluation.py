#%% Relevant packages

import numpy as np
import pandas as pd
from typing import Union

#%% Custom packages

from .cross_section import interpolate_cross_section

#%% Constants

# Physical constants
Z = 1 # Charge of projectile, here protons
Q = 1.6021766208e-19 # elemental charge [C]
N_AVO = 6.022140857e23 # Avogadro's number [#/mol]

# Conversion factors
uA_TO_A = 10**(-6) 
MBARN_TO_CM2 = 10**(-27)
uM_TO_CM = 10**(-4)

#%% Calculation of activity for thin foils (micrometric thickness)

def evaluate_eob_activity(
        nuclide: str,
        energy: Union[float, np.ndarray],
        t_irr: float,
        current: float,
        density: float,
        foil_thickness: float,
        lambda_: float,
        mol_weight: float
        ) -> Union[float, np.ndarray]:
    '''
    Evaluate expected end-of-bombardment (EoB) activity using recommended monitor cross sections

    Parameters
    ----------
    nuclide : str
        Produced radionuclide whose cross section is used as monitor cross section.
    energy : Union[float, np.ndarray]
        Ingoing proton energy (MeV).
    t_irr : float
        Irradiation time.
    current : float
        Proton current (uA).
    density : float
        Target material density (g/cm^3).
    foil_thickness : float
        Foil thickness (um).
    lambda_ : float
        Decay constant of produced radionuclide (s^-1).
    mol_weight : float
        Molecular weight of target material.

    Returns
    -------
    act_eob : Union[float, np.ndarray]
        End-of-bombardment (EoB) activity of produced radionuclide.

    '''
    
    # Get interpolation function of cross section data
    f_xs, f_unc_xs = interpolate_cross_section(nuclide)
    
    # Convert current from micro-Ampere to Ampere
    current_conv = current * uA_TO_A
    
    # Convert foil thickness from micrometer to centimeter
    foil_thickness_conv = foil_thickness * uM_TO_CM
    
    # Evaluate EoB activity
    act_eob = f_xs(energy) * MBARN_TO_CM2 * current_conv * N_AVO * density * foil_thickness_conv * (1-np.exp(-lambda_ * t_irr)) / (Q * Z * mol_weight)
    
    return act_eob
    