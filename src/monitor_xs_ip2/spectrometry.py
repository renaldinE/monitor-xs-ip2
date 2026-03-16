#%% Relevant packages

import numpy as np
from numpy import matmul
from numpy.linalg import inv
import pandas as pd
from typing import Dict, List, Tuple
import re
from datetime import datetime as dt
from os.path import join
from dataclasses import dataclass

#%% Custom packages

from .core.utils import get_lambda_err, time_difference, load_config

#%% Efficiency function

@dataclass(frozen=True)
class CalibrationFit:
    p: int
    b: np.ndarray
    X: np.ndarray
    W: np.ndarray
    MSE: float


def resolve_efficiency_filepath(detector: str, config: dict) -> str:
    if "Lab 109" in detector:
        return join(config["DIR_EFF"], "Eff_Cal_EastHPGe_Lab109.xlsx")
    if "Lab 35" in detector:
        return join(config["DIR_EFF"], "Efficiency_Calibrations_OIPA_Lab35.xlsx")
    if "West HPGe" in detector:
        return join(config["DIR_EFF"], "230724_Eff_Cal_EastHPGe_Lab35.xlsx")
    raise ValueError(f"Unknown detector '{detector}'.")


def load_efficiency_sheet(filepath: str, level: str) -> pd.DataFrame:
    df = pd.read_excel(
        filepath,
        sheet_name=level,
        dtype={"ref_date": str, "datetime_meas": str},
    ).dropna()
    if df.empty:
        raise ValueError(f"No efficiency calibration data found for level '{level}' in '{filepath}'.")
    return df


def compute_decay_constants(df: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
    lambdas, err_lambdas = zip(
        *[get_lambda_err(hf, err_hf, uom) for hf, err_hf, uom in zip(df.half_life, df.err_half_life, df.uom)]
    )
    return np.array(lambdas, dtype=float), np.array(err_lambdas, dtype=float)


def compute_correction_factors(df: pd.DataFrame, lambdas: np.ndarray, err_lambdas: np.ndarray) -> Dict[str, np.ndarray]:
    t_cool = np.array(
        [time_difference(ref_date, dt_meas) for ref_date, dt_meas in zip(df.ref_date, df.datetime_meas)],
        dtype=float,
    )
    t_real = df.t_real.to_numpy(dtype=float)
    t_live = df.t_live.to_numpy(dtype=float)
    if np.any(t_real <= 0) or np.any(t_live <= 0):
        raise ValueError("t_real and t_live must be strictly positive.")

    c_cool = np.exp(lambdas * t_cool)
    err_c_cool = t_cool * c_cool * err_lambdas

    one_minus_exp = 1 - np.exp(-lambdas * t_real)
    c_meas = lambdas / one_minus_exp
    err_c_meas = (1 - np.exp(-lambdas * t_real) * (1 + lambdas * t_real)) / (one_minus_exp ** 2) * err_lambdas

    c_dt = t_real / t_live
    return {
        "c_cool": c_cool,
        "err_c_cool": err_c_cool,
        "c_meas": c_meas,
        "err_c_meas": err_c_meas,
        "c_dt": c_dt,
    }


def compute_activity_and_uncertainty(df: pd.DataFrame, corrections: Dict[str, np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
    iy = df.intensity.to_numpy(dtype=float) / 100
    err_iy = df.err_intensity.to_numpy(dtype=float) / 100
    net_counts = df.net_counts.to_numpy(dtype=float)
    err_net_counts = df.err_net_counts.to_numpy(dtype=float)

    with np.errstate(divide="ignore", invalid="ignore"):
        act = net_counts * corrections["c_cool"] * corrections["c_meas"] * corrections["c_dt"] / iy
        rel_err_sq = np.sum(
            np.square([err_net_counts, corrections["err_c_cool"], corrections["err_c_meas"], err_iy]) /
            np.square([net_counts, corrections["c_cool"], corrections["c_meas"], iy]),
            axis=0,
        )

    act = np.nan_to_num(act, nan=0.0, posinf=0.0, neginf=0.0)
    rel_err_sq = np.nan_to_num(rel_err_sq, nan=0.0, posinf=0.0, neginf=0.0)
    err_act = np.sqrt(rel_err_sq) * act
    return act, err_act


def compute_experimental_efficiency(
    df: pd.DataFrame, act: np.ndarray, err_act: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    ref_activities = df.ref_act.to_numpy(dtype=float)
    err_ref_activities = df.err_ref_act.to_numpy(dtype=float)

    with np.errstate(divide="ignore", invalid="ignore"):
        exp_eff = act / ref_activities
        err_exp_eff = np.sqrt((err_act / act) ** 2 + (err_ref_activities / ref_activities) ** 2) * exp_eff

    exp_eff = np.nan_to_num(exp_eff, nan=0.0, posinf=0.0, neginf=0.0)
    err_exp_eff = np.nan_to_num(err_exp_eff, nan=0.0, posinf=0.0, neginf=0.0)
    return exp_eff, err_exp_eff


def filter_valid_points(
    energies: np.ndarray,
    exp_eff: np.ndarray,
    err_exp_eff: np.ndarray,
    p: int,
    level: str,
    detector: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    if np.any(energies <= 0):
        raise ValueError("All gamma energies must be > 0 to evaluate logarithms.")

    valid = (exp_eff > 0) & (err_exp_eff > 0)
    valid_count = np.count_nonzero(valid)
    if valid_count <= p:
        raise ValueError(
            f"Not enough valid calibration points for level '{level}' and detector '{detector}'. "
            f"Need more than p={p}, found {valid_count}."
        )
    return energies[valid], exp_eff[valid], err_exp_eff[valid]


def fit_weighted_log_polynomial(
    energies: np.ndarray, exp_eff: np.ndarray, err_exp_eff: np.ndarray, p: int
) -> CalibrationFit:
    n = energies.size
    X = np.power(np.log(energies)[:, None], np.arange(p))
    Y = np.log(exp_eff)
    W = np.diag(1 / (err_exp_eff / exp_eff) ** 2)

    xtw = np.matmul(X.T, W)
    xtwx = np.matmul(xtw, X)
    b = np.matmul(np.linalg.inv(xtwx), np.matmul(xtw, Y))

    res = matmul(X, b) - Y
    sse = matmul(res.T, matmul(W, res))
    mse = float(sse / (n - p))

    return CalibrationFit(p=p, b=b, X=X, W=W, MSE=mse)


def efficiency_function_calibration(level, detector, p=5) -> CalibrationFit:
    """Fit detector efficiency coefficients for a given level and detector."""
    config = load_config()
    filepath = resolve_efficiency_filepath(detector, config)
    df = load_efficiency_sheet(filepath, level)
    lambdas, err_lambdas = compute_decay_constants(df)
    corrections = compute_correction_factors(df, lambdas, err_lambdas)
    act, err_act = compute_activity_and_uncertainty(df, corrections)
    exp_eff, err_exp_eff = compute_experimental_efficiency(df, act, err_act)

    energies = df.energy.to_numpy(dtype=float)
    energies, exp_eff, err_exp_eff = filter_valid_points(energies, exp_eff, err_exp_eff, p, level, detector)
    return fit_weighted_log_polynomial(energies, exp_eff, err_exp_eff, p)

#%% Calcualtion of coefficients for each detector and level

_CALIBRATION_LEVELS = {
    'OIPA Lab 109': ['200 cm', '150 cm', '100 cm', '50 cm'],
    'OIPA Lab 35': ['Rille_2', '25Feb_Rille_2', 'Rille_6', '240927_Rille_6', 'Rille_9', 'Rille_17', '240916_Rille_17'],
    'West HPGe - WILA': [],
}

HPGe_calibration_params = {detector: {} for detector in _CALIBRATION_LEVELS}

def _get_calibration_params(level: str, detector: str) -> CalibrationFit:
    if detector not in HPGe_calibration_params:
        raise KeyError(f"Detector '{detector}' is not configured.")
    if level not in _CALIBRATION_LEVELS.get(detector, []):
        raise KeyError(f"Level '{level}' is not configured for detector '{detector}'.")

    if level not in HPGe_calibration_params[detector]:
        HPGe_calibration_params[detector][level] = efficiency_function_calibration(level, detector)
    return HPGe_calibration_params[detector][level]


def efficiency_fun(Ey_s, level, detector):
    '''
    Function that returns the efficiency and its error at a given energy
    or as a function of a gamma energy range

    Parameters
    ----------
    Ey_s : int, float or list
        Gamma energy/ies where to evaluate the detector efficiency.

    Returns
    -------
    eff_m : flot or list of floats
        Efficiency value.
    err_eff_m : float or list of floats
        Error of the efficiency values.

    '''
    
    fit = _get_calibration_params(level, detector)
    p, b, X, W, MSE = fit.p, fit.b, fit.X, fit.W, fit.MSE
    xtwx_inv = inv(matmul(matmul(X.T, W), X))
    
    if isinstance(Ey_s, list) or isinstance(Ey_s, np.ndarray):
        eff_m = []
        err_eff_m = []
        for ey in Ey_s:
            if ey <= 0:
                raise ValueError("Gamma energy must be > 0.")
            X_m = np.zeros(0)
            for i in range(p):
                X_m = np.append(X_m, np.log(ey) ** i)
            
            Y_m = matmul(b, X_m)
            err_Y_m = np.sqrt(MSE * matmul(X_m.T, matmul(xtwx_inv, X_m)))
            
            eff_m.append(float(np.exp(Y_m)))
            err_eff_m.append(float(eff_m[-1] * err_Y_m))
    elif isinstance(Ey_s, float) or isinstance(Ey_s, int) or isinstance(Ey_s, np.int64):
        if Ey_s <= 0:
            raise ValueError("Gamma energy must be > 0.")
        X_m = np.zeros(0)
        for i in range(p):
            X_m = np.append(X_m, np.log(Ey_s) ** i)
        
        Y_m = matmul(b, X_m)
        err_Y_m = np.sqrt(MSE * matmul(X_m.T, matmul(xtwx_inv, X_m)))
        
        eff_m = float(np.exp(Y_m))
        err_eff_m = float(eff_m * err_Y_m)
    else:
        raise TypeError(f"Unsupported Ey_s type: {type(Ey_s)}")
    
    return eff_m, err_eff_m

#%% Report class of measurements

class Report:
    ''' Model the spectrometry report '''
    
    def __init__(self, filepath: str):
        
        self.filepath: str = filepath
        
        self.energy: List[float] = []
        self.net_counts: List[float] = []
        self.err_net_counts: List[float] = []
        self.detector: str = None
    
    def get_report_InterWinner(self):
        
        # Extract the target ID
        self.tar_id = self.filepath.split('\\')[-1].split('.')[0].split('_')[0].lower()
        
        # Detector
        self.detector = 'OIPA Lab 35'
        
        # Data patterns to match
        acq_data_patterns = [r'Acq\.time \(live\):\s+([\d\.]+)\s+s',
                             r'Acq\.time \(real\):\s+([\d\.]+)\s+s',
                             r'Acquisition date:\s+([\d\.:\s]+)',
                             r'\s*(Rille\s\d+)']
        # r'Comment:\s*(?:[^\n]*)\n\s*(Rille\s\d+)'
        
        keys_data_patterns = ['live_time',
                              'real_time',
                              'datetime_meas',
                              'detector_level']
        
        # Check when to start collecting
        start_collecting = False
        
        # Pattern data
        # When there is NETTO-NULL
        pattern_data_1 = r'\|\s*\d+\|[\\\/\|\s]*(\d+\.\d+)\s*\|[\s\d\.\d\*]*\|[\d\.\s\*]+?\|\s*\d+\|\s*(\d+\.\d*)\s*\|[\d\.\s]+?\|\s*(\d+\.\d+).*$'
        # When there is NET only
        pattern_data_2 = r'\|\s*\d+\|[\\\/\|\s]*(\d+\.\d+)\s*\|[\s\d\.\d\*]*\|[\d\.\s\*]+?\|\s*\d+\|\s*(\d+\.\d*)\s*\|\s*(\d+\.\d+).*$'
        
        # Header pattern to check is there is NET only or NET-NUL
        header_pattern = r'\|No\.\|\s+Energy\s+\|[\s+]?FWHM\s+\|[\s+]?FWTM\s+\|\s+GROSS\s+\|\s+NET\s+\|\s+NETTO-NUL\s+\|UNCERT\[\%\]\|\s+EFF\.\[\%\]\s+\|\s+ISOTOPE\s+\|.*$'
        
        # Load report file in binary mode (should be faster)
        with open(self.filepath, 'r', encoding='utf-8', errors='ignore') as file:
            report = file.read()
            has_netto_null = re.search(header_pattern, report, re.MULTILINE | re.DOTALL) is not None
            pattern_data = pattern_data_1 if has_netto_null else pattern_data_2
            
            for line in report.splitlines():
                for pattern, key in zip(acq_data_patterns, keys_data_patterns):
                    match = re.search(pattern, line)
                    if match:
                        value = match.group(1).strip()
                        
                        if 'date' in key:
                            try:
                                parsed_date = dt.strptime(value, '%d.%m.%Y %H:%M:%S')
                                parsed_date = dt.strftime(parsed_date, "%Y-%m-%d %H:%M:%S")
                                setattr(self, key, parsed_date)
                            except ValueError as e:
                                print(f"Error parsing date '{value}': {e}")
                        
                        elif '_time' in key:
                            setattr(self, key, float(value))
                        
                        elif 'detector' in key:
                            setattr(self, key, value.replace(' ', '_'))

                if 'List by energies (with candidate isotopes)' in line:
                    start_collecting = True
                    continue
                elif 'List by energies (with confirmed isotopes)' in line:
                    break
                
                if start_collecting:
                    match = re.search(pattern_data, line)
                    
                    if match:
                        energy, net, uncert = match.groups()
                        self.energy.append(float(energy))
                        self.net_counts.append(float(net))
                        # uncert is in % if you load data from InterWinner report
                        self.err_net_counts.append(float(net) * float(uncert) / 100 )
            
            # Convert to numpy
            self.energy = np.array(self.energy)
            self.net_counts = np.array(self.net_counts)
            self.err_net_counts = np.array(self.err_net_counts)
    
    def get_report_Genie2K(self):
        '''
        Extracts the relevant quantities from the Genie2000 report

        Parameters
        ----------
        filepath : raw str
            File path to the report.txt location.

        Returns
        -------
        None.

        '''
        # Data patterns to match
        acq_data_patterns = [r'Sample ID\s+:\s(\S+)\s+$',
                             r'Detector level\s+:\s(.+?\scm)\s+$',
                             r'Detector\s+:\s(.*?)\s*$',
                             r'Acquisition Time\s+:\s+(\d{2}\.\d{2}\.\d{4} \d{2}:\d{2}:\d{2})\s*$',
                             r'Live Time\s+:\s+([\d.]+)\sseconds$',
                             r'Real Time\s+:\s+([\d.]+)\sseconds$']
        
        keys_data_patterns = ['tar_id',
                              'detector_level',
                              'detector',
                              'datetime_meas',
                              'live_time',
                              'real_time']
        
        header_data = r'\s*Peak\s+ROI\s+ROI\s+Peak\s+Energy\s+FWHM\s+Net Peak\s+Net Area\s+Continuum\s*$'
        
        with open(self.filepath, 'r', encoding='utf-8', errors='ignore') as file:
            report = file.read()
            
            for pattern, key in zip(acq_data_patterns,keys_data_patterns):
                match = re.findall(pattern,report, re.MULTILINE)
                if match:
                    value = match[0]
                    if key == 'datetime_meas':
                        try:
                            parsed_date = dt.strptime(value, '%d.%m.%Y %H:%M:%S')
                            parsed_date = dt.strftime(parsed_date, "%Y-%m-%d %H:%M:%S")
                            setattr(self, key, parsed_date)
                        except ValueError as e:
                            print(f"Error parsing date '{value}': {e}")
                    elif '_time' in key:
                        setattr(self, key, float(value))
                    elif key == 'detector':
                        self.detector = value
                    elif key == 'tar_id':
                        setattr(self, key, value.lower())
                    else:
                        setattr(self, key, value)
                else:
                    raise ValueError(f"Value for '{key}' not found.")
            
            # Store the data for each detected peak
            pattern_data = r'\s*[mM]?\s*\d+\s+\d+-\s*\d+\s+[\d.]+\s+([\d.]+)\s+[\d.]+\s+([\d.E+]+)\s+([\d.]+)\s+[\d.E+]+\s*$'
            
            # Check when to start collecting
            start_collecting = False
            
            for line in report.splitlines():
                
                if re.search(header_data, line):
                    start_collecting = True
                    continue
                
                if start_collecting:
                    match = re.search(pattern_data, line)
                    
                    if match:
                        energy, net, uncert = match.groups()
                        self.energy.append(float(energy))
                        self.net_counts.append(float(net))
                        self.err_net_counts.append(float(uncert))

        # Convert to numpy for consistency with InterWinner parser.
        self.energy = np.array(self.energy)
        self.net_counts = np.array(self.net_counts)
        self.err_net_counts = np.array(self.err_net_counts)
                        
    def __str__(self):
        return f"Report of '{self.filepath}'."
    
    def __repr__(self):
        return f"Report('{self.filepath}')"
