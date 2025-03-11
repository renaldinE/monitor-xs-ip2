#%% Relevant packages

import numpy as np
from numpy import matmul
from numpy.linalg import inv
import pandas as pd
from typing import List
import re
from datetime import datetime as dt
from os.path import join

#%% Custom packages

from .core.utils import get_lambda_err, time_difference, load_config

#%% Efficiency function

def efficiency_function_calibration(level,detector,p=5):
    '''
    Function that fit the efficiency function based on the selected level.
    It returns the relevant quantities to evaluate the efficiency at whichever
    gamma energy the user desires.

    Parameters
    ----------
    filepath : r-string
        File path to the Excel file containing the net counts and all the
        relevant data to perform the fitting.
    level : string
        Detector level. Example: "Rille_17" if using HPGe in Lab 35.

    Returns
    -------
    list
        [p,b,X,W,MSE].
        p : int
            Number of coefficient to estimate.
        b : list of floats
            Numeric values of the estimators of the coefficients of the linear
            regression.
        X : numpy.ndarray
            Matrix (array of arrays) containing the independent variables.
        W : numpy.ndarray
            Matrix (array of arrays) containing the weights.
        MSE : float
            Error of mean square: estimator of the sample variance.

    '''
    
    # Load configuration file
    config = load_config()
    
    if "Lab 109" in detector:
        filepath = join(config["DIR_EFF"],'Eff_Cal_EastHPGe_Lab109.xlsx')
    elif "Lab 35" in detector:
        filepath = join(config["DIR_EFF"],'Efficiency_Calibrations_OIPA_Lab35.xlsx')
    elif "West HPGe" in detector:
        filepath = join(config["DIR_EFF"],'230724_Eff_Cal_EastHPGe_Lab35.xlsx')
    
    # Load the report data of the corresponding detector level
    df = pd.read_excel(filepath,
                       sheet_name=level,
                       dtype={"ref_date": str,
                              "datetime_meas": str}
                       )
    df = df.dropna()
    
    # Reference activities of each radionuclide at the reference date
    ref_activities = df.ref_act.to_numpy() # Reference activities (Bq)
    err_ref_activities = df.err_ref_act.to_numpy() # errors or reference activities (Bq)
    
    # Conversion of the half-lives to decay constants (s^-1)
    lambdas, err_lambdas = zip(*[get_lambda_err(hf, err_hf, uom) for hf, err_hf, uom in zip(df.half_life, df.err_half_life, df.uom)])
    lambdas = np.array(lambdas)
    err_lambdas = np.array(err_lambdas)
    
    # Nuclear data of each radionuclide
    Ey_m = df.energy.to_numpy()
    Iy_m = df.intensity.to_numpy() / 100
    err_Iy_m = df.err_intensity.to_numpy() / 100
    
    # Load the net counts
    net_counts = df.net_counts.to_numpy()
    err_net_counts = df.err_net_counts.to_numpy()
    
    # Calculate the cooling time
    t_cool = [time_difference(ref_date, dt_meas) for ref_date, dt_meas in zip(df.ref_date,df.datetime_meas)]
    
    # Correction factor due to source decay
    C_cool = [np.exp( Lambda * t_c ) for Lambda,t_c in zip(lambdas, t_cool)]
    err_C_cool = [t_c * np.exp( Lambda * t_c ) * err_Lambda for Lambda,err_Lambda,t_c in zip(lambdas, err_lambdas,t_cool)]
    
    # Correction factor due to decay during measurement
    C_meas = lambdas / (1 - np.exp(-lambdas * df.t_real.values))
    err_C_meas = err_C_meas = (1 - np.exp(-lambdas * df.t_real.values) * (1 + lambdas * df.t_real.values)) / (1 - np.exp(-lambdas * df.t_real.values))**2 * err_lambdas
    
    # Correction due to detector dead time
    C_dt = df.t_real.values / df.t_live.values
    
    # Calculate activity for each gamma line
    act = net_counts * C_cool * C_meas * C_dt / Iy_m

    # Calculate err_act using numpy operations
    err_act = np.sqrt(np.sum(np.square([err_net_counts, err_C_cool, err_C_meas, err_Iy_m]) / 
                             np.square([net_counts, C_cool, C_meas, Iy_m]), axis=0)) * act
    
    # Experimental efficiency calculation
    exp_eff = act / ref_activities
    
    # Error in experimental efficiency
    err_exp_eff = np.sqrt( (err_act / act)**2 + (err_ref_activities / ref_activities)**2 ) * exp_eff
    
    # Relative error in experimental efficiency
    errRel_exp_eff = err_exp_eff / exp_eff * 100
    
    # *** Weighted Least Square method  ***
    # to calculate the parameters of the efficiency function 
    #(Chp 11 - Applied Linear Statistical Models, M. H. Kutner)
    # No. of parameters to fit, p: No. of independent variables -> POLYNOMIAL of (p-1)th GRADE
    n = Ey_m.size # No. of points sampled considering the most intense gamma lines of the nuclides selected
    
    # X: matrix of known constants, values of the predictor variable. dim{X} = (n,p)
    Ey_m_log = np.log(Ey_m)  # Compute log(Ey) for all Ey in Ey_m
    X = np.power(Ey_m_log[:, None], np.arange(p))  # Raise log(Ey) to the power of i for each i in range(p)
    Y = np.log(exp_eff) # observations
    
    # Weights matrix made of errors of experimental efficiency
    W = np.diag( 1 / (err_exp_eff / exp_eff)**2 )
    
    # Find the estimators of the regression coefficients
    XtW = np.matmul(X.T, W)
    XtWX = np.matmul(XtW, X)
    XtWX_inv = np.linalg.inv(XtWX)
    XtWY = np.matmul(XtW, Y)
    b = np.matmul(XtWX_inv, XtWY)
    
    ## Estimation of the variance-covariance matrix of the regressoin coefficients
    ## Case of error variances known
    
    res = matmul(X,b) - Y # Residuals: (y_1-y*_i), y*_i the predicted variable using the estimators of the regression coefficients  
    
    # Error sum of squares: 
    SSE = matmul( res.T, matmul(W , res) )
    # Error mean square or Residual mean square: unbiased estimator of sigma^2
    MSE = SSE / (n - p)
    
    # Calculate the weighted mean of Y
    WeightedMean_Y = np.sum(W * Y) / np.sum(np.diag(W))
    
    # Total sum of squares (SSTO): sum of distances between weighted mean and observations
    SSTO = np.sum(np.diag(W) * (Y - WeightedMean_Y)**2)
    
    # Coefficient of Determination: it express the degree of association between the variables x_i and y_i
    SqR_adj = 1-(n-1)/(n-p)*(SSE/SSTO)
    
    return [p,b,X,W,MSE]

#%% Calcualtion of coefficients for each detector and level

Lab109_detector = {
    '200 cm': efficiency_function_calibration('200 cm','OIPA Lab 109'),
    '150 cm': efficiency_function_calibration('150 cm','OIPA Lab 109'),
    '100 cm': efficiency_function_calibration('100 cm','OIPA Lab 109'),
    '50 cm': efficiency_function_calibration('50 cm','OIPA Lab 109')
    }
Lab35_detector = {
    'Rille_2': efficiency_function_calibration('Rille_2','OIPA Lab 35'),
    '25Feb_Rille_2': efficiency_function_calibration('25Feb_Rille_2','OIPA Lab 35'),
    'Rille_6': efficiency_function_calibration('Rille_6','OIPA Lab 35'),
    '240927_Rille_6': efficiency_function_calibration('240927_Rille_6','OIPA Lab 35'),
    'Rille_9': efficiency_function_calibration('Rille_9','OIPA Lab 35'),
    'Rille_17': efficiency_function_calibration('Rille_17','OIPA Lab 35'),
    '240916_Rille_17': efficiency_function_calibration('240916_Rille_17','OIPA Lab 35')
    }
WILA_detector = {}

HPGe_calibration_params = {
    'OIPA Lab 109': Lab109_detector,
    'OIPA Lab 35': Lab35_detector,
    'West HPGe - WILA': WILA_detector
    }

def efficiency_fun(Ey_s,level,detector):
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
    
    [p,b,X,W,MSE] = HPGe_calibration_params[detector][level]
    
    if isinstance(Ey_s, list) or isinstance(Ey_s,np.ndarray):
        eff_m = []
        err_eff_m = []
        for ey in Ey_s:
            X_m = np.zeros(0)
            for i in range(p):
                X_m = np.append( X_m, np.log(ey)**i )
            
            Y_m = matmul(b, X_m)
            err_Y_m = np.sqrt( MSE * matmul( X_m.T , matmul( inv( matmul( matmul(X.T, W), X ) ) , X_m ) ) )
            
            eff_m.append( float( np.exp( Y_m ) ) )
            err_eff_m.append( float( eff_m[-1] * err_Y_m ) )
    elif isinstance(Ey_s, float) or isinstance(Ey_s, int) or isinstance(Ey_s,np.int64):
        X_m = np.zeros(0)
        for i in range(p):
            X_m = np.append( X_m, np.log(Ey_s)**i )
        
        Y_m = matmul(b, X_m)
        err_Y_m = np.sqrt( MSE * matmul( np.transpose(X_m) , matmul( inv( matmul( matmul(np.transpose(X), W), X ) ) , X_m ) ) )
        
        eff_m = float( np.exp( Y_m ) )
        err_eff_m = float( eff_m * err_Y_m )
    
    return eff_m, err_eff_m

#%% Report class of measurements

class Report():
    ''' Model the spectrometry report '''
    
    def __init__(self, filepath: str):
        
        self.filepath: str = filepath
        
        self.energy: List[float] = []
        self.net_counts: List[float] = []
        self.err_net_counts: List[float] = []
        self.detector: str = None
    
    def get_report_InterWinner(self):
        
        # Extract the target ID
        self.tar_id = self.filepath.split('\\')[-1].split('.')[0].split('_')[0]
        
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
        with open(self.filepath, 'r') as file:
            report = file.read()
            
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
                
                if re.search(header_pattern, report, re.MULTILINE | re.DOTALL):
                    pattern_data = pattern_data_1
                else:
                    pattern_data = pattern_data_2
                
                if start_collecting:
                    match = re.search(pattern_data, line)
                    
                    if match:
                        energy, net, uncert = match.groups()
                        self.energy.append(float(energy))
                        self.net_counts.append(float(net))
                        # uncert is in % if you load data from InterWinner report
                        self.err_net_counts.append(float(net) * float(uncert) / 100 )
            
            # Convert to numpy
            self.energy =  np.array(self.energy)
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
        
        with open(self.filepath,'r') as file:
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
                        
    def __str__(self):
        return f"Report of '{self.filepath}'."
    
    def __repr__(self):
        return f"Report('{self.filepath}')"
