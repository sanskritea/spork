import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import json
from datetime import datetime
import pandas as pd


def double_lorentzian(f, A1, f1, gamma1, A2, f2, gamma2, offset):
    """
    Double Lorentzian function for NV-ODMR fitting
    
    Parameters:
    f: frequency array
    A1, A2: amplitudes (negative for dips in ODMR)
    f1, f2: center frequencies of the two peaks
    gamma1, gamma2: linewidths (FWHM) of the two peaks
    offset: baseline offset
    """
    lorentz1 = A1 * (gamma1/2)**2 / ((f - f1)**2 + (gamma1/2)**2)
    lorentz2 = A2 * (gamma2/2)**2 / ((f - f2)**2 + (gamma2/2)**2)
    return offset + lorentz1 + lorentz2


def constrained_double_lorentzian(f, A1, f_center, gamma1, A2, gamma2, offset):
    """
    Constrained double Lorentzian with fixed 3 MHz separation

    Parameters:
    f: frequency array
    A1, A2: amplitudes
    f_center: center frequency between the two peaks
    gamma1, gamma2: linewidths of the two peaks
    offset: baseline offset
    
    The two peaks are automatically positioned at f_center ± 1.5 MHz
    """
    f1 = f_center - 1.5e6  # -1.5 MHz from center
    f2 = f_center + 1.5e6  # +1.5 MHz from center
    
    lorentz1 = A1 * (gamma1/2)**2 / ((f - f1)**2 + (gamma1/2)**2)
    lorentz2 = A2 * (gamma2/2)**2 / ((f - f2)**2 + (gamma2/2)**2)
    return offset + lorentz1 + lorentz2


def podmr_format_fit_results(fit_params, fit_errors, constrain_separation = True):
    """Format fit results for display"""
    if fit_params is None:
        return "No fit results"
    
    # lines = ["Fit Results:"]
    lines = []
    # for param, value in fit_params.items():
    #     error = fit_errors[param]
    #     if 'f' in param:  # Frequency parameters
    #         lines.append(f"{param}: {value/1e6:.3f} ± {error/1e6:.3f} MHz")
    #     elif 'gamma' in param:  # Linewidth parameters
    #         lines.append(f"{param}: {value/1e3:.1f} ± {error/1e3:.1f} kHz")
    #     # else:  # Other parameters
    #     #     lines.append(f"{param}: {value:.4f} ± {error:.4f}")
    
    if constrain_separation:
        f_center = fit_params['f_center']
        lines.append(f"Peak separation: 3.000 MHz (fixed)")
        lines.append(f"f1: {(f_center-1.5e6)/1e6:.3f} MHz")
        lines.append(f"f2: {(f_center+1.5e6)/1e6:.3f} MHz")
    else:
        f1, f2 = fit_params['f1'], fit_params['f2']
        separation = abs(f2 - f1)
        lines.append(f"Peak separation: {separation/1e6:.3f} MHz")
    
    return '\n'.join(lines)

def pulsed_odmr(constrain_separation=True):
# def pulsed_odmr(data, constrain_separation=True): #, filename):
# def pulsed_odmr(freqs_data, MW_ON_data, MW_OFF_data, constrain_separation=True): #, filename):
    """
    Initialize NV-ODMR fitter
    
    Parameters:
    constrain_separation: if True, constrains peak separation to exactly 3 MHz
    filename: json containting P-ODMR data
    """
    # print('setup')
    constrain_separation = constrain_separation
    fit_params = None
    fit_errors = None
    fit_function = None

    frequency = [] #freqs_data
    MW_ON = []
    MW_OFF = []


    # extract data from csv file
    data = pd.read_csv('C://Users/awschlab/Desktop/data/2.26um.csv')
    frequency = data['MW_ON_x'].to_numpy()*1e9
    MW_ON = data['MW_ON_y'].to_numpy()
    MW_OFF = data['MW_OFF_y'].to_numpy()

    # extract data from the json file if not passed
    # if data != False:
    #     data = pd.read_json('C://Users/awschlab/Box/ScanProbeData/250915/Pulsed ODMR/_175449.json')
    #     frequency = np.array(data['datasets']['freqs'])

    #     for f in frequency:
    #     	MW_ON.append(np.mean(data['datasets']['mwCountsDict'][str(f)]))
    #     	MW_OFF.append(np.mean(data['datasets']['noMwCountsDict'][str(f)]))

    # else:
    #     frequency, MW_ON_raw, MW_OFF_raw = data
    #     for f in frequency:
    #         MW_ON.append(np.mean(MW_ON_raw[str(f)]))
    #         MW_OFF.append(np.mean(MW_OFF_raw[str(f)]))

    # data = pd.read_json('C://Users/awschlab/Box/ScanProbeData/250917/Pulsed ODMR/_113251.json')
    # for f in frequency:
    #     MW_ON.append(np.mean(MW_ON_data[f]))
    #     MW_OFF.append(np.mean(MW_OFF_data[f]))

    counts = np.array(MW_ON) - np.array(MW_OFF)

    """
    Estimate initial parameters from the data
    """
    # Find approximate center frequency
    f_center = np.mean(frequency)
    f_center = 2.8685e9

    # Estimate baseline (assume edges are baseline)
    baseline = np.mean([np.mean(counts[:5]), np.mean(counts[:5])])
    baseline = 0
    
    # Find peaks (ODMR typically shows dips, so look for minima)
    inverted_counts = -counts + baseline
    peaks, properties = find_peaks(inverted_counts, height=0.1*np.max(inverted_counts), distance=len(counts)//10)
    
    # Estimate amplitude (depth of dips)
    if len(peaks) >= 2:
        # Use the two highest peaks
        peak_heights = inverted_counts[peaks]
        sorted_indices = np.argsort(peak_heights)[-2:]
        selected_peaks = peaks[sorted_indices]
        amplitude_est = -np.mean(peak_heights[sorted_indices])
    else:
        amplitude_est = -(np.min(counts) - baseline) / 2
    
    # Estimate linewidth (rough estimate)
    linewidth_est = (np.max(frequency) - np.min(frequency)) / 50  # Conservative estimate
    print('amplitude_est ', amplitude_est)
    print('linewidth_est ', linewidth_est)
    print('baseline ', baseline)
    
    p0 = amplitude_est, f_center, linewidth_est, amplitude_est, linewidth_est, baseline
    
    if constrain_separation:
        fit_function = constrained_double_lorentzian
        param_names = ['A1', 'f_center', 'gamma1', 'A2', 'gamma2', 'offset']
    else:
        fit_function = double_lorentzian
        param_names = ['A1', 'f1', 'gamma1', 'A2', 'f2', 'gamma2', 'offset']
    
    try:
        popt, pcov = curve_fit(fit_function, frequency, counts, p0=p0)
        param_errors = np.sqrt(np.diag(pcov))
        
        fit_params = dict(zip(param_names, popt))
        fit_errors = dict(zip(param_names, param_errors))
                    
    except Exception as e:
        print(f"Fitting failed: {e}")
        return None, None


    """
    Plot the data and fit
    """
    # print('plotting')

    if fit_params is None:
        print("No fit available. Run fit() first.")
        return
    
    # if data != False:
        # Generate fitted curve
    f_fine = np.linspace(np.min(frequency), np.max(frequency), 1000)
    if constrain_separation:
        fitted_counts = constrained_double_lorentzian(f_fine, *list(fit_params.values()))
    else:
        fitted_counts = double_lorentzian(f_fine, *list(fit_params.values()))


    fig = plt.figure(figsize=(10, 6))
    plt.plot(frequency/1e6, counts, 'bo', label='Data', alpha=0.7)
    plt.plot(f_fine/1e6, fitted_counts, 'r-', label='Double Lorentzian Fit', linewidth=2)
    
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Normalized Fluorescence')
    plt.title('Fitted Pulsed ODMR')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Add fit parameters as text
    textstr = podmr_format_fit_results(fit_params, fit_errors)
    plt.text(0.02, 0.98, textstr, transform=plt.gca().transAxes, 
            verticalalignment='bottom', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    filetitle = 'C://Users/awschlab/Desktop/data/fits/' + str(datetime.now().strftime("%Y_%m_%d_%H_%M_%S")) + '.png'
    # plt.savefig(filetitle)
    plt.show()

# else:
    f_center = fit_params['f_center']
    f_L = (f_center - 1.5e6) / 1e9
    f_U = (f_center + 1.5e6) / 1e9

    print('f_center ', f_center)
    return f_center, fig

        # return f_L, f_U 


##############################################################################################
##############################################################################################


def constrained_quadruple_lorentzian(f, A1, f_center_1, gamma1, A2, gamma2, A3, f_center_2, gamma3, A4, gamma4, offset):
    """
    Constrained double Lorentzian with fixed 3 MHz separation

    Parameters:
    f: frequency array
    A1, A2: amplitudes
    f_center: center frequency between the two peaks
    gamma1, gamma2: linewidths of the two peaks
    offset: baseline offset
    
    The two peaks are automatically positioned at f_center ± 1.5 MHz
    """
    f1 = f_center_1 - 1.5e6  # -1.5 MHz from center
    f2 = f_center_1 + 1.5e6  # +1.5 MHz from center

    f3 = f_center_2 - 1.5e6  # -1.5 MHz from center
    f4 = f_center_2 + 1.5e6  # +1.5 MHz from center
    
    lorentz1 = A1 * (gamma1/2)**2 / ((f - f1)**2 + (gamma1/2)**2)
    lorentz2 = A2 * (gamma2/2)**2 / ((f - f2)**2 + (gamma2/2)**2)
    lorentz3 = A3 * (gamma3/2)**2 / ((f - f3)**2 + (gamma3/2)**2)
    lorentz4 = A4 * (gamma4/2)**2 / ((f - f4)**2 + (gamma4/2)**2)

    return offset + lorentz1 + lorentz2 + lorentz3 + lorentz4


def four_podmr_format_fit_results(fit_params, fit_errors, constrain_separation = True):
    """Format fit results for display"""
    if fit_params is None:
        return "No fit results"
    
    # lines = ["Fit Results:"]
    lines = []
    # for param, value in fit_params.items():
    #     error = fit_errors[param]
    #     if 'f' in param:  # Frequency parameters
    #         lines.append(f"{param}: {value/1e6:.3f} ± {error/1e6:.3f} MHz")
    #     elif 'gamma' in param:  # Linewidth parameters
    #         lines.append(f"{param}: {value/1e3:.1f} ± {error/1e3:.1f} kHz")
    #     # else:  # Other parameters
    #     #     lines.append(f"{param}: {value:.4f} ± {error:.4f}")
    
    if constrain_separation:
        f_center_1 = fit_params['f_center_1']
        lines.append(f"Peak separation: 3.000 MHz (fixed)")
        lines.append(f"f1: {(f_center_1-1.5e6)/1e6:.3f} MHz")
        lines.append(f"f2: {(f_center_1+1.5e6)/1e6:.3f} MHz")

        f_center_2 = fit_params['f_center_2']
        lines.append(f"Peak separation: 3.000 MHz (fixed)")
        lines.append(f"f3: {(f_center_2-1.5e6)/1e6:.3f} MHz")
        lines.append(f"f4: {(f_center_2+1.5e6)/1e6:.3f} MHz")
    else:
        f1, f2 = fit_params['f1'], fit_params['f2']
        separation = abs(f2 - f1)
        lines.append(f"Peak separation: {separation/1e6:.3f} MHz")
    
    return '\n'.join(lines)


def four_pulsed_odmr(constrain_separation=True): #, filename):
# def pulsed_odmr_four_peaks(freqs_data, MW_ON_data, MW_OFF_data, constrain_separation=True): #, filename):
    """
    Initialize NV-ODMR fitter
    
    Parameters:
    constrain_separation: if True, constrains peak separation to exactly 3 MHz
    filename: json containting P-ODMR data
    """
    # print('setup')
    constrain_separation = constrain_separation
    fit_params = None
    fit_errors = None
    fit_function = None

    # frequency = freqs_data
    MW_ON = []
    MW_OFF = []


    # extract data from csv file
    data = pd.read_csv('C://Users/awschlab/Desktop/data/700nm.csv')
    frequency = data['MW_ON_x'].to_numpy()*1e9
    MW_ON = data['MW_ON_y'].to_numpy()
    MW_OFF = data['MW_OFF_y'].to_numpy()


    # # extract data from the json file if not passed
    # if data != False:
    #     data = pd.read_json('C://Users/awschlab/Box/ScanProbeData/250915/Pulsed ODMR/_175449.json')
    #     frequency = np.array(data['datasets']['freqs'])

    #     for f in frequency:
    #       MW_ON.append(np.mean(data['datasets']['mwCountsDict'][str(f)]))
    #       MW_OFF.append(np.mean(data['datasets']['noMwCountsDict'][str(f)]))

    # else:
    #     frequency, MW_ON_raw, MW_OFF_raw = data
    #     for f in frequency:
    #         MW_ON.append(np.mean(MW_ON_raw[str(f)]))
    #         MW_OFF.append(np.mean(MW_OFF_raw[str(f)]))


    # # data = pd.read_json('C://Users/awschlab/Box/ScanProbeData/250917/Pulsed ODMR/_113251.json')
    # for f in frequency:
    #     MW_ON.append(np.mean(MW_ON_data[f]))
    #     MW_OFF.append(np.mean(MW_OFF_data[f]))

    counts = np.array(MW_ON) - np.array(MW_OFF)
    # print('counts ', counts)

    """
    Estimate initial parameters from the data
    """
    # Find approximate center frequency
    f_center = np.mean(frequency)
    # print('f_center ', f_center)
    
    # Estimate baseline (assume edges are baseline)
    baseline = np.mean([np.mean(counts[:5]), np.mean(counts[-5:])])
    
    # Find peaks (ODMR typically shows dips, so look for minima)
    inverted_counts = -counts + baseline
    peaks, properties = find_peaks(inverted_counts, height=0.1*np.max(inverted_counts), distance=len(counts)//10)
    
    # Estimate amplitude (depth of dips)
    if len(peaks) >= 4:
        # Use the two highest peaks
        peak_heights = inverted_counts[peaks]
        sorted_indices = np.argsort(peak_heights)[-4:]
        selected_peaks = peaks[sorted_indices]
        amplitude_est = -np.mean(peak_heights[sorted_indices])
    else:
        amplitude_est = -(np.min(counts) - baseline) / 2
    
    # Estimate linewidth (rough estimate)
    linewidth_est = (np.max(frequency) - np.min(frequency)) / 50  # Conservative estimate

    f_center_1 = 2.8685e9 #f_center - 0.001e9
    f_center_2 = 2.8735e9 #f_center + 0.001e9

    p0 = amplitude_est, f_center_1, linewidth_est, amplitude_est, linewidth_est, amplitude_est, f_center_2, linewidth_est, amplitude_est, linewidth_est, baseline
    
    if constrain_separation:
        fit_function = constrained_quadruple_lorentzian
        param_names = ['A1', 'f_center_1', 'gamma1', 'A2', 'gamma2', 'A3', 'f_center_2', 'gamma3', 'A4', 'gamma4', 'offset']
    else:
        fit_function = double_lorentzian
        param_names = ['A1', 'f_center_1', 'gamma1', 'A2', 'gamma2', 'A3', 'f_center_2', 'gamma3', 'A4', 'gamma4', 'offset']
    
    try:
        popt, pcov = curve_fit(fit_function, frequency, counts, p0=p0, maxfev = 5000)
        param_errors = np.sqrt(np.diag(pcov))
        
        fit_params = dict(zip(param_names, popt))
        fit_errors = dict(zip(param_names, param_errors))
                    
    except Exception as e:
        print(f"Fitting failed: {e}")
        return None, None


    """
    Plot the data and fit
    """
    # print('plotting')

    if fit_params is None:
        print("No fit available. Run fit() first.")
        return
    
    # if data != False:
        # Generate fitted curve
    f_fine = np.linspace(np.min(frequency), np.max(frequency), 1000)
    if constrain_separation:
        fitted_counts = constrained_quadruple_lorentzian(f_fine, *list(fit_params.values()))
    else:
        fitted_counts = double_lorentzian(f_fine, *list(fit_params.values()))

    print('f1 ', fit_params['f_center_1'] - 1.5e6)
    fig = plt.figure(figsize=(10, 6))
    plt.plot(frequency/1e6, counts, 'bo', label='Data', alpha=0.7)
    plt.plot(f_fine/1e6, fitted_counts, 'r-', label='Quadruple Lorentzian Fit', linewidth=2)
    
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Normalized Fluorescence')
    plt.title('Fitted Pulsed ODMR')
    plt.legend()
    plt.grid(True, alpha=0.3)

    
    # Add fit parameters as text
    textstr = four_podmr_format_fit_results(fit_params, fit_errors)
    plt.text(0.02, 0.98, textstr, transform=plt.gca().transAxes, 
            verticalalignment='bottom', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    filetitle = 'C://Users/awschlab/Desktop/data/fits/' + str(datetime.now().strftime("%Y_%m_%d_%H_%M_%S")) + '.png'
    # plt.savefig(filetitle)
    plt.show()

# else:
    f_center_1 = fit_params['f_center_1']
    f_1 = (f_center_1 - 1.5e6) / 1e9
    f_2 = (f_center_1 + 1.5e6) / 1e9

    f_center_2 = fit_params['f_center_2']
    f_3 = (f_center_2 - 1.5e6) / 1e9
    f_4 = (f_center_2 + 1.5e6) / 1e9

    print('f_center_1 ', f_center_1)

    return f_center_1, fig

##############################################################################################
##############################################################################################


def rabi_fitting(t, a, b, c):
    return a + b * np.cos(np.pi * c * t)


def rabi_format_fit_results(popt):
    """Format fit results for display"""
    
    # lines = ["Fit Results:"]
    pi_time = 1 / popt[2]
    lines = []
    lines.append(f"Pi pulse duartion: {pi_time:.3f} ns")
    
    return '\n'.join(lines)


def rabi(mw_times_data, MW_ON_data, MW_OFF_data):
# def rabi():

    # data = pd.read_csv('C://Users/awschlab/Desktop/data/rabi_20250923.csv')
    # times = data['MW_ON_x'].to_numpy()
    # MW_ON = data['MW_ON_y'].to_numpy()
    # MW_OFF = data['MW_OFF_y'].to_numpy()

    # data = pd.read_json('C://Users/awschlab/Box/ScanProbeData/250917/Rabi/_105335.json')
    # times = np.array(data['datasets']['mw_times'])
    # MW_ON = []
    # MW_OFF = []
    # for t in times:
    #     MW_ON.append(np.mean(data['datasets']['mwCountsDict'][str(t)]))
    #     MW_OFF.append(np.mean(data['datasets']['noMwCountsDict'][str(t)]))

    times = mw_times_data
    MW_ON = []
    MW_OFF = []
    for t in times:
    	MW_ON.append(np.mean(MW_ON_data[t]))
    	MW_OFF.append(np.mean(MW_OFF_data[t]))

    counts = (np.array(MW_ON) - np.array(MW_OFF))
    counts = (counts - ((np.max(counts) + np.min(counts)) / 2)) / ((np.max(counts) - np.min(counts)) / 2)
    # print('counts ', counts)
    # counts = MW_ON
    a_est = 0.01
    b_est =  counts[0] #1
    c_est = 0 # 1 / 80

    c_est_temp = []
    for i in range(2, len(times) - 2):
        if counts[i] < 0 and counts[i+2] > counts[i] and counts[i+1] > counts[i] and counts[i-1] > counts[i] and counts[i-2] > counts[i]:
            c_est_temp.append(times[i])

    print('c_est_temp ', c_est_temp)

    c_est = 1 / c_est_temp[0]
    print('c_est ', c_est)
    p0 = a_est, b_est, c_est 

    popt, pcov = curve_fit(rabi_fitting, times, counts, p0)
    pi_time = 1 / popt[2]
    print('pi time ', pi_time)

    tau = np.linspace(np.min(times), np.max(times), 1000)
    fit = rabi_fitting(tau, popt[0], popt[1], popt[2]) #, popt[3])
    fig = plt.figure(figsize=(10, 6))
    # print('popt[0] ', popt[0])
    # print('popt[1] ', popt[1])
    # print('popt[2] ', popt[2])
    plt.plot(times, counts, 'o-', color = 'black', label='data')
    plt.plot(tau, fit, 'r-', label='fit')
    plt.xlabel('tau (ns)')
    plt.ylabel('Counts')
    plt.title('Fitted Rabi')
    plt.legend()
    plt.grid()
    plt.show()

    textstr = rabi_format_fit_results(popt)
    plt.text(0.02, 0.98, textstr, transform=plt.gca().transAxes, 
            verticalalignment='bottom', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # filetitle = 'C://Users/awschlab/Desktop/data/fits/' + str(datetime.now().strftime("%Y_%m_%d_%H_%M_%S")) + '.png'
    # plt.savefig(filetitle)

    # return pi_time, fig




##############################################################################################
##############################################################################################


# # Example usage and test data generation
# if __name__ == "__main__":

# 	# pulsed_odmr()
# 	rabi()
    