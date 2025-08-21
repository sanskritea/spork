import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import json
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


def pulsed_odmr(constrain_separation=True, data): #, filename):
    """
    Initialize NV-ODMR fitter
    
    Parameters:
    constrain_separation: if True, constrains peak separation to exactly 3 MHz
    filename: json containting P-ODMR data
    """
    print('setup')
    constrain_separation = constrain_separation
    fit_params = None
    fit_errors = None
    fit_function = None

    frequency = []
    MW_ON = []
    MW_OFF = []

    # extract data from the json file if not passed
    if data != False:
        data = pd.read_json('C://Users/awschlab/Desktop/data/250803/Pulsed ODMR/82G_300nm_Pulsed ODMR231406final.json')
        frequency = np.array(data['datasets']['freqs'])

        for f in frequency:
        	MW_ON.append(np.mean(data['datasets']['mwCountsDict'][str(f)]))
        	MW_OFF.append(np.mean(data['datasets']['noMwCountsDict'][str(f)]))

    else:
        frequency, MW_ON_raw, MW_OFF_raw = data
        for f in frequency:
        MW_ON.append(np.mean(MW_ON_raw[str(f)]))
        MW_OFF.append(np.mean(MW_OFF_raw[str(f)]))

    counts = np.array(MW_ON) - np.array(MW_OFF)

    """
    Estimate initial parameters from the data
    """
    # Find approximate center frequency
    f_center = np.mean(frequency)
    
    # Estimate baseline (assume edges are baseline)
    baseline = np.mean([np.mean(counts[:5]), np.mean(counts[-5:])])
    
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
    
    # if constrain_separation:
    #     return [amplitude_est, f_center, linewidth_est, amplitude_est, linewidth_est, baseline]
    # else:
    #     f1_est = f_center - 1.5e6
    #     f2_est = f_center + 1.5e6
    #     # return [amplitude_est, f1_est, linewidth_est, amplitude_est, f2_est, linewidth_est, baseline]

    print('fitting')

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
    print('plotting')

    if fit_params is None:
        print("No fit available. Run fit() first.")
        return
    
    if data != False:
        # Generate fitted curve
        f_fine = np.linspace(np.min(frequency), np.max(frequency), 1000)
        if constrain_separation:
            fitted_counts = constrained_double_lorentzian(f_fine, *list(fit_params.values()))
        else:
            fitted_counts = double_lorentzian(f_fine, *list(fit_params.values()))


        plt.figure(figsize=(10, 6))
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
        plt.show()

    else:
        f_center = fit_params['f_center']
        f_L = (f_center - 1.5e6) / 1e9
        f_U = (f_center + 1.5e6) / 1e9
        
        return f_L, f_U 


##############################################################################################
##############################################################################################

def rabi_fitting(t, a, b, c):
    return a + b * np.cos(np.pi * c * t)


def rabi():

    data = pd.read_json('C://Users/awschlab/Desktop/data/250803/Rabi/H_82G_300nm_Rabi230918final.json')
    times = np.array(data['datasets']['mw_times'])
    MW_ON = []
    MW_OFF = []

    for t in times:
    	MW_ON.append(np.mean(data['datasets']['mwCountsDict'][str(t)]))
    	MW_OFF.append(np.mean(data['datasets']['noMwCountsDict'][str(t)]))

    counts = (np.array(MW_ON) - np.array(MW_OFF))

    a_est = (np.max(counts) - np.min(counts)) / 2
    b_est = counts[0]
    c_est = 1 / times[np.argmin(counts)]
    # d_est = 1
    p0 = a_est, b_est, c_est #, d_est

    popt, pcov = curve_fit(rabi_fitting, times, counts, p0)
    print('pi time ', 1 / popt[2])

    tau = np.linspace(np.min(times), np.max(times), 1000)
    fit = rabi_fitting(tau, popt[0], popt[1], popt[2])#, popt[3])

    plt.plot(times, counts, 'o-', color = 'black', label='data')
    plt.plot(tau, fit, 'r-', label='fit')
    plt.xlabel('tau (ns)')
    plt.ylabel('Counts')
    plt.title('Fitted Rabi')
    plt.legend()
    plt.grid()
    plt.show()




##############################################################################################
##############################################################################################


# # Example usage and test data generation
# if __name__ == "__main__":

# 	# pulsed_odmr()
# 	rabi()
    