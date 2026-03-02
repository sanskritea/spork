import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import json
from datetime import datetime
import pandas as pd

import json
import csv
import scipy as scp
import pandas as pd
import ast
import numpy as np
import matplotlib.pyplot as plt
import emcee
import datetime
from scipy.signal import find_peaks


################################################################################################################################
############# Exponential
################################################################################################################################

def expfit(x, a, b, c):
    return ((a * np.exp(-x * b)) + c)

################################################################################################################################
################################################################################################################################


################################################################################################################################
############# BAYESIAN T1
################################################################################################################################

def calculate_M(tau_plus, tau_minus, gamma_plus, gamma_minus):
    
    G = (gamma_plus**2 + gamma_minus**2 - gamma_plus * gamma_minus) ** 0.5
    beta_plus = gamma_plus + gamma_minus + G
    beta_minus = gamma_plus + gamma_minus - G

    M_tilde_plus = (G + gamma_plus) * np.exp(-tau_plus * beta_plus) + (
        G - gamma_plus
    ) * np.exp(-tau_plus * beta_minus)
    M_tilde_minus = (G + gamma_minus) * np.exp(-tau_minus * beta_plus) + (
        G - gamma_minus
    ) * np.exp(-tau_minus * beta_minus)
    
    M_tilde_plus = M_tilde_plus / (2 * G)
    M_tilde_minus = M_tilde_minus / (2 * G)
    
    return M_tilde_plus, M_tilde_minus


def calculate_single_M(tau, gamma):
    
    G = gamma
    beta = (2 * gamma) + G

    M_tilde = (G + gamma) * np.exp(-tau * beta) / (2 * G)
    
    return M_tilde

################################################################################################################################
################################################################################################################################



################################################################################################################################
############# TWO PEAK LORENTIZAN
################################################################################################################################

def constrained_double_lorentzian(f, A1, f_center, gamma1, A2, gamma2, offset):
    f1 = f_center - 1.5e6  # -1.5 MHz from center
    f2 = f_center + 1.5e6  # +1.5 MHz from center
    
    lorentz1 = A1 * (gamma1/2)**2 / ((f - f1)**2 + (gamma1/2)**2)
    lorentz2 = A2 * (gamma2/2)**2 / ((f - f2)**2 + (gamma2/2)**2)
    
    return offset + lorentz1 + lorentz2


def double_lorentzian_fit(frequency, counts):
    
    f_center = np.mean(frequency)
    baseline = 0

    inverted_counts = -counts + baseline
    peaks, properties = find_peaks(inverted_counts, height=0.1*np.max(inverted_counts), distance=len(counts)//10)
    peak_heights = inverted_counts[peaks]
    sorted_indices = np.argsort(peak_heights)[-2:]
    selected_peaks = peaks[sorted_indices]
    amplitude_est = -np.mean(peak_heights[sorted_indices])
    linewidth_est = (np.max(frequency) - np.min(frequency)) / 10  # Conservative estimate

    p0 = amplitude_est, f_center, linewidth_est, amplitude_est, linewidth_est, baseline
    popt, pcov = scp.optimize.curve_fit(constrained_double_lorentzian, frequency, counts, p0=p0)
    param_errors = np.sqrt(np.diag(pcov))
    a1, fc, l1, a2, l2, b = p0
    
    freq_list = 1e9*np.array(np.linspace(2.77, 2.87, 2001))
    fitted_line = constrained_double_lorentzian(freq_list, a1, fc, l1, a2, l2, b)
    
    return fitted_line, freq_list

################################################################################################################################
################################################################################################################################


################################################################################################################################
############# SINGLE PEAK LORENTIZAN
################################################################################################################################

def single_lorentzian(f, A, f_center, gamma, offset):
    
    lorentz = A * (gamma/2)**2 / ((f - f_center)**2 + (gamma/2)**2)
    
    return offset + lorentz


def single_lorentzian_fit(frequency, mw_on_counts, mw_off_counts):

    freqs = frequency
    MW_ON = []
    MW_OFF = []
    for f in freqs:
        MW_ON.append(np.mean(mw_on_counts[f]))
        MW_OFF.append(np.mean(mw_off_counts[f]))
    norm_diff = -(np.array(MW_OFF) - np.array(MW_ON)) / np.mean(MW_OFF)

    f_center = np.mean(frequency)
    baseline = 0

    inverted_counts = -norm_diff + baseline
    peaks, properties = find_peaks(inverted_counts, height=0.1*np.max(inverted_counts), distance=len(norm_diff)//10)
    peak_heights = inverted_counts[peaks]
    sorted_indices = np.argsort(peak_heights)[-1:]
    selected_peaks = peaks[sorted_indices]
    amplitude_est = -np.mean(peak_heights[sorted_indices])
    linewidth_est = (np.max(frequency) - np.min(frequency)) / 5  # Conservative estimate

    p0 = amplitude_est, f_center, linewidth_est, baseline
    popt, pcov = scp.optimize.curve_fit(single_lorentzian, frequency, norm_diff, p0=p0)
    param_errors = np.sqrt(np.diag(pcov))
    a, fc, l, b = popt
    
    freq_list = np.array(np.linspace(frequency[0], frequency[-1], 501))
    fit = single_lorentzian(freq_list, a, fc, l, b)

    fc = fc*1e-9  # Convert to GHz
    param_errors[1] = param_errors[1]*1e-9  # Convert
    print(f"Fitted center frequency: {fc:.3f} ± {param_errors[1]:.3f} GHz")

    plt.figure(figsize = (6,5))
    plt.title("Pulsed ODMR")
    plt.plot(freqs*1e-9, norm_diff, 'o', color = 'black', label = 'data')
    plt.plot(freq_list*1e-9, fit, ls = '--', color = 'black', label = 'fit')
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Normalized PL (AU)')
    plt.grid(alpha=0.2)
    plt.legend()
    plt.tight_layout()
    # plt.show()
    
    return fc, fit, freq_list

################################################################################################################################
################################################################################################################################


################################################################################################################################
############# FOUR PEAK LORENTIZAN
################################################################################################################################

# Sum of four lorentzians L1 L2 L3 L4 (from left to right), L1 and L2 separated by 3MHz, L3 and L4 separated by 3MHz, and L2 and L3possibly almost overlapped; for fitting pulsed ODMR data with hyperfine splitting and possible strain splitting
def constrained_four_lorentzian(f, A1, f_center, f_center_gap, gamma1, A2, gamma2, A3, gamma3, A4, gamma4, offset):
    f_center_global = f_center
    f1 = f_center_global - f_center_gap/2
    f2 = f1 - 3e6  #  3MHz from left center
    f3 = f_center_global + f_center_gap/2
    f4 = f3 + 3e6  # 3MHz from right center

    lorentz1 = A1 * (gamma1/2)**2 / ((f - f1)**2 + (gamma1/2)**2)
    lorentz2 = A2 * (gamma2/2)**2 / ((f - f2)**2 + (gamma2/2)**2)
    lorentz3 = A3 * (gamma3/2)**2 / ((f - f3)**2 + (gamma3/2)**2)
    lorentz4 = A4 * (gamma4/2)**2 / ((f - f4)**2 + (gamma4/2)**2)

    return offset + lorentz1 + lorentz2 + lorentz3 + lorentz4

def constrained_four_lorentzian_fit(frequency, counts):
    f_center = np.mean(frequency)
    baseline = 0
    f_center_gap_est = 1e6  # 1 MHz gap between the two centers

    inverted_counts = -counts + baseline
    peaks, properties = find_peaks(inverted_counts, height=0.1*np.max(inverted_counts), distance=len(counts)//10)
    peak_heights = inverted_counts[peaks]
    sorted_indices = np.argsort(peak_heights)[-4:]
    selected_peaks = peaks[sorted_indices]
    amplitude_est = -np.mean(peak_heights[sorted_indices])
    linewidth_est = (np.max(frequency) - np.min(frequency)) / 10  # Conservative estimate

    p0 = [amplitude_est, f_center, f_center_gap_est, linewidth_est, amplitude_est, linewidth_est, amplitude_est, linewidth_est, amplitude_est, linewidth_est, baseline]
    
    popt, pcov = scp.optimize.curve_fit(constrained_four_lorentzian, frequency, counts, p0=p0)
    
    a1, fc, fc_gap, l1, a2, l2, a3, l3, a4, l4, b = popt
    
    freq_list = np.array(np.linspace(frequency[0], frequency[-1], 501))
    fitted_line = constrained_four_lorentzian(freq_list, a1, fc, fc_gap, l1, a2, l2, a3, l3, a4, l4, b)
    
    return fitted_line, freq_list, a1, fc, fc_gap, l1, a2, l2, a3, l3, a4, l4, b

################################################################################################################################
################################################################################################################################


################################################################################################################################
############# Rabi
################################################################################################################################

# def rabi(t, A, freq, phase, offset, decay):
    
#     return A * np.exp(-decay * t) * np.cos(2 * np.pi * freq * t + phase) + offset


# def rabi_fit(times, mw_on_counts, mw_off_counts):

#     times = np.array(times[1:])
#     MW_ON = []
#     MW_OFF = []

#     for t in times:
#         MW_ON.append(np.mean(mw_on_counts[t]))
#         MW_OFF.append(np.mean(mw_off_counts[t]))
#     MW_ON, MW_OFF = np.array(MW_ON), np.array(MW_OFF)
#     counts = -(np.array(MW_OFF) - np.array(MW_ON)) / np.mean(MW_OFF)
#     # counts = MW_ON
    
#     # Or estimate from peak counting
#     from scipy.signal import find_peaks
#     peaks, _ = find_peaks(counts)
#     if len(peaks) > 1:
#         avg_period = np.mean(np.diff(times[peaks]))
#         freq_estimate = 1 / avg_period
#     else:
#         freq_estimate = dominant_freq

#     p0 = [
#         (counts.max() - counts.min()) / 2,
#         freq_estimate,  # Use estimated frequency instead of 1.0
#         0,
#         counts.mean(),
#         0.0001
#     ]

#     popt, pcov = scp.optimize.curve_fit(rabi, times, counts, p0=p0)
#     perr = np.sqrt(np.diag(pcov))
    
#     A, freq, phase, offset, decay = popt
# #     print(f"Rabi Frequency: {freq:.3f} ± {perr[1]:.3f} MHz")
#     print(f"Pi-time (π-pulse): {0.5/freq:.3f} ns")
# #     print(f"Decay rate: {decay:.3f} ± {perr[4]:.3f} μs^-1")
# #     print(f"T2*: {1/decay:.3f} μs")
    
#     tau = np.linspace(times.min(), times.max(), 1001)
#     fit = rabi(tau, A, freq, phase, offset, decay)

#     plt.figure(figsize=(6, 5))
#     plt.plot(times, counts, 'o', color = 'black', label = 'data')
#     plt.plot(tau, fit, ls = '--', color = 'black', label = 'fit')
#     plt.xlabel('MW Time (ns)')
#     plt.ylabel('Counts')
#     plt.title('Rabi Oscillation Fit')
#     plt.grid(alpha=0.2)
#     plt.legend()
#     plt.show()
    
#     pi_time = np.round(0.5/freq, 0)
    
#     return pi_time, fit, tau

def rabi(t, A, freq, phase, offset, decay):
    return A * np.exp(-decay * t) * np.cos(2 * np.pi * freq * t + phase) + offset

def rabi_fit(times, mw_on_counts, mw_off_counts):
    """
    Robust Rabi oscillation fitter with better initial parameter estimation.
    """
    times = np.array(times[1:])
    MW_ON = []
    MW_OFF = []

    for t in times:
        MW_ON.append(np.mean(mw_on_counts[t]))
        MW_OFF.append(np.mean(mw_off_counts[t]))
    MW_ON, MW_OFF = np.array(MW_ON), np.array(MW_OFF)
    
    # Correct normalization (remove negative sign)
    counts = (np.array(MW_OFF) - np.array(MW_ON)) / np.mean(MW_OFF)
    
    # === BETTER FREQUENCY ESTIMATION ===
    # Method 1: Find first minimum for pi-time estimate
    first_half = min(15, len(counts)//3)  # Look in first third of data
    first_min_idx = np.argmin(counts[:first_half])
    pi_time_estimate = times[first_min_idx]
    freq_from_pi = 0.5 / pi_time_estimate if pi_time_estimate > 0 else 0.001
    
    # Method 2: FFT-based frequency estimation
    from scipy.fft import fft, fftfreq
    dt = times[1] - times[0]  # Assuming uniform spacing
    N = len(counts)
    yf = np.abs(fft(counts - counts.mean()))
    xf = fftfreq(N, dt)[:N//2]  # Frequency in 1/ns
    # Find dominant frequency (skip DC component)
    dominant_idx = np.argmax(yf[1:N//2]) + 1
    freq_from_fft = xf[dominant_idx]
    
    # Method 3: Peak detection (but be more careful)
    peaks, peak_props = find_peaks(counts, prominence=0.01)  # Add prominence threshold
    if len(peaks) > 1:
        # Use median instead of mean to be robust to outliers
        periods = np.diff(times[peaks])
        median_period = np.median(periods)
        freq_from_peaks = 1 / median_period
    else:
        freq_from_peaks = None
    
    # # Choose the most reliable estimate
    # print(f"Frequency estimates:")
    # print(f"  From first minimum: {freq_from_pi:.6f} /ns ({freq_from_pi*1000:.3f} MHz)")
    # print(f"  From FFT: {freq_from_fft:.6f} /ns ({freq_from_fft*1000:.3f} MHz)")
    # if freq_from_peaks:
    #     print(f"  From peaks: {freq_from_peaks:.6f} /ns ({freq_from_peaks*1000:.3f} MHz)")
    
    # Use FFT as primary, first-minimum as backup
    freq_estimate = freq_from_fft if freq_from_fft > 0 else freq_from_pi
    
    # === BETTER INITIAL GUESS ===
    amplitude_estimate = (counts.max() - counts.min()) / 2
    offset_estimate = counts.mean()
    
    # CRITICAL: Decay should be MUCH smaller!
    # For T2* ~ 1-10 μs, decay ~ 0.001-0.0001 per ns
    decay_estimate = 0.0001  # Conservative guess
    
    # Phase: estimate from whether we start at max or min
    if counts[0] > offset_estimate:
        phase_estimate = 0  # Starting near maximum
    else:
        phase_estimate = np.pi  # Starting near minimum
    
    p0 = [
        amplitude_estimate,
        freq_estimate,
        phase_estimate,
        offset_estimate,
        decay_estimate
    ]
    
    # print(f"\nInitial guess (p0):")
    # print(f"  A = {p0[0]:.6f}")
    # print(f"  freq = {p0[1]:.6f} /ns ({p0[1]*1000:.3f} MHz)")
    # print(f"  phase = {p0[2]:.6f} rad")
    # print(f"  offset = {p0[3]:.6f}")
    # print(f"  decay = {p0[4]:.6f} /ns")
    
    # === ADD BOUNDS TO CONSTRAIN THE FIT ===
    # This prevents the fitter from going crazy
    bounds = (
        [0,           # A must be positive
         0,           # freq must be positive
         -2*np.pi,    # phase can be anything
         -np.inf,     # offset can be anything
         0],          # decay must be positive
        [1,           # A shouldn't exceed normalized range
         0.05,        # freq shouldn't be > 50 MHz (period < 20 ns is unlikely)
         2*np.pi,     # phase upper bound
         np.inf,      # offset upper bound
         0.01]        # decay < 0.01 /ns means T2* > 100 ns
    )
    
    # === TRY FITTING WITH ROBUST SETTINGS ===
    try:
        popt, pcov = scp.optimize.curve_fit(
            rabi, times, counts, p0=p0,
            bounds=bounds,
            maxfev=20000,  # More iterations
            ftol=1e-8,     # Tighter tolerance
            xtol=1e-8
        )
        perr = np.sqrt(np.diag(pcov))
        
        A, freq, phase, offset, decay = popt
        
        # print(f"\n{'='*50}")
        # print(f"FIT RESULTS:")
        # print(f"{'='*50}")
        # print(f"Amplitude: {A:.6f} ± {perr[0]:.6f}")
        print(f"Rabi Frequency: {freq*1000:.3f} ± {perr[1]*1000:.3f} MHz")
        # print(f"Phase: {phase:.3f} ± {perr[2]:.3f} rad")
        # print(f"Offset: {offset:.6f} ± {perr[3]:.6f}")
        # print(f"Decay: {decay:.6f} ± {perr[4]:.6f} /ns")
        print(f"\nPi-time (π-pulse): {0.5/freq:.1f} ns")
        if decay > 0:
            print(f"T2*: {1/(decay*1000):.3f} μs")
        
    except RuntimeError as e:
        print(f"\n✗ Fit failed: {e}")
        print("Try adjusting bounds or collecting better data")
        return None, None, None
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return None, None, None
    
    # === PLOT RESULTS ===
    tau = np.linspace(times.min(), times.max(), 1001)
    fit = rabi(tau, *popt)
    initial_guess_curve = rabi(tau, *p0)

    # fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # # Left: Initial guess
    # axes[0].plot(times, counts, 'o', color='blue', markersize=5, label='Data')
    # axes[0].plot(tau, initial_guess_curve, '--', color='red', linewidth=2, label='Initial guess')
    # axes[0].set_xlabel('MW Time (ns)')
    # axes[0].set_ylabel('Normalized Contrast')
    # axes[0].set_title('Initial Guess')
    # axes[0].grid(alpha=0.2)
    # axes[0].legend()
    
    # # Right: Fitted result
    # axes[1].plot(times, counts, 'o', color='black', markersize=5, label='Data')
    # axes[1].plot(tau, fit, '-', color='red', linewidth=2, label='Fit')
    # axes[1].set_xlabel('MW Time (ns)')
    # axes[1].set_ylabel('Normalized Contrast')
    # axes[1].set_title(f'Rabi Fit (π-time = {0.5/freq:.1f} ns)')
    # axes[1].grid(alpha=0.2)
    # axes[1].legend()
    
    # plt.tight_layout()
    # plt.show()
    
    pi_time = np.round(0.5/freq, 0)
    
    return pi_time, fit, tau

################################################################################################################################
################################################################################################################################


################################################################################################################################
############# T2 STAR 
################################################################################################################################

def t2_star(t, A, T2star, offset):
    return A * np.exp(- (t / T2star) ** 2) + offset

def t2_star_fit(times, counts):
    
    A_est = counts.max() - counts.min()
    offset_est = counts.min()
    T2star_est = (times.max() - times.min()) / 5  # Rough estimate

    p0 = [A_est, T2star_est, offset_est]

    popt, pcov = scp.optimize.curve_fit(t2_star, times, counts, p0=p0)
    perr = np.sqrt(np.diag(pcov))
    
    A, T2star, offset = popt
    T2star_err = perr[1]
    T2star = T2star

    tau = np.linspace(times.min(), times.max(), 1001)
    fit = t2_star(tau, A, T2star, offset)

    T2star_err = np.round(perr[1]/1000, 3)  #   Convert to μs
    T2star = np.round(T2star/1000, 3)  # Convert to μs
    print('T2*:', T2star, '±', T2star_err, 'μs')
    
    return fit, tau, T2star

################################################################################################################################
################################################################################################################################

    

