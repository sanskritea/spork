import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg, interpolate, constants
import time

# Material and condition parameters
def calculate_ms(delta_h=1):
    return ((2870/2.8) - (4-2)/2) / 2 - delta_h

# Core parameters
M0 = 1716 #calculate_ms()  # Oe
L = 3.0  # μm
DD = 5.4e-9 * 1e8  # Oe μm^-2
gamma = 2.8e-3  # GHz/G
omega_M = gamma * M0  # GHz
print(f"M0 = {M0} Oe")

def F(q, n):
    return 2 * (1 - (-1)**n * np.exp(-q)) / q

def P(q, n, m):
    kronecker = 1.0 if n == m else 0.0
    sqrt_term = np.sqrt((1 + (1.0 if n == 0 else 0)) * (1 + (1.0 if m == 0 else 0)))
    first_term = q**2 / (q**2 + n**2 * np.pi**2) * kronecker
    second_term = 1/sqrt_term * q**4 / ((q**2 + n**2 * np.pi**2) * (q**2 + m**2 * np.pi**2)) * F(q, n) * (1 + (-1)**(n+m)) / 2
    return first_term - second_term

def Q(q, n, m):
    sqrt_term = np.sqrt((1 + (1.0 if n == 0 else 0)) * (1 + (1.0 if m == 0 else 0)))
    first_factor = q**2 / (q**2 + m**2 * np.pi**2)
    second_factor = (m**2 / (m**2 - n**2 + (1 + (-1)**(n+m)) / 2) * 2 / q - 
                    q**2 / (2 * (q**2 + n**2 * np.pi**2)) * F(q, n))
    third_factor = 1 / sqrt_term * (1 - (-1)**(n+m)) / 2
    return first_factor * second_factor * third_factor

def Omega(omega_H, q, n):
    return (omega_H + (gamma * DD) / L**2 * (q**2 + n**2 * np.pi**2)) / omega_M

def H(omega_H, q, phi_k, n, m):
    kronecker = 1.0 if n == m else 0.0
    first_term = Omega(omega_H, q, n) * np.array([[1, 0], [0, 1]]) * kronecker
    second_term = 0.5 * np.array([[1, 1], [1, 1]]) * kronecker
    third_term = 0.5 * np.array([[1 - np.sin(phi_k)**2, 1 + np.sin(phi_k)**2], 
                                [1 + np.sin(phi_k)**2, 1 - np.sin(phi_k)**2]]) * P(q, n, m)
    fourth_term = 0.5 * np.array([[0, -4], [4, 0]]) * np.sin(phi_k) * Q(q, n, m)
    return first_term + second_term - third_term - fourth_term

def generate_HBdG(omega_H, q, phi_k, Nmax):
    HBdG = np.zeros((2*Nmax, 2*Nmax), dtype=complex)
    for m in range(Nmax):
        for n in range(Nmax):
            result = H(omega_H, q, phi_k, n, m)
            HBdG[n, m] = result[0, 0]
            HBdG[n, m+Nmax] = result[0, 1]
            HBdG[n+Nmax, m] = result[1, 0]
            HBdG[n+Nmax, m+Nmax] = result[1, 1]
    return HBdG

def paraunitary_diag(H_matrix):
    """Perform paraunitary diagonalization"""
    # Similar to CholeskyDecomposition in Mathematica
    K = linalg.cholesky(H_matrix)
    
    # Create sigma3 matrix
    dim = H_matrix.shape[0]
    sigma3 = np.diag([(-1)**np.floor((2*n-1)/dim) for n in range(1, dim+1)])
    
    # Calculate W matrix
    W = K @ sigma3 @ np.conjugate(K).T
    
    # Find eigenvalues and eigenvectors
    eval, evec = linalg.eigh(W)
    evec = np.array([v/np.linalg.norm(v) for v in evec.T]).T
    
    # Create permutation
    preperm = np.concatenate([np.arange(dim//2, dim), np.arange(dim//2, 0, -1) - 1])
    ordering = np.argsort(eval)
    permutation = np.argsort(preperm)[ordering]
    
    # Apply permutation
    eval = eval[permutation]
    evec = evec[:, permutation]
    U = evec.T
    
    # Create diagonal matrix
    Hdiag = sigma3 @ np.diag(eval)
    
    # Calculate T matrix
    T = np.linalg.inv(K) @ U @ np.sqrt(np.abs(Hdiag))
    
    # Extract blocks
    Tpp = T[:dim//2, :dim//2]
    Tnp = T[dim//2:, :dim//2]
    Tpn = T[:dim//2, dim//2:]
    Tnn = T[dim//2:, dim//2:]
    
    # Phase correction
    phase_array_p = np.exp(1j * np.angle(np.diag(Tpp)))
    phase_array_n = np.exp(1j * np.angle(np.diag(Tnn)))
    V = np.diag(np.concatenate([np.conjugate(phase_array_p), np.conjugate(phase_array_n)]))
    T = T @ V
    
    # Recompute blocks after phase correction
    Tpp = T[:dim//2, :dim//2]
    Tnp = T[dim//2:, :dim//2]
    Tpn = T[:dim//2, dim//2:]
    Tnn = T[dim//2:, dim//2:]
    
    return eval[:dim//2], Tpp, Tnp, Tpn, Tnn, T

def f(q, n, hover_L):
    return (-1)**n / np.sqrt(2 * (1 + (1.0 if n == 0 else 0))) * q**2 / (q**2 + n**2 * np.pi**2) * np.exp(-q*hover_L) * (1 - (-1)**n * np.exp(-q))

def multi_para_diag(hNV_array, omega_H, qtable, phi_table, Nmax):
    """Execute the paraunitary diagonalization and return interpolation functions"""
    print("Starting multi_para_diag...")
    start_time = time.time()
    
    Numq = len(qtable)
    Num_phi = len(phi_table)
    NumhNV = len(hNV_array)
    
    omega_BdG_table = np.zeros((Numq, 2*Num_phi, Nmax))
    coupling_plus_table = np.zeros((Numq, 2*Num_phi, NumhNV, Nmax), dtype=complex)
    coupling_minus_table = np.zeros((Numq, 2*Num_phi, NumhNV, Nmax), dtype=complex)
    coupling_z_table = np.zeros((Numq, 2*Num_phi, NumhNV, Nmax), dtype=complex)
    
    for countq in range(Numq):
        if countq % 10 == 0:
            print(f"Progress: {countq/Numq*100:.1f}% (q={countq+1}/{Numq})")
            print(f"Time elapsed: {time.time() - start_time:.1f} seconds")
        
        for count_phi in range(Num_phi):
            q = qtable[countq]
            phi_k = phi_table[count_phi]
            
            HBdG = generate_HBdG(omega_H, q, phi_k, Nmax)
            result = paraunitary_diag((HBdG + np.conjugate(HBdG.T))/2)
            
            omega_BdG_table[countq, count_phi] = result[0]
            omega_BdG_table[countq, count_phi + Num_phi] = result[0]
            
            Tpp = result[1]
            Tnp = result[2]
            Tpn = result[3]
            Tnn = result[4]
            
            # Calculate gamma values
            gamma_plus = (1 + np.sin(phi_k))/2 * (Tpp + Tnp + np.sin(phi_k)*(Tpp - Tnp))
            gamma_plus_mirror = (1 + np.sin(phi_k + np.pi))/2 * np.conjugate((Tnn + Tpn + np.sin(phi_k + np.pi)*(Tnn - Tpn)))
            
            gamma_minus = (1 - np.sin(phi_k))/2 * (Tpp + Tnp + np.sin(phi_k)*(Tpp - Tnp))
            gamma_minus_mirror = (1 - np.sin(phi_k + np.pi))/2 * np.conjugate((Tnn + Tpn + np.sin(phi_k + np.pi)*(Tnn - Tpn)))
            
            gamma_z = -1j * np.cos(phi_k)/2 * (Tpp + Tnp + np.sin(phi_k)*(Tpp - Tnp))
            gamma_z_mirror = -1j * np.cos(phi_k + np.pi)/2 * np.conjugate((Tnn + Tpn + np.sin(phi_k + np.pi)*(Tnn - Tpn)))
            
            # Calculate fbar values
            for tt in range(NumhNV):
                fbar_array = np.array([f(q, nn, hNV_array[tt]/L) for nn in range(Nmax)])
                
                vplus = fbar_array @ gamma_plus
                vplus_mirror = fbar_array @ gamma_plus_mirror
                vminus = fbar_array @ gamma_minus
                vminus_mirror = fbar_array @ gamma_minus_mirror
                vz = fbar_array @ gamma_z
                vz_mirror = fbar_array @ gamma_z_mirror
                
                coupling_plus_table[countq, count_phi, tt] = vplus
                coupling_plus_table[countq, count_phi + Num_phi, tt] = vplus_mirror
                coupling_minus_table[countq, count_phi, tt] = vminus
                coupling_minus_table[countq, count_phi + Num_phi, tt] = vminus_mirror
                coupling_z_table[countq, count_phi, tt] = vz
                coupling_z_table[countq, count_phi + Num_phi, tt] = vz_mirror
    
    # Create interpolation functions
    print("Creating interpolation functions...")
    phi_table_mod = np.array([phi_table[0] + i*np.pi for i in range(2)])
    
    int_omega_BdG_table = []
    int_coupling_plus_table = []
    int_coupling_minus_table = []
    int_coupling_z_table = []
    
    # For simplicity, we'll use a simple approach for interpolation
    # This is a simplified version of the original 2D interpolation
    
    return {
        'qtable': qtable,
        'omega_min': np.min(omega_BdG_table),
        'phi_table_mod': phi_table_mod,
        'omega_BdG_table': omega_BdG_table,
        'coupling_plus_table': coupling_plus_table,
        'coupling_minus_table': coupling_minus_table,
        'coupling_z_table': coupling_z_table
    }

def NBose(omega):
    """Calculate Bose-Einstein distribution"""
    h_plank = constants.h  # J*s
    mu0 = constants.mu_0  # H/m
    kB = constants.k  # J/K
    temperature = 300  # K
    return (1e-9 * kB * temperature / h_plank) / omega

def gamma_values_UL(eta_small, NQest, N_capital_phi, H0, omega_BdG_table, coupling_plus_table, coupling_minus_table, num_hNV):
    """Calculate Gamma values"""
    DNV = 2.87
    omega_target_L = DNV - H0 * gamma
    omega_target_U = DNV + H0 * gamma
    
    # Simplified implementation
    print("Calculating Gamma values...")
    DOS_L = 0
    DOS_U = 0
    gamma_L_Hz_array = np.zeros(num_hNV)
    gamma_U_Hz_array = np.zeros(num_hNV)
    
    # Use available data without interpolation in this simplified version
    numq, num_phi, _ = omega_BdG_table.shape
    
    for s in range(Nmax):  # Limit to 5 for simplicity
        print(f"Processing mode s={s+1}")
        
        for q_idx in range(numq):
            for phi_idx in range(num_phi):
                q = qtable[q_idx]
                if q_idx > 0:
                    delta_q = qtable[q_idx] - qtable[q_idx-1]
                else:
                    delta_q = qtable[1] - qtable[0]
                
                phi = phi_table[phi_idx % len(phi_table)]
                delta_phi = np.pi / len(phi_table)
                
                omega = omega_BdG_table[q_idx, phi_idx, s]
                
                # Calculate DOS contributions
                lorentzian_L = (eta_small/np.pi)/(eta_small**2 + (omega_M*omega - omega_target_L)**2)
                lorentzian_U = (eta_small/np.pi)/(eta_small**2 + (omega_M*omega - omega_target_U)**2)
                
                DOS_L += (1/L**2) * (delta_phi/(2*np.pi)**2) * delta_q * q * lorentzian_L
                DOS_U += (1/L**2) * (delta_phi/(2*np.pi)**2) * delta_q * q * lorentzian_U
                
                # Calculate Gamma contributions
                for tt in range(num_hNV):
                    coupling_plus = coupling_plus_table[q_idx, phi_idx, tt, s]
                    coupling_minus = coupling_minus_table[q_idx, phi_idx, tt, s]
                    
                    bose_factor = 2*NBose(omega_M*omega) + 1
                    
                    gamma_L_Hz_array[tt] += (2*np.pi)**2 * omega_M * omega_d * (delta_phi/(2*np.pi)**2) * \
                                           bose_factor * delta_q * q * (abs(coupling_plus)**2) * lorentzian_L
                    
                    gamma_U_Hz_array[tt] += (2*np.pi)**2 * omega_M * omega_d * (delta_phi/(2*np.pi)**2) * \
                                           bose_factor * delta_q * q * (abs(coupling_minus)**2) * lorentzian_U
    
    return DOS_L, DOS_U, gamma_L_Hz_array, gamma_U_Hz_array

def gamma_1_from_H0(H0):
    """Calculate Gamma1 values from H0"""
    print(f"Calculating for H0 = {H0} Oe")
    omega_H = gamma * H0  # GHz
    
    # Run multi_para_diag
    result = multi_para_diag(hNV_array, omega_H, qtable, phi_table, Nmax)
    
    # Extract results
    omega_BdG_table = result['omega_BdG_table']
    coupling_plus_table = result['coupling_plus_table']
    coupling_minus_table = result['coupling_minus_table']
    
    # Calculate Gamma values
    eta_small = 0.003
    NQ = 2*200  # Reduced for performance
    N_capital_phi = 2*360  # Reduced for performance
    
    DOS_L, DOS_U, gamma_L_Hz_array, gamma_U_Hz_array = gamma_values_UL(
        eta_small, NQ, N_capital_phi, H0, omega_BdG_table, 
        coupling_plus_table, coupling_minus_table, len(hNV_array)
    )
    
    print(f"Gamma(omega=omega_L) = {gamma_L_Hz_array} Hz")
    print(f"Gamma(omega=omega_U) = {gamma_U_Hz_array} Hz")
    
    return DOS_L, DOS_U, gamma_L_Hz_array, gamma_U_Hz_array

# Setup constants for calculations
h_plank = constants.h  # J*s
mu0 = constants.mu_0  # H/m
omega_d = (h_plank * mu0 * (gamma*1e9)**2 / (L*1e-6)**3) * 1e8  # Hz

# Parameters for calculation
Num_phi = 2*90  # Reduced from 2*90 for performance
NumQ = 2*100   # Reduced from 2*100 for performance
Del_phi = np.pi / Num_phi
phi_table = np.linspace(0, np.pi - Del_phi, Num_phi)

Qmax = 50 * L  # till 50 rad/μm
# Create log-spaced qtable
qtable = np.concatenate([[1e-6], np.logspace(np.log10(L/1000), np.log10(Qmax), NumQ-1)])

Fmax = 5
Nmax = int(np.ceil(L/np.pi * np.sqrt(Fmax/(gamma * DD))))  # till 5GHz
print(f"Nmax is {Nmax}")

# NV array positions in μm
hNV_array = np.array([0.4])
print(f"NV array positions: {hNV_array} μm")

# Field array
H0_array = np.array([80, 81.5, 82, 82.5, 83, 85])

# Run calculation for each field value
DOS_L_array = np.zeros_like(H0_array, dtype=float)
DOS_U_array = np.zeros_like(H0_array, dtype=float)
gamma_L_Hz_array = np.zeros((len(hNV_array), len(H0_array)), dtype=float)
gamma_U_Hz_array = np.zeros((len(hNV_array), len(H0_array)), dtype=float)

for i in range(len(H0_array)):
   result_gamma_UL = gamma_1_from_H0(H0_array[i])
   DOS_L_array[i] = result_gamma_UL[0]
   DOS_U_array[i] = result_gamma_UL[1]
   gamma_L_Hz_array[:, i] = result_gamma_UL[2]
   gamma_U_Hz_array[:, i] = result_gamma_UL[3]

# Plot results
plt.figure(figsize=(10, 5))
plt.plot(H0_array, gamma_L_Hz_array[0]*1000, 'o-', label='γL x 1000')
plt.plot(H0_array, gamma_U_Hz_array[0]*1000, 'o-', label='γU x 1000')
plt.xlabel('Field (G)')
plt.ylabel('γL, γU x 1,000 (Hz)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Save results
#saving_matrix = np.column_stack([np.concatenate([[0], H0_array]), 
                               #np.row_stack([hNV_array, gamma_L_Hz_array])])
print("Calculation complete")