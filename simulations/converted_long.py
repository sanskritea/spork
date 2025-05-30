import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from scipy.linalg import cholesky
import time

# -------------------------- Parameters ---------------------------------
# Calculate Ms solution value (equivalent to MsSol in Mathematica)
delta_H = 1
Ms_sol = ((2870/2.8) - 82 - delta_H) * 2.8 - 4 + 2
print(f"Ms_sol = {Ms_sol}")

# Material & condition parameters
M0 = Ms_sol  # Oe
L = 3  # μm
DD = 5.4e-9 * 1e8  # Oe μm^-2
gamma = 2.8e-3  # GHz/G
omega_M = gamma * M0  # GHz

# Physical constants
h_plank = 6.626e-34  # J*s
mu0 = 4 * np.pi * 1e-7  # H/m
omega_d = (h_plank * mu0 * (gamma * 1e9)**2 / (L * 1e-6)**3) * 1e8  # Hz
kB = 1.381e-23  # J/K
temperature = 300  # K

def N_bose(omega):
    """Bose-Einstein distribution function approximation."""
    return (1e-9 * kB * temperature / h_plank) / omega

# -------------------------- Utility Functions --------------------------
def F(q, n):
    return 2 * (1 - (-1)**n * np.exp(-q)) / q

def P(q, n, m):
    term1 = (q**2 / (q**2 + n**2 * np.pi**2)) * (n == m)
    
    normalization = 1 / np.sqrt((1 + (n == 0)) * (1 + (m == 0)))
    term2 = normalization * (q**4 / ((q**2 + n**2 * np.pi**2) * (q**2 + m**2 * np.pi**2))) * F(q, n) * (1 + (-1)**(n+m)) / 2
    
    return term1 - term2

def Q(q, n, m):
    if n == m:
        return 0
        
    normalization = 1 / np.sqrt((1 + (n == 0)) * (1 + (m == 0)))
    
    term = q**2 / (q**2 + m**2 * np.pi**2) * (
        m**2 / (m**2 - n**2 + (1 + (-1)**(n+m)) / 2) * 2 / q - 
        q**2 / (2 * (q**2 + n**2 * np.pi**2)) * F(q, n)
    )
    
    return term * normalization * (1 - (-1)**(n+m)) / 2

def Omega(omega_H, q, n):
    return (omega_H + (gamma * DD) / L**2 * (q**2 + n**2 * np.pi**2)) / omega_M

def H(omega_H, q, phi_k, n, m):
    """Generate Hamiltonian element."""
    if n == m:
        # Diagonal term
        term1 = np.array([[Omega(omega_H, q, n), 0], [0, Omega(omega_H, q, n)]])
        term2 = 0.5 * np.array([[1, 1], [1, 1]])
        diag_part = term1 + term2
        return diag_part
    else:
        # Off-diagonal term
        sin_phi_k = np.sin(phi_k)
        sin_sq = sin_phi_k**2
        
        term1 = -0.5 * np.array([[1-sin_sq, 1+sin_sq], [1+sin_sq, 1-sin_sq]]) * P(q, n, m)
        term2 = -0.5 * np.array([[0, -4], [4, 0]]) * sin_phi_k * Q(q, n, m)
        
        return term1 + term2

def f(q, n, hover_L):
    """Coupling function."""
    return (-1)**n / np.sqrt(2 * (1 + (n == 0))) * q**2 / (q**2 + n**2 * np.pi**2) * np.exp(-q * hover_L) * (1 - (-1)**n * np.exp(-q))

def generate_H_BdG(omega_H, q, phi_k, N_max):
    """Generate Bogoliubov-de Gennes Hamiltonian."""
    H_BdG = np.zeros((2*N_max, 2*N_max), dtype=complex)
    
    for m in range(N_max):
        for n in range(N_max):
            result = H(omega_H, q, phi_k, n, m)
            
            H_BdG[n, m] = result[0, 0]
            H_BdG[n, m+N_max] = result[0, 1]
            H_BdG[n+N_max, m] = result[1, 0]
            H_BdG[n+N_max, m+N_max] = result[1, 1]
    
    return H_BdG

def para_unitary_diag(H):
    """Perform paraunitary diagonalization."""
    # Ensure H is Hermitian for numerical stability
    H = (H + H.conj().T) / 2
    
    # Cholesky decomposition
    try:
        K = cholesky(H)
    except np.linalg.LinAlgError:
        # Add small diagonal term if not positive definite
        H = H + 1e-10 * np.eye(H.shape[0])
        K = cholesky(H)
    
    Dim = H.shape[0]
    
    # Create sigma_3 matrix
    sigma_3 = np.diag([(-1)**int(np.floor((2*n-1)/Dim)) for n in range(1, Dim+1)])
    
    # Calculate W matrix
    W = K @ sigma_3 @ K.conj().T
    
    # Eigendecomposition
    eval, evec = np.linalg.eigh(W)
    
    # Normalize eigenvectors
    for i in range(len(evec)):
        evec[:, i] = evec[:, i] / np.sqrt(np.sum(np.abs(evec[:, i])**2))
    
    # Create permutation that gives (+small, -small, +mid, -mid, +large, -large)
    preperm = list(range(Dim//2, Dim)) + list(range(Dim//2, 0, -1))
    ordering = np.argsort(eval)
    
    # Apply permutation
    permutation = [preperm[i] for i in ordering]
    eval = eval[permutation]
    evec = evec[:, permutation]
    
    # Construct U matrix
    U = evec.T
    
    # Compute H_diag
    H_diag = sigma_3 @ np.diag(eval)
    
    # Compute T matrix
    T = np.linalg.inv(K) @ U @ np.sqrt(np.abs(H_diag))
    
    # Extract blocks of T
    Tpp = T[:Dim//2, :Dim//2]
    Tnn = T[Dim//2:, Dim//2:]
    
    # Phase correction
    phase_array_p = np.exp(1j * np.angle(np.diag(Tpp)))
    phase_array_n = np.exp(1j * np.angle(np.diag(Tnn)))
    
    V = np.diag(np.concatenate([np.conjugate(phase_array_p), np.conjugate(phase_array_n)]))
    T = T @ V
    
    # Extract blocks after phase correction
    Tpp = T[:Dim//2, :Dim//2]
    Tnp = T[Dim//2:, :Dim//2]
    Tpn = T[:Dim//2, Dim//2:]
    Tnn = T[Dim//2:, Dim//2:]
    
    return [eval[:Dim//2], Tpp, Tnp, Tpn, Tnn, T]

def multi_para_diag(hNV_array, omega_H, q_table, phi_table, N_max):
    """Perform paraunitary diagonalization for multiple parameters."""
    Num_q = len(q_table)
    Num_phi = len(phi_table)
    Num_hNV = len(hNV_array)
    
    # Calculate measure for integration
    measure = np.mean(np.diff(q_table)) * np.mean(np.diff(phi_table)) / (2 * np.pi)**2
    
    # Initialize arrays for results
    Q_table = np.tile(np.reshape(q_table, (Num_q, 1, 1)), (1, 2*Num_phi, N_max))
    
    # Create phi table with extended range (phi and phi+pi)
    phi_extended = []
    for n in range(Num_q):
        for m in range(2*Num_phi):
            phi_val = phi_table[m % Num_phi] + np.pi * (m // Num_phi)
            phi_extended.append([n, m, phi_val])
    phi_extended = np.array(phi_extended)
    
    # Initialize result arrays
    omega_BdG_table = np.zeros((Num_q, 2*Num_phi, N_max))
    coupling_plus_table = np.zeros((Num_q, 2*Num_phi, Num_hNV, N_max), dtype=complex)
    coupling_minus_table = np.zeros((Num_q, 2*Num_phi, Num_hNV, N_max), dtype=complex)
    coupling_z_table = np.zeros((Num_q, 2*Num_phi, Num_hNV, N_max), dtype=complex)
    
    print("Starting paraunitary diagonalization...")
    start_time = time.time()
    
    # Loop over q and phi values
    for count_q in range(Num_q):
        for count_phi in range(Num_phi):
            if count_q % 10 == 0 and count_phi == 0:
                elapsed = time.time() - start_time
                progress = (count_q * Num_phi + count_phi) / (Num_q * Num_phi)
                est_remaining = elapsed / progress - elapsed if progress > 0 else 0
                print(f"Progress: {100*progress:.1f}% - Est. remaining: {est_remaining/60:.1f} min")
            
            q = q_table[count_q]
            phi_k = phi_table[count_phi]
            
            # Generate and diagonalize Hamiltonian
            H_BdG = generate_H_BdG(omega_H, q, phi_k, N_max)
            result = para_unitary_diag((H_BdG + H_BdG.conj().T) / 2)
            
            # Store eigenvalues
            omega_BdG_table[count_q, count_phi, :] = result[0]
            omega_BdG_table[count_q, count_phi + Num_phi, :] = result[0]
            
            # Extract transformation matrices
            Tpp, Tnp, Tpn, Tnn = result[1], result[2], result[3], result[4]
            
            # Calculate coupling coefficients
            gamma_plus = (1 + np.sin(phi_k))/2 * (Tpp + Tnp + np.sin(phi_k) * (Tpp - Tnp))
            gamma_plus_mirror = (1 + np.sin(phi_k + np.pi))/2 * np.conjugate(Tnn + Tpn + np.sin(phi_k + np.pi) * (Tnn - Tpn))
            
            gamma_minus = (1 - np.sin(phi_k))/2 * (Tpp + Tnp + np.sin(phi_k) * (Tpp - Tnp))
            gamma_minus_mirror = (1 - np.sin(phi_k + np.pi))/2 * np.conjugate(Tnn + Tpn + np.sin(phi_k + np.pi) * (Tnn - Tpn))
            
            gamma_z = -1j * np.cos(phi_k)/2 * (Tpp + Tnp + np.sin(phi_k) * (Tpp - Tnp))
            gamma_z_mirror = -1j * np.cos(phi_k + np.pi)/2 * np.conjugate(Tnn + Tpn + np.sin(phi_k + np.pi) * (Tnn - Tpn))
            
            # Calculate coupling for each NV height
            for tt, h_nv in enumerate(hNV_array):
                fbar_array = np.array([f(q, nn, h_nv/L) for nn in range(N_max)])
                
                v_plus = fbar_array @ gamma_plus
                v_plus_mirror = fbar_array @ gamma_plus_mirror
                v_minus = fbar_array @ gamma_minus
                v_minus_mirror = fbar_array @ gamma_minus_mirror
                v_z = fbar_array @ gamma_z
                v_z_mirror = fbar_array @ gamma_z_mirror
                
                coupling_plus_table[count_q, count_phi, tt, :] = v_plus
                coupling_plus_table[count_q, count_phi + Num_phi, tt, :] = v_plus_mirror
                coupling_minus_table[count_q, count_phi, tt, :] = v_minus
                coupling_minus_table[count_q, count_phi + Num_phi, tt, :] = v_minus_mirror
                coupling_z_table[count_q, count_phi, tt, :] = v_z
                coupling_z_table[count_q, count_phi + Num_phi, tt, :] = v_z_mirror
    
    print(f"Paraunitary diagonalization completed in {time.time() - start_time:.1f} seconds")
    
    # Create interpolation functions
    print("Creating interpolation functions...")
    
    # Prepare data for interpolation - extending phi range for better interpolation
    Q_table_extended = np.zeros((Num_q, 2*Num_phi+1, N_max))
    phi_table_extended = np.zeros((Num_q, 2*Num_phi+1, N_max))
    omega_BdG_table_extended = np.zeros((Num_q, 2*Num_phi+1, N_max))
    coupling_plus_table_extended = np.zeros((Num_q, 2*Num_phi+1, Num_hNV, N_max), dtype=complex)
    coupling_minus_table_extended = np.zeros((Num_q, 2*Num_phi+1, Num_hNV, N_max), dtype=complex)
    coupling_z_table_extended = np.zeros((Num_q, 2*Num_phi+1, Num_hNV, N_max), dtype=complex)
    
    # Fill extended tables
    for count_q in range(Num_q):
        for count_phi in range(2*Num_phi):
            omega_BdG_table_extended[count_q, count_phi, :] = omega_BdG_table[count_q, count_phi, :]
            coupling_plus_table_extended[count_q, count_phi, :, :] = coupling_plus_table[count_q, count_phi, :, :]
            coupling_minus_table_extended[count_q, count_phi, :, :] = coupling_minus_table[count_q, count_phi, :, :]
            coupling_z_table_extended[count_q, count_phi, :, :] = coupling_z_table[count_q, count_phi, :, :]
        
        # Add wrap-around point
        omega_BdG_table_extended[count_q, 2*Num_phi, :] = omega_BdG_table[count_q, 0, :]
        coupling_plus_table_extended[count_q, 2*Num_phi, :, :] = coupling_plus_table[count_q, 0, :, :]
        coupling_minus_table_extended[count_q, 2*Num_phi, :, :] = coupling_minus_table[count_q, 0, :, :]
        coupling_z_table_extended[count_q, 2*Num_phi, :, :] = coupling_z_table[count_q, 0, :, :]
    
    # Create interpolation grid
    q_grid = q_table
    phi_grid = np.linspace(0, 2*np.pi, 2*Num_phi+1)
    
    # Create interpolation functions for each mode
    omega_BdG_interp = []
    coupling_plus_interp = []
    coupling_minus_interp = []
    coupling_z_interp = []
    
    for s in range(N_max):
        # Create interpolator for omega
        points = (q_grid, phi_grid)
        values = omega_BdG_table_extended[:, :, s].real  # Take real part for safety
        omega_BdG_interp.append(RegularGridInterpolator(points, values, method='linear'))
        
        # Create interpolators for couplings
        coupling_plus_mode = []
        coupling_minus_mode = []
        coupling_z_mode = []
        
        for tt in range(Num_hNV):
            # Real and imaginary parts need separate interpolation
            cp_real = coupling_plus_table_extended[:, :, tt, s].real
            cp_imag = coupling_plus_table_extended[:, :, tt, s].imag
            cm_real = coupling_minus_table_extended[:, :, tt, s].real
            cm_imag = coupling_minus_table_extended[:, :, tt, s].imag
            cz_real = coupling_z_table_extended[:, :, tt, s].real
            cz_imag = coupling_z_table_extended[:, :, tt, s].imag
            
            coupling_plus_mode.append((
                RegularGridInterpolator(points, cp_real, method='linear'),
                RegularGridInterpolator(points, cp_imag, method='linear')
            ))
            coupling_minus_mode.append((
                RegularGridInterpolator(points, cm_real, method='linear'),
                RegularGridInterpolator(points, cm_imag, method='linear')
            ))
            coupling_z_mode.append((
                RegularGridInterpolator(points, cz_real, method='linear'),
                RegularGridInterpolator(points, cz_imag, method='linear')
            ))
        
        coupling_plus_interp.append(coupling_plus_mode)
        coupling_minus_interp.append(coupling_minus_mode)
        coupling_z_interp.append(coupling_z_mode)
    
    # Return results
    phi_table_mod = [phi_extended[n * 2*Num_phi, 2] for n in range(Num_phi * 2)]
    
    return [q_table, np.min(omega_BdG_table), phi_table_mod, 
            omega_BdG_interp, coupling_plus_interp, coupling_minus_interp, coupling_z_interp]

def gamma_values_UL(eta_small, NQ_est, N_Capital_Phi, H0, omega_BdG_interp, coupling_plus_interp, coupling_minus_interp, Num_hNV):
    """Calculate Gamma values for upper and lower transitions."""
    DNV = 2.87
    omega_target_L = DNV - H0 * gamma
    omega_target_U = DNV + H0 * gamma
    
    # Define q grid with higher density near middle
    Q_middle = 5 * L
    NQ_half = round(NQ_est / 2)
    
    # Logarithmic spacing for first half
    q_first_half = np.logspace(np.log10(L/1000), np.log10(Q_middle), NQ_half)
    q_table = np.concatenate([[1e-6], q_first_half])
    
    # Linear spacing for second half
    last_delta_q = q_table[-1] - q_table[-2]
    q_second_half = np.arange(q_table[-1] + last_delta_q, 50*L, last_delta_q)
    q_table = np.concatenate([q_table, q_second_half])
    
    delta_q = np.diff(q_table)
    NQ = len(delta_q)
    delta_Capital_Phi = 2 * np.pi / N_Capital_Phi
    
    # Initialize result arrays
    Gamma_L_Hz_array = np.zeros(Num_hNV)
    Gamma_U_Hz_array = np.zeros(Num_hNV)
    DOS_L = 0
    DOS_U = 0
    
    print("Starting Gamma calculations...")
    start_time = time.time()
    
    # Loop over modes
    N_max = len(omega_BdG_interp)
    for s in range(N_max):
        if s % 2 == 0:
            elapsed = time.time() - start_time
            progress = s / N_max
            est_remaining = elapsed / progress - elapsed if progress > 0 else 0
            print(f"Mode {s}/{N_max} - {100*progress:.1f}% complete - Est. remaining: {est_remaining/60:.1f} min")
        
        # Create flattened tables for calculations
        flat_q_table = np.array([q_table[i1] for i1 in range(NQ) for i2 in range(N_Capital_Phi)])
        flat_delta_q_table = np.array([delta_q[i1] for i1 in range(NQ) for i2 in range(N_Capital_Phi)])
        flat_phi_table = np.array([2*np.pi*(i2)/N_Capital_Phi for i1 in range(NQ) for i2 in range(N_Capital_Phi)])
        
        # Interpolate omega and coupling values
        flat_omega_BdG_values = np.zeros(len(flat_q_table))
        flat_coupling_plus_values = np.zeros((Num_hNV, len(flat_q_table)), dtype=complex)
        flat_coupling_minus_values = np.zeros((Num_hNV, len(flat_q_table)), dtype=complex)
        
        for ii, (q, phi) in enumerate(zip(flat_q_table, flat_phi_table)):
            points = np.array([q, phi])
            flat_omega_BdG_values[ii] = omega_BdG_interp[s](points)
            
            for tt in range(Num_hNV):
                real = coupling_plus_interp[s][tt][0](points)
                imag = coupling_plus_interp[s][tt][1](points)
                flat_coupling_plus_values[tt, ii] = real + 1j * imag
                
                real = coupling_minus_interp[s][tt][0](points)
                imag = coupling_minus_interp[s][tt][1](points)
                flat_coupling_minus_values[tt, ii] = real + 1j * imag
        
        # Calculate DOS contributions
        lorentzian_L = (eta_small/np.pi) / (eta_small**2 + (omega_M * flat_omega_BdG_values - omega_target_L)**2)
        lorentzian_U = (eta_small/np.pi) / (eta_small**2 + (omega_M * flat_omega_BdG_values - omega_target_U)**2)
        
        DOS_L += (1/L**2) * (delta_Capital_Phi/(2*np.pi)**2) * np.sum(flat_delta_q_table * flat_q_table * lorentzian_L)
        DOS_U += (1/L**2) * (delta_Capital_Phi/(2*np.pi)**2) * np.sum(flat_delta_q_table * flat_q_table * lorentzian_U)
        
        # Calculate Gamma contributions
        for tt in range(Num_hNV):
            # Lower transition
            coupling_term_L = np.abs(flat_coupling_plus_values[tt, :])**2
            bose_term = 2 * N_bose(omega_M * flat_omega_BdG_values) + 1
            Gamma_L_Hz_array[tt] += (2*np.pi)**2 * omega_M * omega_d * (delta_Capital_Phi/(2*np.pi)**2) * np.sum(
                bose_term * flat_delta_q_table * flat_q_table * coupling_term_L * lorentzian_L
            )
            
            # Upper transition
            coupling_term_U = np.abs(flat_coupling_minus_values[tt, :])**2
            Gamma_U_Hz_array[tt] += (2*np.pi)**2 * omega_M * omega_d * (delta_Capital_Phi/(2*np.pi)**2) * np.sum(
                bose_term * flat_delta_q_table * flat_q_table * coupling_term_U * lorentzian_U
            )
    
    print(f"Gamma calculations completed in {time.time() - start_time:.1f} seconds")
    return [DOS_L, DOS_U, Gamma_L_Hz_array, Gamma_U_Hz_array]

def gamma_1_from_H0(H0, hNV_array, q_table, phi_table, N_max):
    """Calculate Gamma1 for a given magnetic field H0."""
    omega_H = gamma * H0  # GHz
    
    print(f"Starting calculations for H0 = {H0} Oe")
    result = multi_para_diag(hNV_array, omega_H, q_table, phi_table, N_max)
    
    [q_table_temp, omega_min, phi_table_mod, 
     omega_BdG_interp, coupling_plus_interp, coupling_minus_interp, coupling_z_interp] = result
    
    eta_small = 0.003
    NQ = 2 * 2
    N_Capital_Phi = 2 * 360
    
    [DOS_L, DOS_U, Gamma_L_Hz_array, Gamma_U_Hz_array] = gamma_values_UL(
        eta_small, NQ, N_Capital_Phi, H0, 
        omega_BdG_interp, coupling_plus_interp, coupling_minus_interp, len(hNV_array)
    )
    
    print(f"DOS(ω=ωL) = {DOS_L:.6e} 1/GHz μm²")
    print(f"DOS(ω=ωU) = {DOS_U:.6e} 1/GHz μm²")
    print(f"Γ(ω=ωL) = {Gamma_L_Hz_array} Hz")
    print(f"Γ(ω=ωU) = {Gamma_U_Hz_array} Hz")
    
    return [DOS_L, DOS_U, Gamma_L_Hz_array, Gamma_U_Hz_array]

# -------------------------- Main Program --------------------------
def main():
    # Setup parameters similar to the Mathematica code
    Num_phi = 2 * 90  # Reduced from original for speed
    Num_Q = 2 * 1    # Reduced from original for speed
    
    Del_phi = np.pi / Num_phi
    phi_table = np.linspace(0, np.pi - Del_phi, Num_phi)
    
    Q_max = 50 * L  # Up to 50 rad/μm
    
    # Create logarithmically spaced q values
    def f_space(min_val, max_val, steps, f=np.log):
        return np.exp(np.linspace(np.log(min_val), np.log(max_val), steps))
    
    q_table = np.concatenate([[1e-6], f_space(L/1000, Q_max, Num_Q)])
    
    F_max = 5
    N_max = int(np.ceil(L/np.pi * np.sqrt(F_max/(gamma * DD))))
    print(f"N_max = {N_max}")
    
    # Define NV heights 
    hNV_array = np.array([0.4])#, 0.5, 0.6, 0.7])  # positions of NV in μm
    
    # # Test for a single magnetic field value
    # H0 = 82  # Oe
    # result = gamma_1_from_H0(H0, hNV_array, q_table, phi_table, N_max)
    
    # For full field sweep, uncomment and run:
    
    H0_array = np.array([80, 81.5, 82, 82.5, 83, 85])
    DOS_L_array = np.zeros_like(H0_array)
    DOS_U_array = np.zeros_like(H0_array)
    Gamma_L_Hz_array = np.zeros((len(hNV_array), len(H0_array)))
    Gamma_U_Hz_array = np.zeros((len(hNV_array), len(H0_array)))
    
    for ii, H0 in enumerate(H0_array):
        print(f"Processing H0 = {H0} Oe ({ii+1}/{len(H0_array)})")
        result_gamma_UL = gamma_1_from_H0(H0, hNV_array, q_table, phi_table, N_max)
        
        DOS_L_array[ii] = result_gamma_UL[0]
        DOS_U_array[ii] = result_gamma_UL[1]
        Gamma_L_Hz_array[:, ii] = result_gamma_UL[2]
        Gamma_U_Hz_array[:, ii] = result_gamma_UL[3]

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

if __name__ == "__main__":
    main()