import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from scipy.interpolate import griddata, RegularGridInterpolator
from tqdm import tqdm
import multiprocessing as mp
from functools import partial
import numba
from numba import jit, prange

# Material and condition parameters
M0 = 1716  # Oe (MsSol in original)
thickness = 3  # μm
DD = 5.4e-9 * 1e8  # Oe μm^-2
gamma = 2.8e-3  # GHz/G


# Constants
h_planck = 6.626e-34  # J*s
mu_0 = 4 * np.pi * 1e-7  # H/m
k_B = 1.381e-23  # J/K
temperature = 300  # K
omega_M = gamma * M0  # GHz

# Precomputed constants
CONST_FACTOR = 2 * np.pi**2 * omega_M * h_planck * mu_0 * (gamma * 1e9)**2   # Factor for gamma calculation
# CONST_FACTOR = mu_0 * (gamma * 1e9)**2 # Factor for gamma calculation


@jit(nopython=True)
def F(q, n):
    """F function for dipolar exchange spin waves - Optimized with numba"""
    return 2 * (1 - (-1)**n * np.exp(-q)) / q


@jit(nopython=True)
def P(q, n, m):
    """P function for dipolar exchange spin waves - Optimized with numba"""
    # Avoid division by zero and precompute common expressions
    q_squared = q**2
    n_pi_squared = n**2 * np.pi**2
    m_pi_squared = m**2 * np.pi**2
    
    term1 = q_squared / (q_squared + n_pi_squared) * (1 if n == m else 0)
    
    factor = 1 / np.sqrt((1 + (1 if n == 0 else 0)) * (1 + (1 if m == 0 else 0)))
    term2 = q_squared**2 / ((q_squared + n_pi_squared) * (q_squared + m_pi_squared)) * F(q, n) * (1 + (-1)**(n + m)) / 2
    
    return term1 - factor * term2


@jit(nopython=True)
def Q(q, n, m):
    """Q function for dipolar exchange spin waves - Optimized with numba"""
    # Avoid unnecessary calculations
    if m == n and m != 0:  # Avoiding division by zero
        return 0
        
    # Precompute common values
    q_squared = q**2
    m_pi_squared = m**2 * np.pi**2
    n_pi_squared = n**2 * np.pi**2
    nm_parity = (-1)**(n+m)
    
    factor = 1 / np.sqrt((1 + (1 if n == 0 else 0)) * (1 + (1 if m == 0 else 0)))
    
    if m**2 == n**2 and nm_parity == -1:
        # Handle special case to avoid division by zero
        term1 = 0
    else:
        term1 = q_squared / (q_squared + m_pi_squared) * (m**2 / (m**2 - n**2 + (1 + nm_parity) / 2) * 2 / q)
    
    term2 = q_squared / (q_squared + m_pi_squared) * (q_squared / (2 * (q_squared + n_pi_squared)) * F(q, n))
    
    return (term1 - term2) * factor * (1 - nm_parity) / 2


@jit(nopython=True)
def Omega(omega_H, q, n):
    """Omega function for spin waves - Optimized with numba"""
    return (omega_H + (gamma * DD) / thickness**2 * (q**2 + n**2 * np.pi**2)) / omega_M


def H(omega_H, q, phi_k, n, m):
    """Hamiltonian matrix for spin waves"""
    omega_term = Omega(omega_H, q, n) * np.array([[1, 0], [0, 1]])
    identity_term = 0.5 * np.array([[1, 1], [1, 1]])
    
    base_term = (omega_term + identity_term) * (1 if n == m else 0)
    
    sin_phi_squared = np.sin(phi_k)**2
    P_matrix = 0.5 * np.array([
        [1 - sin_phi_squared, 1 + sin_phi_squared],
        [1 + sin_phi_squared, 1 - sin_phi_squared]
    ])
    
    Q_matrix = 0.5 * np.array([
        [0, -4],
        [4, 0]
    ])
    
    P_term = P_matrix * P(q, n, m)
    Q_term = Q_matrix * np.sin(phi_k) * Q(q, n, m)
    
    return base_term - P_term - Q_term


@jit(nopython=True)
def f(q, n, NV_height_by_thickness):
    """f function for coupling calculations - Optimized with numba"""
    n_parity = (-1)**n
    return n_parity / np.sqrt(2 * (1 + (1 if n == 0 else 0))) * \
           q**2 / (q**2 + n**2 * np.pi**2) * \
           np.exp(-q * NV_height_by_thickness) * \
           (1 - n_parity * np.exp(-q))


def generate_H_BdG(omega_H, q, phi_k, Nmax):
    """Generate Bogoliubov-de Gennes Hamiltonian matrix"""
    H_BdG = np.zeros((2*Nmax, 2*Nmax), dtype=complex)
    
    # Vectorized approach - generate all H matrices at once
    H_matrices = np.array([[H(omega_H, q, phi_k, n, m) for m in range(Nmax)] for n in range(Nmax)])
    
    # Fill the Hamiltonian more efficiently
    for m in range(Nmax):
        for n in range(Nmax):
            H_BdG[n, m] = H_matrices[n, m][0, 0]
            H_BdG[n, m + Nmax] = H_matrices[n, m][0, 1]
            H_BdG[n + Nmax, m] = H_matrices[n, m][1, 0]
            H_BdG[n + Nmax, m + Nmax] = H_matrices[n, m][1, 1]
    
    return H_BdG


def cholesky_decomposition(H):
    """Thin wrapper for Cholesky decomposition, using scipy's optimized version"""
    return linalg.cholesky(H)


def para_unitary_diag(H):
    """
    Paraunitary diagonalization of a Hamiltonian matrix - Optimized
    """
    # Convert H to NumPy array if not already and ensure it's complex
    H = np.array(H, dtype=complex)
    
    # Get dimensions
    Dim = H.shape[0]
    half_dim = Dim // 2
    
    # Get Cholesky decomposition
    K = cholesky_decomposition(H)
    
    # Create sigma3 matrix - use precomputation for the diagonal
    sigma3_diag = np.array([(-1)**np.floor((2*n-1)/Dim) for n in range(1, Dim+1)])
    sigma3 = np.diag(sigma3_diag)
    
    # Calculate W more efficiently
    W = K @ sigma3 @ np.conjugate(K).T
    
    # Eigensystem
    evals, evec = linalg.eig(W)
    
    # Normalize eigenvectors (vectorized)
    norms = np.sqrt(np.sum(np.abs(evec)**2, axis=0))
    evec = evec / norms
    
    # Create permutation
    preperm = list(range(half_dim, Dim)) + list(range(half_dim, 0, -1))
    ordering = np.argsort(evals)
    permutation = [preperm[i] for i in ordering]
    
    # Apply permutation
    evals = evals[permutation]
    evec = evec[:, permutation]
    
    # Calculate U and other matrices
    U = evec.T
    Hdiag = sigma3 @ np.diag(evals)
    T = np.linalg.inv(K) @ U @ np.sqrt(Hdiag)
    
    # Extract submatrices
    Tpp = T[:half_dim, :half_dim]        # Upper left
    Tnn = T[half_dim:, half_dim:]        # Lower right
    
    # Calculate phase arrays
    PhaseArrayP = np.exp(1j * np.angle(np.diag(Tpp)))
    PhaseArrayN = np.exp(1j * np.angle(np.diag(Tnn)))
    
    # Create V matrix
    V = np.diag(np.concatenate([np.conjugate(PhaseArrayP), np.conjugate(PhaseArrayN)]))
    
    # Update T
    T = T @ V
    
    # Extract final submatrices
    Tpp = T[:half_dim, :half_dim]      # Upper left
    Tnp = T[half_dim:, :half_dim]      # Lower left
    Tpn = T[:half_dim, half_dim:]      # Upper right
    Tnn = T[half_dim:, half_dim:]      # Lower right
    
    return (evals[:half_dim], Tpp, Tnp, Tpn, Tnn, T)


def process_phi_slice(args):
    """Function to process one phi slice for parallelization"""
    countq, count_phi, q, phi_k, omega_H, Nmax, NV_height_by_thickness = args
    
    # Generate and diagonalize Hamiltonian
    H_BdG = generate_H_BdG(omega_H, q, phi_k, Nmax)
    result = para_unitary_diag((H_BdG + H_BdG.conj().T)/2)
    
    # Extract results
    omega_BdG = result[0]
    Tpp, Tnp, Tpn, Tnn = result[1], result[2], result[3], result[4]
    
    # Calculate gamma values
    sin_phi = np.sin(phi_k)
    sin_phi_pi = np.sin(phi_k + np.pi)
    
    gamma_plus = (1 + sin_phi)/2 * (Tpp + Tnp + sin_phi * (Tpp - Tnp))
    gamma_plus_mirror = (1 + sin_phi_pi)/2 * np.conjugate(Tnn + Tpn + sin_phi_pi * (Tnn - Tpn))
    
    gamma_minus = (1 - sin_phi)/2 * (Tpp + Tnp + sin_phi * (Tpp - Tnp))
    gamma_minus_mirror = (1 - sin_phi_pi)/2 * np.conjugate(Tnn + Tpn + sin_phi_pi * (Tnn - Tpn))
    
    gamma_z = -1j * np.cos(phi_k)/2 * (Tpp + Tnp + sin_phi * (Tpp - Tnp))
    gamma_z_mirror = -1j * np.cos(phi_k + np.pi)/2 * np.conjugate(Tnn + Tpn + sin_phi_pi * (Tnn - Tpn))
    
    # Calculate coupling values using vector operations
    fbar_array = np.array([f(q, nn, NV_height_by_thickness) for nn in range(Nmax)])
    
    vplus_array = np.dot(fbar_array, gamma_plus)
    vplus_mirror_array = np.dot(fbar_array, gamma_plus_mirror)
    vminus_array = np.dot(fbar_array, gamma_minus)
    vminus_mirror_array = np.dot(fbar_array, gamma_minus_mirror)
    vz_array = np.dot(fbar_array, gamma_z)
    vz_mirror_array = np.dot(fbar_array, gamma_z_mirror)
    
    return (countq, count_phi, omega_BdG, 
            vplus_array, vplus_mirror_array, 
            vminus_array, vminus_mirror_array,
            vz_array, vz_mirror_array)


# ---- Important change: We move the interpolation outside the multiprocessing parts ----

class InterpolationResult:
    """Class to hold the interpolation data and compute interpolated values"""
    def __init__(self, qtable, phi_table_mod, omega_BdG_table, coupling_plus_table, 
                 coupling_minus_table, coupling_Z_table, Nmax):
        self.qtable = qtable
        self.phi_table_mod = phi_table_mod
        self.omega_BdG_table = omega_BdG_table
        self.coupling_plus_table = coupling_plus_table
        self.coupling_minus_table = coupling_minus_table
        self.coupling_Z_table = coupling_Z_table
        self.Nmax = Nmax
        
        # Create interpolation functions
        self.omega_interp = []
        self.coupling_plus_interp = []
        self.coupling_minus_interp = []
        self.coupling_Z_interp = []
        
        for s in range(Nmax):
            points = (qtable, phi_table_mod)
            
            # Extract the data for this mode
            omega_values = omega_BdG_table[:, :, s]
            coupling_plus_values = coupling_plus_table[:, :, s].real
            coupling_minus_values = coupling_minus_table[:, :, s].real
            coupling_Z_values = coupling_Z_table[:, :, s].real
            
            # Create interpolation functions
            self.omega_interp.append(
                RegularGridInterpolator(points, omega_values, method='linear', bounds_error=False, fill_value=None)
            )
            self.coupling_plus_interp.append(
                RegularGridInterpolator(points, coupling_plus_values, method='linear', bounds_error=False, fill_value=None)
            )
            self.coupling_minus_interp.append(
                RegularGridInterpolator(points, coupling_minus_values, method='linear', bounds_error=False, fill_value=None)
            )
            self.coupling_Z_interp.append(
                RegularGridInterpolator(points, coupling_Z_values, method='linear', bounds_error=False, fill_value=None)
            )
        
        self.min_omega = np.min(omega_BdG_table)
    
    def get_omega(self, mode, q, phi):
        """Get interpolated omega value for given mode, q and phi"""
        return self.omega_interp[mode-1](np.array([q, phi]))
    
    def get_coupling_plus(self, mode, q, phi):
        """Get interpolated coupling_plus value for given mode, q and phi"""
        return self.coupling_plus_interp[mode-1](np.array([q, phi]))
    
    def get_coupling_minus(self, mode, q, phi):
        """Get interpolated coupling_minus value for given mode, q and phi"""
        return self.coupling_minus_interp[mode-1](np.array([q, phi]))
    
    def get_coupling_Z(self, mode, q, phi):
        """Get interpolated coupling_Z value for given mode, q and phi"""
        return self.coupling_Z_interp[mode-1](np.array([q, phi]))


def multi_para_diag(NV_height, omega_H, qtable, phi_table, Nmax):
    """Execute paraunitary diagonalization with parallelization"""
    # Initialize dimensions
    Numq = len(qtable)
    Num_phi = len(phi_table)
    NV_height_by_thickness = NV_height / thickness
    
    # Initialize tables
    omega_BdG_table = np.zeros((Numq, 2*Num_phi, Nmax))
    coupling_plus_table = np.zeros((Numq, 2*Num_phi, Nmax), dtype=complex)
    coupling_minus_table = np.zeros((Numq, 2*Num_phi, Nmax), dtype=complex)
    coupling_Z_table = np.zeros((Numq, 2*Num_phi, Nmax), dtype=complex)
    
    # Prepare arguments for parallel processing
    arg_list = []
    for countq in range(Numq):
        for count_phi in range(Num_phi):
            q = qtable[countq]
            phi_k = phi_table[count_phi]
            arg_list.append((countq, count_phi, q, phi_k, omega_H, Nmax, NV_height_by_thickness))
    
    # Determine number of processes
    num_cores = mp.cpu_count()
    
    # Run parallel computation
    print(f"Starting parallel computation using {num_cores} cores")
    with mp.Pool(processes=num_cores) as pool:
        results = list(tqdm(pool.imap(process_phi_slice, arg_list), total=len(arg_list), desc="MultiParaDiag"))
    
    # Collect results
    for result in results:
        countq, count_phi, omega, vplus, vplus_mirror, vminus, vminus_mirror, vz, vz_mirror = result
        
        # Store results
        omega_BdG_table[countq, count_phi] = omega
        omega_BdG_table[countq, count_phi + Num_phi] = omega
        
        coupling_plus_table[countq, count_phi] = vplus
        coupling_plus_table[countq, count_phi + Num_phi] = vplus_mirror
        coupling_minus_table[countq, count_phi] = vminus
        coupling_minus_table[countq, count_phi + Num_phi] = vminus_mirror
        coupling_Z_table[countq, count_phi] = vz
        coupling_Z_table[countq, count_phi + Num_phi] = vz_mirror
    
    # Create modified phi table for interpolation
    phi_table_mod = np.concatenate([phi_table, phi_table + np.pi])
    
    # Create interpolation result object
    interp_result = InterpolationResult(
        qtable, phi_table_mod, omega_BdG_table,
        coupling_plus_table, coupling_minus_table, coupling_Z_table, Nmax
    )
    
    return interp_result


@jit(nopython=True)
def nbose(omega):
    """Bose-Einstein distribution function - Optimized with numba."""
    epsilon = 1e-9
    # return (1e-9 * k_B * temperature / h_planck) / max(omega, epsilon)
    return 1 / (np.exp(((h_planck / (2 * np.pi)) * 1e9*max(omega, epsilon)) / (k_B * temperature)) - 1)


# Process a batch of q, phi points for a given mode
def process_mode_batch(mode, q_batch, phi_batch, delta_q_batch, 
                      eta_small, omega_target_l, omega_target_u, interp_result):
    """Process a batch of points for a single mode"""
    # Initialize contributions for this batch
    gamma_l_contribution = 0
    gamma_u_contribution = 0
    
    # Get interpolated values for this batch
    batch_size = len(q_batch)
    omega_values = np.zeros(batch_size)
    coupling_plus_values = np.zeros(batch_size)
    coupling_minus_values = np.zeros(batch_size)
    
    for i in range(batch_size):
        omega_values[i] = interp_result.get_omega(mode, q_batch[i], phi_batch[i])
        coupling_plus_values[i] = interp_result.get_coupling_plus(mode, q_batch[i], phi_batch[i])
        coupling_minus_values[i] = interp_result.get_coupling_minus(mode, q_batch[i], phi_batch[i])
    
    # Calculate contributions
    for i in range(batch_size):

        # Bose factor
        nbose_factor = 2 * nbose(omega_M * omega_values[i]) + 1
        
        # Denominators
        denom_l = eta_small**2 + (omega_M * omega_values[i] - omega_target_l)**2
        denom_u = eta_small**2 + (omega_M * omega_values[i] - omega_target_u)**2
        
        # Accumulate contributions
        gamma_l_contribution += nbose_factor * delta_q_batch[i] * q_batch[i] * \
                               (coupling_plus_values[i]**2) * (eta_small/np.pi) / denom_l
        
        gamma_u_contribution += nbose_factor * delta_q_batch[i] * q_batch[i] * \
                               (coupling_minus_values[i]**2) * (eta_small/np.pi) / denom_u
    
    return gamma_l_contribution, gamma_u_contribution


def process_mode(args):
    """Process a single mode - for parallel execution"""
    mode, q_batches, phi_batches, delta_q_batches, eta_small, \
    omega_target_l, omega_target_u, interp_result = args
    
    gamma_l_contribution = 0
    gamma_u_contribution = 0
    
    # Process each batch
    num_batches = len(q_batches)
    for b in range(num_batches):
        l_contrib, u_contrib = process_mode_batch(
            mode, q_batches[b], phi_batches[b], delta_q_batches[b],
            eta_small, omega_target_l, omega_target_u, interp_result
        )
        gamma_l_contribution += l_contrib
        gamma_u_contribution += u_contrib
    
    return gamma_l_contribution, gamma_u_contribution


def gamma_values_ul(eta_small, nq_est, n_phi, h0, interp_result, width, length, nv_height):
    """
    Calculate Gamma+ and Gamma- values using parallel processing.
    """
    dnv = 2.87
    omega_target_l = (dnv - h0 * gamma)  
    omega_target_u = (dnv + h0 * gamma)
    
    # Calculate omega_dwl (characteristic frequency) - constant factor moved to global
    volume_factor = 1.0 / ((thickness * 1e-6) * (width * 1e-6) * (length * 1e-6))
    omega_dwl = CONST_FACTOR * volume_factor
    
    # Set up q-space grid more efficiently
    q_middle = 5 * thickness
    nq_half = round(nq_est / 2)
    q_table = np.concatenate(([1e-6], fspace(thickness/1000, q_middle, nq_half)))
    
    # Calculate remaining q values more efficiently
    last_delta_q = q_table[-1] - q_table[-2]
    q_max = 50 * thickness
    q_extension = np.arange(q_table[-1] + last_delta_q, q_max + last_delta_q, last_delta_q)
    q_table = np.concatenate((q_table, q_extension))

    # Calculate differentials
    delta_q = np.diff(q_table)
    nq = len(delta_q)
    delta_phi = 2 * np.pi / n_phi
    
    # Create flattened tables for vectorized calculations
    flat_q_table = np.tile(q_table[:-1], n_phi)
    flat_delta_q_table = np.tile(delta_q, n_phi)
    flat_phi_table = np.repeat(np.linspace(0, 2*np.pi*(1-1/n_phi), n_phi), nq)
    
    # Split the flattened tables into batches for efficient processing
    batch_size = 1000
    num_points = len(flat_q_table)
    num_batches = (num_points + batch_size - 1) // batch_size
    
    # Create batches
    q_batches = []
    phi_batches = []
    delta_q_batches = []
    
    for b in range(num_batches):
        start_idx = b * batch_size
        end_idx = min((b + 1) * batch_size, num_points)
        
        q_batches.append(flat_q_table[start_idx:end_idx])
        phi_batches.append(flat_phi_table[start_idx:end_idx])
        delta_q_batches.append(flat_delta_q_table[start_idx:end_idx])
    
    # Calculate maximum mode number
    fmax = 5
    nmax = int(np.ceil(thickness / np.pi * np.sqrt(fmax / (gamma * DD))))
    
    print(f"Processing {nmax} modes in parallel...")
    
    # Set up parallel processing
    num_cores = mp.cpu_count()
    
    # Create arguments for each mode
    mode_args = []
    for mode in range(1, nmax+1):
        mode_args.append((
            mode, q_batches, phi_batches, delta_q_batches,
            eta_small, omega_target_l, omega_target_u, interp_result
        ))
    
    # Run parallel computation for each mode
    with mp.Pool(processes=num_cores) as pool:
        results = list(tqdm(pool.imap(process_mode, mode_args), 
                          total=nmax, desc="Processing magnon modes"))
    
    # Aggregate results
    gamma_l_hz_array = 0
    gamma_u_hz_array = 0
    
    for gamma_l_contribution, gamma_u_contribution in results:
        # Combine results with scaling factors
        common_factor = (2 * np.pi)**2 * omega_M * omega_dwl * (delta_phi / (2 * np.pi)**2)
        gamma_l_hz_array += gamma_l_contribution * common_factor
        gamma_u_hz_array += gamma_u_contribution * common_factor
    
    return gamma_l_hz_array, gamma_u_hz_array


def fspace(min_val, max_val, steps, f=np.log):
    """Create logarithmically spaced values."""
    inv_f = lambda x: np.exp(x)  # Inverse function for log
    return inv_f(np.linspace(f(min_val), f(max_val), steps))


def gamma1_from_h0_width_length(h0, width, length, nv_height, multi_para_diag_func=None):
    """
    Calculate Gamma1 for given H0, width, length, and NV height with caching.
    """
    # Use provided function or default to our implementation
    if multi_para_diag_func is None:
        multi_para_diag_func = multi_para_diag
        
    omega_h = gamma * h0  # GHz
    
    # More efficient setup for calculation
    num_phi = 2 * 90  # Can be reduced for faster computation if needed
    num_q = 2 * 10    # Can be reduced for faster computation if needed
    
    del_phi = np.pi / num_phi
    phi_table = np.arange(0, np.pi, del_phi)
    
    q_max = 50 * thickness
    q_table = np.concatenate(([1e-6], fspace(thickness/1000, q_max, num_q)))
    
    # Calculate maximum mode based on desired frequency
    fmax = 5
    nmax = int(np.ceil(thickness / np.pi * np.sqrt(fmax / (gamma * DD))))
    
    # Call the multi-parameter diagonalization function
    print(f"Starting diagonalization for h0={h0}, width={width}, length={length}, height={nv_height}")
    interp_result = multi_para_diag_func(nv_height, omega_h, q_table, phi_table, nmax)
    
    # Parameters for gamma calculation
    eta_small = 0.003  # Can be adjusted for performance/accuracy tradeoff
    nq = 2 * 10        # Can be reduced for faster computation if needed
    n_phi = 2 * 360    # Can be reduced for faster computation if needed
    
    # Calculate the Gamma values
    print(f"Calculating gamma values...")
    gamma_l_hz_array, gamma_u_hz_array = gamma_values_ul(
        eta_small, nq, n_phi, h0, interp_result, 
        width, length, nv_height
    )
    
    print(f"Γ(ω=ωL)= {gamma_l_hz_array:.2f} Hz")
    print(f"Γ(ω=ωU)= {gamma_u_hz_array:.2f} Hz")
    
    return gamma_l_hz_array, gamma_u_hz_array



# Example usage
if __name__ == "__main__":
    # Parameter arrays - can use fewer points for faster testing
    width_array = [3]  # μm
    length_array = [3]  # μm
    nv_height_array = [0.4, 0.5, 0.6, 0.7]  # μm
    
    # Create a 3D array of T1 values
    t1_upper_table = np.zeros((len(width_array), len(length_array), len(nv_height_array)))
    
    # Calculate T1 values for each combination
    for i, w in enumerate(width_array):
        for j, l in enumerate(length_array):
            for k, h in enumerate(nv_height_array):
                # Calculate Gamma values for H0 = 82 Oersted
                gamma_l, _ = gamma1_from_h0_width_length(
                    82, w, l, h
                )
                # Convert to T1 (in microseconds)
                t1_upper_table[i, j, k] = 1000000 / (gamma_l)
    
    print("T1 Upper Table (μs):")
    print(t1_upper_table)
    

    # Plot 1/T1 vs hNV
    figtitle = 'Decay rate vs h_NV for 82G and d=w=l=3um'
    plt.title(figtitle)
    plt.plot(nv_height_array, t1_upper_table[0][0])
    plt.xlabel('NV-YIG distance (um)')
    plt.ylabel('1/T1 (Hz)')
    plt.tight_layout()
    plt.savefig(figtitle)
    plt.show()


    # # Visualize as a 3D plot
    # try:
    #     from mpl_toolkits.mplot3d import Axes3D
    #     from matplotlib import cm
        
    #     # Create figure
    #     fig = plt.figure(figsize=(12, 10))
    #     ax = fig.add_subplot(111, projection='3d')
        
    #     # Create a colormap of the data
    #     for i in range(len(width_array)):
    #         for j in range(len(length_array)):
    #             for k in range(len(nv_height_array)):
    #                 size = t1_upper_table[i, j, k] / np.max(t1_upper_table) * 200
    #                 color = cm.rainbow(t1_upper_table[i, j, k] / np.max(t1_upper_table))
    #                 ax.scatter(width_array[i], length_array[j], nv_height_array[k], 
    #                           c=[color], s=size, alpha=0.7)
        
    #     # Add a color bar
    #     norm = plt.Normalize(np.min(t1_upper_table), np.max(t1_upper_table))
    #     sm = plt.cm.ScalarMappable(cmap=cm.rainbow, norm=norm)
    #     sm.set_array([])
    #     cbar = plt.colorbar(sm)
    #     cbar.set_label('T1 (μs)')
        
    #     # Label the axes
    #     ax.set_xlabel('Waveguide width (μm)')
    #     ax.set_ylabel('Waveguide length (μm)')
    #     ax.set_zlabel('NV-YIG separation (μm)')
        
    #     plt.title('T1 Values vs Device Parameters')
    #     plt.tight_layout()
    #     plt.show()
        
    # except ImportError:
    #     print("3D plotting requires mpl_toolkits and matplotlib. Skipping visualization.")
 


