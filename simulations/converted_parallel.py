import numpy as np
import scipy.linalg as la
from scipy import interpolate
from numba import jit, prange
import multiprocessing as mp
from tqdm import tqdm
import matplotlib.pyplot as plt
from functools import partial
import time

# ----------------------- Constants and parameters -----------------------
# Material & condition parameters
def initialize_parameters(debug=False):
    # Calculate Ms value similar to MsSol in the original code
    delta_h = 1
    ms_sol = ((2870/2.8) - (1716 + 4 - 2)/2) / 2 - 82 - delta_h
    
    # Set global parameters
    params = {
        'M0': 1716.0,  # Oe (value from MsSol calculation)
        'L': 3.0,      # μm
        'DD': 5.4e-9 * 1e8,  # Oe μm^-2
        'gamma': 2.8e-3,  # GHz/G
        'omega_M': None,  # Will be calculated as gamma*M0
        
        # Physical constants
        'h_plank': 6.626e-34,  # J*s
        'mu0': 4*np.pi*1e-7,  # H/m
        'kB': 1.381e-23,  # J/K
        'temperature': 300,  # K
        
        # Calculation parameters - reduced for debug mode
        'num_phi': 2*90, #if not debug else 2*30,
        'num_q': 2*10, # if not debug else 2*50,
        'qmax': None,  # Will be set to 50*L or 20*L in debug
        'fmax': 5,     # GHz
        'nmax': None,  # Will be calculated or reduced in debug
    }
    
    # Calculate dependent parameters
    params['omega_M'] = params['gamma'] * params['M0']
    params['qmax'] = 50 * params['L'] if not debug else 20 * params['L']
    
    if debug:
        # Use smaller nmax for debugging
        params['nmax'] = min(5, int(np.ceil(params['L']/np.pi * np.sqrt(params['fmax']/(params['gamma'] * params['DD'])))))
    else:
        params['nmax'] = int(np.ceil(params['L']/np.pi * np.sqrt(params['fmax']/(params['gamma'] * params['DD']))))
    
    # Calculate omega_d
    params['omega_d'] = (params['h_plank'] * params['mu0'] * 
                        (params['gamma']*1e9)**2 / 
                        (params['L']*1e-6)**3) * 1e8  # Hz
    
    return params

# ------------------------ Helper functions ------------------------
@jit(nopython=True)
def n_bose(omega, temperature, h_plank, kB):
    """Calculate Bose-Einstein distribution"""
    # Using high temperature approximation as in the original code
    return (1e-9 * kB * temperature / h_plank) / omega

@jit(nopython=True)
def f_function(q, n):
    """Calculate F(q,n) function"""
    if q < 1e-10:  # Handle q→0 limit
        return 0.0
    return 2 * (1 - (-1)**n * np.exp(-q)) / q

@jit(nopython=True)
def p_function(q, n, m):
    """Calculate P(q,n,m) function"""
    if q < 1e-10:  # Handle q→0 limit
        if n == m:
            return 1.0
        return 0.0
    
    result = 0.0
    if n == m:
        result = q**2 / (q**2 + n**2 * np.pi**2)
    
    normalization = 1.0 / np.sqrt((1 + (1 if n == 0 else 0)) * (1 + (1 if m == 0 else 0)))
    
    f_val = f_function(q, n)
    common_factor = q**4 / ((q**2 + n**2 * np.pi**2) * (q**2 + m**2 * np.pi**2))
    parity_factor = (1 + (-1)**(n+m)) / 2
    
    result -= normalization * common_factor * f_val * parity_factor
    
    return result

@jit(nopython=True)
def q_function(q, n, m):
    """Calculate Q(q,n,m) function"""
    if q < 1e-10 or n == m:  # Handle q→0 limit or n=m case
        return 0.0
    
    if (n + m) % 2 == 0:  # Check if n+m is even
        return 0.0
    
    normalization = 1.0 / np.sqrt((1 + (1 if n == 0 else 0)) * (1 + (1 if m == 0 else 0)))
    
    # Calculate the term with careful handling of division
    if m**2 == n**2:
        return 0.0
    
    term1 = q**2 / (q**2 + m**2 * np.pi**2)
    term2 = m**2 / (m**2 - n**2 + (1 + (-1)**(n+m)) / 2)
    term3 = 2/q - q**2 / (2 * (q**2 + n**2 * np.pi**2)) * f_function(q, n)
    
    parity_factor = (1 - (-1)**(n+m)) / 2
    
    return term1 * term2 * term3 * normalization * parity_factor

@jit(nopython=True)
def omega_function(omega_h, q, n, gamma_dd, L, omega_M):
    """Calculate Ω(ωH,q,n) function"""
    return (omega_h + (gamma_dd / L**2) * (q**2 + n**2 * np.pi**2)) / omega_M

# Worker function for multiprocessing, defined outside the main function
def worker_for_multi_para_diag(task, h_nv_array, omega_h, nmax, params):
    q_idx, phi_idx, q, phi_k = task
    result = calculate_for_single_q_phi(q, phi_k, h_nv_array, omega_h, nmax, params)
    return q_idx, phi_idx, result

@jit(nopython=True)
def h_matrix(omega_h, q, phi_k, n, m, gamma_dd, L, omega_M):
    """Calculate H(ωH,q,ϕk,n,m) 2x2 matrix"""
    if n == m:
        # Diagonal part
        omega_val = omega_function(omega_h, q, n, gamma_dd, L, omega_M)
        h00 = omega_val + 0.5
        h01 = 0.5
        h10 = 0.5
        h11 = omega_val + 0.5
    else:
        h00 = h11 = h01 = h10 = 0.0
    
    # Add the P and Q terms
    sin_phi_k = np.sin(phi_k)
    sin_phi_k_sq = sin_phi_k**2
    
    p_val = p_function(q, n, m)
    q_val = q_function(q, n, m)
    
    h00 -= 0.5 * (1 - sin_phi_k_sq) * p_val
    h11 -= 0.5 * (1 - sin_phi_k_sq) * p_val
    h01 -= 0.5 * (1 + sin_phi_k_sq) * p_val - 2 * sin_phi_k * q_val
    h10 -= 0.5 * (1 + sin_phi_k_sq) * p_val + 2 * sin_phi_k * q_val
    
    return np.array([[h00, h01], [h10, h11]])

@jit(nopython=True)
def generate_h_bdg(omega_h, q, phi_k, nmax, gamma_dd, L, omega_M):
    """Generate the Bogoliubov-de Gennes Hamiltonian matrix"""
    dim = 2 * nmax
    h_bdg = np.zeros((dim, dim), dtype=np.complex128)
    
    for m in range(nmax):
        for n in range(nmax):
            h_2x2 = h_matrix(omega_h, q, phi_k, n, m, gamma_dd, L, omega_M)
            h_bdg[n, m] = h_2x2[0, 0]
            h_bdg[n, m + nmax] = h_2x2[0, 1]
            h_bdg[n + nmax, m] = h_2x2[1, 0]
            h_bdg[n + nmax, m + nmax] = h_2x2[1, 1]
    
    return h_bdg

@jit(nopython=True)
def f_coupling(q, n, h_over_L):
    """Calculate f(q,n,h/L) coupling function"""
    if q < 1e-10:  # Handle q→0 limit
        if n == 0:
            return 1.0
        return 0.0
    
    normalization = (-1)**n / np.sqrt(2 * (1 + (1 if n == 0 else 0)))
    term = q**2 / (q**2 + n**2 * np.pi**2) * np.exp(-q * h_over_L) * (1 - (-1)**n * np.exp(-q))
    
    return normalization * term

# ------------------------ Para-unitary diagonalization ------------------------
# Fix 1: Add stability to para_unitary_diag function
def para_unitary_diag(h):
    """Perform para-unitary diagonalization of a Hamiltonian matrix"""
    dim = h.shape[0]
    half_dim = dim // 2
    
    # Ensure h is Hermitian
    h = (h + h.conj().T) / 2
    
    # Add stability to the matrix
    epsilon = 1e-8
    h = h + np.eye(dim) * epsilon
    
    # Cholesky decomposition with robust error handling
    try:
        k = la.cholesky(h)
    except np.linalg.LinAlgError:
        # If matrix is not positive definite, add larger offset to diagonal
        print("Warning: Matrix not positive definite, adding offset")
        h_modified = h + np.eye(dim) * 1e-6
        try:
            k = la.cholesky(h_modified)
        except np.linalg.LinAlgError:
            # If still failing, use eigendecomposition approach
            print("Warning: Cholesky still failing, using eigendecomposition")
            evalss, evecs = la.eigh(h)
            # Ensure all eigenvalues are positive
            evalss = np.maximum(evalss, 1e-10)
            # Reconstruct h with positive eigenvalues
            h_modified = evecs @ np.diag(evalss) @ evecs.conj().T
            k = la.cholesky(h_modified)
    
    # Create sigma_3 matrix
    sigma_3 = np.diag([(-1)**int(np.floor((2*n-1)/dim)) for n in range(1, dim+1)])
    
    # Calculate W matrix
    w = k @ sigma_3 @ k.conj().T
    
    # Compute eigensystem
    try:
        evals, evec = la.eigh(w)
    except np.linalg.LinAlgError:
        print("Warning: Eigendecomposition failed, adding stability")
        w = (w + w.conj().T)/2  # Ensure hermitian
        w = w + np.eye(dim) * 1e-6  # Add stability
        evals, evec = la.eigh(w)
    
    # Check for NaN values in eigenvalues/vectors
    if np.any(np.isnan(evals)) or np.any(np.isnan(evec)):
        print("Warning: NaN values in eigendecomposition, using more stable approach")
        # Alternative approach
        w = (w + w.conj().T)/2  # Ensure hermitian
        w = w + np.eye(dim) * 1e-4  # Add stability
        evals, evec = la.eigh(w)
    
    # Normalize eigenvectors with stability
    for i in range(dim):
        norm = np.abs(evec[:, i].conj() @ evec[:, i])
        if norm > 1e-10:  # Avoid division by very small numbers
            evec[:, i] = evec[:, i] / np.sqrt(norm)
        else:
            print(f"Warning: Near-zero norm for eigenvector {i}")
            evec[:, i] = np.zeros_like(evec[:, i])
            evec[i, i] = 1.0
    
    # Create permutation to order eigenvalues as (+small, -small, +mid, -mid, ...)
    preperm = np.concatenate([np.arange(half_dim, dim), np.arange(half_dim-1, -1, -1)])
    ordering = np.argsort(evals)
    permutation = np.empty_like(preperm)
    for i, idx in enumerate(preperm):
        permutation[i] = ordering[idx]
    
    # Apply permutation
    evals = evals[permutation]
    evec = evec[:, permutation]
    
    # Create U matrix
    u = evec.T
    
    # Create diagonal Hamiltonian
    h_diag = sigma_3 @ np.diag(evals)
    
    # Create transformation matrix T with stability
    try:
        # Use pseudo-inverse for more stability
        k_inv = np.linalg.pinv(k, rcond=1e-10)
        sqrt_h_diag = np.zeros_like(h_diag)
        np.fill_diagonal(sqrt_h_diag, np.sqrt(np.maximum(np.abs(np.diag(h_diag)), 1e-10)))
        t = k_inv @ u @ sqrt_h_diag
    except np.linalg.LinAlgError:
        print("Warning: Matrix inversion failed, using different approach")
        # Add stability to k before inversion
        k_mod = k + np.eye(dim) * 1e-6
        k_inv = np.linalg.pinv(k_mod, rcond=1e-8)
        sqrt_h_diag = np.zeros_like(h_diag)
        np.fill_diagonal(sqrt_h_diag, np.sqrt(np.maximum(np.abs(np.diag(h_diag)), 1e-10)))
        t = k_inv @ u @ sqrt_h_diag
    
    # Check for NaN values in T
    if np.any(np.isnan(t)):
        print("Warning: NaN values in transformation matrix")
        # Replace NaNs with zeros
        t = np.nan_to_num(t)
    
    # Fix phases with stability
    tpp = t[:half_dim, :half_dim]
    tnn = t[half_dim:, half_dim:]
    
    # Safe angle calculation
    def safe_angle(z):
        if np.abs(z) < 1e-10:
            return 0.0
        return np.angle(z)
    
    phase_array_p = np.array([np.exp(1j * safe_angle(tpp[i, i])) for i in range(tpp.shape[0])])
    phase_array_n = np.array([np.exp(1j * safe_angle(tnn[i, i])) for i in range(tnn.shape[0])])
    
    v = np.diag(np.concatenate([np.conjugate(phase_array_p), np.conjugate(phase_array_n)]))
    t = t @ v
    
    # Get submatrices
    tpp = t[:half_dim, :half_dim]
    tnp = t[half_dim:, :half_dim]
    tpn = t[:half_dim, half_dim:]
    tnn = t[half_dim:, half_dim:]
    
    # Final check for NaNs
    if (np.any(np.isnan(evals[:half_dim])) or np.any(np.isnan(tpp)) or 
        np.any(np.isnan(tnp)) or np.any(np.isnan(tpn)) or np.any(np.isnan(tnn))):
        print("Warning: NaN values in final output")
        # Replace NaNs with zeros
        evals = np.nan_to_num(evals)
        tpp = np.nan_to_num(tpp)
        tnp = np.nan_to_num(tnp)
        tpn = np.nan_to_num(tpn)
        tnn = np.nan_to_num(tnn)
    
    return evals[:half_dim], tpp, tnp, tpn, tnn, t

# ------------------------ Main calculation functions ------------------------
def calculate_for_single_q_phi(q, phi_k, h_nv_array, omega_h, nmax, params):
    """Calculate for a single (q, phi) point"""
    h_bdg = generate_h_bdg(omega_h, q, phi_k, nmax, params['DD'], params['L'], params['omega_M'])
    
    # Para-unitary diagonalization
    evalss, tpp, tnp, tpn, tnn, _ = para_unitary_diag(h_bdg)
    
    # Calculate coupling parameters
    num_hNV = len(h_nv_array)
    omega_bdg = evalss
    
    # For phi_k
    sin_phi = np.sin(phi_k)
    gamma_plus = (1 + sin_phi)/2 * (tpp + tnp + sin_phi * (tpp - tnp))
    gamma_minus = (1 - sin_phi)/2 * (tpp + tnp + sin_phi * (tpp - tnp))
    gamma_z = -1j * np.cos(phi_k)/2 * (tpp + tnp + sin_phi * (tpp - tnp))
    
    # For phi_k + pi
    sin_phi_pi = np.sin(phi_k + np.pi)
    gamma_plus_mirror = (1 + sin_phi_pi)/2 * np.conjugate(tnn + tpn + sin_phi_pi * (tnn - tpn))
    gamma_minus_mirror = (1 - sin_phi_pi)/2 * np.conjugate(tnn + tpn + sin_phi_pi * (tnn - tpn))
    gamma_z_mirror = -1j * np.cos(phi_k + np.pi)/2 * np.conjugate(tnn + tpn + sin_phi_pi * (tnn - tpn))
    
    # Calculate coupling for all NV heights
    coupling_plus = np.zeros((num_hNV, nmax), dtype=np.complex128)
    coupling_plus_mirror = np.zeros((num_hNV, nmax), dtype=np.complex128)
    coupling_minus = np.zeros((num_hNV, nmax), dtype=np.complex128)
    coupling_minus_mirror = np.zeros((num_hNV, nmax), dtype=np.complex128)
    coupling_z = np.zeros((num_hNV, nmax), dtype=np.complex128)
    coupling_z_mirror = np.zeros((num_hNV, nmax), dtype=np.complex128)
    
    for tt in range(num_hNV):
        f_bar_array = np.array([f_coupling(q, nn, h_nv_array[tt]/params['L']) for nn in range(nmax)])
        
        for nn in range(nmax):
            coupling_plus[tt, nn] = np.sum(f_bar_array * gamma_plus[:, nn])
            coupling_plus_mirror[tt, nn] = np.sum(f_bar_array * gamma_plus_mirror[:, nn])
            coupling_minus[tt, nn] = np.sum(f_bar_array * gamma_minus[:, nn])
            coupling_minus_mirror[tt, nn] = np.sum(f_bar_array * gamma_minus_mirror[:, nn])
            coupling_z[tt, nn] = np.sum(f_bar_array * gamma_z[:, nn])
            coupling_z_mirror[tt, nn] = np.sum(f_bar_array * gamma_z_mirror[:, nn])
    
    return (omega_bdg, 
            coupling_plus, coupling_plus_mirror,
            coupling_minus, coupling_minus_mirror,
            coupling_z, coupling_z_mirror)

def multi_para_diag(h_nv_array, omega_h, q_table, phi_table, nmax, params):
    """Perform para-unitary diagonalization for multiple q and phi values with multiprocessing"""
    num_q = len(q_table)
    num_phi = len(phi_table)
    num_hNV = len(h_nv_array)
    
    # Prepare output arrays
    omega_bdg_table = np.zeros((num_q, 2*num_phi, nmax))
    coupling_plus_table = np.zeros((num_q, 2*num_phi, num_hNV, nmax), dtype=np.complex128)
    coupling_minus_table = np.zeros((num_q, 2*num_phi, num_hNV, nmax), dtype=np.complex128)
    coupling_z_table = np.zeros((num_q, 2*num_phi, num_hNV, nmax), dtype=np.complex128)
    
    # Prepare task list for multiprocessing
    tasks = []
    for q_idx in range(num_q):
        for phi_idx in range(num_phi):
            q = q_table[q_idx]
            phi_k = phi_table[phi_idx]
            tasks.append((q_idx, phi_idx, q, phi_k))
    
    # Execute tasks with multiprocessing
    worker_partial = partial(worker_for_multi_para_diag, h_nv_array=h_nv_array, omega_h=omega_h, nmax=nmax, params=params)
    with mp.Pool(processes=mp.cpu_count()) as pool:
        results = list(tqdm(pool.imap(worker_partial, tasks), total=len(tasks), desc="Para-unitary diagonalization"))
    
    # Process results
    for q_idx, phi_idx, result in results:
        omega_bdg, cp, cp_m, cm, cm_m, cz, cz_m = result
        
        # Store results for phi_k
        omega_bdg_table[q_idx, phi_idx] = omega_bdg
        for tt in range(num_hNV):
            coupling_plus_table[q_idx, phi_idx, tt] = cp[tt]
            coupling_minus_table[q_idx, phi_idx, tt] = cm[tt]
            coupling_z_table[q_idx, phi_idx, tt] = cz[tt]
        
        # Store results for phi_k + pi
        omega_bdg_table[q_idx, phi_idx + num_phi] = omega_bdg
        for tt in range(num_hNV):
            coupling_plus_table[q_idx, phi_idx + num_phi, tt] = cp_m[tt]
            coupling_minus_table[q_idx, phi_idx + num_phi, tt] = cm_m[tt]
            coupling_z_table[q_idx, phi_idx + num_phi, tt] = cz_m[tt]
    
    # Create interpolation functions
    phi_table_mod = np.array([phi_table[0] + i*np.pi for i in range(2)] + [2*np.pi])
    
    int_omega_bdg_table = []
    int_coupling_plus_table = []
    int_coupling_minus_table = []
    int_coupling_z_table = []
    
    # Create points for interpolation
    points = []
    for q_idx in range(num_q):
        for phi_idx in range(2*num_phi):
            points.append([q_table[q_idx], (phi_idx % num_phi) * (np.pi/num_phi)])
    
    points = np.array(points)
    
    # Create interpolations for each mode
    for s in range(nmax):
        # Flatten the data for interpolation
        omega_values = omega_bdg_table[:, :, s].flatten()
        int_omega_bdg = interpolate.LinearNDInterpolator(points, omega_values)
        int_omega_bdg_table.append(int_omega_bdg)
        
        # Create interpolations for coupling parameters
        mode_coupling_plus = []
        mode_coupling_minus = []
        mode_coupling_z = []
        
        for tt in range(num_hNV):
            cp_values = coupling_plus_table[:, :, tt, s].flatten()
            cm_values = coupling_minus_table[:, :, tt, s].flatten()
            cz_values = coupling_z_table[:, :, tt, s].flatten()
            
            int_cp = interpolate.LinearNDInterpolator(points, cp_values)
            int_cm = interpolate.LinearNDInterpolator(points, cm_values)
            int_cz = interpolate.LinearNDInterpolator(points, cz_values)
            
            mode_coupling_plus.append(int_cp)
            mode_coupling_minus.append(int_cm)
            mode_coupling_z.append(int_cz)
        
        int_coupling_plus_table.append(mode_coupling_plus)
        int_coupling_minus_table.append(mode_coupling_minus)
        int_coupling_z_table.append(mode_coupling_z)
    
    return (q_table, np.min(omega_bdg_table), phi_table_mod, 
            int_omega_bdg_table, int_coupling_plus_table, 
            int_coupling_minus_table, int_coupling_z_table)

def gamma_values_ul(eta_small, nq_est, n_phi, h0, int_omega_bdg_table, 
                   int_coupling_plus_table, int_coupling_minus_table, num_hnv, params):
    """Calculate Gamma values for upper and lower transitions with enhanced numerical stability"""
    dnv = 2.87
    omega_target_l = dnv - h0 * params['gamma']
    omega_target_u = dnv + h0 * params['gamma']
    
    # Create q table with more points around the middle
    q_middle = 5 * params['L']
    nq_half = nq_est // 2
    
    # Create q table with logarithmic spacing - ensure minimum value is positive
    min_q = max(params['L']/1000, 1e-6)
    q_low = np.logspace(np.log10(min_q), np.log10(q_middle), nq_half)
    q_table = np.concatenate([[1e-6], q_low])
    
    # Add linear spacing for higher q values
    last_delta_q = max(q_table[-1] - q_table[-2], 1e-6)  # Ensure positive delta
    q_high = np.arange(q_table[-1] + last_delta_q, params['qmax'] + last_delta_q, last_delta_q)
    q_table = np.concatenate([q_table, q_high])
    
    delta_q = np.diff(q_table)
    nq = len(delta_q)
    delta_phi = 2 * np.pi / n_phi
    
    gamma_l_hz_array = np.zeros(num_hnv)
    gamma_u_hz_array = np.zeros(num_hnv)
    dos_l = 0.0
    dos_u = 0.0
    
    # Calculate contribution from each mode
    for s in range(params['nmax']):
        print(f"Processing mode {s+1}/{params['nmax']}")
        
        # Create mesh grid for q and phi
        q_mesh, phi_mesh = np.meshgrid(q_table[:-1], np.linspace(0, 2*np.pi, n_phi, endpoint=False))
        flat_q = q_mesh.flatten()
        flat_phi = phi_mesh.flatten()
        
        # Calculate step sizes for integration
        delta_q_mesh, _ = np.meshgrid(delta_q, np.linspace(0, 2*np.pi, n_phi, endpoint=False))
        flat_delta_q = delta_q_mesh.flatten()
        
        # Evaluate dispersion and coupling at mesh points
        flat_omega_bdg = np.zeros_like(flat_q)
        flat_coupling_plus = np.zeros((num_hnv, len(flat_q)), dtype=np.complex128)
        flat_coupling_minus = np.zeros((num_hnv, len(flat_q)), dtype=np.complex128)
        
        # Evaluate with error handling
        for i in range(len(flat_q)):
            try:
                result = int_omega_bdg_table[s](flat_q[i], flat_phi[i])
                # Check if result is NaN and replace with fallback value
                if np.isnan(result):
                    flat_omega_bdg[i] = 1.0  # Fallback value
                    print(f"Warning: NaN in omega_bdg at q={flat_q[i]}, phi={flat_phi[i]}")
                else:
                    flat_omega_bdg[i] = result
                    
                for tt in range(num_hnv):
                    cp_result = int_coupling_plus_table[s][tt](flat_q[i], flat_phi[i])
                    cm_result = int_coupling_minus_table[s][tt](flat_q[i], flat_phi[i])
                    
                    # Check for NaNs
                    if np.isnan(cp_result) or np.isinf(cp_result):
                        flat_coupling_plus[tt, i] = 0.0
                    else:
                        flat_coupling_plus[tt, i] = cp_result
                        
                    if np.isnan(cm_result) or np.isinf(cm_result):
                        flat_coupling_minus[tt, i] = 0.0
                    else:
                        flat_coupling_minus[tt, i] = cm_result
            except Exception as e:
                print(f"Error at q={flat_q[i]}, phi={flat_phi[i]}: {e}")
                flat_omega_bdg[i] = 1.0  # Fallback value
                for tt in range(num_hnv):
                    flat_coupling_plus[tt, i] = 0.0
                    flat_coupling_minus[tt, i] = 0.0
        
        # Avoid division by zero in Lorentzian
        denominator_l = eta_small**2 + (params['omega_M'] * flat_omega_bdg - omega_target_l)**2
        denominator_u = eta_small**2 + (params['omega_M'] * flat_omega_bdg - omega_target_u)**2
        
        # Ensure denominators are not too small
        denominator_l = np.maximum(denominator_l, 1e-10)
        denominator_u = np.maximum(denominator_u, 1e-10)
        
        # Calculate DOS contributions
        lorentzian_l = (eta_small/np.pi) / denominator_l
        lorentzian_u = (eta_small/np.pi) / denominator_u
        
        # Replace NaNs with zeros
        lorentzian_l = np.nan_to_num(lorentzian_l)
        lorentzian_u = np.nan_to_num(lorentzian_u)
        
        # Integrate DOS
        dos_l += (1/params['L']**2) * (delta_phi/(2*np.pi)**2) * np.sum(flat_delta_q * flat_q * lorentzian_l)
        dos_u += (1/params['L']**2) * (delta_phi/(2*np.pi)**2) * np.sum(flat_delta_q * flat_q * lorentzian_u)
        
        # Calculate occupation number (Bose-Einstein) with protection against invalid values
        n_bose_vals = np.zeros_like(flat_omega_bdg)
        for i, freq in enumerate(flat_omega_bdg):
            if freq > 0:  # Protect against non-positive frequencies
                n_bose_vals[i] = n_bose(params['omega_M'] * freq, params['temperature'], 
                                       params['h_plank'], params['kB'])
            else:
                n_bose_vals[i] = 0.0
        
        # Replace NaNs and infinities
        n_bose_vals = np.nan_to_num(n_bose_vals)
        
        # Calculate coupling contributions
        for tt in range(num_hnv):
            # Lower transition - use abs to avoid complex issues
            coupling_term_l = (2*n_bose_vals + 1) * flat_delta_q * flat_q * np.abs(flat_coupling_plus[tt])**2 * lorentzian_l
            gamma_l_hz_array[tt] += (2*np.pi)**2 * params['omega_M'] * params['omega_d'] * (delta_phi/(2*np.pi)**2) * np.sum(coupling_term_l)
            
            # Upper transition
            coupling_term_u = (2*n_bose_vals + 1) * flat_delta_q * flat_q * np.abs(flat_coupling_minus[tt])**2 * lorentzian_u
            gamma_u_hz_array[tt] += (2*np.pi)**2 * params['omega_M'] * params['omega_d'] * (delta_phi/(2*np.pi)**2) * np.sum(coupling_term_u)
    
    # Final check for NaNs
    gamma_l_hz_array = np.nan_to_num(gamma_l_hz_array)
    gamma_u_hz_array = np.nan_to_num(gamma_u_hz_array)
    dos_l = np.nan_to_num(dos_l)
    dos_u = np.nan_to_num(dos_u)
    
    return dos_l, dos_u, gamma_l_hz_array, gamma_u_hz_array

def gamma1_from_h0(h0, h_nv_array, params):
    """Calculate Gamma1 for a given external field H0 with debug output"""
    print(f"\n--- Calculating for H0 = {h0} Oe ---")
    omega_h = params['gamma'] * h0  # GHz
    
    # Create q and phi tables
    del_phi = np.pi / params['num_phi']
    phi_table = np.arange(0, np.pi, del_phi)
    
    # Create q table with logarithmic spacing
    q_table = np.concatenate([[1e-6], np.logspace(np.log10(params['L']/1000), np.log10(params['qmax']), params['num_q'])])
    
    print(f"q range: {q_table[0]} to {q_table[-1]}")
    print(f"phi range: {phi_table[0]} to {phi_table[-1]}")
    print(f"omega_h = {omega_h} GHz")
    print(f"nmax = {params['nmax']}")
    
    # Calculate dispersion and coupling
    print("Starting multi_para_diag calculation...")
    result = multi_para_diag(h_nv_array, omega_h, q_table, phi_table, params['nmax'], params)
    
    q_table_temp, omega_min, phi_table_mod, int_omega_bdg_table, int_coupling_plus_table, int_coupling_minus_table, int_coupling_z_table = result
    
    print(f"multi_para_diag completed. omega_min = {omega_min}")
    
    # Check for NaNs in interpolation functions
    test_q = q_table[len(q_table)//2]
    test_phi = phi_table[len(phi_table)//2]
    print(f"Testing interpolation at q={test_q}, phi={test_phi}:")
    print(f"  omega_bdg = {int_omega_bdg_table[0](test_q, test_phi)}")
    print(f"  coupling_plus = {int_coupling_plus_table[0][0](test_q, test_phi)}")
    
    # Calculate Gamma values
    print("Starting gamma_values_ul calculation...")
    eta_small = 0.003
    nq = 2 * 10  # Reduced for testing
    n_phi = 2 * 360  # Reduced for testing
    
    dos_l, dos_u, gamma_l_hz_array, gamma_u_hz_array = gamma_values_ul(
        eta_small, nq, n_phi, h0, int_omega_bdg_table, 
        int_coupling_plus_table, int_coupling_minus_table, 
        len(h_nv_array), params
    )
    
    print(f"DOS(ω=ωL) = {dos_l:.6f} 1/GHz μm²")
    print(f"DOS(ω=ωU) = {dos_u:.6f} 1/GHz μm²")
    print(f"Γ(ω=ωL) = {gamma_l_hz_array} Hz")
    print(f"Γ(ω=ωU) = {gamma_u_hz_array} Hz")
    
    return dos_l, dos_u, gamma_l_hz_array, gamma_u_hz_array

def run_calculation(h0_array, h_nv_array):
    """Run calculations for an array of H0 values"""
    params = initialize_parameters()
    
    # Initialize result arrays
    n_iterations = len(h0_array)
    dos_l_array = np.zeros(n_iterations)
    dos_u_array = np.zeros(n_iterations)
    gamma_l_hz_array = np.zeros((len(h_nv_array), n_iterations))
    gamma_u_hz_array = np.zeros((len(h_nv_array), n_iterations))
    
    # Process each H0 value, potentially in parallel
    for ii, h0 in enumerate(h0_array):
        print(f"Processing {ii+1}/{n_iterations}: H0 = {h0} Oe")
        dos_l, dos_u, gamma_l_hz, gamma_u_hz = gamma1_from_h0(h0, h_nv_array, params)
        
        dos_l_array[ii] = dos_l
        dos_u_array[ii] = dos_u
        gamma_l_hz_array[:, ii] = gamma_l_hz
        gamma_u_hz_array[:, ii] = gamma_u_hz
    
    return dos_l_array, dos_u_array, gamma_l_hz_array, gamma_u_hz_array

def plot_results(h0_array, h_nv_array, dos_l_array, dos_u_array, gamma_l_hz_array, gamma_u_hz_array):
    """Plot results of the calculation"""
    # Create figure with multiple subplots
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Spin Wave Relaxation Analysis', fontsize=16)
    
    # Plot DOS for lower and upper transitions
    axs[0, 0].semilogy(h0_array, dos_l_array, 'b-', label='DOS Lower')
    axs[0, 0].semilogy(h0_array, dos_u_array, 'r-', label='DOS Upper')
    axs[0, 0].set_xlabel('H0 (Oe)')
    axs[0, 0].set_ylabel('DOS (1/GHz μm²)')
    axs[0, 0].set_title('Density of States')
    axs[0, 0].legend()
    axs[0, 0].grid(True, which='both', linestyle='--', alpha=0.6)
    
    # Plot Γ1 for lower transition
    for i, h_nv in enumerate(h_nv_array):
        axs[0, 1].semilogy(h0_array, gamma_l_hz_array[i], label=f'h_NV = {h_nv:.2f} μm')
    axs[0, 1].set_xlabel('H0 (Oe)')
    axs[0, 1].set_ylabel('Γ1 Lower (Hz)')
    axs[0, 1].set_title('Relaxation Rate - Lower Transition')
    axs[0, 1].legend()
    axs[0, 1].grid(True, which='both', linestyle='--', alpha=0.6)
    
    # Plot Γ1 for upper transition
    for i, h_nv in enumerate(h_nv_array):
        axs[1, 0].semilogy(h0_array, gamma_u_hz_array[i], label=f'h_NV = {h_nv:.2f} μm')
    axs[1, 0].set_xlabel('H0 (Oe)')
    axs[1, 0].set_ylabel('Γ1 Upper (Hz)')
    axs[1, 0].set_title('Relaxation Rate - Upper Transition')
    axs[1, 0].legend()
    axs[1, 0].grid(True, which='both', linestyle='--', alpha=0.6)
    
    # Plot ratio of upper/lower relaxation rates
    for i, h_nv in enumerate(h_nv_array):
        ratio = gamma_u_hz_array[i] / gamma_l_hz_array[i]
        axs[1, 1].plot(h0_array, ratio, label=f'h_NV = {h_nv:.2f} μm')
    axs[1, 1].set_xlabel('H0 (Oe)')
    axs[1, 1].set_ylabel('Γ1(Upper) / Γ1(Lower)')
    axs[1, 1].set_title('Ratio of Relaxation Rates')
    axs[1, 1].legend()
    axs[1, 1].grid(True, which='both', linestyle='--', alpha=0.6)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    return fig

def save_results(h0_array, h_nv_array, dos_l_array, dos_u_array, gamma_l_hz_array, gamma_u_hz_array, filename='spin_wave_results'):
    """Save calculation results to a file"""
    np.savez(
        filename,
        h0_array=h0_array,
        h_nv_array=h_nv_array,
        dos_l_array=dos_l_array,
        dos_u_array=dos_u_array,
        gamma_l_hz_array=gamma_l_hz_array,
        gamma_u_hz_array=gamma_u_hz_array
    )
    print(f"Results saved to {filename}.npz")

def load_results(filename='spin_wave_results'):
    """Load previously saved calculation results"""
    data = np.load(f"{filename}.npz")
    return (
        data['h0_array'], 
        data['h_nv_array'], 
        data['dos_l_array'], 
        data['dos_u_array'], 
        data['gamma_l_hz_array'], 
        data['gamma_u_hz_array']
    )

def main():
    """Main function to execute the calculation pipeline"""
    # Initialize parameters in debug mode (reduced complexity)
    params = initialize_parameters(debug=True)
    print("Parameters initialized (DEBUG MODE):")
    for key, value in params.items():
        print(f"  {key}: {value}")
    
    # Define simplified test case
    h0_array = np.array([82])  # Just one H0 value
    h_nv_array = np.array([0.4]) * params['L']  # Just one NV height
    print(f"H0 value: {h0_array[0]} Oe")
    print(f"NV height: {h_nv_array[0]} μm")
    
    # Set multiprocessing start method
    mp.set_start_method('spawn', force=True)
    
    # Run calculation with extensive debug output
    print("Running calculation with debug output...")
    start_time = time.time()
    
    # Calculate directly without using run_calculation
    dos_l, dos_u, gamma_l_hz, gamma_u_hz = gamma1_from_h0(h0_array[0], h_nv_array, params)
    
    # Store in arrays
    dos_l_array = np.array([dos_l])
    dos_u_array = np.array([dos_u])
    gamma_l_hz_array = np.array([gamma_l_hz])
    gamma_u_hz_array = np.array([gamma_u_hz])
    
    end_time = time.time()
    print(f"Calculation completed in {end_time - start_time:.2f} seconds")
    
    # Print results
    print("\nFinal Results:")
    print(f"DOS Lower: {dos_l_array}")
    print(f"DOS Upper: {dos_u_array}")
    print(f"Gamma Lower: {gamma_l_hz_array}")
    print(f"Gamma Upper: {gamma_u_hz_array}")
    
    return h0_array, h_nv_array, dos_l_array, dos_u_array, gamma_l_hz_array, gamma_u_hz_array

if __name__ == "__main__":
    main()

