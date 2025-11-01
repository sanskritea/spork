import numpy as np
from scipy.linalg import cholesky, eig
from scipy.integrate import dblquad, quad
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import warnings
import time
import os

warnings.filterwarnings("ignore")

#=============================================================================
# SECTION 1: CONSTANTS & PARAMETERS
#=============================================================================

# Material parameters (from Mathematica PDF)
M0 = 2458  # Oe
DD = 5.39e-9 * 1e8  # Oe μm^-2
gamma = 2.8e-3  # GHz/G
omega_M = gamma * M0  # GHz

# Bar geometry (from Mathematica PDF - CRITICAL for correct results)
d_bar = 0.005  # μm (5 nm) - thickness
w_bar = 0.03   # μm (30 nm) - width
l_bar = 3.0    # μm - length

# Physical constants
h_plank = 6.626e-34  # J*s
mu_0 = 4 * np.pi * 1e-7  # H/m
k_B = 1.381e-23  # J/K
temperature = 300  # K

# CORRECTED unit conversion factor (uses actual bar volume, not L^3)
omega_dwl = (h_plank * mu_0 * (gamma * 1e9)**2 / 
             ((d_bar * 1e-6) * (w_bar * 1e-6) * (l_bar * 1e-6))) * 1e8  # Hz

print("="*60)
print("INITIALIZATION")
print("="*60)
print(f"Bar geometry: {d_bar*1000:.1f} nm × {w_bar*1000:.1f} nm × {l_bar*1000:.0f} nm")
print(f"ω_M = {omega_M:.4f} GHz")
print(f"ω_dwl = {omega_dwl:.2f} Hz (Mathematica: 1450.66 Hz)")
print(f"Temperature = {temperature} K")
print("="*60 + "\n")

#=============================================================================
# SECTION 2: LINEAR ALGEBRA UTILITIES
#=============================================================================

def scipy_cholesky_decomposition(A):
    """Cholesky decomposition with regularization"""
    try:
        eigenvals = np.linalg.eigvals(A)
        if np.min(eigenvals) <= 0:
            reg_param = 1e-10
            A_reg = A + reg_param * np.eye(A.shape[0])
        else:
            A_reg = A
        
        L = cholesky(A_reg, lower=False)
        return L
    
    except np.linalg.LinAlgError as e:
        print(f"Cholesky decomposition failed: {e}")
        eigenvals, eigenvecs = np.linalg.eigh(A)
        eigenvals = np.maximum(eigenvals, 1e-10)
        A_reconstructed = eigenvecs @ np.diag(eigenvals) @ eigenvecs.T
        return cholesky(A_reconstructed, lower=False)


def mathematica_eigensystem(matrix, tolerance=1e-10):
    """Replicate Mathematica's Eigensystem behavior"""
    matrix = np.array(matrix, dtype=float)
    
    eigenvals, eigenvecs = eig(matrix)
    
    if np.allclose(eigenvals.imag, 0, atol=tolerance):
        eigenvals = eigenvals.real
    if np.allclose(eigenvecs.imag, 0, atol=tolerance):
        eigenvecs = eigenvecs.real
    
    for i in range(eigenvecs.shape[1]):
        norm = np.linalg.norm(eigenvecs[:, i])
        if norm > tolerance:
            eigenvecs[:, i] = eigenvecs[:, i] / norm
    
    for i in range(eigenvecs.shape[1]):
        vec = eigenvecs[:, i]
        first_nonzero_idx = None
        for j in range(len(vec)):
            if abs(vec[j]) > tolerance:
                first_nonzero_idx = j
                break
        
        if first_nonzero_idx is not None:
            if vec[first_nonzero_idx] < 0:
                eigenvecs[:, i] = -eigenvecs[:, i]
    
    def sort_key(i):
        val = eigenvals[i]
        return (-abs(val), -val)
    
    indices = sorted(range(len(eigenvals)), key=sort_key)
    sorted_eigenvals = eigenvals[indices]
    sorted_eigenvecs = eigenvecs[:, indices]
    
    return sorted_eigenvals, sorted_eigenvecs


def paraunitary_diag(H):
    """Paraunitary diagonalization for BdG Hamiltonian"""
    K = scipy_cholesky_decomposition(H)
    dim = H.shape[0]
    
    sigma3 = np.diag([(-1) ** np.floor((2 * n - 1) / dim) for n in range(1, dim + 1)])
    W = K @ sigma3 @ K.conj().T
    
    eval_w, evec_w = mathematica_eigensystem(W)
    evec_w = evec_w.T
    evec_w = evec_w / np.linalg.norm(evec_w, axis=0)
    
    preperm = np.concatenate(
        [np.arange(dim // 2, dim), np.arange(dim // 2 - 1, -1, -1)]
    )
    ordering = np.argsort(eval_w)
    permutation = ordering[preperm]
    
    eval_w = eval_w[permutation]
    evec_w = evec_w[permutation]
    
    U = evec_w.T
    H_diag = sigma3 * np.diag(eval_w)
    
    T = np.linalg.inv(K) @ U @ np.sqrt(np.abs(H_diag))
    
    dim_half = dim // 2
    Tpp = T[:dim_half, :dim_half]
    Tnn = T[dim_half:, dim_half:]
    
    phase_array_p = np.exp(1j * np.angle(np.diag(Tpp)))
    phase_array_n = np.exp(1j * np.angle(np.diag(Tnn)))
    
    V = np.diag(np.concatenate([phase_array_p.conj(), phase_array_n.conj()]))
    T = T @ V
    
    Tpp = T[:dim_half, :dim_half]
    Tnp = T[dim_half:, :dim_half]  # NOTE: This is "Thp" in Mathematica
    Tpn = T[:dim_half, dim_half:]  # NOTE: This is "Tph" in Mathematica
    Tnn = T[dim_half:, dim_half:]
    
    return eval_w[:dim_half], Tpp, Tnp, Tpn, Tnn, T

#=============================================================================
# SECTION 3: DISCRETE BAR HAMILTONIAN
#=============================================================================

def FF(p1, p2, zdif, zave):
    """Helper for z-direction dipolar integrals"""
    if p1 == p2:
        if p1 == 0:
            return zave
        else:
            return (0.5 * zave * np.cos(p1 * np.pi * zdif) + 
                    np.sin(2 * p1 * np.pi * zave) / (4 * p1 * np.pi))
    else:
        num1 = np.sin((p1 + p2) * np.pi * zave + 0.5 * (p1 - p2) * np.pi * zdif)
        num2 = np.sin((p2 - p1) * np.pi * zave - 0.5 * (p1 + p2) * np.pi * zdif)
        return (num1 / (p1 + p2) - num2 / (p2 - p1)) / (2 * np.pi)


def Integrand1(ydif, zdif, p1, p2, w_norm, d_norm):
    """Dipolar integrand for HXX contribution"""
    term1 = 1.0 / (w_norm**2 * ydif**2 + zdif**2 + 1e-12)
    term2 = 1.0 / (d_norm**2 + w_norm**2 * ydif**2 + zdif**2 + 1e-12)
    FF_term = FF(p1, p2, zdif, 1 - abs(zdif)/2) - FF(p1, p2, zdif, abs(zdif)/2)
    return (term1 - term2) * FF_term * (1 - ydif)


def Integrand2(xdif, zdif, p1, p2, w_norm, d_norm):
    """Dipolar integrand for HYY contribution"""
    term1 = 1.0 / (d_norm**2 * xdif**2 + zdif**2 + 1e-12)
    term2 = 1.0 / (d_norm**2 * xdif**2 + w_norm**2 + zdif**2 + 1e-12)
    FF_term = FF(p1, p2, zdif, 1 - abs(zdif)/2) - FF(p1, p2, zdif, abs(zdif)/2)
    return (term1 - term2) * FF_term * (1 - xdif)


def generate_H_BdG_discrete_bar(omega_H, N_max, d_bar, w_bar, l_bar, use_cache=True):
    """
    Generate discrete bar BdG Hamiltonian - FULLY CORRECTED normalization
    """
    cache_file = f'H_BdG_cache_N{N_max}_H{omega_H:.3f}.npz'
    
    if use_cache and os.path.exists(cache_file):
        print(f"Loading cached Hamiltonian from {cache_file}")
        data = np.load(cache_file)
        return data['H_BdG']
    
    H_BdG = np.zeros((2*N_max, 2*N_max), dtype=float)
    d_norm = d_bar / l_bar
    w_norm = w_bar / l_bar
    
    print(f"Building discrete bar Hamiltonian (N_max={N_max})...")
    start_time = time.time()
    
    for p1 in range(N_max):
        for p2 in range(p1, N_max):
            
            if (p1 + p2) % 2 == 1:
                H_XX, H_YY = 0.0, 0.0
            else:
                # Compute H_XX
                try:
                    result_XX, _ = dblquad(
                        lambda zdif, ydif: Integrand1(ydif, zdif, p1, p2, w_norm, d_norm),
                        0, 1, -1, 1,
                        epsabs=1e-6, epsrel=1e-6
                    )
                    prefactor_XX = (2 * w_bar / (np.pi * d_bar)) * (1 + (p1==0)) * (1 + (p2==0))
                    H_XX = prefactor_XX * result_XX
                except:
                    H_XX = 0.0
                
                # Compute H_YY
                try:
                    result_YY, _ = dblquad(
                        lambda zdif, xdif: Integrand2(xdif, zdif, p1, p2, w_norm, d_norm),
                        0, 1, -1, 1,
                        epsabs=1e-6, epsrel=1e-6
                    )
                    prefactor_YY = (2 * d_bar / (np.pi * w_bar)) * (1 + (p1==0)) * (1 + (p2==0))
                    H_YY = prefactor_YY * result_YY
                except:
                    H_YY = 0.0
            
            n, m = p1, p2

            # ✅ FINAL CORRECTION: Dipolar integrals are ALREADY dimensionless!
            result_11 = (H_XX + H_YY) / 2  # NO omega_M division!
            result_12 = (H_XX - H_YY) / 2  # NO omega_M division!

            if p1 == p2:
                # Only the external terms need normalization
                exchange = gamma * DD * (p1 * np.pi / l_bar)**2 / omega_M  
                result_11 += omega_H / omega_M + exchange

            H_BdG[n, m] = result_11
            H_BdG[n, m + N_max] = result_12
            H_BdG[n + N_max, m] = result_12
            H_BdG[n + N_max, m + N_max] = result_11
            
            if p1 != p2:
                H_BdG[m, n] = result_11
                H_BdG[m, n + N_max] = result_12
                H_BdG[m + N_max, n] = result_12
                H_BdG[m + N_max, n + N_max] = result_11
        
        if (p1 + 1) % max(1, N_max // 5) == 0:
            elapsed = time.time() - start_time
            progress = (p1 + 1) / N_max
            eta = elapsed / progress * (1 - progress)
            print(f"  Progress: {progress*100:.1f}% | Elapsed: {elapsed/60:.1f}min | ETA: {eta/60:.1f}min")
    
    total_time = time.time() - start_time
    print(f"  Hamiltonian built in {total_time/60:.1f} minutes")
    
    if use_cache:
        np.savez(cache_file, H_BdG=H_BdG)
        print(f"  Cached to {cache_file}")
    
    return H_BdG


def add_demagnetization_corrections(H_BdG, N_max, d_bar, w_bar, l_bar, use_cache=True):
    """
    Add demagnetization corrections with caching
    """
    cache_file = f'Demag_cache_N{N_max}.npz'
    
    if use_cache and os.path.exists(cache_file):
        print("Loading cached demagnetization corrections...")
        data = np.load(cache_file)
        return H_BdG + data['H_demag']
    
    print("Computing demagnetization corrections...")
    start_time = time.time()
    
    d_norm = d_bar / l_bar
    w_norm = w_bar / l_bar
    
    def HdZMean_integrand(xdif, ydif, CoordZ):
        r1_sq = d_norm**2 * xdif**2 + w_norm**2 * ydif**2 + (CoordZ - 1)**2
        term1 = (CoordZ - 1) / (r1_sq**(3/2) + 1e-12)
        
        r0_sq = d_norm**2 * xdif**2 + w_norm**2 * ydif**2 + (CoordZ - 0)**2
        term2 = (CoordZ - 0) / (r0_sq**(3/2) + 1e-12)
        
        return (1 - xdif) * (1 - ydif) * (term1 - term2)
    
    def HdZMean(CoordZ):
        prefactor = (w_bar * d_bar) / (np.pi * l_bar**2)
        result, _ = dblquad(
            lambda ydif, xdif: HdZMean_integrand(xdif, ydif, CoordZ),
            0, 1, 0, 1,
            epsabs=1e-6, epsrel=1e-6
        )
        return prefactor * result
    
    NumEval = max(20, int(50 * l_bar / d_bar))
    z_list = np.linspace(0.01, 0.99, NumEval)
    HdZMean_values = np.array([HdZMean(z) for z in z_list])
    
    z_list_full = np.concatenate([[0], z_list, [1]])
    HdZMean_full = np.concatenate([[-0.5], HdZMean_values, [-0.5]])
    IntHdZMean = interp1d(z_list_full, HdZMean_full, kind='linear')
    
    Demag_matrix = np.zeros((N_max, N_max), dtype=float)
    
    for p1 in range(N_max):
        for p2 in range(p1, N_max):
            if (p1 + p2) % 2 == 1:
                Demag_matrix[p1, p2] = 0.0
            else:
                def demag_integrand(z):
                    norm1 = np.sqrt(2 / (1 + (p1 == 0)))
                    norm2 = np.sqrt(2 / (1 + (p2 == 0)))
                    return (IntHdZMean(z) * 
                            np.cos(p1 * np.pi * z) * 
                            np.cos(p2 * np.pi * z) * 
                            norm1 * norm2)
                
                result, _ = quad(demag_integrand, 0, 0.5, 
                                epsabs=1e-6, epsrel=1e-6)
                Demag_matrix[p1, p2] = 2 * result
                
                if p1 != p2:
                    Demag_matrix[p2, p1] = Demag_matrix[p1, p2]
    
    H_demag = np.zeros_like(H_BdG)
    
    for p1 in range(N_max):
        for p2 in range(N_max):
            H_demag[p1, p2] = Demag_matrix[p1, p2]
            H_demag[p1 + N_max, p2 + N_max] = Demag_matrix[p2, p1]
    
    total_time = time.time() - start_time
    print(f"  Demagnetization computed in {total_time/60:.1f} minutes")
    print(f"  Max demag correction: {np.max(np.abs(Demag_matrix)):.6f}")
    
    if use_cache:
        np.savez(cache_file, H_demag=H_demag)
        print(f"  Cached to {cache_file}")
    
    return H_BdG + H_demag

#=============================================================================
# SECTION 4: NV-MAGNON COUPLING
#=============================================================================

def calculate_coupling_at_position(coord_x, coord_y, coord_z, N_max, d_bar, w_bar, l_bar, Tpp, Tnp, mode_idx):
    """
    CORRECTED: Integration order matches Mathematica
    """
    from scipy.integrate import dblquad
    
    d_norm = d_bar / l_bar
    w_norm = w_bar / l_bar
    
    print(f"  NV normalized position: ({coord_x:.2f}, {coord_y:.2f}, {coord_z:.6f})")
    
    HalfGamma_mm_array = np.zeros(N_max, dtype=complex)
    HalfGamma_pm_array = np.zeros(N_max, dtype=complex)
    
    for p in range(N_max):
        # Mathematica normalization: sqrt(2/(1+delta_p0))
        norm_factor = np.sqrt(2 / (1 + (p == 0)))
        
        # ========================================
        # Gamma_XX: CORRECTED integration order
        # Mathematica: NIntegrate[..., {z, 0, 1}, {y, 0, 1}]
        # scipy dblquad: ∫∫ f(y,z) dy dz where z is outer, y is inner
        # ========================================
        def integrand_XX(y, z):  # NOTE: y is first arg (inner integral)
            r_top = np.sqrt((d_norm * (coord_x - 1))**2 + (w_norm * (coord_y - y))**2 + (coord_z - z)**2)
            r_bot = np.sqrt((d_norm * (coord_x - 0))**2 + (w_norm * (coord_y - y))**2 + (coord_z - z)**2)
            
            field = (coord_x - 1) / (r_top**3 + 1e-16) - (coord_x - 0) / (r_bot**3 + 1e-16)
            mode = np.cos(p * np.pi * z) * norm_factor
            
            return field * mode
        
        # CORRECTED: z from 0 to 1 (outer), y from 0 to 1 (inner)
        result_XX, _ = dblquad(integrand_XX, 0, 1, lambda z: 0, lambda z: 1, epsabs=1e-8, epsrel=1e-6)
        Gamma_XX = (d_bar * w_bar / (4 * np.pi * l_bar**2)) * result_XX
        
        # ========================================
        # Gamma_XY: CORRECTED
        # Mathematica: NIntegrate[..., {z, 0, 1}, {x, 0, 1}]
        # ========================================
        def integrand_XY(x, z):
            r_right = np.sqrt((d_norm * (coord_x - x))**2 + (w_norm * (coord_y - 1))**2 + (coord_z - z)**2)
            r_left = np.sqrt((d_norm * (coord_x - x))**2 + (w_norm * (coord_y - 0))**2 + (coord_z - z)**2)
            
            field = (coord_x - x) / (r_right**3 + 1e-16) - (coord_x - x) / (r_left**3 + 1e-16)
            mode = np.cos(p * np.pi * z) * norm_factor
            
            return field * mode
        
        result_XY, _ = dblquad(integrand_XY, 0, 1, lambda z: 0, lambda z: 1, epsabs=1e-8, epsrel=1e-6)
        Gamma_XY = (d_bar**2 / (4 * np.pi * l_bar**2)) * result_XY
        
        # ========================================
        # Gamma_YX: CORRECTED
        # ========================================
        def integrand_YX(y, z):
            r_top = np.sqrt((d_norm * (coord_x - 1))**2 + (w_norm * (coord_y - y))**2 + (coord_z - z)**2)
            r_bot = np.sqrt((d_norm * (coord_x - 0))**2 + (w_norm * (coord_y - y))**2 + (coord_z - z)**2)
            
            field = (coord_y - y) / (r_top**3 + 1e-16) - (coord_y - y) / (r_bot**3 + 1e-16)
            mode = np.cos(p * np.pi * z) * norm_factor
            
            return field * mode
        
        result_YX, _ = dblquad(integrand_YX, 0, 1, lambda z: 0, lambda z: 1, epsabs=1e-8, epsrel=1e-6)
        Gamma_YX = (w_bar**2 / (4 * np.pi * l_bar**2)) * result_YX
        
        # ========================================
        # Gamma_YY: CORRECTED
        # ========================================
        def integrand_YY(x, z):
            r_right = np.sqrt((d_norm * (coord_x - x))**2 + (w_norm * (coord_y - 1))**2 + (coord_z - z)**2)
            r_left = np.sqrt((d_norm * (coord_x - x))**2 + (w_norm * (coord_y - 0))**2 + (coord_z - z)**2)
            
            field = (coord_y - 1) / (r_right**3 + 1e-16) - (coord_y - 0) / (r_left**3 + 1e-16)
            mode = np.cos(p * np.pi * z) * norm_factor
            
            return field * mode
        
        result_YY, _ = dblquad(integrand_YY, 0, 1, lambda z: 0, lambda z: 1, epsabs=1e-8, epsrel=1e-6)
        Gamma_YY = (d_bar * w_bar / (4 * np.pi * l_bar**2)) * result_YY
        
        # Combine
        HalfGamma_mm_array[p] = (Gamma_XX - Gamma_YY - 1j * (Gamma_XY + Gamma_YX)) / 2
        HalfGamma_pm_array[p] = (Gamma_XX + Gamma_YY - 1j * (Gamma_XY - Gamma_YX)) / 2
        
        if p <= 2 or p == 5:
            print(f"    p={p}: ΓXX={Gamma_XX:.6e}, ΓYY={Gamma_YY:.6e}")
            print(f"          |HalfΓ_mm|={np.abs(HalfGamma_mm_array[p]):.6e}")
    
    # Transform
    coupling = (HalfGamma_mm_array @ Tpp[:, mode_idx] + 
                HalfGamma_pm_array @ Tnp[:, mode_idx])
    
    print(f"  Mode {mode_idx}: |g|_dimensionless = {np.abs(coupling):.6e}")
    
    return coupling


def n_bose_vectorized(omega):
    """Bose-Einstein distribution (high temperature approximation)"""
    return (1e-9 * k_B * temperature / h_plank) / omega

#=============================================================================
# SECTION 5: T1 CALCULATIONS (CORRECTED)
#=============================================================================

def calculate_all_mode_coupling(h_NV, H0, N_max, d_bar, w_bar, l_bar, position='center', Tpp=None, Tnp=None):
    """
    Calculate coupling to all modes at specific NV position
    
    Parameters:
    -----------
    h_NV : float, NV height above bar (μm)
    H0 : float, magnetic field (Oe)
    N_max : int, number of modes
    position : str or tuple
        - 'center': y=w_bar/2, z=l_bar/2
        - tuple (y, z): specific position in μm
    Tpp, Tnp : arrays, transformation matrices (optional, will compute if None)
    
    Returns:
    --------
    coupling_all_modes : array, coupling to each mode
    eigenfreqs_GHz : array, mode frequencies (GHz)
    Tpp, Tnp : transformation matrices (for reuse)
    """
    # Build and diagonalize if not provided
    if Tpp is None or Tnp is None:
        omega_H = gamma * H0
        H_BdG = generate_H_BdG_discrete_bar(omega_H, N_max, d_bar, w_bar, l_bar)
        H_BdG = add_demagnetization_corrections(H_BdG, N_max, d_bar, w_bar, l_bar)
        
        eigenfreqs, Tpp, Tnp, Tpn, Tnn, T = paraunitary_diag((H_BdG + H_BdG.conj().T) / 2)
        eigenfreqs_GHz = eigenfreqs * omega_M
    else:
        # If matrices provided, we need eigenfreqs
        # This branch should include eigenfreqs as input too
        raise ValueError("If Tpp/Tnp provided, must also provide eigenfreqs_GHz")
    
    # Determine NV position
    if position == 'center':
        coord_x = (d_bar + h_NV) / d_bar
        coord_y = 0.5
        coord_z = 0.5
    else:
        y_pos, z_pos = position
        coord_x = (d_bar + h_NV) / d_bar
        coord_y = (y_pos + w_bar/2) / w_bar
        coord_z = (z_pos + l_bar/2) / l_bar
    
    # Calculate coupling to all modes
    coupling_all_modes = np.zeros(N_max, dtype=complex)
    
    print(f"Calculating coupling for all {N_max} modes at position {position}...")
    start = time.time()
    
    for mode_idx in range(N_max):
        coupling_all_modes[mode_idx] = calculate_coupling_at_position(
            coord_x, coord_y, coord_z,
            N_max, d_bar, w_bar, l_bar,
            Tpp, Tnp, mode_idx
        )
        
        if (mode_idx + 1) % 10 == 0:
            elapsed = time.time() - start
            rate = (mode_idx + 1) / elapsed
            eta = (N_max - mode_idx - 1) / rate
            print(f"  Mode {mode_idx + 1}/{N_max} | Rate: {rate:.1f} modes/s | ETA: {eta/60:.1f} min")
    
    return coupling_all_modes, eigenfreqs_GHz, Tpp, Tnp


def calculate_T1_finite_bar_at_position(h_NV, H0, N_max, position='center', 
                                        coupling_all_modes=None, eigenfreqs_GHz=None,
                                        Tpp=None, Tnp=None):
    """Calculate T1 for finite bar at specific position - CORRECTED"""
    
    # Get coupling if not provided
    if coupling_all_modes is None or eigenfreqs_GHz is None:
        coupling_all_modes, eigenfreqs_GHz, Tpp, Tnp = calculate_all_mode_coupling(
            h_NV, H0, N_max, d_bar, w_bar, l_bar, position
        )
    
    # NV transition frequencies
    DNV = 2.87  # GHz
    omega_target_L = DNV - H0 * gamma
    omega_target_U = DNV + H0 * gamma
    
    # Lorentzian linewidth
    eta_small = 0.003  # GHz
    
    print(f"\nNV transition frequencies:")
    print(f"  ω_L = {omega_target_L:.4f} GHz (lower transition)")
    print(f"  ω_U = {omega_target_U:.4f} GHz (upper transition)")
    print(f"  Linewidth η = {eta_small:.3f} GHz")
    
    # Initialize relaxation rates
    gamma_L = 0.0
    gamma_U = 0.0
    
    print("\nCalculating T1 contributions from each mode...")
    contributions_L = []
    contributions_U = []
    
    for mode_idx in range(N_max):
        omega_mode = eigenfreqs_GHz[mode_idx]
        
        # Lorentzian lineshape at NV frequency
        lorentzian_L = (eta_small / np.pi) / (eta_small**2 + (omega_mode - omega_target_L)**2)
        lorentzian_U = (eta_small / np.pi) / (eta_small**2 + (omega_mode - omega_target_U)**2)
        
        # Bose factor
        bose = 2 * n_bose_vectorized(omega_mode) + 1
        
        # ✅ CORRECTED coupling strength formula (matches Mathematica exactly)
        g_coupling_kHz = np.abs(coupling) * np.sqrt(omega_dwl * 1e-3 * omega_M * 1e6)
        g_coupling_Hz = g_coupling_kHz * 1e3  # Convert to Hz for T1 calculation
        
        # Fermi's Golden Rule contribution
        contrib_L = g_coupling_Hz**2 * bose * lorentzian_L
        contrib_U = g_coupling_Hz**2 * bose * lorentzian_U
        
        contributions_L.append(contrib_L)
        contributions_U.append(contrib_U)
        
        gamma_L += contrib_L
        gamma_U += contrib_U
    
    # Show breakdown
    print(f"\n{'='*70}")
    print(f"Top 10 contributing modes to Gamma_L:")
    print(f"{'='*70}")
    print(f"{'Mode':<6} {'Freq (GHz)':<12} {'|g| (Hz)':<12} {'Lorentzian':<12} {'Contrib (Hz)':<12} {'%':<6}")
    print(f"{'-'*70}")
    
    top_indices_L = np.argsort(contributions_L)[::-1][:10]
    for idx in top_indices_L:
        g_Hz = np.abs(coupling_all_modes[idx]) * np.sqrt(omega_dwl * omega_M * 1e6)
        lorentz_L = (eta_small / np.pi) / (eta_small**2 + (eigenfreqs_GHz[idx] - omega_target_L)**2)
        print(f"{idx:<6} {eigenfreqs_GHz[idx]:<12.5f} {g_Hz:<12.2e} {lorentz_L:<12.2e} "
              f"{contributions_L[idx]:<12.2e} {contributions_L[idx]/gamma_L*100:<6.1f}")
    
    print(f"{'='*70}")
    print(f"Total Gamma_L: {gamma_L:.2e} Hz")
    print(f"Total Gamma_U: {gamma_U:.2e} Hz")
    print(f"T1_L: {1/gamma_L:.2e} s")
    print(f"T1_U: {1/gamma_U:.2e} s")
    print(f"{'='*70}\n")
    
    return gamma_L, gamma_U, eigenfreqs_GHz, coupling_all_modes, Tpp, Tnp


def calculate_T1_spatial_scan(h_NV, H0, N_max, y_positions, z_positions, 
                               Tpp=None, Tnp=None, eigenfreqs_GHz=None,
                               save_every=1, checkpoint_file=None):
    """
    Scan T1 over Y-Z spatial grid
    
    With checkpointing for long calculations
    """
    num_y = len(y_positions)
    num_z = len(z_positions)
    
    gamma_L_map = np.zeros((num_y, num_z))
    gamma_U_map = np.zeros((num_y, num_z))
    coupling_map = np.zeros((num_y, num_z, N_max), dtype=complex)
    
    # Setup checkpoint
    if checkpoint_file is None:
        checkpoint_file = f'T1_scan_checkpoint_H{H0:.1f}Oe_N{N_max}.npz'
    
    # Try to load existing progress
    start_iy = 0
    if os.path.exists(checkpoint_file):
        try:
            data = np.load(checkpoint_file)
            gamma_L_map = data['gamma_L_map']
            gamma_U_map = data['gamma_U_map']
            coupling_map = data['coupling_map']
            start_iy = int(data['last_iy']) + 1
            print(f"✓ Resuming from checkpoint: row {start_iy}/{num_y}")
        except:
            print(f"⚠ Could not load checkpoint, starting from beginning")
    
    # Get transformation matrices if not provided
    if Tpp is None or Tnp is None or eigenfreqs_GHz is None:
        print("Computing transformation matrices...")
        omega_H = gamma * H0
        H_BdG = generate_H_BdG_discrete_bar(omega_H, N_max, d_bar, w_bar, l_bar)
        H_BdG = add_demagnetization_corrections(H_BdG, N_max, d_bar, w_bar, l_bar)
        eigenfreqs, Tpp, Tnp, Tpn, Tnn, T = paraunitary_diag((H_BdG + H_BdG.conj().T) / 2)
        eigenfreqs_GHz = eigenfreqs * omega_M
    
    print(f"\n{'='*70}")
    print(f"Starting spatial T1 scan: {num_y}×{num_z} = {num_y*num_z} positions")
    print(f"Estimated time: {num_y*num_z*N_max*0.1/3600:.1f} hours")
    print(f"{'='*70}\n")
    
    overall_start = time.time()
    
    for iy in range(start_iy, num_y):
        row_start = time.time()
        y_pos = y_positions[iy]
        
        for iz, z_pos in enumerate(z_positions):
            
            # Calculate coupling at this position
            coord_x = (d_bar + h_NV) / d_bar
            coord_y = (y_pos + w_bar/2) / w_bar
            coord_z = (z_pos + l_bar/2) / l_bar
            
            coupling_all_modes = np.zeros(N_max, dtype=complex)
            for mode_idx in range(N_max):
                coupling_all_modes[mode_idx] = calculate_coupling_at_position(
                    coord_x, coord_y, coord_z,
                    N_max, d_bar, w_bar, l_bar,
                    Tpp, Tnp, mode_idx
                )
            
            coupling_map[iy, iz, :] = coupling_all_modes
            
            # Calculate T1 at this position
            gamma_L, gamma_U, _, _, _, _ = calculate_T1_finite_bar_at_position(
                h_NV, H0, N_max, position=(y_pos, z_pos),
                coupling_all_modes=coupling_all_modes,
                eigenfreqs_GHz=eigenfreqs_GHz,
                Tpp=Tpp, Tnp=Tnp
            )
            
            gamma_L_map[iy, iz] = gamma_L
            gamma_U_map[iy, iz] = gamma_U
        
        # Row complete - show progress
        row_time = time.time() - row_start
        elapsed_total = time.time() - overall_start
        rows_done = iy - start_iy + 1
        rows_remaining = num_y - iy - 1
        eta = elapsed_total / rows_done * rows_remaining if rows_done > 0 else 0
        
        print(f"\n{'='*70}")
        print(f"Row {iy+1}/{num_y} complete (y={y_pos*1000:.1f} nm)")
        print(f"  Row time: {row_time/60:.1f} min")
        print(f"  Total elapsed: {elapsed_total/60:.1f} min ({elapsed_total/3600:.2f} hr)")
        print(f"  ETA: {eta/60:.1f} min ({eta/3600:.2f} hr)")
        print(f"  Avg Gamma_L this row: {np.mean(gamma_L_map[iy, :]):.2f} Hz")
        print(f"  T1_L range: [{1/np.max(gamma_L_map[iy, :]):.2e}, {1/np.min(gamma_L_map[iy, :]):.2e}] s")
        print(f"{'='*70}\n")
        
        # Save checkpoint
        if (iy + 1 - start_iy) % save_every == 0 or iy == num_y - 1:
            np.savez_compressed(checkpoint_file,
                               gamma_L_map=gamma_L_map,
                               gamma_U_map=gamma_U_map,
                               coupling_map=coupling_map,
                               last_iy=iy,
                               y_positions=y_positions,
                               z_positions=z_positions,
                               H0=H0, h_NV=h_NV, N_max=N_max)
            print(f"*** Checkpoint saved: {checkpoint_file} ***\n")
    
    # Clean up checkpoint when done
    if os.path.exists(checkpoint_file):
        # Save final result with different name
        final_file = f'T1_spatial_final_H{H0:.1f}Oe_N{N_max}.npz'
        np.savez_compressed(final_file,
                           gamma_L_map=gamma_L_map,
                           gamma_U_map=gamma_U_map,
                           coupling_map=coupling_map,
                           y_positions=y_positions,
                           z_positions=z_positions,
                           H0=H0, h_NV=h_NV, N_max=N_max)
        os.remove(checkpoint_file)
        print(f"✓ Final results saved to {final_file}")
    
    return gamma_L_map, gamma_U_map, coupling_map

#=============================================================================
# SECTION 6: VERIFICATION AGAINST MATHEMATICA
#=============================================================================

def verify_against_mathematica(H0, h_NV, N_max):
    """
    Verify results against Mathematica PDF benchmarks
    
    From PDF page 12-13:
    - Mode 5 frequency: 2.77613 GHz
    - Mode 5 dimensionless coupling: 0.163658
    - Mode 5 coupling strength: 517.119 kHz
    - ω_dwl: 1450.66 Hz
    """
    print("\n" + "="*70)
    print("VERIFICATION AGAINST MATHEMATICA PDF")
    print("="*70)
    
    # Build Hamiltonian
    omega_H = gamma * H0
    H_BdG = generate_H_BdG_discrete_bar(omega_H, N_max, d_bar, w_bar, l_bar)
    H_BdG = add_demagnetization_corrections(H_BdG, N_max, d_bar, w_bar, l_bar)
    
    # Diagonalize
    eigenfreqs, Tpp, Tnp, Tpn, Tnn, T = paraunitary_diag((H_BdG + H_BdG.conj().T) / 2)
    eigenfreqs_GHz = eigenfreqs * omega_M

#     # In verify_against_mathematica or main test:
# # CRITICAL: Use Mathematica's normalized coordinates!
# coord_x = (d_bar + h_NV) / d_bar  # = 2.0 (NV at 2× thickness height)
# coord_y = w_bar / w_bar            # = 1.0 (NV at the edge) 
# coord_z = 0.4 / l_bar              # = 0.133 (NV at 0.4 μm along bar)

# # These are now NORMALIZED coordinates matching Mathematica exactly
    
    # Calculate coupling at Mathematica position
    coord_x = (d_bar + h_NV) / d_bar  # = 2.0
    coord_y = 1.0  # w/w = 1.0
    coord_z = 0.4 / l_bar  # = 0.133333
    
    print(f"\nPosition: CoordX={coord_x:.1f}, CoordY={coord_y:.1f}, CoordZ={coord_z:.6f}")
    print(f"Physical: x={coord_x*d_bar*1000:.1f} nm, y={coord_y*w_bar*1000:.1f} nm, z={coord_z*l_bar*1000:.1f} nm")
    
    coupling_mode_5 = calculate_coupling_at_position(
        coord_x, coord_y, coord_z,
        N_max, d_bar, w_bar, l_bar,
        Tpp, Tnp, 5
    )
    
    # Convert to kHz
    coupling_kHz = np.abs(coupling_mode_5) * np.sqrt(omega_dwl * 1e-3 * omega_M * 1e9)
    
    # Expected values from Mathematica
    expected_freq = 2.77613
    expected_coupling_dimensionless = 0.163658
    expected_coupling_kHz = 517.119
    expected_omega_dwl = 1450.66
    
    print(f"\n{'='*70}")
    print(f"{'Parameter':<30} {'Python':<20} {'Mathematica':<20} {'Error %':<10}")
    print(f"{'-'*70}")
    
    # ω_dwl
    omega_dwl_error = abs(omega_dwl - expected_omega_dwl) / expected_omega_dwl * 100
    print(f"{'ω_dwl (Hz)':<30} {omega_dwl:<20.2f} {expected_omega_dwl:<20.2f} {omega_dwl_error:<10.3f}")
    
    # Mode 5 frequency
    freq_error = abs(eigenfreqs_GHz[5] - expected_freq) / expected_freq * 100
    print(f"{'Mode 5 frequency (GHz)':<30} {eigenfreqs_GHz[5]:<20.5f} {expected_freq:<20.5f} {freq_error:<10.3f}")
    
    # Dimensionless coupling
    coup_dim_error = abs(np.abs(coupling_mode_5) - expected_coupling_dimensionless) / expected_coupling_dimensionless * 100
    print(f"{'Dimensionless coupling':<30} {np.abs(coupling_mode_5):<20.6f} {expected_coupling_dimensionless:<20.6f} {coup_dim_error:<10.3f}")
    
    # Coupling in kHz
    coup_error = abs(coupling_kHz - expected_coupling_kHz) / expected_coupling_kHz * 100
    print(f"{'Coupling strength (kHz)':<30} {coupling_kHz:<20.3f} {expected_coupling_kHz:<20.3f} {coup_error:<10.3f}")
    
    print(f"{'='*70}\n")
    
    # Pass/fail
    tolerance_freq = 0.01  # 0.01%
    tolerance_coupling = 1.0  # 1%
    tolerance_omega = 0.1  # 0.1%
    
    all_pass = (omega_dwl_error < tolerance_omega and 
                freq_error < tolerance_freq and 
                coup_error < tolerance_coupling)
    
    if all_pass:
        print("✓✓✓ VERIFICATION PASSED ✓✓✓")
        print("Python results match Mathematica within tolerances!")
    else:
        print("✗✗✗ VERIFICATION FAILED ✗✗✗")
        print("Results differ from Mathematica - check implementation!")
    
    print("="*70 + "\n")
    
    return all_pass, eigenfreqs_GHz, Tpp, Tnp

#=============================================================================
# SECTION 7: PLOTTING UTILITIES
#=============================================================================

def plot_mode_profiles(eigenfreqs_GHz, Tpp, Tnp, N_max, modes_to_plot=[0, 1, 2, 5]):
    """Plot spatial profiles of selected modes"""
    fig, axes = plt.subplots(len(modes_to_plot), 1, figsize=(10, 3*len(modes_to_plot)))
    
    if len(modes_to_plot) == 1:
        axes = [axes]
    
    z_array = np.linspace(0, l_bar, 300)
    
    for idx, mode_n in enumerate(modes_to_plot):
        # Build mode profile: sum over basis functions
        profile = np.zeros_like(z_array)
        for n1 in range(N_max):
            normalization = np.sqrt(2 / (1 + (n1 == 0)))
            basis = np.cos((n1) * np.pi * z_array / l_bar)
            weight = Tpp[n1, mode_n] + Tnp[n1, mode_n]
            profile += normalization * basis * weight
        
        axes[idx].plot(z_array * 1000, np.real(profile), 'b-', linewidth=2)
        axes[idx].set_xlabel('z (nm)', fontsize=12)
        axes[idx].set_ylabel('$m_x$ (arb.)', fontsize=12)
        axes[idx].set_title(f'Mode {mode_n}: f = {eigenfreqs_GHz[mode_n]:.4f} GHz', fontsize=14)
        axes[idx].grid(True, alpha=0.3)
        axes[idx].axhline(0, color='k', linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    return fig


def plot_spatial_T1_map(gamma_L_map, y_positions, z_positions, H0, h_NV):
    """Plot 2D spatial map of T1"""
    Y, Z = np.meshgrid(y_positions * 1000, z_positions * 1000, indexing='ij')
    T1_L_map = 1 / gamma_L_map
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Gamma_L map
    im1 = ax1.contourf(Z, Y, gamma_L_map, levels=20, cmap='viridis')
    ax1.set_xlabel('Z Position (nm)', fontsize=12)
    ax1.set_ylabel('Y Position (nm)', fontsize=12)
    ax1.set_title(f'$\Gamma_L$ (Hz) at H={H0:.1f} Oe, h={h_NV*1000:.0f} nm', fontsize=14)
    plt.colorbar(im1, ax=ax1, label='$\Gamma_L$ (Hz)')
    
    # T1_L map (log scale)
    im2 = ax2.contourf(Z, Y, np.log10(T1_L_map), levels=20, cmap='plasma')
    ax2.set_xlabel('Z Position (nm)', fontsize=12)
    ax2.set_ylabel('Y Position (nm)', fontsize=12)
    ax2.set_title(f'$T_1$ (s, log scale) at H={H0:.1f} Oe, h={h_NV*1000:.0f} nm', fontsize=14)
    cbar = plt.colorbar(im2, ax=ax2, label='$\log_{10}(T_1)$ (s)')
    
    plt.tight_layout()
    return fig


def plot_coupling_vs_frequency(eigenfreqs_GHz, coupling_all_modes, H0, h_NV, position):
    """Plot coupling strength vs mode frequency"""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    coupling_kHz = np.abs(coupling_all_modes) * np.sqrt(omega_dwl * 1e-3 * omega_M * 1e9)
    
    ax.scatter(eigenfreqs_GHz, coupling_kHz, s=50, alpha=0.6)
    ax.set_xlabel('Mode Frequency (GHz)', fontsize=12)
    ax.set_ylabel('Coupling Strength (kHz)', fontsize=12)
    ax.set_title(f'NV-Magnon Coupling vs Frequency\nH={H0:.1f} Oe, h={h_NV*1000:.0f} nm, pos={position}', 
                 fontsize=14)
    ax.grid(True, alpha=0.3)
    
    # Mark special modes
    target_mode = 5
    ax.scatter(eigenfreqs_GHz[target_mode], coupling_kHz[target_mode], 
               s=200, c='red', marker='*', label=f'Mode {target_mode}', zorder=10)
    ax.legend(fontsize=12)
    
    plt.tight_layout()
    return fig



def diagnose_coupling_calculation(h_NV, H0, N_max, d_bar, w_bar, l_bar):
    """Diagnostic version with detailed printouts"""
    
    print("\n" + "="*70)
    print("DIAGNOSTIC: Full coupling calculation")
    print("="*70)
    
    # Step 1: Build Hamiltonian
    print("\nStep 1: Building Hamiltonian...")
    omega_H = gamma * H0
    H_BdG = generate_H_BdG_discrete_bar(omega_H, N_max, d_bar, w_bar, l_bar, use_cache=True)
    H_BdG = add_demagnetization_corrections(H_BdG, N_max, d_bar, w_bar, l_bar, use_cache=True)
    
    print(f"  H_BdG shape: {H_BdG.shape}")
    print(f"  H_BdG[0,0] = {H_BdG[0,0]:.6e}")
    print(f"  H_BdG[0,N_max] = {H_BdG[0,N_max]:.6e}")
    print(f"  Max |H_BdG| = {np.max(np.abs(H_BdG)):.6e}")
    
    # Step 2: Diagonalize
    print("\nStep 2: Diagonalizing...")
    eigenfreqs, Tpp, Tnp, Tpn, Tnn, T = paraunitary_diag((H_BdG + H_BdG.conj().T) / 2)
    eigenfreqs_GHz = eigenfreqs * omega_M
    
    print(f"  Eigenfreqs (first 10 modes, GHz):")
    for i in range(min(10, len(eigenfreqs_GHz))):
        print(f"    Mode {i}: {eigenfreqs_GHz[i]:.6f} GHz")
    
    print(f"\n  Tpp shape: {Tpp.shape}")
    print(f"  Tpp (first 5x5 block):")
    print(Tpp[:5, :5])
    
    print(f"\n  Tnp shape: {Tnp.shape}")
    print(f"  Tnp (first 5x5 block):")
    print(Tnp[:5, :5])
    
    print(f"\n  |Tpp| range: [{np.min(np.abs(Tpp)):.6e}, {np.max(np.abs(Tpp)):.6e}]")
    print(f"  |Tnp| range: [{np.min(np.abs(Tnp)):.6e}, {np.max(np.abs(Tnp)):.6e}]")
    
    # Step 3: Calculate coupling for mode 5
    print("\nStep 3: Calculating coupling for mode 5...")
    
    coord_x = (d_bar + h_NV) / d_bar  # = 2.0
    coord_y = w_bar / w_bar            # = 1.0
    coord_z = 0.4 / l_bar              # = 0.133333
    
    print(f"  NV position (normalized): ({coord_x:.2f}, {coord_y:.2f}, {coord_z:.6f})")
    
    d_norm = d_bar / l_bar
    w_norm = w_bar / l_bar
    
    # Calculate Gamma arrays for first few modes
    HalfGamma_mm_array = np.zeros(N_max, dtype=complex)
    HalfGamma_pm_array = np.zeros(N_max, dtype=complex)
    
    for p in range(min(10, N_max)):  # Just first 10 for diagnosis
        norm = np.sqrt(2 / (1 + (p == 0)))
        
        # Gamma_XX
        def integrand_XX(y, z):
            r_top_sq = (d_norm * (coord_x - 1.0))**2 + (w_norm * (coord_y - y))**2 + (coord_z - z)**2
            r_bot_sq = (d_norm * (coord_x - 0.0))**2 + (w_norm * (coord_y - y))**2 + (coord_z - z)**2
            field = (coord_x - 1.0) / r_top_sq**(3/2) - (coord_x - 0.0) / r_bot_sq**(3/2)
            mode = np.cos(p * np.pi * z) * norm
            return field * mode
        
        result_XX, _ = dblquad(integrand_XX, 0, 1, 0, 1, epsabs=1e-8, epsrel=1e-6)
        Gamma_XX = (d_bar * w_bar / (4 * np.pi * l_bar**2)) * result_XX
        
        # Gamma_YY (simplified, assuming symmetric)
        def integrand_YY(x, z):
            r_right_sq = (d_norm * (coord_x - x))**2 + (w_norm * (coord_y - 1.0))**2 + (coord_z - z)**2
            r_left_sq = (d_norm * (coord_x - x))**2 + (w_norm * (coord_y - 0.0))**2 + (coord_z - z)**2
            field = (coord_y - 1.0) / r_right_sq**(3/2) - (coord_y - 0.0) / r_left_sq**(3/2)
            mode = np.cos(p * np.pi * z) * norm
            return field * mode
        
        result_YY, _ = dblquad(integrand_YY, 0, 1, 0, 1, epsabs=1e-8, epsrel=1e-6)
        Gamma_YY = (d_bar * w_bar / (4 * np.pi * l_bar**2)) * result_YY
        
        # For simplicity, assume Gamma_XY and Gamma_YX are smaller
        Gamma_XY = 0.0
        Gamma_YX = 0.0
        
        HalfGamma_mm_array[p] = (Gamma_XX - Gamma_YY - 1j * (Gamma_XY + Gamma_YX)) / 2
        HalfGamma_pm_array[p] = (Gamma_XX + Gamma_YY - 1j * (Gamma_XY - Gamma_YX)) / 2
        
        print(f"  p={p}: ΓXX={Gamma_XX:.6e}, ΓYY={Gamma_YY:.6e}")
        print(f"       HalfΓ_mm={HalfGamma_mm_array[p]:.6e}, HalfΓ_pm={HalfGamma_pm_array[p]:.6e}")
    
    # Step 4: Transform to mode 5
    print("\nStep 4: Transform to mode 5...")
    mode_idx = 5
    
    print(f"  Tpp[:, 5] (first 10 elements):")
    print(Tpp[:10, mode_idx])
    
    print(f"  Tnp[:, 5] (first 10 elements):")
    print(Tnp[:10, mode_idx])
    
    # Just use first 10 modes we calculated
    coupling_partial = (HalfGamma_mm_array[:10] @ Tpp[:10, mode_idx] + 
                       HalfGamma_pm_array[:10] @ Tnp[:10, mode_idx])
    
    print(f"\n  Coupling (partial, first 10 modes): {coupling_partial:.6e}")
    print(f"  |Coupling|: {np.abs(coupling_partial):.6e}")
    
    # Step 5: Check unit conversion
    print("\nStep 5: Unit conversion check...")
    
    h_plank = 6.626e-34
    mu_0 = 4 * np.pi * 1e-7
    omega_dwl = (h_plank * mu_0 * (gamma * 1e9)**2 / 
                 ((d_bar * 1e-6) * (w_bar * 1e-6) * (l_bar * 1e-6))) * 1e8
    
    print(f"  ω_dwl = {omega_dwl:.2f} Hz (expected: 1450.66 Hz)")
    print(f"  ω_M = {omega_M:.6f} GHz (expected: 6.88 GHz)")
    
    # Mathematica formula: g_kHz = |coupling| * sqrt(ω_dwl * 1e-3 * ω_M * 1e6)
    g_Hz_v1 = np.abs(coupling_partial) * np.sqrt(omega_dwl * 1e-3 * omega_M * 1e6)
    g_Hz_v2 = np.abs(coupling_partial) * np.sqrt(omega_dwl * omega_M * 1e6)
    g_Hz_v3 = np.abs(coupling_partial) * np.sqrt(omega_dwl * 1e-3 * omega_M * 1e9)
    
    print(f"\n  g (Hz) with formula v1 (ω_dwl * 1e-3 * ω_M * 1e6): {g_Hz_v1:.2e}")
    print(f"  g (Hz) with formula v2 (ω_dwl * ω_M * 1e6): {g_Hz_v2:.2e}")
    print(f"  g (Hz) with formula v3 (ω_dwl * 1e-3 * ω_M * 1e9): {g_Hz_v3:.2e}")
    print(f"\n  Expected from Mathematica: ~517 kHz = 5.17e5 Hz")
    
    return eigenfreqs_GHz, Tpp, Tnp




#=============================================================================
# MAIN EXECUTION
#=============================================================================

if __name__ == "__main__":

    # FORCE DELETE OLD CACHE FILES
    import glob
    for cache_file in glob.glob('H_BdG_cache_*.npz'):
        try:
            os.remove(cache_file)
            print(f"Deleted old cache: {cache_file}")
        except:
            pass
    
    for cache_file in glob.glob('Demag_cache_*.npz'):
        try:
            os.remove(cache_file)
            print(f"Deleted old cache: {cache_file}")
        except:
            pass
    
    print("\n" + "#"*70)
    print("# DIAGNOSTIC TEST")
    print("#"*70 + "\n")

    # Run the diagnostic FIRST
    print("\n" + "#"*70)
    print("# DIAGNOSTIC TEST")
    print("#"*70 + "\n")
    
    H0_test = 36.027
    h_NV_test = 0.005
    N_max_test = 50
    
    eigenfreqs_GHz, Tpp, Tnp = diagnose_coupling_calculation(
        h_NV_test, H0_test, N_max_test, d_bar, w_bar, l_bar
    )

    print('eigenfreqs_GHz ', eigenfreqs_GHz)
    
    print('Tpp ', Tpp)

    print('Tnp ', Tnp)




    # =================================================================
    # STEP 1: VERIFICATION TEST
    # =================================================================
    print("\n" + "#"*70)
    print("# STEP 1: VERIFICATION AGAINST MATHEMATICA")
    print("#"*70 + "\n")
    
    H0_verify = 36.027  # Oe (from PDF)
    h_NV_verify = 0.005  # μm (5 nm)
    N_max_verify = 50
    
    verification_passed, eigenfreqs_GHz, Tpp, Tnp = verify_against_mathematica(
        H0_verify, h_NV_verify, N_max_verify
    )
    
    if not verification_passed:
        print("⚠ WARNING: Verification failed! Check implementation before proceeding.")
        proceed = input("Continue anyway? (y/n): ")
        if proceed.lower() != 'y':
            exit()
    
    # =================================================================
    # STEP 2: PLOT MODE PROFILES
    # =================================================================
    print("\n" + "#"*70)
    print("# STEP 2: PLOTTING MODE PROFILES")
    print("#"*70 + "\n")
    
    fig_modes = plot_mode_profiles(eigenfreqs_GHz, Tpp, Tnp, N_max_verify, modes_to_plot=[0, 1, 2, 5])
    plt.savefig('mode_profiles.png', dpi=300, bbox_inches='tight')
    print("✓ Mode profiles saved to mode_profiles.png\n")
    
    # =================================================================
    # STEP 3: SINGLE POSITION T1 CALCULATION
    # =================================================================
    print("\n" + "#"*70)
    print("# STEP 3: T1 AT CENTER POSITION")
    print("#"*70 + "\n")
    
    gamma_L, gamma_U, eigenfreqs, coupling, Tpp, Tnp = calculate_T1_finite_bar_at_position(
        h_NV_verify, H0_verify, N_max_verify, position='center',
        coupling_all_modes=None, eigenfreqs_GHz=eigenfreqs_GHz,
        Tpp=Tpp, Tnp=Tnp
    )
    
    # Plot coupling vs frequency
    fig_coupling = plot_coupling_vs_frequency(eigenfreqs, coupling, H0_verify, h_NV_verify, 'center')
    plt.savefig('coupling_vs_frequency.png', dpi=300, bbox_inches='tight')
    print("✓ Coupling plot saved to coupling_vs_frequency.png\n")
    
    # =================================================================
    # STEP 4: SPATIAL T1 SCAN (OPTIONAL - VERY SLOW!)
    # =================================================================
    print("\n" + "#"*70)
    print("# STEP 4: SPATIAL T1 SCAN (OPTIONAL)")
    print("#"*70 + "\n")
    
    do_spatial_scan = input("Perform spatial T1 scan? This will take HOURS! (y/n): ")
    
    if do_spatial_scan.lower() == 'y':
        # Small grid for testing (increase for publication quality)
        num_y_points = int(input("Number of Y points (recommend 11 for test, 121 for final): "))
        num_z_points = int(input("Number of Z points (recommend 11 for test, 121 for final): "))
        
        y_positions = np.linspace(-w_bar, w_bar, num_y_points)
        z_positions = np.linspace(-l_bar * 3/5, l_bar * 3/5, num_z_points)
        
        estimated_time = num_y_points * num_z_points * N_max_verify * 0.1 / 3600
        print(f"\nEstimated time: {estimated_time:.1f} hours")
        confirm = input("Continue? (y/n): ")
        
        if confirm.lower() == 'y':
            gamma_L_map, gamma_U_map, coupling_map = calculate_T1_spatial_scan(
                h_NV_verify, H0_verify, N_max_verify,
                y_positions, z_positions,
                Tpp=Tpp, Tnp=Tnp, eigenfreqs_GHz=eigenfreqs_GHz,
                save_every=1
            )
            
            # Plot results
            fig_spatial = plot_spatial_T1_map(gamma_L_map, y_positions, z_positions, 
                                             H0_verify, h_NV_verify)
            plt.savefig(f'T1_spatial_map_H{H0_verify:.1f}Oe.png', dpi=300, bbox_inches='tight')
            print(f"✓ Spatial map saved to T1_spatial_map_H{H0_verify:.1f}Oe.png\n")
    
    # =================================================================
    # DONE
    # =================================================================
    print("\n" + "#"*70)
    print("# ALL CALCULATIONS COMPLETE")
    print("#"*70 + "\n")
    
    print("Summary of results:")
    print(f"  Mode 5 frequency: {eigenfreqs_GHz[5]:.5f} GHz")
    print(f"  Mode 5 coupling: {np.abs(coupling[5]) * np.sqrt(omega_dwl * 1e-3 * omega_M * 1e9):.3f} kHz")
    print(f"  Gamma_L (center): {gamma_L:.2f} Hz")
    print(f"  T1_L (center): {1/gamma_L:.2e} s")
    print(f"\nAll figures saved to current directory.")
    
    plt.show()