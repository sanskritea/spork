import numpy as np
from scipy.linalg import cholesky, eig
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import warnings
warnings.filterwarnings('ignore')

# Constants and initial calculations
delta_h = 1
# Solve for MsSol: ((2870/2.8) - XX/2)/2 == 82 + delta_h
# Rearranging: XX = 2*((2870/2.8) - 2*(82 + delta_h))
xx_solution = 2*((2870/2.8) - 2*(82 + delta_h))
ms_sol = xx_solution - 4 + 2
print('ms_sol ', ms_sol)

# Material & condition parameters
M0 = ms_sol  # Oe
L = 3  # μm
DD = 5.4e-9 * 1e8  # Oe μm^-2
gamma = 2.8e-3  # GHz/G
omega_M = gamma * M0  # GHz


# from wolframclient.evaluation import WolframLanguageSession
# from wolframclient.language import wl

# # Mathematica Cholesky decomposition
def mathematica_cholesky_decomposition(A):
    """
    Method 8: Direct interface to Mathematica (if available)
    """
    
    try:
        from wolframclient.evaluation import WolframLanguageSession
        from wolframclient.language import wl
        
        # Start Mathematica session
        session = WolframLanguageSession()
        
        # Convert matrix to Mathematica format
        A_list = A.tolist()
        
        # Call Mathematica's CholeskyDecomposition
        result = session.evaluate(wl.CholeskyDecomposition(A_list))
        
        # Convert result back to numpy
        L = np.array(result, dtype=float)
        
        session.terminate()
        
        return L
        
    except ImportError:
        return None, "wolframclient not available", None
    except Exception as e:
        return None, f"mathematica_error: {str(e)}", None


def scipy_cholesky_decomposition(A):
    """
    Replace Mathematica Cholesky decomposition with SciPy version
    """
    try:
        # Ensure matrix is positive definite
        # Add small regularization if needed
        eigenvals = np.linalg.eigvals(A)
        if np.min(eigenvals) <= 0:
            reg_param = 1e-10
            A_reg = A + reg_param * np.eye(A.shape[0])
        else:
            A_reg = A
            
        # Compute Cholesky decomposition
        L = cholesky(A_reg, lower=False)
        return L
        
    except np.linalg.LinAlgError as e:
        print(f"Cholesky decomposition failed: {e}")
        # Fallback to eigenvalue decomposition
        eigenvals, eigenvecs = np.linalg.eigh(A)
        # Ensure positive eigenvalues
        eigenvals = np.maximum(eigenvals, 1e-10)
        # Reconstruct matrix and compute Cholesky
        A_reconstructed = eigenvecs @ np.diag(eigenvals) @ eigenvecs.T
        return cholesky(A_reconstructed, lower=False)


# Mathematica sorting of eigenvalues
def simple_mathematica_order(values, positive_first=True):

    values = np.array(values)
    def sort_key(x):
        # Primary: decreasing absolute value (negative for descending sort)
        abs_val = -abs(x)
        
        # Secondary: sign preference
        if positive_first:
            # Positive first: use negative of sign (so positive gets -1, negative gets 1)
            sign_preference = -np.sign(x) if x != 0 else 0
        else:
            # Negative first: use sign directly
            sign_preference = np.sign(x) if x != 0 else 0
            
        return (abs_val, sign_preference)
    
    # Sort using the custom key
    sorted_indices = sorted(range(len(values)), key=lambda i: sort_key(values[i]))
    
    return values[sorted_indices]


# Rounding off eigenvalues to 6 significant figures to match Mathematica
def round_array_sig_figs_vectorized(arr, sig_figs=6):
    """
    Vectorized version using numpy operations (faster for large arrays).
    """
    arr = np.array(arr, dtype=float)
    
    # Handle zeros separately
    zero_mask = (arr == 0)
    
    # For non-zero elements
    non_zero = arr[~zero_mask]
    if len(non_zero) > 0:
        # Calculate magnitude for each element
        magnitudes = np.floor(np.log10(np.abs(non_zero))).astype(int)
        
        # Calculate decimal places needed
        decimal_places = sig_figs - magnitudes - 1
        
        # Round each element
        rounded_non_zero = np.array([
            round(val, dec_places) 
            for val, dec_places in zip(non_zero, decimal_places)
        ])
        
        # Put back into original array
        result = arr.copy()
        result[~zero_mask] = rounded_non_zero
        result[zero_mask] = 0.0
    else:
        result = arr.copy()
    
    return result


def mathematica_eigensystem(matrix, tolerance=1e-10):
    """
    Replicate Mathematica's Eigensystem behavior exactly
    """
    matrix = np.array(matrix, dtype=float)
    
    # Get eigenvalues and eigenvectors from scipy
    eigenvals, eigenvecs = eig(matrix)
    
    # Convert to real if they're essentially real
    if np.allclose(eigenvals.imag, 0, atol=tolerance):
        eigenvals = eigenvals.real
    if np.allclose(eigenvecs.imag, 0, atol=tolerance):
        eigenvecs = eigenvecs.real
    
    # Normalize eigenvectors (scipy usually returns normalized, but let's be sure)
    for i in range(eigenvecs.shape[1]):
        norm = np.linalg.norm(eigenvecs[:, i])
        if norm > tolerance:
            eigenvecs[:, i] = eigenvecs[:, i] / norm
    
    # Apply Mathematica's sign convention
    # For each eigenvector, find the first non-zero element and make it positive
    for i in range(eigenvecs.shape[1]):
        vec = eigenvecs[:, i]
        # Find first non-zero element
        first_nonzero_idx = None
        for j in range(len(vec)):
            if abs(vec[j]) > tolerance:
                first_nonzero_idx = j
                break
        
        if first_nonzero_idx is not None:
            # If the first non-zero element is negative, flip the entire vector
            if vec[first_nonzero_idx] < 0:
                eigenvecs[:, i] = -eigenvecs[:, i]
    
    # Sort eigenvalues and eigenvectors by eigenvalue magnitude (descending)
    # Then by sign (positive first for same magnitude)
    def sort_key(i):
        val = eigenvals[i]
        return (-abs(val), -val)  # Negative signs for descending order
    
    indices = sorted(range(len(eigenvals)), key=sort_key)
    
    sorted_eigenvals = eigenvals[indices]
    sorted_eigenvecs = eigenvecs[:, indices]
    
    return sorted_eigenvals, sorted_eigenvecs





def paraunitary_diag(H):
    """Paraunitary diagonalization function"""

    K = scipy_cholesky_decomposition(H)
    # print('K ', K)
    dim = H.shape[0]

    # Create sigma3 matrix
    sigma3 = np.diag([(-1)**np.floor((2*n-1)/dim) for n in range(1, dim+1)])
    # print('sigma3 ', sigma3)
    W = K @ sigma3 @ K.conj().T
    # print('W ', W)
    
    # Eigendecomposition
    eval_w, evec_w = mathematica_eigensystem(W)
    evec_w = evec_w.T
    evec_w = evec_w / np.linalg.norm(evec_w, axis=0)

    # Permutation for ordering
    preperm = np.concatenate([
        np.arange(dim//2, dim),
        np.arange(dim//2-1, -1, -1)
    ])
    # print('preperm ', preperm)
    ordering = np.argsort(eval_w)
    # print('ordering ', ordering)
    permutation = ordering[preperm]
    # print('normal permutation ', permutation)
    
    eval_w = eval_w[permutation]
    evec_w = evec_w[permutation]
    # print('Post perm eval_w ', eval_w)
    # print('Post perm evec_w ', evec_w)
    
    U = evec_w.T
    # print('U ', U)
    H_diag = sigma3 * np.diag(eval_w)
    # print('H_diag ', H_diag)
    
    T = np.linalg.inv(K) @ U @ np.sqrt(np.abs(H_diag))
    # print('T ', T)
    
    # Phase correction
    dim_half = dim // 2
    Tpp = T[:dim_half, :dim_half]
    Tnn = T[dim_half:, dim_half:]
    
    phase_array_p = np.exp(1j * np.angle(np.diag(Tpp)))
    phase_array_n = np.exp(1j * np.angle(np.diag(Tnn)))
    
    V = np.diag(np.concatenate([phase_array_p.conj(), phase_array_n.conj()]))
    T = T @ V
    
    # Extract blocks
    Tpp = T[:dim_half, :dim_half]
    Tnp = T[dim_half:, :dim_half]
    Tpn = T[:dim_half, dim_half:]
    Tnn = T[dim_half:, dim_half:]
    
    return eval_w[:dim_half], Tpp, Tnp, Tpn, Tnn, T


def F_func(q, n):
    """F function for dipolar calculations"""
    if q == 0:
        return 2 * n
    return 2 * (1 - (-1)**n * np.exp(-q)) / q


def P_func(q, n, m):
    """P function for dipolar calculations"""
    kronecker_nm = 1 if n == m else 0
    kronecker_n0 = 1 if n == 0 else 0
    kronecker_m0 = 1 if m == 0 else 0
    
    term1 = q**2 / (q**2 + n**2 * np.pi**2) * kronecker_nm
    
    normalization = 1 / np.sqrt((1 + kronecker_n0) * (1 + kronecker_m0))
    term2 = (q**4 / ((q**2 + n**2 * np.pi**2) * (q**2 + m**2 * np.pi**2)) * 
             F_func(q, n) * (1 + (-1)**(n + m)) / 2 * normalization)
    
    return term1 - term2


def Q_func(q, n, m):
    """Q function for dipolar calculations"""
    kronecker_n0 = 1 if n == 0 else 0
    kronecker_m0 = 1 if m == 0 else 0
    
    if n == m:
        return 0
        
    normalization = 1 / np.sqrt((1 + kronecker_n0) * (1 + kronecker_m0))
    
    term1 = (q**2 / (q**2 + m**2 * np.pi**2) * 
             (m**2 / (m**2 - n**2 + (1 + (-1)**(n + m))/2) * 2/q - 
              q**2 / (2 * (q**2 + n**2 * np.pi**2)) * F_func(q, n)))
    
    return term1 * normalization * (1 - (-1)**(n + m)) / 2


def capital_omega(omega_H, q, n):
    """Capital Omega function"""
    return (omega_H + (gamma * DD) / L**2 * (q**2 + n**2 * np.pi**2)) / omega_M


def H_matrix(omega_H, q, phi_k, n, m):
    """H matrix calculation"""

    omega_val = capital_omega(omega_H, q, n)
    kronecker_nm = 1 if n == m else 0
    
    H_nm = np.zeros((2, 2), dtype=complex)
    
    # Diagonal terms
    if n == m:
        H_nm += omega_val * np.eye(2)
        H_nm += 0.5 * np.ones((2, 2))
    
    # P terms
    P_val = P_func(q, n, m)
    P_matrix = np.array([[1 - np.sin(phi_k)**2, 1 + np.sin(phi_k)**2],
                         [1 + np.sin(phi_k)**2, 1 - np.sin(phi_k)**2]])
    H_nm -= 0.5 * P_matrix * P_val
    
    # Q terms
    Q_val = Q_func(q, n, m)
    Q_matrix = np.array([[0, -4], [4, 0]])
    H_nm -= 0.5 * Q_matrix * np.sin(phi_k) * Q_val
    
    return H_nm


def generate_H_BdG(omega_H, q, phi_k, N_max):
    """Generate BdG Hamiltonian"""
    H_BdG = np.zeros((2*N_max, 2*N_max), dtype=complex)
    
    for m in range(N_max):
        for n in range(N_max):
            result = H_matrix(omega_H, q, phi_k, n, m)
            H_BdG[n, m] = result[0, 0]
            H_BdG[n, m + N_max] = result[0, 1]
            H_BdG[n + N_max, m] = result[1, 0]
            H_BdG[n + N_max, m + N_max] = result[1, 1]
    
    return H_BdG


def f_func(q, n, h_over_L):
    """f function for coupling calculations"""
    kronecker_n0 = 1 if n == 0 else 0
    normalization = 1 / np.sqrt(2 * (1 + kronecker_n0))
    
    if q == 0:
        return 0
        
    return ((-1)**n * normalization * q**2 / (q**2 + n**2 * np.pi**2) * 
            np.exp(-q * h_over_L) * (1 - (-1)**n * np.exp(-q)))


from quadinterp import QuadraticInterpolator2D

def multi_para_diag(h_NV_array, omega_H, q_table, phi_table, N_max):
    """Multi paraunitary diagonalization"""
    # print('In multi_para_diag')
    
    # make everything a numpy array
    q_table = np.array(q_table)
    phi_table = np.array(phi_table)
    h_NV_array = np.array(h_NV_array)

    # Get dimensions
    num_q = len(q_table)
    num_phi = len(phi_table)
    num_h_NV = len(h_NV_array)
    
    # Calculate measure
    q_diff = np.mean(np.diff(q_table))
    phi_diff = np.mean(np.diff(phi_table))
    measure = q_diff * phi_diff / (2 * np.pi)**2
    
    # Initialize tables - size for phi is 2*num_phi since we get phi+pi results
    Q_table = np.tile(np.array(q_table[:, np.newaxis, np.newaxis]), (1, 2*num_phi, N_max))
    
    # Create Phi table with phi and phi+pi values
    Phi_table = np.zeros((num_q, 2*num_phi, N_max))
    for n in range(num_q):
        for m in range(2*num_phi):
            phi_idx = (m) % num_phi
            pi_factor = (m) // num_phi
            Phi_table[n, m, :] = phi_table[phi_idx] + np.pi * pi_factor
    
    # Initialize result tables
    omega_BdG_table = np.zeros((num_q, 2*num_phi, N_max))
    coupling_plus_table = np.zeros((num_q, 2*num_phi, num_h_NV, N_max), dtype=complex)
    coupling_minus_table = np.zeros((num_q, 2*num_phi, num_h_NV, N_max), dtype=complex)
    coupling_z_table = np.zeros((num_q, 2*num_phi, num_h_NV, N_max), dtype=complex)
    
    # Main calculation loop with progress monitoring
    print("Starting MultiParaDiag calculation...")
    for count_q in range(num_q):
        progress = (count_q + 1) / num_q * 100
        print(f"Progress: {progress:.1f}% MultiParaDiag")
        
        for count_phi in range(num_phi):
            q = Q_table[count_q, count_phi, 0]
            phi_k = Phi_table[count_q, count_phi, 0]
            
            
            # Generate HBdG (this function needs to be implemented separately)
            # print('Generating HBdG')
            H_BdG = generate_H_BdG(omega_H, q, phi_k, N_max)
            
            # Para-diagonalization
            # print('Para-unitary diagonalization')
            result = paraunitary_diag((H_BdG + H_BdG.conj().T) / 2)
            
            # Extract results
            omega_BdG_table[count_q, count_phi] = result[0]
            omega_BdG_table[count_q, count_phi + num_phi] = omega_BdG_table[count_q, count_phi]
            
            T_pp = result[1]
            T_np = result[2]
            T_pn = result[3]
            T_nn = result[4]
            
            # Calculate gamma values
            gamma_plus = ((1 + np.sin(phi_k)) / 2 * 
                         (T_pp + T_np + np.sin(phi_k) * (T_pp - T_np)))
            gamma_plus_mirror = ((1 + np.sin(phi_k + np.pi)) / 2 * 
                               np.conj(T_nn + T_pn + np.sin(phi_k + np.pi) * (T_nn - T_pn)))
            
            gamma_minus = ((1 - np.sin(phi_k)) / 2 * 
                          (T_pp + T_np + np.sin(phi_k) * (T_pp - T_np)))
            gamma_minus_mirror = ((1 - np.sin(phi_k + np.pi)) / 2 * 
                                np.conj(T_nn + T_pn + np.sin(phi_k + np.pi) * (T_nn - T_pn)))
            
            gamma_z = (-1j * np.cos(phi_k) / 2 * 
                      (T_pp + T_np + np.sin(phi_k) * (T_pp - T_np)))
            gamma_z_mirror = (-1j * np.cos(phi_k + np.pi) / 2 * 
                            np.conj(T_nn + T_pn + np.sin(phi_k + np.pi) * (T_nn - T_pn)))
            
            # Calculate f_bar arrays
            f_bar_array = np.zeros((num_h_NV, N_max), dtype=complex)
            for tt in range(num_h_NV):
                for nn in range(N_max):
                    f_bar_array[tt, nn] = f_func(q, nn, h_NV_array[tt] / L)
            
            # Calculate v arrays
            v_plus_array = np.array([f_bar_array[tt] @ gamma_plus for tt in range(num_h_NV)])
            v_plus_mirror_array = np.array([f_bar_array[tt] @ gamma_plus_mirror for tt in range(num_h_NV)])
            v_minus_array = np.array([f_bar_array[tt] @ gamma_minus for tt in range(num_h_NV)])
            v_minus_mirror_array = np.array([f_bar_array[tt] @ gamma_minus_mirror for tt in range(num_h_NV)])
            v_z_array = np.array([f_bar_array[tt] @ gamma_z for tt in range(num_h_NV)])
            v_z_mirror_array = np.array([f_bar_array[tt] @ gamma_z_mirror for tt in range(num_h_NV)])
            
            # Store coupling values
            coupling_plus_table[count_q, count_phi] = v_plus_array
            coupling_plus_table[count_q, count_phi + num_phi] = v_plus_mirror_array
            coupling_minus_table[count_q, count_phi] = v_minus_array
            coupling_minus_table[count_q, count_phi + num_phi] = v_minus_mirror_array
            coupling_z_table[count_q, count_phi] = v_z_array
            coupling_z_table[count_q, count_phi + num_phi] = v_z_mirror_array
    
    # Prepare extended tables for interpolation (append overlapping point)
    Q_table_extend = np.zeros((num_q, 2*num_phi + 1, N_max))
    Phi_table_extend = np.zeros((num_q, 2*num_phi + 1, N_max))
    omega_BdG_table_extend = np.zeros((num_q, 2*num_phi + 1, N_max))
    coupling_plus_table_extend = np.zeros((num_q, 2*num_phi + 1, num_h_NV, N_max), dtype=complex)
    coupling_minus_table_extend = np.zeros((num_q, 2*num_phi + 1, num_h_NV, N_max), dtype=complex)
    coupling_z_table_extend = np.zeros((num_q, 2*num_phi + 1, num_h_NV, N_max), dtype=complex)
    
    # Copy existing data and add periodic boundary
    for count_q in range(num_q):
        # Copy existing data
        Q_table_extend[count_q, :2*num_phi] = Q_table[count_q]
        Phi_table_extend[count_q, :2*num_phi] = Phi_table[count_q]
        omega_BdG_table_extend[count_q, :2*num_phi] = omega_BdG_table[count_q]
        coupling_plus_table_extend[count_q, :2*num_phi] = coupling_plus_table[count_q]
        coupling_minus_table_extend[count_q, :2*num_phi] = coupling_minus_table[count_q]
        coupling_z_table_extend[count_q, :2*num_phi] = coupling_z_table[count_q]
        
        # Add periodic boundary (last point = first point)
        Q_table_extend[count_q, 2*num_phi] = Q_table[count_q, 0]
        Phi_table_extend[count_q, 2*num_phi] = Phi_table[count_q, 0]
        omega_BdG_table_extend[count_q, 2*num_phi] = omega_BdG_table[count_q, 0]
        coupling_plus_table_extend[count_q, 2*num_phi] = coupling_plus_table[count_q, 0]
        coupling_minus_table_extend[count_q, 2*num_phi] = coupling_minus_table[count_q, 0]
        coupling_z_table_extend[count_q, 2*num_phi] = coupling_z_table[count_q, 0]
    
    # # Create interpolation functions
    # print("Creating interpolation functions...")
    
    # # Prepare coordinate arrays for interpolation
    # def create_interpolation_function(data_table):
    #     """Create interpolation function for given data table"""
    #     interpolation_functions = []
        
    #     for s in range(N_max):
    #         # Extract coordinates and values for this s index
    #         q_coords = Q_table_extend[:, :, s].flatten()
    #         phi_coords = Phi_table_extend[:, :, s].flatten()
            
    #         if len(np.array(data_table).shape) == 3:  # omega table
    #             values = data_table[:, :, s].flatten()
    #             # Use griddata for irregular grid interpolation
    #             def interp_func(q_new, phi_new, q_c=q_coords, phi_c=phi_coords, v=values):
    #                 points = np.column_stack([q_c, phi_c])
    #                 return griddata(points, v, (q_new, phi_new), method='cubic')
    #         else:  # coupling tables (4D)
    #             interp_funcs_tt = []
    #             for tt in range(num_h_NV):
    #                 values = data_table[:, :, tt, s].flatten()
    #                 def interp_func(q_new, phi_new, q_c=q_coords, phi_c=phi_coords, v=values):
    #                     points = np.column_stack([q_c, phi_c])
    #                     return griddata(points, v, (q_new, phi_new), method='cubic')
    #                 interp_funcs_tt.append(interp_func)
    #             interpolation_functions.append(interp_funcs_tt)
    #             continue
                
    #         interpolation_functions.append(interp_func)
        
    #     return interpolation_functions

    # Create interpolation functions using QuadraticInterpolator2D
    print("Creating QuadraticInterpolator2D interpolation functions...")
    
    def create_quadratic_interpolation_function(data_table):
        """Create QuadraticInterpolator2D function for given data table"""
        interpolation_functions = []
        
        # Get coordinate arrays - need to handle both q and phi properly
        q_coords_raw = Q_table_extend[:, 0, 0]  # Extract q coordinates
        phi_coords_raw = Phi_table_extend[0, :, 0]  # Extract phi coordinates
        
        # Ensure q coordinates are unique and sorted
        q_unique, q_unique_indices = np.unique(q_coords_raw, return_index=True)
        q_coords_1d = q_unique
        
        # Sort phi coordinates and corresponding data to ensure ascending order
        phi_sort_indices = np.argsort(phi_coords_raw)
        phi_coords_1d = phi_coords_raw[phi_sort_indices]
        
        # Remove duplicate phi values if any
        phi_unique, phi_unique_indices = np.unique(phi_coords_1d, return_index=True)
        phi_coords_1d = phi_unique
        phi_final_indices = phi_sort_indices[phi_unique_indices]
        
        print(f"Debug: q_coords shape: {q_coords_1d.shape}, phi_coords shape: {phi_coords_1d.shape}")
        print(f"Debug: q_coords range: [{q_coords_1d[0]:.3f}, {q_coords_1d[-1]:.3f}]")
        print(f"Debug: phi_coords range: [{phi_coords_1d[0]:.3f}, {phi_coords_1d[-1]:.3f}]")
        
        for s in range(N_max):
            if len(np.array(data_table).shape) == 3:  # omega table
                # Select unique q and sorted unique phi data
                values_selected = data_table[np.ix_(q_unique_indices, phi_final_indices, [s])]
                values_2d = values_selected[:, :, 0].T  # Shape: (len(phi), len(q))
                
                print(f"Debug: values_2d shape for omega: {values_2d.shape}")
                
                # Create interpolator
                try:
                    interp = QuadraticInterpolator2D(q_coords_1d, phi_coords_1d, values_2d, method='spline')
                except Exception as e:
                    print(f"Error creating interpolator for omega at s={s}: {e}")
                    # Fallback to linear interpolation if quadratic fails
                    from scipy.interpolate import RectBivariateSpline
                    interp = RectBivariateSpline(phi_coords_1d, q_coords_1d, values_2d, kx=1, ky=1, s=0)
                
                # Create wrapper that handles periodicity
                def periodic_interp_func(q_new, phi_new, interpolator=interp, 
                                       phi_min=phi_coords_1d[0], phi_max=phi_coords_1d[-1],
                                       q_min=q_coords_1d[0], q_max=q_coords_1d[-1]):
                    # Map phi to the interpolation range using periodicity
                    phi_period = 2 * np.pi
                    phi_mapped = np.mod(phi_new - phi_min, phi_period) + phi_min
                    phi_mapped = np.clip(phi_mapped, phi_min, phi_max)
                    q_clipped = np.clip(q_new, q_min, q_max)
                    
                    if hasattr(interpolator, '__call__'):
                        return interpolator(q_clipped, phi_mapped)
                    else:  # RectBivariateSpline
                        return interpolator(phi_mapped, q_clipped)
                
                interpolation_functions.append(periodic_interp_func)
                
            else:  # coupling tables (4D)
                interp_funcs_tt = []
                for tt in range(num_h_NV):
                    # Select unique q and sorted unique phi data
                    values_selected = data_table[np.ix_(q_unique_indices, phi_final_indices, [tt], [s])]
                    values_2d = values_selected[:, :, 0, 0].T  # Shape: (len(phi), len(q))
                    
                    if np.iscomplexobj(values_2d):
                        # Create separate interpolators for real and imaginary parts
                        try:
                            real_interp = QuadraticInterpolator2D(q_coords_1d, phi_coords_1d, 
                                                                values_2d.real, method='spline')
                            imag_interp = QuadraticInterpolator2D(q_coords_1d, phi_coords_1d, 
                                                                values_2d.imag, method='spline')
                        except Exception as e:
                            print(f"Error creating complex interpolator at s={s}, tt={tt}: {e}")
                            # Fallback
                            from scipy.interpolate import RectBivariateSpline
                            real_interp = RectBivariateSpline(phi_coords_1d, q_coords_1d, values_2d.real, kx=1, ky=1, s=0)
                            imag_interp = RectBivariateSpline(phi_coords_1d, q_coords_1d, values_2d.imag, kx=1, ky=1, s=0)
                        
                        # Create combined interpolator function with periodicity
                        def complex_interp_func(q_new, phi_new, r_interp=real_interp, i_interp=imag_interp, 
                                              phi_min=phi_coords_1d[0], phi_max=phi_coords_1d[-1],
                                              q_min=q_coords_1d[0], q_max=q_coords_1d[-1]):
                            phi_period = 2 * np.pi
                            phi_mapped = np.mod(phi_new - phi_min, phi_period) + phi_min
                            phi_mapped = np.clip(phi_mapped, phi_min, phi_max)
                            q_clipped = np.clip(q_new, q_min, q_max)
                            
                            if hasattr(r_interp, '__call__'):
                                real_part = r_interp(q_clipped, phi_mapped)
                                imag_part = i_interp(q_clipped, phi_mapped)
                            else:  # RectBivariateSpline
                                real_part = r_interp(phi_mapped, q_clipped)
                                imag_part = i_interp(phi_mapped, q_clipped)
                            
                            return real_part + 1j * imag_part
                        
                        interp_funcs_tt.append(complex_interp_func)
                    else:
                        # Real data
                        try:
                            interp = QuadraticInterpolator2D(q_coords_1d, phi_coords_1d, 
                                                           values_2d, method='spline')
                        except Exception as e:
                            print(f"Error creating real interpolator at s={s}, tt={tt}: {e}")
                            # Fallback
                            from scipy.interpolate import RectBivariateSpline
                            interp = RectBivariateSpline(phi_coords_1d, q_coords_1d, values_2d, kx=1, ky=1, s=0)
                        
                        def real_interp_func(q_new, phi_new, interpolator=interp, 
                                           phi_min=phi_coords_1d[0], phi_max=phi_coords_1d[-1],
                                           q_min=q_coords_1d[0], q_max=q_coords_1d[-1]):
                            phi_period = 2 * np.pi
                            phi_mapped = np.mod(phi_new - phi_min, phi_period) + phi_min
                            phi_mapped = np.clip(phi_mapped, phi_min, phi_max)
                            q_clipped = np.clip(q_new, q_min, q_max)
                            
                            if hasattr(interpolator, '__call__'):
                                return interpolator(q_clipped, phi_mapped)
                            else:  # RectBivariateSpline
                                return interpolator(phi_mapped, q_clipped)
                        
                        interp_funcs_tt.append(real_interp_func)
                
                interpolation_functions.append(interp_funcs_tt)
        
        return interpolation_functions

    # Create interpolation tables
    Int_omega_BdG_table = create_quadratic_interpolation_function(omega_BdG_table_extend)
    Int_coupling_plus_table = create_quadratic_interpolation_function(coupling_plus_table_extend)
    Int_coupling_minus_table = create_quadratic_interpolation_function(coupling_minus_table_extend)
    Int_coupling_z_table = create_quadratic_interpolation_function(coupling_z_table_extend)
    
    # # Create interpolation tables
    # Int_omega_BdG_table = create_interpolation_function(omega_BdG_table_extend)
    # Int_coupling_plus_table = create_interpolation_function(coupling_plus_table_extend)
    # Int_coupling_minus_table = create_interpolation_function(coupling_minus_table_extend)
    # Int_coupling_z_table = create_interpolation_function(coupling_z_table_extend)
    
    return (Int_omega_BdG_table, Int_coupling_plus_table, 
            Int_coupling_minus_table, Int_coupling_z_table)


# Temperature and Bose-Einstein distribution
h_plank = 6.626e-34  # J*s
mu_0 = 4*np.pi*1e-7  # H/m
k_B = 1.381e-23  # J/K
temperature = 300  # K

def n_bose(omega):
    """Bose-Einstein distribution (high temperature approximation)"""
    return (1e-9 * k_B * temperature / h_plank) / omega


def f_space(min_val, max_val, steps, f_func=np.log):
    f_min = f_func(min_val)
    f_max = f_func(max_val)
    f_range = np.linspace(f_min, f_max, steps)
    return np.exp(f_range)  # InverseFunction of Log


def gamma_values_UL(eta_small, NQest, N_phi, H0, Int_omega_BdG_table, 
                   Int_coupling_plus_table, Int_coupling_minus_table, num_h_NV):
    print('In gamma_values_UL')
    """Corrected Gamma values calculation"""
    DNV = 2.87  # GHz
    omega_target_L = DNV - H0 * gamma  # GHz
    omega_target_U = DNV + H0 * gamma  # GHz
    
    # Convert to Hz for proper units
    omega_target_L_Hz = omega_target_L # GHz
    omega_target_U_Hz = omega_target_U # GHz
    omega_M_Hz = omega_M # GHz
    
    # Create q and phi grids
    q_max = 50 * L
    q_middle = 5 * L
    NQ_half = round(NQest/ 2)
    def f_space(min_val, max_val, steps):
        """Match Mathematica's fSpace function exactly"""
        f_min = np.log(min_val)
        f_max = np.log(max_val)
        f_range = np.linspace(f_min, f_max, steps)
        return np.exp(f_range)
    Q_table = np.concatenate([[1e-6], f_space(L/1000, q_middle, NQ_half)])

    # Updated Q_table to match MMA
    # Match Mathematica's extension method
    last_delta_q = np.diff(Q_table)[-1]
    additional_points = []
    current_q = max(Q_table) + last_delta_q
    while current_q <= q_max:
        additional_points.append(current_q)
        current_q += last_delta_q
    
    if additional_points:
        Q_table = np.concatenate([Q_table, additional_points])
    
    delta_q = np.diff(Q_table)
    NQ = len(delta_q)
    delta_phi = 2 * np.pi / N_phi
    
    gamma_L_Hz_array = np.zeros(num_h_NV)
    gamma_U_Hz_array = np.zeros(num_h_NV)
    DOS_L = 0
    DOS_U = 0


    print('np.shape(Int_omega_BdG_table) ', np.shape(Int_omega_BdG_table))

    N_max_local = np.shape(Int_omega_BdG_table)[0]
    
    for s in range(N_max_local):
        progress = (s + 1) / N_max_local * 100
        print(f"performing s={s+1} calculation: {progress:.1f}% Gamma1Calc")
        
        # Create flat tables
        # flat_Q_table = np.repeat(Q_table[:-1], N_phi)  # Exclude last element to match delta_Q length
        flat_Q_table = []
        for ii1 in range(NQ):
            for ii2 in range(N_phi):
                flat_Q_table.append(Q_table[ii1])
        flat_delta_Q_table = np.repeat(delta_q, N_phi)
        flat_phi_table = np.tile(np.linspace(0, 2*np.pi*(1-1/N_phi), N_phi), NQ)

        
        length_flat_table = len(flat_Q_table)
        flat_omega_BdG_table = np.zeros(length_flat_table)
        flat_coupling_plus_table_array = np.zeros((num_h_NV, length_flat_table), dtype=complex)
        flat_coupling_minus_table_array = np.zeros((num_h_NV, length_flat_table), dtype=complex)
        
        # Fill the flat tables using interpolation functions
        for ii in range(length_flat_table):
            flat_omega_BdG_table[ii] = Int_omega_BdG_table[s](
                flat_Q_table[ii], flat_phi_table[ii])
            
            for tt in range(num_h_NV):
                flat_coupling_plus_table_array[tt, ii] = Int_coupling_plus_table[s][tt](
                    flat_Q_table[ii], flat_phi_table[ii])
                flat_coupling_minus_table_array[tt, ii] = Int_coupling_minus_table[s][tt](
                    flat_Q_table[ii], flat_phi_table[ii])
        
        # Calculate DOS contributions
        prefactor_DOS = (1/L**2) * (delta_phi/(2*np.pi)**2)
        
        for jj in range(length_flat_table):
            # Match Mathematica's Lorentzian exactly
            lorentzian_L = (eta_small/np.pi) / (eta_small**2 + 
                          (omega_M * flat_omega_BdG_table[jj] - omega_target_L)**2)
            lorentzian_U = (eta_small/np.pi) / (eta_small**2 + 
                          (omega_M * flat_omega_BdG_table[jj] - omega_target_U)**2)
            
            DOS_L += prefactor_DOS * flat_delta_Q_table[jj] * flat_Q_table[jj] * lorentzian_L
            DOS_U += prefactor_DOS * flat_delta_Q_table[jj] * flat_Q_table[jj] * lorentzian_U
        
        # Calculate Gamma contributions
        omega_d = (h_plank * mu_0 * (gamma*1e9)**2 / (L*1e-6)**3) * 1e8 
        prefactor_gamma = (2*np.pi)**2 * omega_M * omega_d * (delta_phi/(2*np.pi)**2)
        
        for jj in range(length_flat_table):
            # Bose factor - ensure this matches Mathematica exactly
            bose_factor = 2 * n_bose(omega_M * flat_omega_BdG_table[jj]) + 1
            
            lorentzian_L = (eta_small/np.pi) / (eta_small**2 + 
                          (omega_M * flat_omega_BdG_table[jj] - omega_target_L)**2)
            lorentzian_U = (eta_small/np.pi) / (eta_small**2 + 
                          (omega_M * flat_omega_BdG_table[jj] - omega_target_U)**2)
            
            # Gamma L contributions (uses coupling_plus)
            coupling_plus_squared = np.abs(flat_coupling_plus_table_array[:, jj])**2
            gamma_L_Hz_array += (prefactor_gamma * bose_factor * 
                               flat_delta_Q_table[jj] * flat_Q_table[jj] * 
                               coupling_plus_squared * lorentzian_L)
            
            # Gamma U contributions (uses coupling_minus)
            coupling_minus_squared = np.abs(flat_coupling_minus_table_array[:, jj])**2
            gamma_U_Hz_array += (prefactor_gamma * bose_factor * 
                               flat_delta_Q_table[jj] * flat_Q_table[jj] * 
                               coupling_minus_squared * lorentzian_U)
    
    return DOS_L, DOS_U, gamma_L_Hz_array, gamma_U_Hz_array


def gamma_1_from_H0(H0):
    print('In gamma_1_from_H0')
    """Calculate Gamma_1 for given H0"""
    omega_H = gamma * H0
    
    # Setup calculation parameters
    num_phi = 2 * 9
    num_q = 2 * 2
    del_phi = np.pi / num_phi
    phi_table = np.arange(0, np.pi, del_phi)
    
    q_max = 50 * L
    def f_space(min_val, max_val, steps):
        f_min = np.log(min_val)
        f_max = np.log(max_val)
        f_range = np.linspace(f_min, f_max, steps)
        return np.exp(f_range)
    q_table = np.concatenate([[1e-6], f_space(L/1000, q_max, num_q)])
    
    F_max = 5
    N_max = int(np.ceil(L/np.pi * np.sqrt(F_max/(gamma * DD))))
    # print(f"N_max is {N_max}")
    N_max = 5
    
    h_NV_array = np.array([0.4]) #, 0.5, 0.6, 0.7])  # positions in μm
    
    # Perform multi-paraunitary diagonalization
    result = multi_para_diag(h_NV_array, omega_H, q_table, phi_table, N_max)
    (Int_omega_BdG_table, Int_coupling_plus_table, Int_coupling_minus_table, Int_coupling_z_table) = result
    
    # Calculate Gamma values
    eta_small = 0.003
    NQ_calc = 2 * 2
    N_phi_calc = 2 * 3
    
    DOS_L, DOS_U, gamma_L_Hz_array, gamma_U_Hz_array = gamma_values_UL(
        eta_small, NQ_calc, N_phi_calc, H0, Int_omega_BdG_table,
        Int_coupling_plus_table, Int_coupling_minus_table, len(h_NV_array))
    
    print(f"Γ(ω=ωL)={gamma_L_Hz_array} Hz")
    print(f"Γ(ω=ωU)={gamma_U_Hz_array} Hz")
    
    return DOS_L, DOS_U, gamma_L_Hz_array, gamma_U_Hz_array


########################################################################################
########################################################################################


# Example usage
if __name__ == "__main__":

    # Compare Python results with MMA results for a single mode and angle

    # # H_BdG_test = generate_H_BdG(gamma * 82, 1e-6, np.pi, 3)
    # print('Checked Hamiltonian generation (generate_H_BdG), matches MMA')
    # # H_hermitian_test = (H_BdG_test + H_BdG_test.conj().T) / 2
    # # H = H_hermitian_test
    # # test_result = paraunitary_diag(H)
    # print('Checked Para-unitary diagonalization (paraunitary_diag), matches MMA')
    # # test_result = multi_para_diag([0.4], gamma * 82, [1e-6], [np.pi], 2)
    # print('Checked Multi-parameter Para-unitary diagonalization (multi_para_diag), matches MMA')

    # print('Checking Gamma calculation (gamma_values_UL)')
    # omega, cplus, cminus, cz = multi_para_diag([0.4], gamma * 82, [1e-6], [np.pi], 4)
    # gamma_test_result = gamma_values_UL(0.003, 4, 2, 82, np.array(omega), np.array(cplus), np.array(cminus), 1)
    # print('gamma_test_result ', gamma_test_result)


   # Test with a single H0 value
    H0_test = 82
    print(f"Calculating for H0 = {H0_test} Oe")
    
    try:
        DOS_L, DOS_U, gamma_L, gamma_U = gamma_1_from_H0(H0_test)
        print(f"DOS_L: {DOS_L}")
        print(f"DOS_U: {DOS_U}")
        print(f"gamma_L: {gamma_L}")
        print(f"gamma_U: {gamma_U}")
    except Exception as e:
        print(f"Error in calculation: {e}")
    
   #  For plotting multiple H0 values (uncomment to use)
   #  H0_array = np.array([76, 77, 78, 79, 80, 81, 81.5, 82, 82.5, 83, 84, 85])
   #  results would be stored and plotted here