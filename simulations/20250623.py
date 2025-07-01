import numpy as np
from scipy.linalg import cholesky, eig, eigh

# from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.interpolate import LinearNDInterpolator
import warnings
import cProfile
import pstats
import time as time

warnings.filterwarnings("ignore")

# Constants and initial calculations
delta_h = 1
xx_solution = 2 * ((2870 / 2.8) - 2 * (82 + delta_h))
ms_sol = xx_solution - 4 + 2
print("Ms_sol ", ms_sol)

# Material & condition parameters
M0 = ms_sol  # Oe
L = 3  # μm
DD = 5.4e-9 * 1e8  # Oe μm^-2
gamma = 2.8e-3  # GHz/G
omega_M = gamma * M0  # GHz


# Scipy Cholesky decomposition
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


def mathematica_eigensystem(matrix, tolerance=1e-10):
    """
    Replicate Mathematica's Eigensystem behavior exactly
    """
    matrix = np.array(matrix, dtype=float)

    # Get eigenvalues and eigenvectors from scipy
    eigenvals, eigenvecs = eig(matrix)

    # # Convert to real if they're essentially real
    # if np.allclose(eigenvals.imag, 0, atol=tolerance): # CHANGE THIS TO SPEEDUP
    #     eigenvals = eigenvals.real
    # if np.allclose(eigenvecs.imag, 0, atol=tolerance): # CHANGE THIS TO SPEEDUP
    #     eigenvecs = eigenvecs.real

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
    dim = H.shape[0]

    # Create sigma3 matrix
    sigma3 = np.diag([(-1) ** np.floor((2 * n - 1) / dim) for n in range(1, dim + 1)])
    W = K @ sigma3 @ K.conj().T

    # Eigendecomposition
    eval_w, evec_w = mathematica_eigensystem(W)
    evec_w = evec_w.T
    evec_w = evec_w / np.linalg.norm(evec_w, axis=0)

    # Permutation for ordering
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


def F_func_vectorized(q, n):
    """Vectorized F function for dipolar calculations"""
    q = np.asarray(q)
    n = np.asarray(n)

    # Handle q=0 case
    zero_mask = q == 0
    result = np.zeros_like(q, dtype=float)

    # For q=0: return 2*n
    if np.isscalar(n):
        result[zero_mask] = 2 * n
    else:
        result[zero_mask] = 2 * n[zero_mask] if zero_mask.any() else 0

    # For q!=0: return 2 * (1 - (-1)^n * exp(-q)) / q
    nonzero_mask = ~zero_mask
    if nonzero_mask.any():
        q_nz = q[nonzero_mask] if q.shape else q
        n_nz = n[nonzero_mask] if hasattr(n, "shape") and n.shape else n
        result[nonzero_mask] = 2 * (1 - (-1) ** n_nz * np.exp(-q_nz)) / q_nz

    return result


def P_func_vectorized(q, n, m):
    """Vectorized P function for dipolar calculations"""
    q = np.asarray(q)
    n = np.asarray(n)
    m = np.asarray(m)

    # Broadcast to common shape
    q, n, m = np.broadcast_arrays(q, n, m)

    kronecker_nm = (n == m).astype(float)
    kronecker_n0 = (n == 0).astype(float)
    kronecker_m0 = (m == 0).astype(float)

    term1 = q**2 / (q**2 + n**2 * np.pi**2) * kronecker_nm

    normalization = 1 / np.sqrt((1 + kronecker_n0) * (1 + kronecker_m0))
    F_vals = F_func_vectorized(q, n)

    term2 = (
        q**4
        / ((q**2 + n**2 * np.pi**2) * (q**2 + m**2 * np.pi**2))
        * F_vals
        * (1 + (-1) ** (n + m))
        / 2
        * normalization
    )

    return term1 - term2


def Q_func_vectorized(q, n, m):
    """Vectorized Q function for dipolar calculations"""
    q = np.asarray(q)
    n = np.asarray(n)
    m = np.asarray(m)

    # Broadcast to common shape
    q, n, m = np.broadcast_arrays(q, n, m)

    kronecker_n0 = (n == 0).astype(float)
    kronecker_m0 = (m == 0).astype(float)

    # Return zeros where n == m
    result = np.zeros_like(q, dtype=float)
    mask = n != m

    if not mask.any():
        return result

    # Apply mask for non-equal n, m
    q_masked = q[mask]
    n_masked = n[mask]
    m_masked = m[mask]

    normalization = 1 / np.sqrt((1 + kronecker_n0[mask]) * (1 + kronecker_m0[mask]))
    F_vals = F_func_vectorized(q_masked, n_masked)

    term1 = (
        q_masked**2
        / (q_masked**2 + m_masked**2 * np.pi**2)
        * (
            m_masked**2
            / (m_masked**2 - n_masked**2 + (1 + (-1) ** (n_masked + m_masked)) / 2)
            * 2
            / q_masked
            - q_masked**2 / (2 * (q_masked**2 + n_masked**2 * np.pi**2)) * F_vals
        )
    )

    result[mask] = term1 * normalization * (1 - (-1) ** (n_masked + m_masked)) / 2
    return result


def capital_omega_vectorized(omega_H, q, n):
    """Vectorized Capital Omega function"""
    q = np.asarray(q)
    n = np.asarray(n)

    # Broadcast arrays
    q, n = np.broadcast_arrays(q, n)

    return (omega_H + (gamma * DD) / L**2 * (q**2 + n**2 * np.pi**2)) / omega_M


def H_matrix_vectorized(omega_H, q, phi_k, n_array, m_array):
    """Vectorized H matrix calculation for all n,m pairs"""

    # Create meshgrids for all n,m combinations
    N, M = np.meshgrid(n_array, m_array, indexing="ij")
    N_flat = N.flatten()
    M_flat = M.flatten()

    # Vectorized calculations
    omega_vals = capital_omega_vectorized(omega_H, q, N_flat)
    P_vals = P_func_vectorized(q, N_flat, M_flat)
    Q_vals = Q_func_vectorized(q, N_flat, M_flat)

    # Create all H matrices at once
    num_pairs = len(N_flat)
    H_matrices = np.zeros((num_pairs, 2, 2), dtype=complex)

    # Diagonal terms (where n == m)
    diag_mask = N_flat == M_flat
    if diag_mask.any():
        H_matrices[diag_mask, :, :] += (
            omega_vals[diag_mask, np.newaxis, np.newaxis] * np.eye(2)[np.newaxis, :, :]
        )
        H_matrices[diag_mask, :, :] += 0.5 * np.ones((2, 2))[np.newaxis, :, :]

    # P terms
    sin_phi = np.sin(phi_k)
    P_matrix = np.array(
        [[1 - sin_phi**2, 1 + sin_phi**2], [1 + sin_phi**2, 1 - sin_phi**2]]
    )
    H_matrices -= 0.5 * P_vals[:, np.newaxis, np.newaxis] * P_matrix[np.newaxis, :, :]

    # Q terms
    Q_matrix = np.array([[0, -4], [4, 0]])
    H_matrices -= (
        0.5 * Q_vals[:, np.newaxis, np.newaxis] * sin_phi * Q_matrix[np.newaxis, :, :]
    )

    # Reshape back to (N_max, N_max, 2, 2)
    return H_matrices.reshape(len(n_array), len(m_array), 2, 2)


def generate_H_BdG_vectorized(omega_H, q, phi_k, N_max):
    """Vectorized BdG Hamiltonian generation"""
    n_array = np.arange(N_max)
    m_array = np.arange(N_max)

    # Get all H matrices at once
    H_matrices = H_matrix_vectorized(omega_H, q, phi_k, n_array, m_array)

    # Assemble the full BdG matrix
    H_BdG = np.zeros((2 * N_max, 2 * N_max), dtype=complex)

    for n in range(N_max):
        for m in range(N_max):
            H_nm = H_matrices[n, m]
            H_BdG[n, m] = H_nm[0, 0]
            H_BdG[n, m + N_max] = H_nm[0, 1]
            H_BdG[n + N_max, m] = H_nm[1, 0]
            H_BdG[n + N_max, m + N_max] = H_nm[1, 1]

    return H_BdG


def f_func_vectorized(q_array, n_array, h_over_L_array):
    """Vectorized f function for coupling calculations"""
    # Create 3D arrays for broadcasting: (h_NV, q, n)
    Q, N, H = np.meshgrid(q_array, n_array, h_over_L_array, indexing="ij")
    Q = Q.transpose(2, 0, 1)  # Shape: (h_NV, q, n)
    N = N.transpose(2, 0, 1)
    H = H.transpose(2, 0, 1)

    kronecker_n0 = (N == 0).astype(float)
    normalization = 1 / np.sqrt(2 * (1 + kronecker_n0))

    # Handle q=0 case
    zero_mask = Q == 0
    result = np.zeros_like(Q, dtype=complex)

    # For q!=0
    nonzero_mask = ~zero_mask
    if nonzero_mask.any():
        Q_nz = Q[nonzero_mask]
        N_nz = N[nonzero_mask]
        H_nz = H[nonzero_mask]
        norm_nz = normalization[nonzero_mask]

        result[nonzero_mask] = (
            (-1) ** N_nz
            * norm_nz
            * Q_nz**2
            / (Q_nz**2 + N_nz**2 * np.pi**2)
            * np.exp(-Q_nz * H_nz)
            * (1 - (-1) ** N_nz * np.exp(-Q_nz))
        )

    return result


def multi_para_diag_vectorized(h_NV_array, omega_H, q_table, phi_table, N_max):
    """Highly vectorized multi paraunitary diagonalization"""

    q_table = np.array(q_table)
    phi_table = np.array(phi_table)
    h_NV_array = np.array(h_NV_array)

    num_q = len(q_table)
    num_phi = len(phi_table)
    num_h_NV = len(h_NV_array)

    # Pre-compute f_bar arrays for all combinations - VECTORIZED
    print("Pre-computing f_bar arrays (vectorized)...")
    n_array = np.arange(N_max)
    h_over_L_array = h_NV_array / L

    # Vectorized f_bar computation: shape (num_h_NV, num_q, N_max)
    f_bar_all = f_func_vectorized(q_table, n_array, h_over_L_array)

    # Initialize result arrays
    omega_BdG_all = np.zeros((num_q, 2 * num_phi, N_max))
    coupling_plus_all = np.zeros((num_q, 2 * num_phi, num_h_NV, N_max), dtype=complex)
    coupling_minus_all = np.zeros((num_q, 2 * num_phi, num_h_NV, N_max), dtype=complex)
    coupling_z_all = np.zeros((num_q, 2 * num_phi, num_h_NV, N_max), dtype=complex)

    print("Starting vectorized MultiParaDiag calculation...")

    # Pre-compute trigonometric values
    sin_phi = np.sin(phi_table)
    cos_phi = np.cos(phi_table)
    sin_phi_mirror = np.sin(phi_table + np.pi)
    cos_phi_mirror = np.cos(phi_table + np.pi)

    # Main calculation loop
    for count_q in range(num_q):
        progress = (count_q + 1) / num_q * 100
        print(f"Progress: {progress:.1f}% MultiParaDiag (vectorized)")

        q_current = q_table[count_q]

        for count_phi, phi_k in enumerate(phi_table):
            # Generate H_BdG matrix (this part is still hard to vectorize due to matrix operations)
            H_BdG = generate_H_BdG_vectorized(omega_H, q_current, phi_k, N_max)

            # Para-diagonalization
            result = paraunitary_diag((H_BdG + H_BdG.conj().T) / 2)

            # Store omega results for both phi and phi+pi
            omega_BdG_all[count_q, count_phi] = result[0]
            omega_BdG_all[count_q, count_phi + num_phi] = result[0]

            T_pp, T_np, T_pn, T_nn = result[1], result[2], result[3], result[4]

            # Faster gamma calculations using pre-computed trig values
            sin_phi_val = sin_phi[count_phi]
            cos_phi_val = cos_phi[count_phi]
            sin_phi_mirror_val = sin_phi_mirror[count_phi]
            cos_phi_mirror_val = cos_phi_mirror[count_phi]

            # Calculate gamma values for phi
            base_term = T_pp + T_np + sin_phi_val * (T_pp - T_np)
            gamma_plus = (1 + sin_phi_val) / 2 * base_term
            gamma_minus = (1 - sin_phi_val) / 2 * base_term
            gamma_z = -1j * cos_phi_val / 2 * base_term

            # Calculate gamma values for phi+pi (mirror)
            base_term_mirror = T_nn + T_pn + sin_phi_mirror_val * (T_nn - T_pn)
            gamma_plus_mirror = (1 + sin_phi_mirror_val) / 2 * np.conj(base_term_mirror)
            gamma_minus_mirror = (
                (1 - sin_phi_mirror_val) / 2 * np.conj(base_term_mirror)
            )
            gamma_z_mirror = -1j * cos_phi_mirror_val / 2 * np.conj(base_term_mirror)

            # Vectorized coupling calculations using pre-computed f_bar
            f_bar_current = f_bar_all[:, count_q, :]  # Shape: (num_h_NV, N_max)

            # Matrix multiplication for all h_NV at once
            v_plus_array = f_bar_current @ gamma_plus
            v_minus_array = f_bar_current @ gamma_minus
            v_z_array = f_bar_current @ gamma_z

            v_plus_mirror_array = f_bar_current @ gamma_plus_mirror
            v_minus_mirror_array = f_bar_current @ gamma_minus_mirror
            v_z_mirror_array = f_bar_current @ gamma_z_mirror

            # Store results
            coupling_plus_all[count_q, count_phi] = v_plus_array
            coupling_plus_all[count_q, count_phi + num_phi] = v_plus_mirror_array

            coupling_minus_all[count_q, count_phi] = v_minus_array
            coupling_minus_all[count_q, count_phi + num_phi] = v_minus_mirror_array

            coupling_z_all[count_q, count_phi] = v_z_array
            coupling_z_all[count_q, count_phi + num_phi] = v_z_mirror_array

    # Create extended tables with periodic boundaries (vectorized)
    print("Creating extended tables (vectorized)...")

    # Pre-allocate extended arrays
    omega_BdG_extended = np.zeros((num_q, 2 * num_phi + 1, N_max))
    coupling_plus_extended = np.zeros(
        (num_q, 2 * num_phi + 1, num_h_NV, N_max), dtype=complex
    )
    coupling_minus_extended = np.zeros(
        (num_q, 2 * num_phi + 1, num_h_NV, N_max), dtype=complex
    )
    coupling_z_extended = np.zeros(
        (num_q, 2 * num_phi + 1, num_h_NV, N_max), dtype=complex
    )

    # Vectorized copy operations
    omega_BdG_extended[:, : 2 * num_phi] = omega_BdG_all
    omega_BdG_extended[:, 2 * num_phi] = omega_BdG_all[:, 0]

    coupling_plus_extended[:, : 2 * num_phi] = coupling_plus_all
    coupling_plus_extended[:, 2 * num_phi] = coupling_plus_all[:, 0]

    coupling_minus_extended[:, : 2 * num_phi] = coupling_minus_all
    coupling_minus_extended[:, 2 * num_phi] = coupling_minus_all[:, 0]

    coupling_z_extended[:, : 2 * num_phi] = coupling_z_all
    coupling_z_extended[:, 2 * num_phi] = coupling_z_all[:, 0]

    # Create coordinate grids for interpolation (vectorized)
    print("Creating interpolation functions (vectorized)...")

    Q_extended = np.tile(
        q_table[:, np.newaxis, np.newaxis], (1, 2 * num_phi + 1, N_max)
    )

    # Vectorized Phi coordinate creation
    Phi_extended = np.zeros((num_q, 2 * num_phi + 1, N_max))
    phi_coords = np.concatenate([phi_table, phi_table + np.pi, [phi_table[0]]])
    Phi_extended[:, :, :] = phi_coords[np.newaxis, :, np.newaxis]

    def create_interpolation_function_vectorized(data_table):
        """Optimized vectorized interpolation function creation"""
        interpolation_functions = []

        for s in range(N_max):
            xcoords = Q_extended[:, :, s].flatten()
            ycoords = Phi_extended[:, :, s].flatten()

            if len(data_table.shape) == 3:  # omega table
                zcoords = data_table[:, :, s].flatten()
                func = LinearNDInterpolator(
                    np.column_stack([xcoords, ycoords]), zcoords
                )
                interpolation_functions.append(func)
            else:  # coupling tables
                interp_functions_tt = []
                for tt in range(num_h_NV):
                    zcoords = data_table[:, :, tt, s].flatten()
                    func = LinearNDInterpolator(
                        np.column_stack([xcoords, ycoords]), zcoords
                    )
                    interp_functions_tt.append(func)
                interpolation_functions.append(interp_functions_tt)

        return interpolation_functions

    # Create all interpolation functions
    Int_omega_BdG_table = create_interpolation_function_vectorized(omega_BdG_extended)
    Int_coupling_plus_table = create_interpolation_function_vectorized(
        coupling_plus_extended
    )
    Int_coupling_minus_table = create_interpolation_function_vectorized(
        coupling_minus_extended
    )
    Int_coupling_z_table = create_interpolation_function_vectorized(coupling_z_extended)

    return (
        Int_omega_BdG_table,
        Int_coupling_plus_table,
        Int_coupling_minus_table,
        Int_coupling_z_table,
    )


# Temperature and Bose-Einstein distribution
h_plank = 6.626e-34  # J*s
mu_0 = 4 * np.pi * 1e-7  # H/m
k_B = 1.381e-23  # J/K
temperature = 300  # K


def n_bose_vectorized(omega):
    """Bose-Einstein distribution (high temperature approximation)"""
    return (1e-9 * k_B * temperature / h_plank) / omega


def gamma_values_UL_vectorized(
    eta_small,
    NQest,
    N_phi,
    H0,
    Int_omega_BdG_table,
    Int_coupling_plus_table,
    Int_coupling_minus_table,
    num_h_NV,
):
    """Fully vectorized Gamma values calculation"""
    print("In gamma_values_UL_vectorized")

    DNV = 2.87  # GHz
    omega_target_L = DNV - H0 * gamma  # GHz
    omega_target_U = DNV + H0 * gamma  # GHz

    # Create q and phi grids (vectorized)
    q_max = 50 * L
    q_middle = 5 * L
    NQ_half = round(NQest / 2)

    def f_space(min_val, max_val, steps):
        """Vectorized logarithmic spacing"""
        f_range = np.linspace(np.log(min_val), np.log(max_val), steps)
        return np.exp(f_range)

    Q_table = np.concatenate([[1e-6], f_space(L / 1000, q_middle, NQ_half)])

    # Vectorized Q_table extension
    last_delta_q = np.diff(Q_table)[-1]
    n_additional = int((q_max - Q_table[-1]) / last_delta_q)
    if n_additional > 0:
        additional_points = Q_table[-1] + last_delta_q * np.arange(1, n_additional + 1)
        Q_table = np.concatenate([Q_table, additional_points])

    delta_q = np.diff(Q_table)
    NQ = len(delta_q)
    delta_phi = 2 * np.pi / N_phi

    # Pre-allocate arrays
    gamma_L_Hz_array = np.zeros(num_h_NV)
    gamma_U_Hz_array = np.zeros(num_h_NV)
    DOS_L = 0
    DOS_U = 0

    N_max_local = len(Int_omega_BdG_table)

    # Create meshgrids for all calculations at once
    q_mesh, phi_mesh = np.meshgrid(
        Q_table[:-1], np.linspace(0, 2 * np.pi * (1 - 1 / N_phi), N_phi), indexing="ij"
    )
    delta_q_mesh = np.broadcast_to(delta_q[:, np.newaxis], q_mesh.shape)

    # Flatten for vectorized operations
    flat_Q = q_mesh.flatten()
    flat_phi = phi_mesh.flatten()
    flat_delta_Q = delta_q_mesh.flatten()

    length_flat_table = len(flat_Q)

    for s in range(N_max_local):
        progress = (s + 1) / N_max_local * 100
        print(f"performing s={s+1} calculation: {progress:.1f}% Gamma1Calc")

        # Vectorized function evaluations
        flat_omega_BdG = Int_omega_BdG_table[s](flat_Q, flat_phi)

        # Pre-allocate coupling arrays
        flat_coupling_plus = np.zeros((num_h_NV, length_flat_table), dtype=complex)
        flat_coupling_minus = np.zeros((num_h_NV, length_flat_table), dtype=complex)

        # Vectorized coupling calculations
        for tt in range(num_h_NV):
            flat_coupling_plus[tt] = Int_coupling_plus_table[s][tt](flat_Q, flat_phi)
            flat_coupling_minus[tt] = Int_coupling_minus_table[s][tt](flat_Q, flat_phi)

        # Vectorized DOS calculations
        prefactor_DOS = (1 / L**2) * (delta_phi / (2 * np.pi) ** 2)

        # All Lorentzian calculations at once
        omega_diff_L = omega_M * flat_omega_BdG - omega_target_L
        omega_diff_U = omega_M * flat_omega_BdG - omega_target_U

        lorentzian_L = (eta_small / np.pi) / (eta_small**2 + omega_diff_L**2)
        lorentzian_U = (eta_small / np.pi) / (eta_small**2 + omega_diff_U**2)

        # Vectorized DOS accumulation
        dos_weight = prefactor_DOS * flat_delta_Q * flat_Q
        DOS_L += np.sum(dos_weight * lorentzian_L)
        DOS_U += np.sum(dos_weight * lorentzian_U)

        # Vectorized Gamma calculations
        omega_d = (h_plank * mu_0 * (gamma * 1e9) ** 2 / (L * 1e-6) ** 3) * 1e8
        prefactor_gamma = (
            (2 * np.pi) ** 2 * omega_M * omega_d * (delta_phi / (2 * np.pi) ** 2)
        )

        # Vectorized Bose factor calculation
        bose_factor = 2 * n_bose_vectorized(omega_M * flat_omega_BdG) + 1

        # Common prefactor for all gamma calculations
        gamma_weight = prefactor_gamma * bose_factor * flat_delta_Q * flat_Q

        # Vectorized gamma accumulation using broadcasting
        coupling_plus_squared = (
            np.abs(flat_coupling_plus) ** 2
        )  # shape: (num_h_NV, length_flat_table)
        coupling_minus_squared = np.abs(flat_coupling_minus) ** 2

        gamma_L_Hz_array += np.sum(
            coupling_plus_squared * (gamma_weight * lorentzian_L), axis=1
        )
        gamma_U_Hz_array += np.sum(
            coupling_minus_squared * (gamma_weight * lorentzian_U), axis=1
        )

    return DOS_L, DOS_U, gamma_L_Hz_array, gamma_U_Hz_array


def gamma_1_from_H0_vectorized(H0):
    """Vectorized version of gamma_1_from_H0"""
    print("In gamma_1_from_H0_vectorized")

    omega_H = gamma * H0

    # Setup calculation parameters
    num_phi = 2 * 180
    num_q = 2 * 200
    del_phi = np.pi / num_phi

    # Vectorized phi_table creation
    phi_table = np.linspace(0, np.pi - del_phi, num_phi)

    q_max = 50 * L

    def f_space(min_val, max_val, steps):
        f_range = np.linspace(np.log(min_val), np.log(max_val), steps)
        return np.exp(f_range)

    q_table = np.concatenate([[1e-6], f_space(L / 1000, q_max, num_q)])

    # F_max = 5
    # N_max = int(np.ceil(L / np.pi * np.sqrt(F_max / (gamma * DD))))
    # # N_max = 1
    # print(f"N_max is {N_max}")

    h_NV_array = np.array([0.4])  # Can be extended: [0.4, 0.5, 0.6, 0.7]

    # Perform multi-paraunitary diagonalization
    result = multi_para_diag(h_NV_array, omega_H, q_table, phi_table, N_max)
    (
        Int_omega_BdG_table,
        Int_coupling_plus_table,
        Int_coupling_minus_table,
        Int_coupling_z_table,
    ) = result

    # Calculate Gamma values using vectorized function
    eta_small = 0.003
    NQ_calc = 2 * 200
    N_phi_calc = 2 * 360

    print("Int_omega_BdG_table ", Int_omega_BdG_table)
    print("np.shape(Int_omega_BdG_table) ", np.shape(Int_omega_BdG_table))

    DOS_L, DOS_U, gamma_L_Hz_array, gamma_U_Hz_array = gamma_values_UL_vectorized(
        eta_small,
        NQ_calc,
        N_phi_calc,
        H0,
        Int_omega_BdG_table,
        Int_coupling_plus_table,
        Int_coupling_minus_table,
        len(h_NV_array),
    )

    return DOS_L, DOS_U, gamma_L_Hz_array, gamma_U_Hz_array


def batch_gamma_calculation(H0_array):
    """Vectorized batch calculation for multiple H0 values"""
    print("Starting batch gamma calculation")

    n_H0 = len(H0_array)
    DOSL_list = np.zeros(n_H0)
    DOSU_list = np.zeros(n_H0)
    GammaL_list = np.zeros(n_H0)
    GammaU_list = np.zeros(n_H0)

    # Pre-compute common parameters that don't depend on H0
    eta_small = 0.003
    NQ_calc = 2 * 200  # to solve for the 55 relevant modes
    N_phi_calc = 2 * 360  # to solve for the 55 relevant modes
    num_phi = 2 * 180  # to create interpolation tables
    num_q = 2 * 200  # to create interpolation tables

    # Pre-compute grids
    del_phi = np.pi / num_phi
    phi_table = np.linspace(0, np.pi - del_phi, num_phi)

    q_max = 50 * L

    def f_space(min_val, max_val, steps):
        f_range = np.linspace(np.log(min_val), np.log(max_val), steps)
        return np.exp(f_range)

    q_table = np.concatenate([[1e-6], f_space(L / 1000, q_max, num_q)])
    h_NV_array = np.array([0.4])

    F_max = 5
    # N_max = int(np.ceil(L / np.pi * np.sqrt(F_max / (gamma * DD))))
    N_max = 1
    print(f"N_max is {N_max}")

    for h, H0_test in enumerate(H0_array):
        print(f"Calculating for H0 = {H0_test} Oe ({h+1}/{n_H0})")

        omega_H = gamma * H0_test

        # This is the main bottleneck - diagonalization for each H0
        result = multi_para_diag_vectorized(
            h_NV_array, omega_H, q_table, phi_table, N_max
        )
        (
            Int_omega_BdG_table,
            Int_coupling_plus_table,
            Int_coupling_minus_table,
            Int_coupling_z_table,
        ) = result

        DOSL_list[h], DOSU_list[h], GammaL_list[h], GammaU_list[h] = (
            gamma_values_UL_vectorized(
                eta_small,
                NQ_calc,
                N_phi_calc,
                H0_test,
                Int_omega_BdG_table,
                Int_coupling_plus_table,
                Int_coupling_minus_table,
                len(h_NV_array),
            )
        )

    return DOSL_list, DOSU_list, GammaL_list, GammaU_list


# Additional optimization utilities
def optimize_grid_creation(NQ, N_phi, L):
    """Pre-compute and cache grid arrays for repeated use"""
    q_max = 50 * L
    q_middle = 5 * L
    NQ_half = round(NQ / 2)

    def f_space(min_val, max_val, steps):
        f_range = np.linspace(np.log(min_val), np.log(max_val), steps)
        return np.exp(f_range)

    Q_table = np.concatenate([[1e-6], f_space(L / 1000, q_middle, NQ_half)])

    # Vectorized extension
    last_delta_q = np.diff(Q_table)[-1]
    n_additional = int((q_max - Q_table[-1]) / last_delta_q)
    if n_additional > 0:
        additional_points = Q_table[-1] + last_delta_q * np.arange(1, n_additional + 1)
        Q_table = np.concatenate([Q_table, additional_points])

    delta_q = np.diff(Q_table)

    # Pre-compute meshgrids
    q_mesh, phi_mesh = np.meshgrid(
        Q_table[:-1], np.linspace(0, 2 * np.pi * (1 - 1 / N_phi), N_phi), indexing="ij"
    )
    delta_q_mesh = np.broadcast_to(delta_q[:, np.newaxis], q_mesh.shape)

    return {
        "Q_table": Q_table,
        "delta_q": delta_q,
        "q_mesh": q_mesh,
        "phi_mesh": phi_mesh,
        "delta_q_mesh": delta_q_mesh,
        "flat_Q": q_mesh.flatten(),
        "flat_phi": phi_mesh.flatten(),
        "flat_delta_Q": delta_q_mesh.flatten(),
    }


########################################################################################
# Example usage with vectorization
########################################################################################

if __name__ == "__main__":

    # Profile the main execution
    profiler = cProfile.Profile()
    profiler.enable()

    # Your existing code
    start_time = time.time()
    H0_array = np.array([82])
    DOSL_list, DOSU_list, GammaL_list, GammaU_list = batch_gamma_calculation(H0_array)

    print("GammaL ", GammaL_list)
    print("GammaU ", GammaU_list)
    print("Total time taken with vectorization ", time.time() - start_time, "s")

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(H0_array, GammaL_list, "o-", label="GammaL", linewidth=2, markersize=6)
    plt.plot(H0_array, GammaU_list, "s-", label="GammaU", linewidth=2, markersize=6)
    plt.xlabel("B field (G)", fontsize=12)
    plt.ylabel("Gamma (Hz)", fontsize=12)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    # plt.show()

    # disable profiler
    profiler.disable()

    # Save and display results
    stats = pstats.Stats(profiler)
    stats.sort_stats("cumulative")
    stats.print_stats(20)  # Show top 20 functions

    # Or save to file
    stats.dump_stats("profile_results3.prof")


########################################################################################
########################################################################################
