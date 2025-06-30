import numpy as np
from scipy.linalg import cholesky, eig

# from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.interpolate import LinearNDInterpolator
import warnings
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


def F_func(q, n):
    """F function for dipolar calculations"""
    if q == 0:
        return 2 * n
    return 2 * (1 - (-1) ** n * np.exp(-q)) / q


def P_func(q, n, m):
    """P function for dipolar calculations"""
    kronecker_nm = 1 if n == m else 0
    kronecker_n0 = 1 if n == 0 else 0
    kronecker_m0 = 1 if m == 0 else 0

    term1 = q**2 / (q**2 + n**2 * np.pi**2) * kronecker_nm

    normalization = 1 / np.sqrt((1 + kronecker_n0) * (1 + kronecker_m0))
    term2 = (
        q**4
        / ((q**2 + n**2 * np.pi**2) * (q**2 + m**2 * np.pi**2))
        * F_func(q, n)
        * (1 + (-1) ** (n + m))
        / 2
        * normalization
    )

    return term1 - term2


def Q_func(q, n, m):
    """Q function for dipolar calculations"""

    kronecker_n0 = 1 if n == 0 else 0
    kronecker_m0 = 1 if m == 0 else 0

    if n == m:
        return 0

    normalization = 1 / np.sqrt((1 + kronecker_n0) * (1 + kronecker_m0))

    term1 = (
        q**2
        / (q**2 + m**2 * np.pi**2)
        * (
            m**2 / (m**2 - n**2 + (1 + (-1) ** (n + m)) / 2) * 2 / q
            - q**2 / (2 * (q**2 + n**2 * np.pi**2)) * F_func(q, n)
        )
    )

    return term1 * normalization * (1 - (-1) ** (n + m)) / 2


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
    P_matrix = np.array(
        [
            [1 - np.sin(phi_k) ** 2, 1 + np.sin(phi_k) ** 2],
            [1 + np.sin(phi_k) ** 2, 1 - np.sin(phi_k) ** 2],
        ]
    )
    H_nm -= 0.5 * P_matrix * P_val

    # Q terms
    Q_val = Q_func(q, n, m)
    Q_matrix = np.array([[0, -4], [4, 0]])
    H_nm -= 0.5 * Q_matrix * np.sin(phi_k) * Q_val

    return H_nm


def generate_H_BdG(omega_H, q, phi_k, N_max):
    """Generate BdG Hamiltonian"""
    H_BdG = np.zeros((2 * N_max, 2 * N_max), dtype=complex)

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

    return (
        (-1) ** n
        * normalization
        * q**2
        / (q**2 + n**2 * np.pi**2)
        * np.exp(-q * h_over_L)
        * (1 - (-1) ** n * np.exp(-q))
    )


def multi_para_diag(h_NV_array, omega_H, q_table, phi_table, N_max):
    """Vectorized multi paraunitary diagonalization"""

    q_table = np.array(q_table)
    phi_table = np.array(phi_table)
    h_NV_array = np.array(h_NV_array)

    num_q = len(q_table)
    num_phi = len(phi_table)
    num_h_NV = len(h_NV_array)

    # Pre-compute f_bar arrays for all combinations
    print("Pre-computing f_bar arrays...")
    f_bar_all = np.zeros((num_h_NV, num_q, N_max), dtype=complex)
    for tt in range(num_h_NV):
        for q_idx, q in enumerate(q_table):
            for nn in range(N_max):
                f_bar_all[tt, q_idx, nn] = f_func(q, nn, h_NV_array[tt] / L)

    # Initialize result arrays - note: we'll store both original and mirror data
    omega_BdG_all = np.zeros((num_q, 2 * num_phi, N_max))
    coupling_plus_all = np.zeros((num_q, 2 * num_phi, num_h_NV, N_max), dtype=complex)
    coupling_minus_all = np.zeros((num_q, 2 * num_phi, num_h_NV, N_max), dtype=complex)
    coupling_z_all = np.zeros((num_q, 2 * num_phi, num_h_NV, N_max), dtype=complex)

    print("Starting vectorized MultiParaDiag calculation...")

    # Main calculation loop
    for count_q in range(num_q):
        progress = (count_q + 1) / num_q * 100
        print(f"Progress: {progress:.1f}% MultiParaDiag")

        q_current = q_table[count_q]

        for count_phi, phi_k in enumerate(phi_table):
            # Generate H_BdG matrix
            H_BdG = generate_H_BdG(omega_H, q_current, phi_k, N_max)

            # Para-diagonalization
            result = paraunitary_diag((H_BdG + H_BdG.conj().T) / 2)

            # Store omega results for both phi and phi+pi
            omega_BdG_all[count_q, count_phi] = result[0]
            omega_BdG_all[count_q, count_phi + num_phi] = result[0]  # Same for mirror

            T_pp, T_np, T_pn, T_nn = result[1], result[2], result[3], result[4]

            # Calculate gamma values for phi
            sin_phi = np.sin(phi_k)
            cos_phi = np.cos(phi_k)

            gamma_plus = (1 + sin_phi) / 2 * (T_pp + T_np + sin_phi * (T_pp - T_np))
            gamma_minus = (1 - sin_phi) / 2 * (T_pp + T_np + sin_phi * (T_pp - T_np))
            gamma_z = -1j * cos_phi / 2 * (T_pp + T_np + sin_phi * (T_pp - T_np))

            # Calculate gamma values for phi+pi (mirror)
            phi_k_mirror = phi_k + np.pi
            sin_phi_mirror = np.sin(phi_k_mirror)
            cos_phi_mirror = np.cos(phi_k_mirror)

            gamma_plus_mirror = (
                (1 + sin_phi_mirror)
                / 2
                * np.conj(T_nn + T_pn + sin_phi_mirror * (T_nn - T_pn))
            )
            gamma_minus_mirror = (
                (1 - sin_phi_mirror)
                / 2
                * np.conj(T_nn + T_pn + sin_phi_mirror * (T_nn - T_pn))
            )
            gamma_z_mirror = (
                -1j
                * cos_phi_mirror
                / 2
                * np.conj(T_nn + T_pn + sin_phi_mirror * (T_nn - T_pn))
            )

            # Vectorized coupling calculations using pre-computed f_bar
            f_bar_current = f_bar_all[:, count_q, :]  # Shape: (num_h_NV, N_max)

            # Calculate v arrays for original phi
            v_plus_array = f_bar_current @ gamma_plus
            v_minus_array = f_bar_current @ gamma_minus
            v_z_array = f_bar_current @ gamma_z

            # Calculate v arrays for mirror phi+pi
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

    # Create extended tables with periodic boundaries
    print("Creating extended tables...")

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

    # Copy data and add periodic boundary
    omega_BdG_extended[:, : 2 * num_phi] = omega_BdG_all
    omega_BdG_extended[:, 2 * num_phi] = omega_BdG_all[:, 0]  # Periodic

    coupling_plus_extended[:, : 2 * num_phi] = coupling_plus_all
    coupling_plus_extended[:, 2 * num_phi] = coupling_plus_all[:, 0]

    coupling_minus_extended[:, : 2 * num_phi] = coupling_minus_all
    coupling_minus_extended[:, 2 * num_phi] = coupling_minus_all[:, 0]

    coupling_z_extended[:, : 2 * num_phi] = coupling_z_all
    coupling_z_extended[:, 2 * num_phi] = coupling_z_all[:, 0]

    # Create coordinate grids for interpolation
    print("Creating interpolation functions...")

    # Create Q coordinates (same for all phi)
    Q_extended = np.tile(
        q_table[:, np.newaxis, np.newaxis], (1, 2 * num_phi + 1, N_max)
    )

    # Create Phi coordinates properly
    Phi_extended = np.zeros((num_q, 2 * num_phi + 1, N_max))
    for n in range(num_q):
        for m in range(2 * num_phi + 1):
            if m < num_phi:
                # Original phi values
                Phi_extended[n, m, :] = phi_table[m]
            elif m < 2 * num_phi:
                # Mirror phi values (phi + pi)
                phi_idx = m - num_phi
                Phi_extended[n, m, :] = phi_table[phi_idx] + np.pi
            else:
                # Periodic boundary (same as first point)
                Phi_extended[n, m, :] = phi_table[0]

    def create_interpolation_function_vectorized(data_table):
        """Vectorized interpolation function creation"""
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


def n_bose(omega):
    """Bose-Einstein distribution (high temperature approximation)"""
    return (1e-9 * k_B * temperature / h_plank) / omega


def gamma_values_UL(
    eta_small,
    NQest,
    N_phi,
    H0,
    Int_omega_BdG_table,
    Int_coupling_plus_table,
    Int_coupling_minus_table,
    num_h_NV,
):
    print("In gamma_values_UL")
    """Corrected Gamma values calculation"""
    DNV = 2.87  # GHz
    omega_target_L = DNV - H0 * gamma  # GHz
    omega_target_U = DNV + H0 * gamma  # GHz

    # Convert to Hz for proper units
    omega_target_L_Hz = omega_target_L  # GHz
    omega_target_U_Hz = omega_target_U  # GHz
    omega_M_Hz = omega_M  # GHz

    # Create q and phi grids
    q_max = 50 * L
    q_middle = 5 * L
    NQ_half = round(NQest / 2)

    def f_space(min_val, max_val, steps):
        """Match Mathematica's fSpace function exactly"""
        f_min = np.log(min_val)
        f_max = np.log(max_val)
        f_range = np.linspace(f_min, f_max, steps)
        return np.exp(f_range)

    Q_table = np.concatenate([[1e-6], f_space(L / 1000, q_middle, NQ_half)])

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

    N_max_local = np.shape(Int_omega_BdG_table)[0]

    for s in range(N_max_local):
        progress = (s + 1) / N_max_local * 100
        print(f"performing s={s+1} calculation: {progress:.1f}% Gamma1Calc")

        # Create flat tables
        flat_Q_table = []
        for ii1 in range(NQ):
            for ii2 in range(N_phi):
                flat_Q_table.append(Q_table[ii1])
        flat_delta_Q_table = np.repeat(delta_q, N_phi)
        flat_phi_table = np.tile(np.linspace(0, 2 * np.pi * (1 - 1 / N_phi), N_phi), NQ)

        length_flat_table = len(flat_Q_table)
        flat_omega_BdG_table = np.zeros(length_flat_table)
        flat_coupling_plus_table_array = np.zeros(
            (num_h_NV, length_flat_table), dtype=complex
        )
        flat_coupling_minus_table_array = np.zeros(
            (num_h_NV, length_flat_table), dtype=complex
        )

        # Vectorize all at once
        flat_omega_BdG_table[:] = Int_omega_BdG_table[s](flat_Q_table, flat_phi_table)

        for tt in range(num_h_NV):
            flat_coupling_plus_table_array[tt, :] = Int_coupling_plus_table[s][tt](
                flat_Q_table, flat_phi_table
            )
            flat_coupling_minus_table_array[tt, :] = Int_coupling_minus_table[s][tt](
                flat_Q_table, flat_phi_table
            )

        # Calculate DOS contributions
        prefactor_DOS = (1 / L**2) * (delta_phi / (2 * np.pi) ** 2)

        # Vectorized version:
        lorentzian_L_array = (eta_small / np.pi) / (
            eta_small**2 + (omega_M * flat_omega_BdG_table - omega_target_L) ** 2
        )
        lorentzian_U_array = (eta_small / np.pi) / (
            eta_small**2 + (omega_M * flat_omega_BdG_table - omega_target_U) ** 2
        )

        DOS_L += np.sum(
            prefactor_DOS * flat_delta_Q_table * flat_Q_table * lorentzian_L_array
        )
        DOS_U += np.sum(
            prefactor_DOS * flat_delta_Q_table * flat_Q_table * lorentzian_U_array
        )

        # Calculate Gamma contributions
        omega_d = (h_plank * mu_0 * (gamma * 1e9) ** 2 / (L * 1e-6) ** 3) * 1e8
        prefactor_gamma = (
            (2 * np.pi) ** 2 * omega_M * omega_d * (delta_phi / (2 * np.pi) ** 2)
        )
        bose_factor = 2 * n_bose(omega_M * flat_omega_BdG_table) + 1

        # Calculate all contributions at once
        prefactor_common = (
            prefactor_gamma * bose_factor * flat_delta_Q_table * flat_Q_table
        )

        # Gamma L: shape (num_h_NV, length_flat_table) -> sum over length_flat_table
        coupling_plus_squared = np.abs(flat_coupling_plus_table_array) ** 2
        gamma_L_contributions = coupling_plus_squared * (
            prefactor_common * lorentzian_L_array
        )
        gamma_L_Hz_array += np.sum(gamma_L_contributions, axis=1)

        # Same for gamma U
        coupling_minus_squared = np.abs(flat_coupling_minus_table_array) ** 2
        gamma_U_contributions = coupling_minus_squared * (
            prefactor_common * lorentzian_U_array
        )
        gamma_U_Hz_array += np.sum(gamma_U_contributions, axis=1)

    return DOS_L, DOS_U, gamma_L_Hz_array, gamma_U_Hz_array


def gamma_1_from_H0(H0):
    print("In gamma_1_from_H0")
    """Calculate Gamma_1 for given H0"""
    omega_H = gamma * H0

    # Setup calculation parameters
    num_phi = 2 * 180
    num_q = 2 * 200
    del_phi = np.pi / num_phi
    phi_table = np.arange(0, np.pi, del_phi)

    q_max = 50 * L

    def f_space(min_val, max_val, steps):
        f_min = np.log(min_val)
        f_max = np.log(max_val)
        f_range = np.linspace(f_min, f_max, steps)
        return np.exp(f_range)

    q_table = np.concatenate([[1e-6], f_space(L / 1000, q_max, num_q)])

    F_max = 5
    # N_max = int(np.ceil(L / np.pi * np.sqrt(F_max / (gamma * DD))))
    N_max = 1
    print(f"N_max is {N_max}")

    h_NV_array = np.array([0.4])  # , 0.5, 0.6, 0.7])  # positions in μm

    # Perform multi-paraunitary diagonalization
    result = multi_para_diag(h_NV_array, omega_H, q_table, phi_table, N_max)
    (
        Int_omega_BdG_table,
        Int_coupling_plus_table,
        Int_coupling_minus_table,
        Int_coupling_z_table,
    ) = result

    # Calculate Gamma values
    eta_small = 0.003
    NQ_calc = 2 * 200
    N_phi_calc = 2 * 360

    print("Int_omega_BdG_table ", Int_omega_BdG_table)
    print("np.shape(Int_omega_BdG_table) ", np.shape(Int_omega_BdG_table))

    DOS_L, DOS_U, gamma_L_Hz_array, gamma_U_Hz_array = gamma_values_UL(
        eta_small,
        NQ_calc,
        N_phi_calc,
        H0,
        Int_omega_BdG_table,
        Int_coupling_plus_table,
        Int_coupling_minus_table,
        len(h_NV_array),
    )

    # print(f"Γ(ω=ωL)={gamma_L_Hz_array} Hz")
    # print(f"Γ(ω=ωU)={gamma_U_Hz_array} Hz")

    return DOS_L, DOS_U, gamma_L_Hz_array, gamma_U_Hz_array


########################################################################################
########################################################################################


# Example usage
if __name__ == "__main__":

    start_time = time.time()

    #  For plotting multiple H0 values (uncomment to use)
    H0_array = np.array([82])
    # H0_array = np.array([75, 78, 80, 81, 81.5, 82, 82.5, 83, 84, 86, 90, 95])
    DOSL_list = np.zeros(len(H0_array))
    DOSU_list = np.zeros(len(H0_array))
    GammaL_list = np.zeros(len(H0_array))
    GammaU_list = np.zeros(len(H0_array))

    for h in range(len(H0_array)):

        H0_test = H0_array[h]
        print(f"Calculating for H0 = {H0_test} Oe")
        DOSL_list[h], DOSU_list[h], GammaL_list[h], GammaU_list[h] = gamma_1_from_H0(
            H0_test
        )

    print("GammaL ", GammaL_list)
    print("GammaU ", GammaU_list)
    print("Total time taken after vectorization ", time.time() - start_time, "s")

    plt.plot(H0_array, GammaL_list, label="GammaL")
    plt.plot(H0_array, GammaU_list, label="GammaU")
    plt.xlabel("B field (G)")
    plt.ylabel("Gamma (Hz)")
    plt.legend()
    plt.show()


########################################################################################
"""
"""
