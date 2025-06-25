import numpy as np
import scipy.linalg as la
from scipy.interpolate import RectBivariateSpline
from tqdm.auto import tqdm  # For progress bars
import time
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# Paraunitary Diagonalization
# -----------------------------------------------------------------------------

def para_unitary_diag(H):
    """
    Performs paraunitary diagonalization of a Hamiltonian matrix H.
    """
    dim = H.shape[0]
    try:
        K = la.cholesky(H)
    except np.linalg.LinAlgError:
        raise np.linalg.LinAlgError(
            "Cholesky decomposition failed. The matrix is not positive definite."
        )

    sigma3_diag = [(-1) ** np.floor((2 * (n + 1) - 1) / dim) for n in range(dim)]
    sigma3 = np.diag(sigma3_diag)

    W = K @ sigma3 @ K.conj().T
    evals, evecs = la.eig(W)
    evals, evecs = np.real(evals), np.real_if_close(evecs)  # Eigenvalues of W are real
    evecs = evecs / np.linalg.norm(evecs, axis=0)

    half_dim = dim // 2
    preperm_indices = np.concatenate(
        [np.arange(half_dim, dim), np.arange(half_dim - 1, -1, -1)]
    )

    ordering_indices = np.argsort(evals)
    permutation = preperm_indices[ordering_indices].astype(int)

    evals = evals[permutation]
    evecs = evecs[:, permutation]
    U = evecs

    Hdiag = sigma3 @ np.diag(evals)
    T = la.inv(K) @ U @ la.sqrtm(Hdiag)

    Tpp = T[0:half_dim, 0:half_dim]
    Tnn = T[half_dim:dim, half_dim:dim]

    phase_array_p = np.exp(1j * np.angle(np.diag(Tpp)))
    phase_array_n = np.exp(1j * np.angle(np.diag(Tnn)))

    v_diag = np.concatenate([np.conj(phase_array_p), np.conj(phase_array_n)])
    V = np.diag(v_diag)

    T = T @ V

    Tpp = T[0:half_dim, 0:half_dim]
    Tnp = T[half_dim:dim, 0:half_dim]
    Tpn = T[0:half_dim, half_dim:dim]
    Tnn = T[half_dim:dim, half_dim:dim]

    return evals[0:half_dim], Tpp, Tnp, Tpn, Tnn, T


def calculate_gamma_and_dos(
    eta_small,
    nq_est,
    n_phi,
    H0,
    interp_omega,
    interp_cplus,
    interp_cminus,
    num_hNV,
    Nmax,
):
    """
    Calculates relaxation rates (Gamma) and Density of States (DOS)
    using the previously generated interpolation functions.
    """
    # Constants
    DNV = 2.87  # GHz
    omega_target_L = DNV - H0 * gamma
    omega_target_U = DNV + H0 * gamma

    # Setup integration grid for q
    q_middle = 5 * L
    nq_half = round(nq_est / 2)
    q_table_part1 = log_space(L / 1000, q_middle, nq_half)
    last_delta_q = np.diff(q_table_part1)[-1]
    q_table_part2 = np.arange(q_table_part1[-1] + last_delta_q, Qmax, last_delta_q)
    q_table = np.concatenate(([1e-6], q_table_part1, q_table_part2))

    delta_q = np.diff(q_table)

    # Setup integration grid for phi
    delta_phi = 2 * np.pi / n_phi
    phi_table_integration = np.arange(0, 2 * np.pi, delta_phi)

    # Initialize result arrays
    gamma_L_hz_array = np.zeros(num_hNV)
    gamma_U_hz_array = np.zeros(num_hNV)
    dos_l = 0.0
    dos_u = 0.0

    # Main loop over each magnon mode 's'
    pbar = tqdm(range(Nmax), desc="Calculating Gamma/DOS")
    for s in pbar:
        # Create 2D grids for q and phi for vectorized calculation
        q_grid, phi_grid = np.meshgrid(q_table[:-1], phi_table_integration)
        delta_q_grid, _ = np.meshgrid(delta_q, phi_table_integration)

        # Evaluate interpolators on the grid
        omega_bdg_grid = interp_omega[s](
            q_grid.flatten(), phi_grid.flatten(), grid=False
        ).reshape(q_grid.shape)

        # --- Lorentzian term for DOS and Gamma ---
        lorentzian_L = (eta_small / np.pi) / (
            eta_small**2 + (omega_M * omega_bdg_grid - omega_target_L) ** 2
        )
        lorentzian_U = (eta_small / np.pi) / (
            eta_small**2 + (omega_M * omega_bdg_grid - omega_target_U) ** 2
        )

        # --- Summation for DOS ---
        # This is the integral of the Lorentzian over k-space
        dos_integrand = delta_q_grid * q_grid * lorentzian_L
        dos_l += (1 / L**2) * (delta_phi / (2 * np.pi) ** 2) * np.sum(dos_integrand)

        dos_integrand_U = delta_q_grid * q_grid * lorentzian_U
        dos_u += (1 / L**2) * (delta_phi / (2 * np.pi) ** 2) * np.sum(dos_integrand_U)

        # --- Summation for Gamma ---
        # Bose-Einstein distribution (high-temperature approximation from script)
        n_bose = NBose(omega_M * omega_bdg_grid)

        for tt in range(num_hNV):
            cplus_grid_real = interp_cplus[tt][s](
                q_grid.flatten(), phi_grid.flatten(), grid=False
            ).reshape(q_grid.shape)
            cminus_grid_real = interp_cminus[tt][s](
                q_grid.flatten(), phi_grid.flatten(), grid=False
            ).reshape(q_grid.shape)
            # Assuming complex part was negligible or handled inside interpolator if needed

            # Gamma for lower frequency
            gamma_L_integrand = (
                (2 * n_bose + 1)
                * delta_q_grid
                * q_grid
                * (np.abs(cplus_grid_real)) ** 2
                * lorentzian_L
            )
            gamma_L_hz_array[tt] += np.sum(gamma_L_integrand)

            # Gamma for upper frequency
            gamma_U_integrand = (
                (2 * n_bose + 1)
                * delta_q_grid
                * q_grid
                * (np.abs(cminus_grid_real)) ** 2
                * lorentzian_U
            )
            gamma_U_hz_array[tt] += np.sum(gamma_U_integrand)

    # Final scaling factor
    gamma_factor = omega_M * omega_d * delta_phi
    gamma_L_hz_array *= gamma_factor
    gamma_U_hz_array *= gamma_factor

    return dos_l, dos_u, gamma_L_hz_array, gamma_U_hz_array


def gamma1_from_h0(H0):
    """
    Top-level function to run the full calculation for a given magnetic field H0.
    """
    print(f"\n--- Starting calculation for H0 = {H0} Oe ---")
    omega_H = gamma * H0  # GHz

    # Step 1: Generate interpolation tables
    _, _, _, interp_omega, interp_cplus, interp_cminus, _ = multi_para_diag(
        hNVarray, omega_H, qtable_setup, phitable_setup, Nmax
    )

    # Step 2: Use tables to calculate DOS and Gamma
    # Parameters from original script
    eta_small = 0.003
    NQ_integration = 2 * 200
    NPhi_integration = 2 * 360

    dos_l, dos_u, gamma_l, gamma_u = calculate_gamma_and_dos(
        eta_small,
        NQ_integration,
        NPhi_integration,
        H0,
        interp_omega,
        interp_cplus,
        interp_cminus,
        len(hNVarray),
        Nmax,
    )

    print("\n--- Final Results ---")
    print(f"DOS(ω=ωL) = {dos_l:.4g} 1/(GHz μm^2)")
    print(f"DOS(ω=ωU) = {dos_u:.4g} 1/(GHz μm^2)")
    print(f"Γ(ω=ωL) = {gamma_l} Hz")
    print(f"Γ(ω=ωU) = {gamma_u} Hz")

    return dos_l, dos_u, gamma_l, gamma_u


def F(q, n):
    if abs(q) < 1e-9:
        return 2.0 if n % 2 != 0 else 0.0
    return 2 * (1 - (-1) ** n * np.exp(-q)) / q


def P(q, n, m):
    term1 = (q**2 / (q**2 + n**2 * np.pi**2)) * (1 if n == m else 0)
    norm_factor = np.sqrt((1 + (1 if n == 0 else 0)) * (1 + (1 if m == 0 else 0)))
    f_term = F(q, n)
    parity_term = (1 + (-1) ** (n + m)) / 2
    term2 = (
        (1 / norm_factor)
        * (q**4 / ((q**2 + n**2 * np.pi**2) * (q**2 + m**2 * np.pi**2)))
        * f_term
        * parity_term
    )
    return term1 - term2


def Q(q, n, m):
    parity_term = (1 - (-1) ** (n + m)) / 2
    if parity_term == 0:
        return 0.0
    denom = m**2 - n**2 + parity_term
    if abs(denom) < 1e-9:
        return 0.0
    norm_factor = np.sqrt((1 + (1 if n == 0 else 0)) * (1 + (1 if m == 0 else 0)))
    term1 = q**2 / (q**2 + m**2 * np.pi**2)
    term2 = ((m**2 / denom) * (2 / q)) - ((q**2 / (2 * (q**2 + n**2 * np.pi**2))) * F(q, n))
    return term1 * term2 * (1 / norm_factor) * parity_term


def Omega_func(omega_H, q, n):
    return (omega_H + (((gamma * DD) / L**2) * (q**2 + n**2 * np.pi**2))) / omega_M


def H_func(omega_H, q, phi_k, n, m):
    omega_val = Omega_func(omega_H, q, n)
    kronecker = 1 if n == m else 0
    mat1 = np.array([[1, 0], [0, 1]])
    mat2 = np.array([[1, 1], [1, 1]])
    term1 = (omega_val * mat1 + 0.5 * mat2) * kronecker
    sin_phi_k_sq = np.sin(phi_k) ** 2
    mat3 = np.array(
        [[1 - sin_phi_k_sq, 1 + sin_phi_k_sq], [1 + sin_phi_k_sq, 1 - sin_phi_k_sq]]
    )
    term2 = -0.5 * mat3 * P(q, n, m)
    mat4 = np.array([[0, -4], [4, 0]], dtype=float)
    term3 = -0.5 * mat4 * np.sin(phi_k) * Q(q, n, m)
    return term1 + term2 + term3


def f_func(q, n, hoverL):
    norm = np.sqrt(2 * (1 + (1 if n == 0 else 0)))
    return (
        ((-1) ** n / norm)
        * (q**2 / (q**2 + n**2 * np.pi**2))
        * np.exp(-q * hoverL)
        * (1 - (-1) ** n * np.exp(-q))
    )


def generate_h_bdg(omega_H, q, phi_k, Nmax):
    H_bdg = np.zeros((2 * Nmax, 2 * Nmax), dtype=np.complex128)
    for m in range(1, Nmax + 1):
        for n in range(1, Nmax + 1):
            result = H_func(omega_H, q, phi_k, n - 1, m - 1)
            H_bdg[n - 1, m - 1] = result[0, 0]
            H_bdg[n - 1, m - 1 + Nmax] = result[0, 1]
            H_bdg[n - 1 + Nmax, m - 1] = result[1, 0]
            H_bdg[n - 1 + Nmax, m - 1 + Nmax] = result[1, 1]
    return H_bdg


def multi_para_diag(hNVarray, omega_H, qtable, phitable, Nmax):
    num_q = len(qtable)
    num_phi = len(phitable)
    num_hNV = len(hNVarray)
    phitable_full = np.concatenate([phitable, phitable + np.pi])
    num_phi_full = len(phitable_full)

    omega_bdg_table = np.zeros((num_q, num_phi_full, Nmax))
    coupling_plus_table = np.zeros(
        (num_q, num_phi_full, num_hNV, Nmax), dtype=np.complex128
    )
    coupling_minus_table = np.zeros(
        (num_q, num_phi_full, num_hNV, Nmax), dtype=np.complex128
    )
    coupling_z_table = np.zeros(
        (num_q, num_phi_full, num_hNV, Nmax), dtype=np.complex128
    )

    pbar = tqdm(total=num_q * num_phi, desc="MultiParaDiag")
    for iq, q in enumerate(qtable):
        for iphi, phi_k in enumerate(phitable):
            H_bdg = generate_h_bdg(omega_H, q, phi_k, Nmax)
            H_sym = (H_bdg + H_bdg.conj().T) / 2.0

            try:
                eigvals, Tpp, Tnp, Tpn, Tnn, _ = para_unitary_diag(H_sym)
            except np.linalg.LinAlgError:
                eigvals = np.full(Nmax, np.nan)
                Tpp = Tnp = Tpn = Tnn = np.full((Nmax, Nmax), np.nan)

            omega_bdg_table[iq, iphi, :] = eigvals
            omega_bdg_table[iq, iphi + num_phi, :] = eigvals

            phi_k_plus_pi = phi_k + np.pi
            gamma_plus = (
                ((1 + np.sin(phi_k)) / 2) * (Tpp + Tnp + np.sin(phi_k) * (Tpp - Tnp))
            )
            gamma_plus_mirror = (
                ((1 + np.sin(phi_k_plus_pi))
                / 2)
                * np.conj(Tnn + Tpn + np.sin(phi_k_plus_pi) * (Tnn - Tpn))
            )
            gamma_minus = (
                ((1 - np.sin(phi_k)) / 2) * (Tpp + Tnp + np.sin(phi_k) * (Tpp - Tnp))
            )
            gamma_minus_mirror = (
                ((1 - np.sin(phi_k_plus_pi))
                / 2)
                * np.conj(Tnn + Tpn + np.sin(phi_k_plus_pi) * (Tnn - Tpn))
            )
            gamma_z = (
                -1j * (np.cos(phi_k) / 2) * (Tpp + Tnp + np.sin(phi_k) * (Tpp - Tnp))
            )
            gamma_z_mirror = (
                -1j
                * (np.cos(phi_k_plus_pi)
                / 2)
                * np.conj(Tnn + Tpn + np.sin(phi_k_plus_pi) * (Tnn - Tpn))
            )

            for ih, h_nv in enumerate(hNVarray):
                fbar = np.array([f_func(q, n, h_nv / L) for n in range(Nmax)])
                coupling_plus_table[iq, iphi, ih, :] = fbar @ gamma_plus
                coupling_plus_table[iq, iphi + num_phi, ih, :] = (
                    fbar @ gamma_plus_mirror
                )
                coupling_minus_table[iq, iphi, ih, :] = fbar @ gamma_minus
                coupling_minus_table[iq, iphi + num_phi, ih, :] = (
                    fbar @ gamma_minus_mirror
                )
                coupling_z_table[iq, iphi, ih, :] = fbar @ gamma_z
                coupling_z_table[iq, iphi + num_phi, ih, :] = fbar @ gamma_z_mirror
            pbar.update(1)
    pbar.close()

    int_omega_bdg, int_coupling_plus, int_coupling_minus, int_coupling_z = (
        [],
        [[] for _ in range(num_hNV)],
        [[] for _ in range(num_hNV)],
        [[] for _ in range(num_hNV)],
    )
    phitable_interp = np.append(phitable_full, 2 * np.pi)

    for s in range(Nmax):
        omega_s = np.hstack(
            [omega_bdg_table[:, :, s], omega_bdg_table[:, 0, s][:, np.newaxis]]
        )
        int_omega_bdg.append(
            RectBivariateSpline(
                qtable, phitable_interp, np.nan_to_num(omega_s), kx=2, ky=2
            )
        )

        for tt in range(num_hNV):
            cplus_s = np.hstack(
                [
                    coupling_plus_table[:, :, tt, s],
                    coupling_plus_table[:, 0, tt, s][:, np.newaxis],
                ]
            )
            int_coupling_plus[tt].append(
                RectBivariateSpline(
                    qtable, phitable_interp, np.nan_to_num(cplus_s.real), kx=2, ky=2
                )
            )
            cminus_s = np.hstack(
                [
                    coupling_minus_table[:, :, tt, s],
                    coupling_minus_table[:, 0, tt, s][:, np.newaxis],
                ]
            )
            int_coupling_minus[tt].append(
                RectBivariateSpline(
                    qtable, phitable_interp, np.nan_to_num(cminus_s.real), kx=2, ky=2
                )
            )
            cz_s = np.hstack(
                [
                    coupling_z_table[:, :, tt, s],
                    coupling_z_table[:, 0, tt, s][:, np.newaxis],
                ]
            )
            int_coupling_z[tt].append(
                RectBivariateSpline(
                    qtable, phitable_interp, np.nan_to_num(cz_s.real), kx=2, ky=2
                )
            )

    return (
        qtable,
        np.nanmin(omega_bdg_table),
        phitable_full,
        int_omega_bdg,
        int_coupling_plus,
        int_coupling_minus,
        int_coupling_z,
    )


if __name__ == "__main__":
    MsSol = 1716
    L = 3.0
    DD = 5.4e-9 * 1e8
    gamma = 2.8e-3
    omega_M = gamma * MsSol

    # Calculation Grid Parameters
    NumPhi_setup = 2 * 90
    NumQ_setup = 2 * 100  # Reduced for faster example run
    del_phi = np.pi / NumPhi_setup
    phitable_setup = np.arange(0, np.pi, del_phi)

    Qmax = 50 * L

    def log_space(min_val, max_val, steps):
        return np.logspace(np.log10(min_val), np.log10(max_val), steps)

    qtable_setup = np.concatenate(([1e-6], log_space(L / 1000, Qmax, NumQ_setup)))

    Fmax = 5
    Nmax = int(np.ceil((L / np.pi) * np.sqrt(Fmax / (gamma * DD))))
    print(f"Using Nmax = {Nmax}")

    # hNVarray = np.array([0.4, 0.5, 0.6, 0.7])  # pos of NV in um
    hNVarray = np.array([0.4])  # pos of NV in um

    # Temperature and Bose-Einstein Factor
    hPlank = 6.626e-34
    mu0 = 4 * np.pi * 1e-7
    kB = 1.381e-23
    omega_d = ((hPlank * mu0 * (gamma * 1e9) ** 2) / (L * 1e-6) ** 3) * 1e8  # Hz
    Temperature = 300  # K

    def NBose(omega_ghz):
        # Using high-temperature approximation from Mathematica script
        return (1e-9 * kB * Temperature / hPlank) / omega_ghz

    # --- Run the main calculation for a given H0 ---

    H0_field = np.array([81.5, 82, 82.5])  # Oe

    # Storage arrays
    DOS_L_array = np.zeros(len(H0_field))
    DOS_U_array = np.zeros(len(H0_field))
    gamma_L_Hz_array = np.zeros((len(H0_field)))
    gamma_U_Hz_array = np.zeros((len(H0_field)))

    for i, h in enumerate(H0_field):

        DOS_L, DOS_U, gamma_L_Hz, gamma_U_Hz  = gamma1_from_h0(h)

        DOS_L_array[i] = DOS_L
        DOS_U_array[i] = DOS_U
        gamma_L_Hz_array[i] = gamma_L_Hz
        gamma_U_Hz_array[i] = gamma_U_Hz

    # --- Plot the results ---
    # Plotting results
    fig, ((ax1, ax2)) = plt.subplots(1, 2)
    
    # Plot dissipation rates for first NV height
    ax1.plot(H0_field, gamma_L_Hz_array, 'o-', label='Γ_L')
    ax1.plot(H0_field, gamma_U_Hz_array * 1000, 's-', label='Γ_U × 1000')
    ax1.set_xlabel('Field (G)')
    ax1.set_ylabel('Γ_L, Γ_U × 1000 (Hz)')
    ax1.legend()
    ax1.grid(True)
    
    # Plot density of states
    ax2.plot(H0_field, DOS_L_array / L, 'o-', label='DOS_L')
    ax2.plot(H0_field, DOS_U_array / L, 's-', label='DOS_U')
    ax2.set_xlabel('Field (G)')
    ax2.set_ylabel('2π*DOS (1/GHz μm³)')
    ax2.legend()
    ax2.grid(True)
    
    # # Plot all NV heights for Gamma_L
    # for i, h_NV in enumerate(hNVarray):
    #     ax3.plot(H0_field, gamma_L_Hz_array[i], 'o-', label=f'h_NV = {h_NV} μm')
    # ax3.set_xlabel('Field (G)')
    # ax3.set_ylabel('Γ_L (Hz)')
    # ax3.legend()
    # ax3.grid(True)
    
    # # Calculate and plot ratio
    # ratio = np.zeros_like(H0_field)
    # for i in range(len(H0_field)):
    #     if DOS_L_array[i] > 0:
    #         ratio[i] = (1e-3 * gamma_L_Hz_array[show_h_NV_index, i] * calc.L * 
    #                    1e-6 * DOS_L_array[i] * 2 * calc.N_bose(2.87 - calc.gamma * H0_field[i]))
    
    # ax4.plot(H0_field, ratio, 'o-')
    # ax4.set_xlabel('Field (G)')
    # ax4.set_ylabel('Ratio (kHz² μm³)')
    # ax4.grid(True)
    
    timestr = time.strftime("%Y%m%d_%H%M%S")
    basename = str(__file__) + '_' + timestr 
    figname = basename + '.png'
    plt.savefig(figname)

    plt.tight_layout()
    plt.show()
