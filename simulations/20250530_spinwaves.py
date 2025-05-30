import numpy as np
import scipy.linalg as la
from scipy.interpolate import RegularGridInterpolator, griddata
from scipy.optimize import fsolve
from multiprocessing import Pool, cpu_count
import matplotlib.pyplot as plt
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

class SpinWaveCalculator:
    def __init__(self):
        # Physical constants
        self.h_plank = 6.626e-34  # J*s
        self.mu_0 = 4 * np.pi * 1e-7  # H/m
        self.k_B = 1.381e-23  # J/K
        self.temperature = 300  # K
        
        # Material parameters
        self.delta_H = 1
        self.M0 = self._calculate_Ms()  # Oe
        self.L = 3  # μm
        self.DD = 5.4e-9 * 1e8  # Oe μm^-2
        self.gamma = 2.8e-3  # GHz/G
        self.omega_M = self.gamma * self.M0  # GHz
        
        # Derived parameters
        self.omega_d = (self.h_plank * self.mu_0 * self.gamma * 1e9 / 2 * 
                       (self.L * 1e-6)**3 * 1e8)  # Hz
        
    def _calculate_Ms(self):
        """Calculate saturation magnetization"""
        return ((2870 / 2.8) - 82 / 2) / 2 * 82 + self.delta_H - 4 + 2
    
    def N_bose(self, omega):
        """Bose-Einstein distribution"""
        return 1e-9 * self.k_B * self.temperature / self.h_plank * omega
    
    def paraunitary_diag(self, H):
        """Paraunitary diagonalization of Hamiltonian matrix"""
        try:
            # Cholesky decomposition
            K = la.cholesky(H, lower=True)
            dim = H.shape[0]
            
            # Create sigma_3 matrix
            sigma_3 = np.diag([(-1)**np.floor((2*n-1)/dim) for n in range(1, dim+1)])
            
            # Compute W matrix
            W = K @ sigma_3 @ K.conj().T
            
            # Eigendecomposition
            eigenvals, eigenvecs = la.eig(W)
            
            # Normalize eigenvectors
            eigenvecs = eigenvecs / np.linalg.norm(eigenvecs, axis=0)
            
            # Create permutation for ordering
            preperm = np.concatenate([
                np.arange(dim//2, dim),
                np.arange(dim//2-1, -1, -1)
            ])
            ordering = np.argsort(eigenvals)
            permutation = preperm[ordering]
            
            # Apply permutation
            eigenvals = eigenvals[permutation]
            eigenvecs = eigenvecs[:, permutation]
            
            U = eigenvecs.T
            H_diag = sigma_3 @ np.diag(eigenvals)
            
            # Compute transformation matrix T
            T = la.inv(K) @ U @ la.sqrtm(H_diag)
            
            # Phase correction
            T_pp = T[:dim//2, :dim//2]
            T_nn = T[dim//2:, dim//2:]
            
            phase_array_p = np.exp(1j * np.angle(np.diag(T_pp)))
            phase_array_n = np.exp(1j * np.angle(np.diag(T_nn)))
            
            V = np.diag(np.concatenate([
                phase_array_p.conj(),
                phase_array_n.conj()
            ]))
            
            T = T @ V
            
            # Extract submatrices
            T_pp = T[:dim//2, :dim//2]
            T_np = T[dim//2:, :dim//2]
            T_pn = T[:dim//2, dim//2:]
            T_nn = T[dim//2:, dim//2:]
            
            return eigenvals[:dim//2], T_pp, T_np, T_pn, T_nn, T
            
        except Exception as e:
            print(f"Error in paraunitary diagonalization: {e}")
            return None
    
    def F_function(self, q, n):
        """F function for Hamiltonian construction"""
        if q == 0:
            return 2 if n == 0 else 0
        return 2 * (1 - (-1)**n * np.exp(-q)) / q
    
    def P_function(self, q, n, m):
        """P function for Hamiltonian construction"""
        delta_nm = 1 if n == m else 0
        delta_n0 = 1 if n == 0 else 0
        delta_m0 = 1 if m == 0 else 0
        
        term1 = (q**2 / (q**2 + n**2 * np.pi**2)) * delta_nm
        
        denominator = (1 + delta_n0) * (1 + delta_m0)
        if denominator == 0:
            return 0
            
        term2 = (1/denominator * 
                q**4 / ((q**2 + n**2 * np.pi**2) * (q**2 + m**2 * np.pi**2)) * 
                self.F_function(q, n) * (1 + (-1)**(n+m)) / 2)
        
        return term1 - term2
    
    def Q_function(self, q, n, m):
        """Q function for Hamiltonian construction"""
        if m == 0:
            return 0
            
        delta_n0 = 1 if n == 0 else 0
        delta_m0 = 1 if m == 0 else 0
        
        term1 = (q**2 / (q**2 + m**2 * np.pi**2)) * (m**2 / (m**2 - n**2 + 
                (1 + (-1)**(n+m))/2))
        
        term2 = q - (q**2 / (2 * (q**2 + n**2 * np.pi**2))) * self.F_function(q, n)
        
        denominator = (1 + delta_n0) * (1 + delta_m0)
        if denominator == 0:
            return 0
            
        term3 = (1/denominator * (1 - (-1)**(n+m))/2)
        
        return term1 * term2 * term3
    
    def omega_function(self, omega_H, q, n):
        """Frequency function"""
        return omega_H + (self.gamma * self.DD / self.L**2) * (q**2 + n**2 * np.pi**2) * self.omega_M
    
    def generate_hamiltonian_bdg(self, omega_H, q, phi_k, N_max):
        """Generate Bogoliubov-de Gennes Hamiltonian"""
        H_bdg = np.zeros((2 * N_max, 2 * N_max), dtype=complex)
        
        for m in range(N_max):
            for n in range(N_max):
                # Diagonal terms
                omega_val = self.omega_function(omega_H, q, n)
                identity = np.eye(2)
                ones = np.ones((2, 2))
                
                # Construct H matrix elements
                term1 = omega_val * identity
                term2 = 0.5 * ones * (1 if n == m else 0)
                
                pauli_like = np.array([
                    [1 - np.sin(phi_k)**2, 1 + np.sin(phi_k)**2],
                    [1 + np.sin(phi_k)**2, 1 - np.sin(phi_k)**2]
                ])
                term3 = -0.5 * pauli_like * self.P_function(q, n, m)
                
                off_diag = np.array([[0, -4], [4, 0]])
                term4 = -0.5 * off_diag * np.sin(phi_k) * self.Q_function(q, n, m)
                
                H_element = term1 + term2 + term3 + term4
                
                # Fill BdG matrix
                H_bdg[n, m] = H_element[0, 0]
                H_bdg[n, m + N_max] = H_element[0, 1]
                H_bdg[n + N_max, m] = H_element[1, 0]
                H_bdg[n + N_max, m + N_max] = H_element[1, 1]
                
        return H_bdg
    
    def f_function(self, q, n, h_over_L):
        """f function for coupling calculations"""
        delta_n0 = 1 if n == 0 else 0
        if q == 0:
            return 0
        return ((-1)**n / (2 * (1 + delta_n0)) * 
                (q**2 / (q**2 + n**2 * np.pi**2)) * 
                np.exp(-q * h_over_L) * 
                (1 - (-1)**n * np.exp(-q)))
    
    def calculate_single_point(self, args):
        """Calculate paraunitary diagonalization for a single (q, phi) point"""
        omega_H, q, phi_k, N_max, h_NV_array = args
        
        try:
            # Generate Hamiltonian
            H_bdg = self.generate_hamiltonian_bdg(omega_H, q, phi_k, N_max)
            
            # Make Hermitian
            H_hermitian = (H_bdg + H_bdg.conj().T) / 2
            
            # Paraunitary diagonalization
            result = self.paraunitary_diag(H_hermitian)
            if result is None:
                return None
                
            eigenvals, T_pp, T_np, T_pn, T_nn, T = result
            
            # Coupling calculations
            gamma_plus = ((1 + np.sin(phi_k))/2 * 
                         (T_pp + T_np + np.sin(phi_k) * (T_pp - T_np)))
            gamma_minus = ((1 - np.sin(phi_k))/2 * 
                          (T_pp + T_np + np.sin(phi_k) * (T_pp - T_np)))
            gamma_z = (-1j * np.cos(phi_k)/2 * 
                      (T_pp + T_np + np.sin(phi_k) * (T_pp - T_np)))
            
            # Mirror terms for phi + pi
            phi_mirror = phi_k + np.pi
            gamma_plus_mirror = ((1 + np.sin(phi_mirror))/2 * 
                                (T_nn + T_pn + np.sin(phi_mirror) * (T_nn - T_pn)).conj())
            gamma_minus_mirror = ((1 - np.sin(phi_mirror))/2 * 
                                 (T_nn + T_pn + np.sin(phi_mirror) * (T_nn - T_pn)).conj())
            gamma_z_mirror = (-1j * np.cos(phi_mirror)/2 * 
                             (T_nn + T_pn + np.sin(phi_mirror) * (T_nn - T_pn)).conj())
            
            # Calculate f_bar arrays for different NV heights
            results = []
            for h_NV in h_NV_array:
                f_bar = np.array([self.f_function(q, n, h_NV / self.L) 
                                 for n in range(N_max)])
                
                v_plus = f_bar @ gamma_plus
                v_plus_mirror = f_bar @ gamma_plus_mirror
                v_minus = f_bar @ gamma_minus
                v_minus_mirror = f_bar @ gamma_minus_mirror
                v_z = f_bar @ gamma_z
                v_z_mirror = f_bar @ gamma_z_mirror
                
                results.append({
                    'eigenvals': eigenvals,
                    'v_plus': v_plus,
                    'v_plus_mirror': v_plus_mirror,
                    'v_minus': v_minus,
                    'v_minus_mirror': v_minus_mirror,
                    'v_z': v_z,
                    'v_z_mirror': v_z_mirror
                })
            
            return results
            
        except Exception as e:
            print(f"Error in calculation: {e}")
            return None
    
    def multi_para_diag(self, h_NV_array, omega_H, q_table, phi_table, N_max, n_processes=None):
        """Parallelized paraunitary diagonalization over q and phi grid"""
        if n_processes is None:
            n_processes = min(cpu_count(), 8)  # Limit to 8 processes
        
        print(f"Using {n_processes} processes for parallelization")
        
        num_q = len(q_table)
        num_phi = len(phi_table)
        
        # Prepare arguments for parallel processing
        args_list = []
        for i, q in enumerate(q_table):
            for j, phi in enumerate(phi_table):
                args_list.append((omega_H, q, phi, N_max, h_NV_array))
        
        # Parallel computation
        print("Starting parallel computation...")
        with Pool(n_processes) as pool:
            results = list(tqdm(
                pool.imap(self.calculate_single_point, args_list),
                total=len(args_list),
                desc="Computing paraunitary diagonalization"
            ))
        
        # Organize results
        organized_results = self._organize_results(results, num_q, num_phi, N_max, len(h_NV_array))
        
        return organized_results
    
    def _organize_results(self, results, num_q, num_phi, N_max, num_h_NV):
        """Organize parallel computation results into structured arrays"""
        # Initialize arrays
        omega_bdg_table = np.zeros((num_q, 2 * num_phi, N_max))
        coupling_plus_table = np.zeros((num_q, 2 * num_phi, num_h_NV, N_max), dtype=complex)
        coupling_minus_table = np.zeros((num_q, 2 * num_phi, num_h_NV, N_max), dtype=complex)
        coupling_z_table = np.zeros((num_q, 2 * num_phi, num_h_NV, N_max), dtype=complex)
        
        # Fill arrays
        idx = 0
        for i in range(num_q):
            for j in range(num_phi):
                if results[idx] is not None:
                    for h_idx, result in enumerate(results[idx]):
                        # Regular phi
                        omega_bdg_table[i, j, :] = result['eigenvals']
                        coupling_plus_table[i, j, h_idx, :] = result['v_plus']
                        coupling_minus_table[i, j, h_idx, :] = result['v_minus']
                        coupling_z_table[i, j, h_idx, :] = result['v_z']
                        
                        # Mirror phi (phi + pi)
                        omega_bdg_table[i, j + num_phi, :] = result['eigenvals']
                        coupling_plus_table[i, j + num_phi, h_idx, :] = result['v_plus_mirror']
                        coupling_minus_table[i, j + num_phi, h_idx, :] = result['v_minus_mirror']
                        coupling_z_table[i, j + num_phi, h_idx, :] = result['v_z_mirror']
                
                idx += 1
        
        return {
            'omega_bdg_table': omega_bdg_table,
            'coupling_plus_table': coupling_plus_table,
            'coupling_minus_table': coupling_minus_table,
            'coupling_z_table': coupling_z_table
        }
    
    def create_interpolators(self, q_table, phi_table, results):
        """Create interpolation functions for the computed data"""
        num_q, num_phi_extended, N_max = results['omega_bdg_table'].shape
        num_h_NV = results['coupling_plus_table'].shape[2]
        
        # Create coordinate grids
        q_grid = np.tile(q_table[:, np.newaxis], (1, num_phi_extended))
        phi_extended = np.concatenate([phi_table, phi_table + np.pi])
        phi_grid = np.tile(phi_extended[np.newaxis, :], (len(q_table), 1))
        
        # Create interpolators for each mode and height
        interpolators = {
            'omega': [],
            'coupling_plus': [],
            'coupling_minus': [],
            'coupling_z': []
        }
        
        for s in range(N_max):
            # Omega interpolator
            points = np.column_stack((q_grid.ravel(), phi_grid.ravel()))
            values = results['omega_bdg_table'][:, :, s].ravel()
            
            interpolators['omega'].append(
                lambda q, phi, vals=values, pts=points: 
                griddata(pts, vals, (q, phi), method='linear', fill_value=0)
            )
            
            # Coupling interpolators for each height
            for h_idx in range(num_h_NV):
                for coupling_type in ['plus', 'minus', 'z']:
                    if len(interpolators[f'coupling_{coupling_type}']) <= h_idx:
                        interpolators[f'coupling_{coupling_type}'].append([])
                    
                    values = results[f'coupling_{coupling_type}_table'][:, :, h_idx, s].ravel()
                    
                    interpolators[f'coupling_{coupling_type}'][h_idx].append(
                        lambda q, phi, vals=values, pts=points:
                        griddata(pts, vals, (q, phi), method='linear', fill_value=0+0j)
                    )
        
        return interpolators
    
    def calculate_gamma_values(self, eta_small, N_Q_est, N_phi, H0, interpolators, num_h_NV, N_max):
        """Calculate dissipation rates Gamma_L and Gamma_U"""
        D_NV = 2.87  # GHz
        omega_target_L = D_NV - H0 * self.gamma
        omega_target_U = D_NV + H0 * self.gamma
        
        # Create Q table with better sampling
        Q_middle = 5 * self.L
        N_Q_half = N_Q_est // 2
        q_table_fine = np.logspace(-6, np.log10(Q_middle), N_Q_half)
        q_table_fine = np.concatenate([q_table_fine, 
                                      np.arange(Q_middle, 50*self.L, 
                                               (50*self.L - Q_middle)/N_Q_half)])
        
        delta_q = np.diff(q_table_fine)
        delta_phi = 2 * np.pi / N_phi
        
        # Initialize results
        gamma_L_Hz_array = np.zeros(num_h_NV)
        gamma_U_Hz_array = np.zeros(num_h_NV)
        DOS_L = 0
        DOS_U = 0
        
        # Integration over modes
        for s in range(N_max):
            print(f"Processing mode {s+1}/{N_max}")
            
            for i in range(len(delta_q)):
                for j in range(N_phi):
                    q = q_table_fine[i]
                    phi = 2 * np.pi * j / N_phi
                    
                    # Get interpolated values
                    omega_val = interpolators['omega'][s](q, phi)
                    if np.isnan(omega_val) or omega_val <= 0:
                        continue
                    
                    # Calculate Lorentzian factors
                    lorentz_L = (eta_small / np.pi) / (eta_small**2 + 
                                (self.omega_M * omega_val - omega_target_L)**2)
                    lorentz_U = (eta_small / np.pi) / (eta_small**2 + 
                                (self.omega_M * omega_val - omega_target_U)**2)
                    
                    # Density of states contribution
                    dos_factor = (1 / self.L**2) * delta_phi / (2*np.pi)**2 * delta_q[i] * q
                    DOS_L += dos_factor * lorentz_L
                    DOS_U += dos_factor * lorentz_U
                    
                    # Dissipation rate calculation
                    bose_factor = 2 * self.N_bose(self.omega_M * omega_val) + 1
                    gamma_prefactor = ((2*np.pi)**2 * self.omega_M * self.omega_d * 
                                     delta_phi / (2*np.pi)**2 * bose_factor * 
                                     delta_q[i] * q)
                    
                    # For each NV height
                    for h_idx in range(num_h_NV):
                        # Get coupling values
                        v_plus = interpolators['coupling_plus'][h_idx][s](q, phi)
                        v_minus = interpolators['coupling_minus'][h_idx][s](q, phi)
                        
                        if not (np.isnan(v_plus) or np.isnan(v_minus)):
                            gamma_L_Hz_array[h_idx] += (gamma_prefactor * 
                                                       np.abs(v_plus)**2 * lorentz_L)
                            gamma_U_Hz_array[h_idx] += (gamma_prefactor * 
                                                       np.abs(v_minus)**2 * lorentz_U)
        
        return DOS_L, DOS_U, gamma_L_Hz_array, gamma_U_Hz_array

def main():
    """Main execution function"""
    # Initialize calculator
    calc = SpinWaveCalculator()
    
    # Parameters for calculation
    h_NV_array = np.array([0.4, 0.5, 0.6, 0.7])  # NV heights in μm
    H0_array = np.array([76, 77, 78, 79, 80, 81, 81.5, 82, 82.5, 83, 
                        84, 85, 86, 87, 88, 89, 90, 91, 92, 95])  # Magnetic fields in Oe
    
    # Grid parameters
    num_phi = 180  # Reduced for faster computation
    num_Q = 200    # Reduced for faster computation
    F_max = 5  # GHz
    N_max = int(np.ceil(calc.L / np.pi * F_max / (calc.gamma * calc.DD)))
    
    print(f"N_max: {N_max}")
    
    # Create grids
    del_phi = np.pi / num_phi
    phi_table = np.arange(0, np.pi, del_phi)
    
    # Logarithmic spacing for q
    Q_max = 50 * calc.L
    q_table = np.logspace(-6, np.log10(Q_max), num_Q)
    
    # Storage arrays
    DOS_L_array = np.zeros(len(H0_array))
    DOS_U_array = np.zeros(len(H0_array))
    gamma_L_Hz_array = np.zeros((len(h_NV_array), len(H0_array)))
    gamma_U_Hz_array = np.zeros((len(h_NV_array), len(H0_array)))
    
    # Calculate for each magnetic field
    for i, H0 in enumerate(tqdm(H0_array, desc="Processing magnetic fields")):
        print(f"\nProcessing H0 = {H0} Oe ({i+1}/{len(H0_array)})")
        
        omega_H = calc.gamma * H0
        
        # Parallel computation
        results = calc.multi_para_diag(h_NV_array, omega_H, q_table, phi_table, N_max)
        
        # Create interpolators
        interpolators = calc.create_interpolators(q_table, phi_table, results)
        
        # Calculate dissipation rates
        eta_small = 0.003
        N_Q_est = 400
        N_phi_fine = 720
        
        DOS_L, DOS_U, gamma_L_Hz, gamma_U_Hz = calc.calculate_gamma_values(
            eta_small, N_Q_est, N_phi_fine, H0, interpolators, len(h_NV_array), N_max
        )
        
        # Store results
        DOS_L_array[i] = DOS_L
        DOS_U_array[i] = DOS_U
        gamma_L_Hz_array[:, i] = gamma_L_Hz
        gamma_U_Hz_array[:, i] = gamma_U_Hz
        
        print(f"DOS_L: {DOS_L:.2e}, DOS_U: {DOS_U:.2e}")
        print(f"Gamma_L: {gamma_L_Hz}")
        print(f"Gamma_U: {gamma_U_Hz}")
    
    # Plotting results
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
    
    # Plot dissipation rates for first NV height
    show_h_NV_index = 0
    ax1.plot(H0_array, gamma_L_Hz_array[show_h_NV_index, :], 'o-', label='Γ_L')
    ax1.plot(H0_array, gamma_U_Hz_array[show_h_NV_index, :] * 1000, 's-', label='Γ_U × 1000')
    ax1.set_xlabel('Field (G)')
    ax1.set_ylabel('Γ_L, Γ_U × 1000 (Hz)')
    ax1.legend()
    ax1.grid(True)
    
    # Plot density of states
    ax2.plot(H0_array, DOS_L_array / calc.L, 'o-', label='DOS_L')
    ax2.plot(H0_array, DOS_U_array / calc.L, 's-', label='DOS_U')
    ax2.set_xlabel('Field (G)')
    ax2.set_ylabel('2π*DOS (1/GHz μm³)')
    ax2.legend()
    ax2.grid(True)
    
    # Plot all NV heights for Gamma_L
    for i, h_NV in enumerate(h_NV_array):
        ax3.plot(H0_array, gamma_L_Hz_array[i, :], 'o-', label=f'h_NV = {h_NV} μm')
    ax3.set_xlabel('Field (G)')
    ax3.set_ylabel('Γ_L (Hz)')
    ax3.legend()
    ax3.grid(True)
    
    # Calculate and plot ratio
    ratio = np.zeros_like(H0_array)
    for i in range(len(H0_array)):
        if DOS_L_array[i] > 0:
            ratio[i] = (1e-3 * gamma_L_Hz_array[show_h_NV_index, i] * calc.L * 
                       1e-6 * DOS_L_array[i] * 2 * calc.N_bose(2.87 - calc.gamma * H0_array[i]))
    
    ax4.plot(H0_array, ratio, 'o-')
    ax4.set_xlabel('Field (G)')
    ax4.set_ylabel('Ratio (kHz² μm³)')
    ax4.grid(True)
    
    plt.tight_layout()
    plt.show()
    
    # Save results
    results_dict = {
        'H0_array': H0_array,
        'DOS_L_array': DOS_L_array,
        'DOS_U_array': DOS_U_array,
        'h_NV_array': h_NV_array,
        'gamma_L_Hz_array': gamma_L_Hz_array,
        'gamma_U_Hz_array': gamma_U_Hz_array
    }
    
    np.savez('spin_wave_results.npz', **results_dict)
    print("Results saved to 'spin_wave_results.npz'")
    
    return results_dict

if __name__ == "__main__":
    results = main()
