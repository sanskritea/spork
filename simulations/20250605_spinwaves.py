import numpy as np
import scipy.linalg as la
import scipy.interpolate as interp
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from tqdm import tqdm
import warnings
import time
warnings.filterwarnings('ignore')


###########################################################################
###########################################################################

##      FIXED MATHEMATICA-TO-PYTHON CONVERSION TYPOS (FROM 20250529 FILE)
##      FIXED POSITIVE-DEFINITE TERM ERRORS (20250605)
##      RESULT?

###########################################################################
###########################################################################


class SpinWaveCalculator:
    def __init__(self):
        # Physical constants
        self.h_plank = 6.626e-34  # J*s
        self.mu_0 = 4 * np.pi * 1e-7  # H/m
        self.k_B = 1.381e-23  # J/K
        self.temperature = 300  # K
        
        # Material parameters
        self.delta_H = 1
        self.Ms_sol = self.calculate_Ms_sol()
        self.M0 = self.Ms_sol  # Oe
        self.L = 3  # μm
        self.DD = 5.4e-9 * 1e8  # Oe μm^-2
        self.gamma = 2.8e-3  # GHz/G
        self.omega_M = self.gamma * self.M0  # GHz
        
        # Derived parameters
        self.omega_d = (self.h_plank * self.mu_0 * ((self.gamma * 1e9)**2) / 
                       ((self.L * 1e-6)**3)) * 1e8 # Hz
        
    def calculate_Ms_sol(self):
        """Calculate Ms solution from the given equation"""
        return 1716.0
    
    def para_unitary_diag(self, H):
        """Robust paraunitary diagonalization of Hamiltonian H"""
        dim = H.shape[0]
        
        # Check if matrix is positive definite
        try:
            # Try standard Cholesky decomposition first
            K = la.cholesky(H, lower=True)
        except la.LinAlgError:
            # If not positive definite, use regularization
            print(f"Matrix not positive definite, applying regularization...")
            
            # Method 1: Add small regularization to diagonal
            reg_param = 1e-12
            H_reg = H + reg_param * np.eye(dim)
            
            try:
                K = la.cholesky(H_reg, lower=True)
            except la.LinAlgError:
                # Method 2: Use eigenvalue regularization
                eigenvals, eigenvecs = la.eigh(H)
                
                # Set negative eigenvalues to small positive value
                eigenvals_reg = np.maximum(eigenvals, reg_param)
                H_reg = eigenvecs @ np.diag(eigenvals_reg) @ eigenvecs.T
                
                try:
                    K = la.cholesky(H_reg, lower=True)
                except la.LinAlgError:
                    # Method 3: Use SVD-based approach
                    U, s, Vh = la.svd(H)
                    s_reg = np.maximum(s, reg_param)
                    H_reg = U @ np.diag(s_reg) @ Vh
                    
                    # Force symmetry
                    H_reg = (H_reg + H_reg.T) / 2
                    
                    try:
                        K = la.cholesky(H_reg, lower=True)
                    except la.LinAlgError:
                        # Final fallback: use square root of regularized matrix
                        eigenvals, eigenvecs = la.eigh(H_reg)
                        eigenvals_sqrt = np.sqrt(np.maximum(eigenvals, reg_param))
                        K = eigenvecs @ np.diag(eigenvals_sqrt)
        
        # Create σ3 matrix
        sigma3 = np.diag([(-1)**(np.floor((2*n-1)/dim)) for n in range(1, dim+1)])
        
        # Calculate W matrix
        W = K @ sigma3 @ K.conj().T
        
        # Eigendecomposition with error handling
        try:
            evals, evec = la.eig(W)
        except la.LinAlgError:
            # Fallback to real symmetric eigendecomposition if possible
            if np.allclose(W, W.T.conj()):
                evals, evec = la.eigh(W)
            else:
                # Use SVD as ultimate fallback
                U, s, Vh = la.svd(W)
                evals = s
                evec = U
        
        # Handle complex eigenvalues
        if np.any(np.iscomplex(evals)):
            # Take real part if imaginary part is small
            if np.max(np.abs(np.imag(evals))) < 1e-12:
                evals = np.real(evals)
                evec = np.real(evec)
            else:
                print(f"Warning: Complex eigenvalues detected, max imag part: {np.max(np.abs(np.imag(evals)))}")
        
        # Normalize eigenvectors
        for i in range(evec.shape[1]):
            norm = np.linalg.norm(evec[:, i])
            if norm > 1e-15:
                evec[:, i] = evec[:, i] / norm
        
        # Create permutation for ordering
        preperm = np.concatenate([
            np.arange(dim//2, dim),
            np.arange(dim//2-1, -1, -1)
        ])
        ordering = np.argsort(np.real(evals))  # Sort by real part
        permutation = preperm[ordering]
        
        # Apply permutation
        evals = evals[permutation]
        evec = evec[:, permutation]
        
        U = evec.T
        H_diag = sigma3 @ np.diag(evals)
        
        # Compute T with robust matrix square root
        try:
            sqrt_H_diag = la.sqrtm(H_diag)
        except:
            # Fallback square root computation
            eigenvals_diag, eigenvecs_diag = la.eig(H_diag)
            sqrt_eigenvals = np.sqrt(eigenvals_diag + 0j)  # Handle negative values
            sqrt_H_diag = eigenvecs_diag @ np.diag(sqrt_eigenvals) @ la.inv(eigenvecs_diag)
        
        try:
            K_inv = la.inv(K)
        except la.LinAlgError:
            # Use pseudoinverse if singular
            K_inv = la.pinv(K)
        
        T = K_inv @ U @ sqrt_H_diag
        
        # Phase correction
        Tpp = T[:dim//2, :dim//2]
        Tnn = T[dim//2:, dim//2:]
        
        # Robust phase extraction
        diag_Tpp = np.diag(Tpp)
        diag_Tnn = np.diag(Tnn)
        
        # Handle zero or very small diagonal elements
        phase_array_p = np.exp(1j * np.angle(diag_Tpp + 1e-15))
        phase_array_n = np.exp(1j * np.angle(diag_Tnn + 1e-15))
        
        V = np.diag(np.concatenate([
            phase_array_p.conj(),
            phase_array_n.conj()
        ]))
        
        T = T @ V
        
        # Extract blocks
        Tpp = T[:dim//2, :dim//2]
        Tnp = T[dim//2:, :dim//2] 
        Tpn = T[:dim//2, dim//2:]
        Tnn = T[dim//2:, dim//2:]
        
        return evals[:dim//2], Tpp, Tnp, Tpn, Tnn, T
    
    def F_function(self, q, n):
        """F function for spin wave calculations"""
        if q == 0:
            return 0
        return 2 * (1 - ((-1)**n) * np.exp(-q)) / q
    
    def P_function(self, q, n, m):
        """P function for spin wave calculations"""
        if n == m == 0:
            delta_factor = 1
        else:
            delta_factor = 1 / np.sqrt((1 + (n == 0)) * (1 + (m == 0)))
        
        term1 = (q**2 / (q**2 + n**2 * np.pi**2)) * (n == m)
        
        if q == 0:
            term2 = 0
        else:
            term2 = (-delta_factor * q**4 / ((q**2 + n**2 * np.pi**2) * (q**2 + m**2 * np.pi**2))) * self.F_function(q, n) * ((1 + (-1)**(n+m)) / 2)
        
        return term1 + term2
    
    def Q_function(self, q, n, m):
        """Q function for spin wave calculations"""
        if n == m == 0:
            delta_factor = 1
        else:
            delta_factor = 1 / np.sqrt((1 + (n == 0)) * (1 + (m == 0)))
        
        term1 = (q**2 / (q**2 + m**2 * np.pi**2)) 

        # Add small regularization to avoid division by zero
        denominator = m**2 - n**2 + ((1 + (-1)**(n+m))/2)
        if m == n:
            denominator = 1
        
        if abs(denominator) < 1e-15:
            denominator = 1e-15

        term2 = ((m**2 / denominator) * (2 / (q + 1e-15))) - ((q**2 / (2 * (q**2 + n**2 * np.pi**2))) * self.F_function(q, n))
        
        term3 = delta_factor * (1 - (-1)**(n+m)) / 2
        
        return term1 * term2 * term3
    
    def Omega_function(self, omega_H, q, n):
        """Omega function for frequency calculation"""
        return (omega_H + ((self.gamma * self.DD / (self.L**2)) * (q**2 + n**2 * np.pi**2))) / self.omega_M
    
    def H_matrix(self, omega_H, q, phi_k, n, m):
        """Hamiltonian matrix element"""
        omega_val = self.Omega_function(omega_H, q, n)
        
        term1 = ((omega_val * np.array([[1, 0],[0, 1]])) + (0.5 * np.array([[1, 1], [1, 1]]))) * (n == m)
        
        term2 = (-0.5 * np.array([[1 - np.sin(phi_k)**2, 1 + np.sin(phi_k)**2], [1 + np.sin(phi_k)**2, 1 - np.sin(phi_k)**2]]) * self.P_function(q, n, m))
        
        term3 = (-0.5 * np.array([[0, -4], [4, 0]]) * np.sin(phi_k) * self.Q_function(q, n, m))
        
        return term1 + term2 + term3
    
    def f_function(self, q, n, h_over_L):
        """f function for coupling calculations"""
        factor = ((-1)**n) / np.sqrt(2 * (1 + (n == 0)))
        
        if q == 0:
            return 0
        
        return (factor * (q**2 / (q**2 + n**2 * np.pi**2)) * np.exp(-q * h_over_L) * (1 - (-1)**n * np.exp(-q)))
    
    def generate_H_BdG(self, omega_H, q, phi_k, N_max):
        """Generate the Bogoliubov-de Gennes Hamiltonian"""
        H_BdG = np.zeros((2 * N_max, 2 * N_max), dtype=complex)
        
        for m in range(N_max):
            for n in range(N_max):
                result = self.H_matrix(omega_H, q, phi_k, n, m)
                H_BdG[n, m] = result[0, 0]
                H_BdG[n, m + N_max] = result[0, 1]
                H_BdG[n + N_max, m] = result[1, 0]
                H_BdG[n + N_max, m + N_max] = result[1, 1]
        
        return H_BdG
    
    def N_Bose(self, omega):
        """Bose-Einstein distribution approximation"""
        return (1e-9 * self.k_B * self.temperature / self.h_plank) / omega
    
    def fspace(self, min_val, max_val, steps, func=np.log):
        """Create logarithmically spaced array"""
        if func == np.log:
            return np.exp(np.linspace(func(min_val), func(max_val), steps))
        else:
            return np.linspace(min_val, max_val, steps)
    
    def multi_para_diag(self, h_NV_array, omega_H, q_table, phi_table, N_max):
        """Multiple paraunitary diagonalization calculations"""
        num_q = len(q_table)
        num_phi = len(phi_table)
        num_h_NV = len(h_NV_array)
        
        # Initialize storage arrays
        omega_BdG_table = np.zeros((num_q, 2 * num_phi, N_max))
        coupling_plus_table = np.zeros((num_q, 2 * num_phi, num_h_NV, N_max), dtype=complex)
        coupling_minus_table = np.zeros((num_q, 2 * num_phi, num_h_NV, N_max), dtype=complex)
        coupling_z_table = np.zeros((num_q, 2 * num_phi, num_h_NV, N_max), dtype=complex)
        
        print("Performing paraunitary diagonalization...")
        
        failed_count = 0
        total_count = num_q * num_phi
        
        for count_q in tqdm(range(num_q), desc="Q values"):
            for count_phi in range(num_phi):
                q = q_table[count_q]
                phi_k = phi_table[count_phi]
                
                # Generate Hamiltonian
                H_BdG = self.generate_H_BdG(omega_H, q, phi_k, N_max)
                
                # Symmetrize Hamiltonian more robustly
                H_sym = (H_BdG + H_BdG.conj().T) / 2
                
                # Add small regularization if needed
                if not np.allclose(H_sym, H_sym.conj().T):
                    H_sym = (H_sym + H_sym.conj().T) / 2
                
                # Paraunitary diagonalization
                try:
                    eval_result, Tpp, Tnp, Tpn, Tnn, T = self.para_unitary_diag(H_sym)
                    
                    # Check for valid results
                    if np.any(np.isnan(eval_result)) or np.any(np.isinf(eval_result)):
                        raise ValueError("Invalid eigenvalues (NaN or Inf)")
                    
                    # Store eigenvalues (with mirroring)
                    omega_BdG_table[count_q, count_phi] = np.real(eval_result)
                    omega_BdG_table[count_q, count_phi + num_phi] = np.real(eval_result)
                    
                    # Calculate coupling coefficients
                    gamma_plus = ((1 + np.sin(phi_k)) / 2) * (Tpp + Tnp + np.sin(phi_k) * (Tpp - Tnp))
                    gamma_plus_mirror = (((1 + np.sin(phi_k + np.pi)) / 2) * (Tnn + Tpn + np.sin(phi_k + np.pi) * (Tnn - Tpn)).conj())
                    
                    gamma_minus = ((1 - np.sin(phi_k)) / 2) * (Tpp + Tnp + np.sin(phi_k) * (Tpp - Tnp))
                    gamma_minus_mirror = (((1 - np.sin(phi_k + np.pi)) / 2) * (Tnn + Tpn + np.sin(phi_k + np.pi) * (Tnn - Tpn)).conj())
                    
                    gamma_z = (-1j * np.cos(phi_k) / 2) * (Tpp + Tnp + np.sin(phi_k) * (Tpp - Tnp))
                    gamma_z_mirror = (-1j * np.cos(phi_k + np.pi) / 2) * ((Tnn + Tpn + np.sin(phi_k + np.pi) * (Tnn - Tpn)).conj())
                    
                    # Calculate f_bar arrays
                    for tt in range(num_h_NV):
                        f_bar_array = np.array([self.f_function(q, nn, h_NV_array[tt] / self.L) 
                                               for nn in range(N_max)])
                        
                        coupling_plus_table[count_q, count_phi, tt] = f_bar_array @ gamma_plus
                        coupling_plus_table[count_q, count_phi + num_phi, tt] = f_bar_array @ gamma_plus_mirror
                        
                        coupling_minus_table[count_q, count_phi, tt] = f_bar_array @ gamma_minus
                        coupling_minus_table[count_q, count_phi + num_phi, tt] = f_bar_array @ gamma_minus_mirror
                        
                        coupling_z_table[count_q, count_phi, tt] = f_bar_array @ gamma_z
                        coupling_z_table[count_q, count_phi + num_phi, tt] = f_bar_array @ gamma_z_mirror
                        
                except Exception as e:
                    failed_count += 1
                    if failed_count < 10:  # Only print first few failures
                        print(f"Warning: Diagonalization failed at q={q:.3f}, phi={phi_k:.3f}: {e}")
                    continue
        
        if failed_count > 0:
            print(f"Total failed diagonalizations: {failed_count}/{total_count} ({100*failed_count/total_count:.1f}%)")
        
        # Create interpolation functions
        print("Creating interpolation functions...")
        
        # Extend tables for periodic boundary conditions
        q_extended = np.tile(q_table[:, np.newaxis], (1, 2 * num_phi + 1))
        phi_extended = np.zeros((num_q, 2 * num_phi + 1))
        
        for n in range(num_q):
            for m in range(2 * num_phi + 1):
                phi_val = phi_table[m % num_phi] + np.pi * (m // num_phi)
                phi_extended[n, m] = phi_val
        
        # Create interpolation objects
        int_omega_BdG_table = []
        int_coupling_plus_table = []
        int_coupling_minus_table = []
        int_coupling_z_table = []
        
        for s in range(N_max):
            # Extend omega table
            omega_extended = np.zeros((num_q, 2 * num_phi + 1))
            omega_extended[:, :-1] = omega_BdG_table[:, :, s]
            omega_extended[:, -1] = omega_BdG_table[:, 0, s]  # Periodic BC
            
            points = np.column_stack([q_extended.flatten(), phi_extended.flatten()])
            values = omega_extended.flatten()
            
            # Remove NaN and Inf values for interpolation
            valid_mask = np.isfinite(values)
            if np.sum(valid_mask) > 0:
                int_omega_BdG_table.append(interp.LinearNDInterpolator(points[valid_mask], values[valid_mask]))
            else:
                int_omega_BdG_table.append(lambda q, phi: 0.0)  # Fallback
            
            # Extend coupling tables
            coupling_plus_extended = np.zeros((num_q, 2 * num_phi + 1, num_h_NV), dtype=complex)
            coupling_minus_extended = np.zeros((num_q, 2 * num_phi + 1, num_h_NV), dtype=complex)
            coupling_z_extended = np.zeros((num_q, 2 * num_phi + 1, num_h_NV), dtype=complex)
            
            coupling_plus_extended[:, :-1] = coupling_plus_table[:, :, :, s]
            coupling_minus_extended[:, :-1] = coupling_minus_table[:, :, :, s]
            coupling_z_extended[:, :-1] = coupling_z_table[:, :, :, s]
            
            coupling_plus_extended[:, -1] = coupling_plus_table[:, 0, :, s]
            coupling_minus_extended[:, -1] = coupling_minus_table[:, 0, :, s]  
            coupling_z_extended[:, -1] = coupling_z_table[:, 0, :, s]
            
            int_coupling_plus_s = []
            int_coupling_minus_s = []
            int_coupling_z_s = []
            
            for tt in range(num_h_NV):
                plus_values = coupling_plus_extended[:, :, tt].flatten()
                minus_values = coupling_minus_extended[:, :, tt].flatten()
                z_values = coupling_z_extended[:, :, tt].flatten()
                
                # Remove invalid values for interpolation
                plus_valid = np.isfinite(plus_values)
                minus_valid = np.isfinite(minus_values)
                z_valid = np.isfinite(z_values)
                
                if np.sum(plus_valid) > 0:
                    int_coupling_plus_s.append(
                        interp.LinearNDInterpolator(points[plus_valid], plus_values[plus_valid]))
                else:
                    int_coupling_plus_s.append(lambda q, phi: 0.0)
                    
                if np.sum(minus_valid) > 0:
                    int_coupling_minus_s.append(
                        interp.LinearNDInterpolator(points[minus_valid], minus_values[minus_valid]))
                else:
                    int_coupling_minus_s.append(lambda q, phi: 0.0)
                    
                if np.sum(z_valid) > 0:
                    int_coupling_z_s.append(
                        interp.LinearNDInterpolator(points[z_valid], z_values[z_valid]))
                else:
                    int_coupling_z_s.append(lambda q, phi: 0.0)
            
            int_coupling_plus_table.append(int_coupling_plus_s)
            int_coupling_minus_table.append(int_coupling_minus_s)
            int_coupling_z_table.append(int_coupling_z_s)
        
        return (q_table, np.min(omega_BdG_table[omega_BdG_table > 0]), phi_table, 
                int_omega_BdG_table, int_coupling_plus_table, 
                int_coupling_minus_table, int_coupling_z_table)
    
    def gamma_values_UL(self, eta_small, NQ_est, N_phi, H0, int_omega_BdG_table,
                       int_coupling_plus_table, int_coupling_minus_table, num_h_NV):
        """Calculate Gamma values for upper and lower transitions"""
        D_NV = 2.87  # GHz
        omega_target_L = D_NV - (H0 * self.gamma)
        omega_target_U = D_NV + (H0 * self.gamma)
        
        # Create Q table
        Q_middle = 5 * self.L
        NQ_half = NQ_est // 2
        Q_max = 50 * self.L
        
        Q_table = np.concatenate([
            [1e-6],
            self.fspace(self.L / 1000, Q_middle, NQ_half)
        ])
        
        last_delta_Q = Q_table[-1] - Q_table[-2]
        Q_extended = np.arange(Q_table[-1] + last_delta_Q, Q_max, last_delta_Q)
        Q_table = np.concatenate([Q_table, Q_extended])
        
        delta_Q = np.diff(Q_table)
        NQ = len(delta_Q)
        delta_phi = 2 * np.pi / N_phi
        
        gamma_L_Hz_array = np.zeros(num_h_NV)
        gamma_U_Hz_array = np.zeros(num_h_NV)
        DOS_L = 0
        DOS_U = 0
        
        N_max = len(int_omega_BdG_table)
        
        print("Calculating Gamma values...")
        
        for s in tqdm(range(N_max), desc="Mode calculation"):
            # Create flat tables for integration
            flat_Q_table = np.tile(Q_table[:-1], N_phi)
            flat_delta_Q_table = np.tile(delta_Q, N_phi)
            flat_phi_table = np.repeat(np.linspace(0, 2*np.pi - delta_phi, N_phi), NQ)
            
            # Evaluate interpolation functions
            flat_omega_BdG_table = np.zeros(len(flat_Q_table))
            flat_coupling_plus_array = np.zeros((num_h_NV, len(flat_Q_table)), dtype=complex)
            flat_coupling_minus_array = np.zeros((num_h_NV, len(flat_Q_table)), dtype=complex)
            
            for ii in range(len(flat_Q_table)):
                q_val = flat_Q_table[ii]
                phi_val = flat_phi_table[ii]
                
                omega_val = int_omega_BdG_table[s](q_val, phi_val)
                flat_omega_BdG_table[ii] = omega_val if (omega_val is not None and np.isfinite(omega_val)) else 0
                
                for tt in range(num_h_NV):
                    plus_val = int_coupling_plus_table[s][tt](q_val, phi_val)
                    minus_val = int_coupling_minus_table[s][tt](q_val, phi_val)
                    
                    flat_coupling_plus_array[tt, ii] = plus_val if (plus_val is not None and np.isfinite(plus_val)) else 0
                    flat_coupling_minus_array[tt, ii] = minus_val if (minus_val is not None and np.isfinite(minus_val)) else 0
            
            # Calculate DOS contributions
            dos_factor = (1 / self.L**2) * (delta_phi / (2 * np.pi)**2)
            
            for jj in range(len(flat_Q_table)):
                if flat_omega_BdG_table[jj] <= 0:
                    continue
                    
                lorentzian_L = (eta_small / np.pi) / ((eta_small**2 + 
                              (self.omega_M * flat_omega_BdG_table[jj] - omega_target_L)**2))
                lorentzian_U = (eta_small / np.pi) / ((eta_small**2 + 
                              (self.omega_M * flat_omega_BdG_table[jj] - omega_target_U)**2))
                
                DOS_L += dos_factor * flat_delta_Q_table[jj] * flat_Q_table[jj] * lorentzian_L
                DOS_U += dos_factor * flat_delta_Q_table[jj] * flat_Q_table[jj] * lorentzian_U
                
                # Calculate Gamma contributions
                bose_factor = 2 * self.N_Bose(self.omega_M * flat_omega_BdG_table[jj]) + 1
                gamma_factor = (2 * np.pi)**2 * self.omega_M * self.omega_d * (delta_phi / (2 * np.pi)**2)
                
                for tt in range(num_h_NV):
                    coupling_plus_sq = np.abs(flat_coupling_plus_array[tt, jj])**2
                    coupling_minus_sq = np.abs(flat_coupling_minus_array[tt, jj])**2
                    
                    gamma_L_Hz_array[tt] += (gamma_factor * bose_factor * 
                                           flat_delta_Q_table[jj] * flat_Q_table[jj] * 
                                           coupling_plus_sq * lorentzian_L)
                    
                    gamma_U_Hz_array[tt] += (gamma_factor * bose_factor * 
                                           flat_delta_Q_table[jj] * flat_Q_table[jj] * 
                                           coupling_minus_sq * lorentzian_U)
        
        return DOS_L, DOS_U, gamma_L_Hz_array, gamma_U_Hz_array
    
    def gamma_1_from_H0(self, H0, h_NV_array):
        """Calculate Gamma_1 for given magnetic field H0"""
        omega_H = self.gamma * H0  # GHz
        
        # Set up calculation parameters
        num_phi = 2 * 90
        num_Q = 2 * 100
        del_phi = np.pi / num_phi
        phi_table = np.arange(0, np.pi, del_phi)
        
        Q_max = 50 * self.L
        q_table = np.concatenate([
            [1e-6],
            self.fspace(self.L / 1000, Q_max, num_Q)
        ])
        
        F_max = 5 # in GHz
        N_max = int(np.ceil((self.L / np.pi) * np.sqrt(F_max / (self.gamma * self.DD))))
        
        print(f"Calculating for H0 = {H0} Oe, N_max = {N_max}")
        
        # Perform paraunitary diagonalization
        result = self.multi_para_diag(h_NV_array, omega_H, q_table, phi_table, N_max)
        (q_table_temp, omega_min, phi_table_mod, int_omega_BdG_table,
         int_coupling_plus_table, int_coupling_minus_table, int_coupling_z_table) = result
        
        # Calculate Gamma values
        eta_small = 0.003
        NQ = 2 * 200
        N_phi = 2 * 360
        
        DOS_L, DOS_U, gamma_L_Hz_array, gamma_U_Hz_array = self.gamma_values_UL(
            eta_small, NQ, N_phi, H0, int_omega_BdG_table,
            int_coupling_plus_table, int_coupling_minus_table, len(h_NV_array))
        
        return DOS_L, DOS_U, gamma_L_Hz_array, gamma_U_Hz_array

def main():
    """Main calculation function"""
    calculator = SpinWaveCalculator()
    
    # NV center heights in micrometers
    h_NV_array = [0.4]#, 0.5, 0.6, 0.7]
    
    # Magnetic field array (in Gauss)
    H0_array = [81, 81.5, 82, 82.5, 83]
    # H0_array = [76, 77, 78, 79, 80, 81, 81.5, 82, 82.5, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 95]
    
    # Initialize result arrays
    DOS_L_array = np.zeros(len(H0_array))
    DOS_U_array = np.zeros(len(H0_array))
    gamma_L_Hz_array = np.zeros((len(h_NV_array), len(H0_array)))
    gamma_U_Hz_array = np.zeros((len(h_NV_array), len(H0_array)))
    
    # Main calculation loop
    print(f"Starting calculations for {len(H0_array)} field values...")
    
    for ii, H0 in enumerate(tqdm(H0_array, desc="Field values")):
        try:
            result_gamma_UL = calculator.gamma_1_from_H0(H0, h_NV_array)
            DOS_L_array[ii] = result_gamma_UL[0]
            DOS_U_array[ii] = result_gamma_UL[1]
            gamma_L_Hz_array[:, ii] = result_gamma_UL[2]
            gamma_U_Hz_array[:, ii] = result_gamma_UL[3]
            
            print(f"Completed H0 = {H0} Oe")
            
        except Exception as e:
            print(f"Error at H0 = {H0}: {e}")
            continue
    
    # Plotting results
    plt.figure(figsize=(12, 8))
    
    # Plot Gamma values for first NV height
    plt.subplot(2, 2, 1)
    plt.plot(H0_array, gamma_L_Hz_array[0], 'o-', label='ΓL')
    plt.plot(H0_array, gamma_U_Hz_array[0] * 1000, 's-', label='ΓU × 1000')
    plt.xlabel('Field (G)')
    plt.ylabel('ΓL, ΓU × 1,000 (Hz)')
    plt.legend()
    plt.grid(True)
    plt.title(f'Gamma values for h_NV = {h_NV_array[0]} μm')
    
    # Plot DOS
    plt.subplot(2, 2, 2)
    plt.plot(H0_array, DOS_L_array / calculator.L, 'o-', label='DOS_L')
    plt.plot(H0_array, DOS_U_array / calculator.L, 's-', label='DOS_U')
    plt.xlabel('Field (G)')
    plt.ylabel('2π×DOS (1/GHz μm³)')
    plt.legend()
    plt.grid(True)
    plt.title('Density of States')
    
    # Plot all NV heights for Gamma_L
    plt.subplot(2, 2, 3)
    for i, h_NV in enumerate(h_NV_array):
        plt.plot(H0_array, gamma_L_Hz_array[i], 'o-', label=f'h_NV = {h_NV} μm')
    plt.xlabel('Field (G)')
    plt.ylabel('ΓL (Hz)')
    plt.legend()
    plt.grid(True)
    plt.title('ΓL for different NV heights')
    
    # Calculate and plot ratio
    plt.subplot(2, 2, 4)
    ratio = np.zeros_like(H0_array, dtype=float)
    for ii in range(len(H0_array)):
        if DOS_L_array[ii] > 0:
            ratio[ii] = (1e-3 * gamma_L_Hz_array[0, ii] * calculator.L / 
                        (1e-6 * DOS_L_array[ii] * 2 * calculator.N_Bose(2.87 - calculator.gamma * H0_array[ii])))
    
    plt.plot(H0_array, ratio, 'o-')
    plt.xlabel('Field (G)')
    plt.ylabel('Ratio (kHz² μm³)')
    plt.grid(True)
    plt.title('Ratio calculation')

    plt.tight_layout()
    timestr = time.strftime("%Y%m%d_%H%M%S")
    basename = str(__file__) + '_' + timestr 
    figname = basename + '.png'
    plt.savefig(figname)
    plt.show()
    
    # # Save results
    # results_dict = {
    #     'H0_array': H0_array,
    #     'DOS_L_array': DOS_L_array,
    #     'DOS_U_array': DOS_U_array,
    #     'h_NV_array': h_NV_array,
    #     'gamma_L_Hz_array': gamma_L_Hz_array,
    #     'gamma_U_Hz_array': gamma_U_Hz_array
    # }
    
    # filename = basename + '.npz'
    # np.savez(filename, **results_dict)
    # print("Results saved")
    
    # return results_dict

if __name__ == "__main__":
    results = main()
