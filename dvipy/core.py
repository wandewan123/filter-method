import numpy as np
from . import solver_lib

def energi_eigen(dx, vpot, psi, stencil=3):
    vpot_arr = np.asarray(vpot, dtype=np.float64)
    psi_arr  = np.asarray(psi, dtype=np.float64)
    n        = len(psi_arr)
    
    energy_out, status = solver_lib.eigene_mod.eigene(dx, vpot_arr, psi_arr, n, stencil)
    
    if status != 0:
        raise RuntimeError(f"Gagal menghitung energi. Status Fortran: {status}")
        
    return energy_out

def psi_eigen(dx, vpot, psi_old, E_guess, stencil=3):
    vpot_arr     = np.asarray(vpot, dtype = np.float64)
    psi_old_arr  = np.asarray(psi_old, dtype = np.float64)
    n            = len(psi_old_arr)
    
    if stencil == 3:
        phi_new, status = solver_lib.eigenf.eigenf2(dx, vpot_arr, psi_old_arr, n, E_guess)
    elif stencil == 5:
        phi_new, status = solver_lib.eigenf.eigenf4(dx, vpot_arr, psi_old_arr, n, E_guess)
    elif stencil == 7:
        phi_new, status = solver_lib.eigenf.eigenf6(dx, vpot_arr, psi_old_arr, n, E_guess)
    else:
        raise ValueError("Stencil tidak valid. Gunakan 3, 5, atau 7.")
        
    if status != 0:
        raise RuntimeError(f"Gagal menyelesaikan (H-E)phi. Status Fortran: {status}")
        
    return phi_new

def solve_full(x, vpot, initial_e, psi0, stencil=3, xi=0.5, tol=1e-7, maxiter=1000):
    x_arr     = np.asarray(x, dtype = np.float64)
    vpot_arr  = np.asarray(vpot, dtype = np.float64)
    psi0_arr  = np.asarray(psi0, dtype = np.float64)
    
    energy_final, phi_final, status = solver_lib.main_filter.solver_filter(
        x_arr, vpot_arr, initial_e, psi0_arr, xi, stencil, tol, maxiter
    )
    
    if status != 0:
        if status == 1:
            print(f"Warning: Solver tidak konvergen setelah {maxiter} iterasi!")
        elif status == -1:
            raise MemoryError("Fortran Error: Gagal alokasi memori.")
        elif status == -2:
            raise ValueError("Fortran Error: Stencil tidak valid. Gunakan 3, 5, atau 7.")
        elif status == -4:
            raise ValueError(f"Fortran Error: Jumlah grid terlalu kecil untuk stencil {stencil}.")
        else:
            raise RuntimeError(f"Fortran Error: Kode status {status}")
            
    return energy_final, phi_final, status
