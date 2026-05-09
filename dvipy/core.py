import numpy as np
from . import solver_lib

def tridag(a, b, c, r):
    return solver_lib.solver_module.tridag(a, b, c, r)

def psi_eigen(dx, vpot, psi, energy):
    return solver_lib.solver_module.psi_eigen_fortran(
        dx, vpot, psi, energy
    )

def energi_eigen(dx, vpot, psi):
    return solver_lib.solver_module.energi_eigen_fortran(
        dx, vpot, psi
    )

def solve_full(x, vpot, initial_e, dx, xi=0.95, tol=1e-7):
    return solver_lib.solver_module.main_filter(
        x, vpot, initial_e, dx, xi, tol
    )
