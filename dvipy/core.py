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

def solve_full(x, vpot, initial_e, psi0, stencil=3, xi=0.5, tol=1e-7, maxiter=1000):
    return solver_lib.main_filter.main(
        x, vpot, initial_e, psi0, xi, stencil, tol, maxiter
    )
