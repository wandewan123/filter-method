import numpy as np
try:
    import solver_lib
except ImportError:
    solver_lib = None

def tridag(a, b, c, r):
    n = len(r)
    return solver_lib.solver_module.tridag(a, b, c, r, n)

def psi_eigen(dx, vpot, psi, energy):
    # n otomatis diurus f2py jika intent(in) diatur benar
    return solver_lib.solver_module.psi_eigen_fortran(dx, vpot, psi, energy, len(vpot))

def energi_eigen(dx, vpot, psi):
    return solver_lib.solver_module.energi_eigen_fortran(dx, vpot, psi, len(vpot))

def solve_full(x, vpot, initial_e, dx, xi=0.95, tol=1e-7):
    # Ini fungsi "all-in-one" kencang yang kamu mau
    return solver_lib.solver_module.main_filter(x, len(x), vpot, initial_e, dx, xi, tol)
