import numpy as np
# Kita asumsikan nanti hasil kompilasi f2py bernama 'solver_lib'
try:
    import solver_lib
except ImportError:
    solver_lib = None

def tridag(a, b, c, r):
    """
    Wrapper untuk memanggil subroutine tridag dari Fortran.
    """
    if solver_lib is None:
        raise ImportError("Modul Fortran 'solver_lib' belum dikompilasi. Jalankan f2py dulu.")
    
    n = len(r)
    # Memanggil subroutine dari module 'solver_module' yang ada di Fortran
    u = solver_lib.solver_module.tridag(a, b, c, r, n)
    return u
