from numpy.distutils.core import Extension, setup

ext = Extension(
    name       =    'dvipy.solver_lib', 
    sources    =    [
        'src/solver.f90',    # Level dasar: dipanggil oleh eigenf
        'src/utils.f90',
        'src/matrix.f90',    # Level dasar: dipanggil oleh eigene
        'src/eigenf.f90',    # Level menengah: butuh solver, dipanggil oleh main
        'src/eigene.f90',    # Level menengah: butuh matrix, dipanggil oleh main
        'src/main.f90'       # Level puncak: butuh eigenf dan eigene
    ]
)

setup(
    name                =    'dvipy',
    version             =    '0.1.0',
    packages            =    ['dvipy'],
    ext_modules         =    [ext],
    install_requires    =    ['numpy'],
)
