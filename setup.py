from numpy.distutils.core import Extension, setup

ext = Extension(
    name       =    'dvipy.solver_lib', 
    sources    =    [
        'src/solver.f90',
        'src/utils.f90',
        'src/matrix.f90',
        'src/eigenf.f90',
        'src/eigene.f90',
        'src/main.f90'
    ]
)

setup(
    name                =    'dvipy',
    version             =    '0.1.0',
    packages            =    ['dvipy'],
    ext_modules         =    [ext],
    install_requires    =    ['numpy'],
)
