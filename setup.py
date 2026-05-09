from setuptools import setup
from numpy.distutils.core import Extension, setup

# Menghubungkan ke file Fortran di folder src
ext = Extension(
    name='solver_lib', 
    sources=['src/solver.f90']
)

setup(
    name='dvipy',
    version='0.1.0',
    packages=['dvipy'],
    ext_modules=[ext],
    install_requires=[
        'numpy',
    ],
)
