# Hapus 'from setuptools import setup' yang berlebihan
from numpy.distutils.core import Extension, setup

# Menghubungkan ke file Fortran di folder src
# Urutan file sangat penting: Modul dasar didahulukan!
ext = Extension(
    name='dvipy.solver_lib', # <-- Ditambahkan awalan 'dvipy.'
    sources=[
        'src/matrix.f90',    # Asumsi: ini modul dasar tidak memanggil yang lain
        'src/eigene.f90',    # Asumsi: modul dasar
        'src/eigenf.f90',    # Asumsi: modul dasar
        'src/main.f90',      # Memanggil matrix, eigene, eigenf
        'src/solver.f90'     # (Jika ini adalah file wrapper utamanya)
    ]
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
