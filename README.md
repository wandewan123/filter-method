# FILTER METHOD 1D

This project provides a lightweight scientific-computing framework for solving quantum eigenvalue problems with a hybrid Fortran–Python architecture. The numerical core is written in modern Fortran for performance, while Python provides a simple and flexible interface for experimentation and visualization.

---

# Overview

The library is designed to solve the stationary Schrödinger equation

$\hat{H}\phi = E\phi$

where the Hamiltonian operator is discretized using finite differences.  
Instead of performing direct matrix diagonalization, this method applies an iterative filtering process to converge toward the desired eigenstate.

The approach is particularly attractive for:

- large grid systems,
- sparse Hamiltonians,
- iterative quantum simulations,
- and memory-efficient eigenvalue computations.

---

# Mathematical Background

For a one-dimensional system, the Hamiltonian operator is written as

$\hat{H} = -\frac{\hbar^{2}}{2m}\frac{\partial^{2}}{\partial x^{2}} + V(x)$

By applying a finite-difference discretization scheme, the continuous differential operator is mapped onto a discrete grid.

Depending on the chosen stencil size (3, 5, or 7 points), the kinetic energy term transforms the problem into a banded matrix system (tridiagonal, pentadiagonal, or heptadiagonal, respectively).

Instead of deploying computationally expensive global diagonalization algorithms, the filtering iteration in this framework isolates the target eigenstate by repeatedly solving the shifted linear system:

$\left(\hat{H} - E_{n}^{k}\right)\phi_{n}^{k+1}(x) = \phi_{n}^{k}(x)$

This step is followed by an iterative wavefunction normalization, an energy expectation value evaluation, and a relaxation mixing update until the specified convergence tolerance is achieved.

---

# Numerical Algorithm

The implemented procedure consists of:

1. Initial wavefunction guess
2. Construction of the matrix Hamiltonian
3. Iterative filtering step
4. Wavefunction normalization
5. Energy expectation evaluation
6. Energy mixing update
7. Convergence check

The matrix systems are solved efficiently using the Thomas algorithm implemented directly in Fortran.

---

# How to Use (will be updated)

Since the library contains optimized Fortran compiled cores, it currently requires a manual compilation step. You can easily clone, compile, and run the solver in a **Google Colab** environment or any Linux-based system with the following steps:

### 1. Environment Setup & Cloning
First, navigate to your working directory, clone the repository, and install the required build tools:
```bash
%cd /content
!rm -rf filter-method
!git clone [https://github.com/wandewan123/filter-method.git](https://github.com/wandewan123/filter-method.git)
%cd filter-method
!pip install meson ninja numpy
```
### 2. Compiling the Fortran Core
Compile the high-performance Fortran source files into a Python-binded shared library (.so) using ```f2py```
```bash
!f2py -c -m solver_lib \
    src/utils.f90 \
    src/solver.f90 \
    src/matrix.f90 \
    src/eigenf.f90 \
    src/eigene.f90 \
    src/main.f90

# Move the compiled library into the python interface directory
!mv solver_lib*.so dvipy/
```
### 3. Running the Solver in Python
Once the compilation is finished, you can import and use the high-level Python interface ```solve_full``` directly in your scripts:
```bash
from dvipy import solve_full

# Example usage (adjust parameters based on your physical system)
# energy, wavefunction, status = solve_full(x, vpot, guessE, ...)
```

To help you get started quickly without setting up anything locally, we provide an interactive Google Colab notebook containing complete, runnable examples of the solver:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1sRqObdd0HV1WzeZinW3yt4oJXPpSeLEO?usp=sharing)

Simply click the badge above, make a copy to your drive, and run the cells!

---

# References

This project is based on research involving iterative filtering methods for quantum eigenvalue problems.
If you use this library in academic work, please cite the corresponding research publication.

📄[Filter method without boundary-value condition for simultaneous calculation of eigenfunction and eigenvalue of a stationary Schrödinger equation on a grid](https://doi.org/10.1103/PhysRevE.96.033302)

---

# License

MIT License
