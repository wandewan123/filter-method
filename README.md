# FILTER_METHOD

A fast Fortran-based iterative eigenvalue solver for one-dimensional Schrödinger equations using the filtering method.

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

Using finite-difference discretization, the problem becomes a tridiagonal system.

The filtering iteration employed in this work is based on repeatedly solving

$\left(\hat{H} - E_{n}^{k}\right)\phi_{n}^{k+1}(x) = \phi_{n}^{k}(x)$

followed by normalization and energy updates until convergence is achieved.

The eigenvalue is computed from the expectation value

:contentReference[oaicite:3]{index=3}

with iterative mixing applied for numerical stability.

---

# Numerical Algorithm

The implemented procedure consists of:

1. Initial wavefunction guess
2. Construction of the tridiagonal Hamiltonian
3. Iterative filtering step
4. Wavefunction normalization
5. Energy expectation evaluation
6. Energy mixing update
7. Convergence check

The tridiagonal systems are solved efficiently using the Thomas algorithm implemented directly in Fortran.

---

# References

This project is based on research involving iterative filtering methods for quantum eigenvalue problems.
If you use this library in academic work, please cite the corresponding research publication:

"Filter method without boundary-value condition for simultaneous calculation of eigenfunction and eigenvalue of a stationary Schrödinger equation on a grid."
DOI:10.1103/PhysRevE.96.033302

---

# License

MIT License
