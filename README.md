# NumericalMethodsForSoftMatter

This repository contains the homework assignments and code developed for the Numerical Methods for Soft Matter course at the University of Padova, taught by Prof. Orlandini and Prof. Locatelli during the 2025/2026 academic year. A comprehensive report detailing the theoretical background and results is available in `report_Part1.pdf` and `report_Part2.pdf`. 
In particular, it includes:
(Part1)
- Homework 1: Sampling random numbers from a given arbitrary 1D distribution (using the inverse CDF or rejection method) or from a multi-dimensional distribution.
- Homework 2: Direct Monte Carlo method to numerically estimate integrals in one or more dimensions, including the use of importance sampling to reduce variance
- Homework 3: Theory of Markov chains (aperiodicity, irreducibility, ergodicity, invariant distribution, etc.).
- Homework 4: 2D Ising model at zero field using the Metropolis algorithm with local dynamics. Study of phase transitions and order parameters (magnetization, energy, specific heat, and susceptibility). Diagnostics of the Markov chain near criticality, including estimates of the correlation length, critical slowing down, and integrated correlation time. The full code for the 2D Ising simulation is written in C++ for efficiency and is hosted in a separate dedicated repository, which extends beyond this specific homework with different MCMC engines: [YAIS repo](https://github.com/lorenzorizzi17/YAIS). Have a look at it!
- Homework 5: Simulation of the [Active Ising model](https://arxiv.org/abs/1506.05749) in 2D with periodic boundary conditions, illustrating a non-equilibrium steady-state (NESS) system and out-of-equilibrium phase transitions. Characterization of the three phases in the $(\rho, T)$ plane (gas, liquid, coexistence) and phase portrait. The code for the Active Ising model is hosted in a [dedicated repository](https://github.com/lorenzorizzi17/ActiveIsingModel).
(Part 2)
- Homework 6: Simulation of the (Patchy) Hard Sphere model implementing Cell Lists and Verlet Lists for efficient neighbor searching in MD
- Homework 7: Equations of motion numerical integrators: 1st order Euler, symplectic integrators like Velocity Verlet. LAMMPS software as an efficient MD simulator.
- Homework 8: Molecular dynamics simulation in the canonical ensemble; algorithmic implementation of canonical and non-canonical thermostats (Berendsen, Andersen)
- Homework 9: Multiple Markov chains and Parallel Tempering. Multiple Histograms method to extract thermal averages at an arbitrary $\beta$
- Homework 10: Gillespie algorithm for continuous time Markov processes
- Homework 11: Numerical implementation of Langevin and Brownian processes.

The codebase can be found in the `code/` folder. Instructions on how to compile and build are to be found in the local README of each directory
