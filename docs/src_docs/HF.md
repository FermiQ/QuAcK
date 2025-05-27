# `src/HF` - Hartree-Fock Module

## Overview

The `src/HF` directory contains the Fortran source code for performing various types of Hartree-Fock (HF) calculations. Hartree-Fock is a fundamental *ab initio* method that approximates the many-electron wavefunction as a single Slater determinant. The routines in this module are responsible for the Self-Consistent Field (SCF) procedure to determine the optimal molecular orbitals and the corresponding HF energy.

This module implements several variants of the Hartree-Fock method, including:
*   Restricted Hartree-Fock (RHF): For closed-shell systems where all electrons are paired.
*   Unrestricted Hartree-Fock (UHF): For open-shell systems where alpha and beta electron spins are treated independently.
*   Generalized Hartree-Fock (GHF): A more general formulation that does not assume spin purity and allows mixing of alpha and beta orbitals.
*   Restricted Open-shell Hartree-Fock (ROHF): For open-shell systems, but with restrictions on the spatial part of the orbitals.
*   Complex Restricted Hartree-Fock (cRHF): For calculations involving complex-valued orbitals, often used in the presence of magnetic fields or with certain types of complex absorbing potentials.
*   Hartree-Fock-Bogoliubov (HFB): A method that generalizes HF to treat pairing correlations, particularly important in nuclear physics and some condensed matter systems.

The module also includes routines for:
*   Generating initial guesses for molecular orbitals.
*   Accelerating SCF convergence (e.g., DIIS extrapolation, level shifting).
*   Searching for stable HF solutions by iteratively performing HF calculations and stability analyses.
*   Performing stability analysis of a converged HF solution to ensure it corresponds to a local minimum.
*   Calculating properties like dipole moments from the HF wavefunction.
*   Printing detailed output and results of the HF calculations.

## Key Components

### Main Hartree-Fock Methods
*   **`RHF(dotest, maxSCF, thresh, max_diis, guess_type, level_shift, ... ERHF, eHF, c, P, F)` (`RHF.f90`)**:
    *   Performs a Restricted Hartree-Fock calculation.
    *   Key Inputs: SCF control parameters (`maxSCF`, `thresh`, `max_diis`, `level_shift`), molecular/basis information, integrals (`S`, `T`, `V`, `Hc`, `ERI`), initial guess type.
    *   Key Outputs: RHF energy (`ERHF`), orbital energies (`eHF`), MO coefficients (`c`), density matrix (`P`), Fock matrix (`F`).
*   **`UHF(dotest, maxSCF, thresh, max_diis, guess_type, mix, level_shift, ... EUHF, eHF, c, P, F)` (`UHF.f90`)**:
    *   Performs an Unrestricted Hartree-Fock calculation.
    *   Treats alpha and beta spins separately. Includes a mixing parameter (`mix`) for initial guess in singlet states.
    *   Key Inputs/Outputs: Similar to RHF, but matrices for `eHF`, `c`, `P`, `F` have a spin dimension.
*   **`GHF(dotest, maxSCF, thresh, max_diis, guess_type, mix, level_shift, ... EGHF, eHF, c, P, F)` (`GHF.f90`)**:
    *   Performs a Generalized Hartree-Fock calculation.
    *   Uses a "super" matrix representation doubling the basis size to handle alpha and beta components together.
    *   Key Inputs/Outputs: Similar to RHF/UHF but adapted for the GHF formalism (e.g., `nBas2` for doubled basis size).
*   **`ROHF(...)` (`ROHF.f90`)**:
    *   Implements Restricted Open-shell Hartree-Fock. (Details to be filled in by further reading if a specific file is documented).
*   **`cRHF(...)` (`cRHF.f90`)**:
    *   Implements complex Restricted Hartree-Fock. (Details to be filled in by further reading if a specific file is documented).
*   **`HFB(...)` (`HFB.f90`)**:
    *   Implements Hartree-Fock-Bogoliubov theory. (Details to be filled in by further reading if a specific file is documented).

### SCF Procedure Support
*   **Initial Guess Routines**:
    *   `mo_guess(nBas, nOrb, guess_type, S, Hc, X, c)` (`mo_guess.f90`): Main dispatcher for MO guess.
        *   `guess_type = 0`: Read from existing `c` matrix (restart).
        *   `guess_type = 1`: Core Hamiltonian guess (`core_guess.f90`).
        *   `guess_type = 2`: Huckel guess (`huckel_guess.f90`).
        *   `guess_type = 3`: Random guess.
    *   `core_guess(nBas, nOrb, Hc, X, c)` (`core_guess.f90`): Diagonalizes the core Hamiltonian in the orthogonal basis.
    *   `huckel_guess(...)` (`huckel_guess.f90`): Implements a Huckel-type guess.
    *   `complex_core_guess.f90`, `complex_mo_guess.f90`: Versions for complex HF.
    *   `mix_guess(nBas, nO, mix, c)` (`mix_guess.f90`): Mixes HOMO/LUMO for UHF calculations on singlets to break spin symmetry.
*   **Convergence Acceleration**:
    *   `DIIS_extrapolation(...)` (called by HF routines, likely in `src/utils/DIIS_extrapolation.f90`): Implements the Direct Inversion in the Iterative Subspace method.
    *   `level_shifting(level_shift, nBas, nOrb, nO, S, c, F)` (called by HF routines, likely in `src/utils/level_shifting.f90`): Applies a level shift to virtual orbitals to aid convergence.
*   **Search and Stability Analysis**:
    *   `RHF_search(...)` (`RHF_search.f90`): Iteratively calls an RHF calculation and then performs stability analysis. If an instability is found, it perturbs the MO coefficients based on the eigenvector of the instability and restarts the RHF calculation.
    *   `UHF_search(...)` (`UHF_search.f90`), `GHF_search(...)` (`GHF_search.f90`): Similar search routines for UHF and GHF.
    *   `RHF_stability(nBas, nC, nO, nV, nR, nS, eHF, ERI)` (`RHF_stability.f90`): Performs stability analysis for a converged RHF solution by checking eigenvalues of stability matrices (e.g., Real RHF -> Real RHF, Real RHF -> Complex RHF, Real RHF -> Real UHF).
    *   `UHF_stability(...)` (`UHF_stability.f90`), `GHF_stability(...)` (`GHF_stability.f90`): Similar stability routines for UHF and GHF.
    These routines typically build and diagonalize matrices related to the Linear Response (LR) or Random Phase Approximation (RPA) problem (e.g., using `phRLR_A`, `phRLR_B` from the `src/LR` module).

### Utility and Output Routines
*   `density_matrix(nBas, nO, c, P)` (`density_matrix.f90`): Calculates the density matrix from MO coefficients. (Note: HF routines often calculate P directly).
*   `density(nBas, P, RGRID, rho)` (`density.f90`): Calculates electron density on a grid.
*   `dipole_moment(nBas, P, nNuc, ZNuc, rNuc, dipole_int, dipole)` (`dipole_moment.f90`): Computes the electric dipole moment.
*   `print_RHF(...)`, `print_UHF(...)`, `print_GHF(...)`, etc.: Subroutines for printing detailed results of the respective HF calculations (energies, MO coefficients, orbital energies, dipole moments, etc.).
*   `read_restart_HFB.f90`, `write_restart_HFB.f90`: Specific to HFB restart capabilities.

## Important Variables/Constants

Common variables across many HF routines include:
*   `maxSCF`: Maximum number of SCF iterations.
*   `thresh`: Convergence threshold for the SCF procedure (typically on density or energy change).
*   `max_diis`: Maximum number of vectors for DIIS extrapolation.
*   `level_shift`: Energy shift applied to virtual orbitals.
*   `nBas`: Number of atomic basis functions.
*   `nOrb`: Number of molecular orbitals (often same as `nBas` unless linear dependencies are removed).
*   `nO` (or `nO(nspin)` for UHF/GHF): Number of occupied orbitals (per spin).
*   `S`, `T`, `V`, `Hc`: One-electron integral matrices (Overlap, Kinetic, Nuclear Attraction, Core Hamiltonian).
*   `ERI` (or `ERI_AO`): Two-electron repulsion integrals in AO basis.
*   `X`: Orthogonalization matrix.
*   `EHF` (e.g., `ERHF`, `EUHF`, `EGHF`): Total Hartree-Fock energy (electronic part).
*   `eHF`: Molecular orbital energies.
*   `c`: Molecular orbital coefficients.
*   `P`: Density matrix.
*   `F`: Fock matrix.
*   `ENuc`: Nuclear repulsion energy (usually passed in to be added to `EHF` for total energy).

## Dependencies and Interactions

*   **`parameters.h` (`include/parameters.h`)**: Defines global parameters and constants (e.g., `nspin`, `ncart`).
*   **Integral Generation**: Relies on one- and two-electron integrals, typically provided by `PyDuck.py` (via `pyscf`) and read from files or passed as arguments.
*   **AO to MO Transformation Module (`src/AOtoMO/`)**: Used by search/stability routines (e.g., `AOtoMO_ERI_RHF`) to transform ERIs from AO to MO basis.
*   **Linear Response Module (`src/LR/`)**: Stability analysis routines utilize functions from this module (e.g., `phRLR_A`, `phRLR_B`) to construct and solve parts of the LR eigenvalue problem.
*   **Utility Modules (`src/utils/`)**:
    *   `DIIS_extrapolation.f90`: For DIIS.
    *   `level_shifting.f90`: For level shifting.
    *   `diagonalize_matrix` (likely from `wrap_lapack.f90`): For matrix diagonalization.
    *   `trace_matrix` (likely from `utils.f90`): For matrix traces.
*   **Input Files**:
    *   Molecular geometry, basis set info (indirectly via integrals).
    *   Integral files (e.g., `Ov.dat`, `Nuc.dat`, `Kin.dat`, `ERI.bin`).
*   **Output Files**:
    *   Writes detailed log output to standard output.
    *   May write restart information (e.g., `write_restart_HFB.f90`).
    *   Test values via `dump_test_value` (from `src/test/`) if `dotest` is true.
