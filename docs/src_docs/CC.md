# `src/CC` - Coupled Cluster Module

## Overview

The `src/CC` directory hosts Fortran routines for performing various Coupled Cluster (CC) calculations. Coupled Cluster theory is a high-accuracy quantum chemistry method for treating electron correlation. It parameterizes the correlated wavefunction using an exponential ansatz of cluster operators acting on a reference determinant (usually Hartree-Fock).

This module implements several standard CC methods, Equation-of-Motion CC (EOM-CC) for excited states, and potentially some specialized or corrected CC variants. The core of these methods involves solving systems of non-linear equations for cluster amplitudes (T-amplitudes for ground state, R-amplitudes for EOM states).

Key functionalities include:
*   Standard ground-state CC methods like CCD, CCSD, and potentially up to CCSDT.
*   Spin-adapted (e.g., RHF-based RCC) and possibly spin-generalized (GHF-based GCC) implementations.
*   Routines for calculating the (T) correction in CCSD(T).
*   Equation-of-Motion CC methods for calculating excitation energies (EE), ionization potentials (IP), and electron affinities (EA).
*   Numerous routines (`form_*`) dedicated to constructing the various intermediates, Fock matrix elements dressed with amplitudes, and residual equations required for solving the amplitude equations.
*   Lambda equation solvers (`l*CCD`) for calculating properties and gradients.
*   Calculation of correlation energies.

## Key Components

### Main Coupled Cluster Methods
*   **`RCC(dotest, doCCD, doCCSD, doCCSDT, ...)` (`RCC.f90`)**:
    *   Acts as a driver routine for various Restricted (RHF-based) Coupled Cluster methods.
    *   Calls specific subroutines like `CCD`, `DCD`, `CCSD`, `drCCD`, `rCCD`, `crCCD`, `lCCD`, `pCCD` based on logical flags.
*   **`CCSD(dotest, maxSCF, thresh, max_diis, doCCSDT, ...)` (`CCSD.f90`)**:
    *   Performs a Coupled Cluster Singles and Doubles (CCSD) calculation, typically based on spin-orbitals.
    *   Transforms spatial MO energies and ERIs to spin-orbital form.
    *   Iteratively solves for T1 (single) and T2 (double) excitation amplitudes.
    *   Can optionally call `CCSDT` to compute the perturbative triples (T) correction if `doCCSDT` is true.
    *   Key Inputs: SCF control parameters, MO energies, MO-basis ERIs.
    *   Key Outputs: CCSD correlation energy, T1 and T2 amplitudes.
*   **`CCD(...)` (`CCD.f90`)**:
    *   Performs a Coupled Cluster Doubles (CCD) calculation.
    *   Solves for T2 amplitudes.
*   **`CCSDT(...)` (`CCSDT.f90` and `form_T.f90`)**:
    *   Likely contributes to or performs the full iterative CCSDT or calculates the (T) correction. The `form_T.f90` specifically calculates the perturbative triples correction `EcCCT` using pre-computed `ub` and `ubb` intermediates.
*   **`GCC(...)` (`GCC.f90`, `GCCSD.f90`)**:
    *   Likely implements Generalized Coupled Cluster methods, possibly based on GHF references.
*   **Specialized/Corrected Methods**:
    *   `DCD.f90`: Distinguishable Cluster Doubles.
    *   `drCCD.f90`, `rCCD.f90`, `crCCD.f90`: Variants of CCD, possibly direct ring, ring, and crossed-ring CCD, often aimed at improving efficiency or accuracy for specific systems or breaking bonds.
    *   `lCCD.f90`, `lGCCD.f90`: Lambda equation solvers for CCD and GCCD, necessary for response properties or gradients.
    *   `pCCD.f90`: Pair CCD.

### Equation-of-Motion (EOM) Coupled Cluster
*   **`EE_EOM_CCD_1h1p(nC, nO, nV, nR, eO, eV, OOVV, OVVO, t)` (`EE_EOM_CCD_1h1p.f90`)**:
    *   Calculates Electron Excitation (EE) energies using EOM-CCD with 1-hole-1-particle (1h1p) excitations (similar to Configuration Interaction Singles for the CCD reference).
    *   Constructs and diagonalizes an EOM Hamiltonian using CCD T2 amplitudes and MO-basis ERIs.
*   **`DEA_EOM_CCD_2p.f90`**: Double Electron Affinity EOM-CCD (2-particle states).
*   **`DIP_EOM_CCD_2h.f90`**: Double Ionization Potential EOM-CCD (2-hole states).
*   Other EOM routines might exist for CCSD (e.g., EOM-CCSD for EE, IP, EA).

### Amplitude Equation Construction (`form_*` routines)
This module contains a large number of subroutines named `form_*` which are crucial for CC calculations. These routines construct various intermediates and terms appearing in the CC amplitude and residual equations. Examples include:
*   `form_delta_OV(...)`, `form_delta_OOVV(...)`: Form energy denominators for T1 and T2 amplitudes.
*   `form_tau(..., t1, t2, tau)`: Forms the effective T2 amplitude `tau = t2 + t1*t1`.
*   `form_h(..., t1, tau, hvv, hoo, hvo)`: Forms dressed Fock matrix elements (h_vv, h_oo, h_vo).
*   `form_g(...)`: Forms other intermediates.
*   `form_abh(...)`: Forms specific ERI-like intermediates.
*   `form_r1(...)`, `form_r2(...)`: Compute the residuals for T1 and T2 amplitude equations.
*   `form_EOM_one_body.f90`, `form_EOM_two_body.f90`: Likely form parts of the EOM Hamiltonians.
*   Many others for specific diagrams or terms in various CC methods.

### Correlation Energy Calculation
*   `CCD_correlation_energy(...)`
*   `CCSD_correlation_energy(..., OOVV, tau, EcCCSD)`: Calculates the CCSD correlation energy using the converged T1 (`tau` contains T1 effects) and T2 amplitudes and OOVV ERIs.
*   `MP3_correlation_energy.f90`: Calculates MP3 correlation energy, which can be derived from some CC intermediates.

## Important Variables/Constants

*   `maxSCF`: Maximum iterations for solving amplitude equations.
*   `thresh`: Convergence threshold for amplitude residuals or energy.
*   `nBasin` (often in driver), `nOrb` (in specific methods, refers to MOs): Number of basis functions/molecular orbitals.
*   `nC`, `nO`, `nV`, `nR`: Number of core, occupied, virtual, and removed orbitals, respectively.
*   `eHF` (or `eO`, `eV`): Molecular orbital energies (occupied, virtual).
*   `ERI` (or `OOVV`, `OVVO`, etc.): Two-electron repulsion integrals in the MO basis, often partitioned into occupied-virtual blocks.
*   `t1(nO,nV)`: T1 (single excitation) amplitudes.
*   `t2(nO,nO,nV,nV)`: T2 (double excitation) amplitudes.
*   `EcCCSD`, `EcCCT`: CCSD correlation energy, (T) correlation energy correction.
*   `delta_OV`, `delta_OOVV`: Energy denominators used in updating amplitudes.
*   `r1`, `r2`: Residuals for T1 and T2 equations.

## Dependencies and Interactions

*   **Hartree-Fock Module (`src/HF/`)**: Requires converged HF results (MO energies `eHF`, MO coefficients `cHF` for integral transformation, reference energy `ERHF`).
*   **AO to MO Transformation (`src/AOtoMO/`)**: ERIs must be transformed from AO to MO basis before being used in most CC routines. This is usually handled before calling the main CC drivers.
*   **`parameters.h`**: Defines global parameters.
*   **Utility Modules (`src/utils/`)**:
    *   `DIIS_extrapolation.f90` (potentially): For DIIS convergence of amplitude equations.
    *   Matrix/tensor manipulation routines.
    *   Sorting routines (e.g., `quick_sort` in EOM).
*   **Input Data**:
    *   Converged HF solution (MO energies, ERIs in MO basis).
    *   Control parameters (max iterations, threshold).
*   **Output**:
    *   Correlation energies (e.g., E_CCSD, E_CCSD(T)).
    *   Excitation energies (for EOM methods).
    *   Amplitudes (T1, T2) are usually intermediate but might be printed for debugging.
    *   Log output detailing convergence and results.
