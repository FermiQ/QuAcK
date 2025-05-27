# `src/AOtoMO/Hartree_matrix_AO_basis.f90`

## Overview

This file provides subroutines to compute the Hartree matrix (also known as the Coulomb matrix, J) in the Atomic Orbital (AO) basis. The Hartree matrix represents the classical electrostatic repulsion between electrons. Its elements are calculated as J<sub>μν</sub> = Σ<sub>λσ</sub> P<sub>λσ</sub> (μλ|νσ), where P is the density matrix and (μλ|νσ) are the two-electron repulsion integrals (ERIs).

The file contains two implementations:
1.  `Hartree_matrix_AO_basis`: A straightforward implementation using nested loops for summation.
2.  `Hartree_matrix_AO_basis_hpc`: A high-performance computing (HPC) version that utilizes OpenMP for parallelization and assumes a packed storage format for the ERIs to optimize memory access and potentially reduce computational redundancy.

## Key Components

*   **Subroutine `Hartree_matrix_AO_basis`**
    *   Description: Computes the Hartree matrix in the AO basis using a direct four-index summation. It assumes ERIs are provided in a standard 4D array format (μ,λ,ν,σ).
    *   Type: Subroutine
    *   Key Arguments/Parameters:
        *   `nBas` (Integer, Input): The number of atomic basis functions.
        *   `P(nBas,nBas)` (Double Precision, Input): The density matrix P<sub>λσ</sub>.
        *   `ERI(nBas,nBas,nBas,nBas)` (Double Precision, Input): The two-electron repulsion integrals (μλ|νσ) in chemist's notation.
        *   `H(nBas,nBas)` (Double Precision, Output): The computed Hartree matrix J<sub>μν</sub>.
    *   Return Value: None.

*   **Subroutine `Hartree_matrix_AO_basis_hpc`**
    *   Description: An optimized version for computing the Hartree matrix, designed for high-performance computing environments. It employs OpenMP for parallel execution of loops. This version expects the two-electron integrals (`ERI_chem`) to be in a packed (condensed) format, likely utilizing permutation symmetries of the integrals (e.g., (μν|λσ) where μ≥ν, λ≥σ, and the pair index (μν) ≥ (λσ)). The loops and indexing are structured to work with this packed storage.
    *   Type: Subroutine
    *   Key Arguments/Parameters:
        *   `nBas` (Integer, Input): The number of atomic basis functions.
        *   `ERI_size` (Integer\*8, Input): The total number of elements in the packed `ERI_chem` array.
        *   `P(nBas,nBas)` (Double Precision, Input): The density matrix P<sub>λσ</sub>.
        *   `ERI_chem(ERI_size)` (Double Precision, Input): The two-electron repulsion integrals in a packed (condensed) format.
        *   `H(nBas,nBas)` (Double Precision, Output): The computed Hartree matrix J<sub>μν</sub>.
    *   Return Value: None.

## Important Variables/Constants

*   **Inside `Hartree_matrix_AO_basis`:**
    *   `mu, nu, la, si`
        *   Description: Integer loop indices iterating over the atomic basis functions to perform the summation.
        *   Type: Integer
        *   Scope: Local to `Hartree_matrix_AO_basis`.

*   **Inside `Hartree_matrix_AO_basis_hpc`:**
    *   `mu, nu, la, si, nBas8`
        *   Description: Loop indices and the number of basis functions, typed as 8-byte integers for potentially large loop bounds or compatibility with packed index calculations.
        *   Type: Integer\*8
        *   Scope: Local, some are private to OpenMP threads.
    *   `nunu, lala, nula, lasi, numu, ...`
        *   Description: Various 8-byte integer variables used to calculate and store the packed indices for accessing elements within the `ERI_chem` array. These calculations typically involve formulas like `idx = I*(I-1)/2 + J` for unique pairs.
        *   Type: Integer\*8
        *   Scope: Local, some are private to OpenMP threads.
    *   `!$OMP ...` directives
        *   Description: OpenMP compiler directives used to specify parallel regions, loop parallelization, and variable scoping (private, shared) for multi-threaded execution.
        *   Type: Compiler Directive
        *   Scope: Applies to structured blocks of code.

## Usage Examples

```fortran
! Example for Hartree_matrix_AO_basis
IMPLICIT NONE
INCLUDE 'parameters.h' ! Assuming this might define default kinds or constants
INTEGER, PARAMETER :: MAX_BASIS = 100 ! Example
INTEGER :: n_basis_actual
DOUBLE PRECISION :: density_matrix(MAX_BASIS,MAX_BASIS), eri_full(MAX_BASIS,MAX_BASIS,MAX_BASIS,MAX_BASIS)
DOUBLE PRECISION :: hartree_matrix(MAX_BASIS,MAX_BASIS)

! ... (Initialize n_basis_actual, density_matrix, and eri_full for relevant dimensions) ...
n_basis_actual = 50 ! Example actual size

CALL Hartree_matrix_AO_basis(n_basis_actual, density_matrix, eri_full, hartree_matrix)

! hartree_matrix now holds the Coulomb matrix.

! ----------------------------------------------------

! Example for Hartree_matrix_AO_basis_hpc
IMPLICIT NONE
INTEGER :: n_basis_hpc
INTEGER*8 :: eri_packed_array_size
DOUBLE PRECISION, ALLOCATABLE :: density_hpc(:,:), eri_packed(:), hartree_hpc(:,:)

n_basis_hpc = 50 ! Example actual size
eri_packed_array_size = INT(n_basis_hpc * (n_basis_hpc+1) / 2, KIND=8)
eri_packed_array_size = INT(eri_packed_array_size * (eri_packed_array_size+1) / 2, KIND=8)

ALLOCATE(density_hpc(n_basis_hpc, n_basis_hpc))
ALLOCATE(eri_packed(eri_packed_array_size))
ALLOCATE(hartree_hpc(n_basis_hpc, n_basis_hpc))

! ... (Initialize density_hpc and eri_packed) ...

CALL Hartree_matrix_AO_basis_hpc(n_basis_hpc, eri_packed_array_size, density_hpc, eri_packed, hartree_hpc)

! hartree_hpc now holds the Coulomb matrix.
! ... (DEALLOCATE arrays) ...
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   `parameters.h`: Included in `Hartree_matrix_AO_basis`. This file likely contains globally used parameters or type definitions.
*   **External Libraries:**
    *   OpenMP: The `Hartree_matrix_AO_basis_hpc` subroutine uses OpenMP directives (`!$OMP`) for parallel processing, indicating a dependency on an OpenMP-enabled Fortran compiler and runtime library.
    *   Bitwise Intrinsic: The `shiftr` intrinsic function (bitwise right shift) is used in `Hartree_matrix_AO_basis_hpc` for index calculations, which is a standard Fortran intrinsic.
*   **Interactions:**
    *   These subroutines are critical components of Self-Consistent Field (SCF) methods (like Hartree-Fock and Density Functional Theory).
    *   They are called during the construction of the Fock (or Kohn-Sham) matrix, providing the electron-electron Coulomb repulsion term.
    *   The computed Hartree matrix is then used in the subsequent steps of the SCF cycle, typically involving diagonalization of the Fock matrix to obtain new molecular orbitals and energies.

---
*This documentation was auto-generated by a script.*
