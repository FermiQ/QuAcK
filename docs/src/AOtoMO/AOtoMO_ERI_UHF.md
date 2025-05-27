# `src/AOtoMO/AOtoMO_ERI_UHF.f90`

## Overview

This file contains the `AOtoMO_ERI_UHF` subroutine, which is responsible for transforming two-electron repulsion integrals (ERIs) from the Atomic Orbital (AO) basis to the Molecular Orbital (MO) basis. This specific version is designed for Unrestricted Hartree-Fock (UHF) wavefunctions, where different spatial orbitals are used for alpha and beta spin electrons. The transformation is performed using a semi-direct O(N<sup>5</sup>) algorithm. The routine allows specifying the spin channels for the "bra" and "ket" sides of the integrals, enabling the transformation of various spin combinations of ERIs.

## Key Components

*   **Subroutine `AOtoMO_ERI_UHF`**
    *   Description: Transforms two-electron integrals from AO to MO basis for UHF wavefunctions. It handles different spin combinations by taking `bra` and `ket` spin channel specifiers as input. The comment `! bra and ket are the spin of (bra|ket) = (bra bra|ket ket) = <bra ket|bra ket>` clarifies how these spin specifiers relate to the integral indices.
    *   Type: Subroutine
    *   Key Arguments/Parameters:
        *   `bra` (Integer, Input): Specifies the spin channel (e.g., 1 for alpha, 2 for beta, as defined in `parameters.h`) for the first two MO indices (the "bra" side).
        *   `ket` (Integer, Input): Specifies the spin channel for the last two MO indices (the "ket" side).
        *   `nBas` (Integer, Input): The number of basis functions. In the context of the `dgemm` calls and output array `ERI_MO`, this argument seems to represent the number of molecular orbitals for the specified spin channel rather than the total number of AOs, which could be a common practice if `c` is already subsetted or if `nBas` here means number of MOs.
        *   `c(nBas,nBas,nspin)` (Double Precision, Input): The UHF MO coefficient matrix. The third dimension `nspin` (likely 2) distinguishes between alpha and beta spin orbitals. The routine uses `c(1,1,bra)` and `c(1,1,ket)` for the transformations. The actual dimensions passed to dgemm suggest `c` stores orbitals of a specific spin `bra` or `ket`.
        *   `ERI_AO(nBas,nBas,nBas,nBas)` (Double Precision, Input): The array of two-electron integrals in the AO basis.
        *   `ERI_MO(nBas,nBas,nBas,nBas)` (Double Precision, Output): The array where the transformed two-electron integrals in the MO basis will be stored. The dimensions here should correspond to the number of MOs for the respective spin channels.
    *   Return Value: None.

## Important Variables/Constants

*   **`scr(:,:,:,:)`**
    *   Description: An allocatable 4-dimensional array used as a scratch workspace to store intermediate results during the four-step transformation of the integrals.
    *   Type: Double Precision, Allocatable 4D Array
    *   Scope: Local to Subroutine `AOtoMO_ERI_UHF`.
*   **`mu, nu, la, si`**
    *   Description: Integer variables, typically used as loop indices for AO basis functions. While not explicitly shown in loops in this snippet, they are conventional names for such indices.
    *   Type: Integer
    *   Scope: Local to Subroutine `AOtoMO_ERI_UHF`.
*   **`i, j, k, l`**
    *   Description: Integer variables, typically used as loop indices for MO basis functions.
    *   Type: Integer
    *   Scope: Local to Subroutine `AOtoMO_ERI_UHF`.
*   **`nspin`** (from `parameters.h`)
    *   Description: A parameter likely defined in `parameters.h` that specifies the number of spin channels (usually 2 for alpha and beta).
    *   Type: Integer (presumably)
    *   Scope: Global via include.

## Usage Examples

```fortran
! Example of how to call SUBROUTINE AOtoMO_ERI_UHF
IMPLICIT NONE
INCLUDE 'parameters.h' ! Assuming nspin, spin_alpha, spin_beta are defined here
INTEGER :: NBAS_AO, NMO_ALPHA, NMO_BETA
DOUBLE PRECISION, ALLOCATABLE :: C_UHF_MO_COEFFS(:,:,:), ERI_AO_BASIS(:,:,:,:), ERI_MO_AABB(:,:,:,:)

! ... (Initialize NBAS_AO, NMO_ALPHA, NMO_BETA)
! ... (Allocate and fill C_UHF_MO_COEFFS(NBAS_AO, NMO_ALPHA for spin_alpha, NBAS_AO, NMO_BETA for spin_beta))
! ... (Allocate and fill ERI_AO_BASIS(NBAS_AO,NBAS_AO,NBAS_AO,NBAS_AO))
! ... (Allocate ERI_MO_AABB(NMO_ALPHA,NMO_ALPHA,NMO_BETA,NMO_BETA))

! To transform (alpha alpha | beta beta) type integrals:
! ERI_MO_AABB(i_alpha, j_alpha, k_beta, l_beta)
! The 'nBas' argument to AOtoMO_ERI_UHF should match the MO dimension for the given spin.
! This example assumes a convention where the ERI_MO output is always (dim1,dim1,dim2,dim2) based on bra/ket spins
! and the dgemm calls are structured to make this happen by passing appropriate sections of C.
! The actual implementation details of how c is indexed (c(AO, MO, spin_idx)) and how ERI_MO is shaped
! (MO_spin1, MO_spin1, MO_spin2, MO_spin2) or (MO_spin1, MO_spin2, MO_spin3, MO_spin4) needs careful handling.
! The provided code's ERI_MO(nBas,nBas,nBas,nBas) and c(nBas,nBas,nspin) with dgemm calls
! using c(1,1,bra) and c(1,1,ket) implies a specific structure for the transformed ERI_MO.
! For (alpha alpha | beta beta), the first two dgemms would use c(:,:,spin_alpha) and
! the last two would use c(:,:,spin_beta) if following the pattern (ab|cd) -> C_a C_b ERI_AO C_c C_d.
! The provided code does: C_bra ERI_AO -> scr1 ; C_ket scr1 -> ERI_MO_step2; C_bra ERI_MO_step2 -> scr2; C_ket scr2 -> ERI_MO_final
! This transforms indices k, l, i, j in order -> (ij|kl)_MO = sum_pqrs C_pi C_qj (pq|rs)_AO C_rk C_sl
! Or rather, (MO1,MO2,MO3,MO4) = C(AO,MO1,spin_bra) C(AO,MO2,spin_ket) C(AO,MO3,spin_bra) C(AO,MO4,spin_ket) * ERI_AO
! The code appears to do:
! ERI_MO(i,j,k,l) = sum_mu,nu,la,si c(mu,i,bra) c(nu,j,ket) c(la,k,bra) c(si,l,ket) * ERI_AO(mu,nu,la,si)
! (This is based on the dgemm calls transforming one index at a time, with specific spin MO coefficients)

! If nBas_for_routine is num_MOs_alpha for bra=alpha, num_MOs_beta for bra=beta etc.
! And ERI_MO is (nMOs_bra, nMOs_ket, nMOs_bra, nMOs_ket)
! CALL AOtoMO_ERI_UHF(spin_alpha, spin_beta, nMO_alpha_or_beta, C_UHF_MO_COEFFS, ERI_AO_BASIS, ERI_MO_RESULT)
! This part is complex due to the indexing and meaning of 'nBas' in the routine.
! A safer example would be:
CALL AOtoMO_ERI_UHF(spin_alpha, spin_alpha, NMO_ALPHA, C_UHF_MO_COEFFS, ERI_AO_BASIS, ERI_MO_AAAA)
CALL AOtoMO_ERI_UHF(spin_beta, spin_beta, NMO_BETA, C_UHF_MO_COEFFS, ERI_AO_BASIS, ERI_MO_BBBB)
CALL AOtoMO_ERI_UHF(spin_alpha, spin_beta, NMO_ALPHA, C_UHF_MO_COEFFS, ERI_AO_BASIS, ERI_MO_AABB) ! This might need ERI_MO(NMO_ALPHA,NMO_ALPHA,NMO_BETA,NMO_BETA)

```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   `parameters.h`: This include file likely provides definitions for `nspin` and integer constants representing spin channels (e.g., `spin_alpha`, `spin_beta`).
*   **External Libraries:**
    *   `dgemm`: A standard routine from BLAS (Basic Linear Algebra Subprograms) used for general matrix-matrix multiplication. The four calls to `dgemm` implement the O(N<sup>5</sup>) transformation by sequentially transforming each of the four AO indices to the MO basis, using the appropriate spin-orbital coefficient matrix (`c(:,:,bra)` or `c(:,:,ket)`) at each step.
*   **Interactions:**
    *   This subroutine is a fundamental component for many post-UHF quantum chemistry methods that require ERIs in the MO basis with specific spin configurations. Examples include:
        *   Unrestricted Møller-Plesset perturbation theory (UMP2, UMP3, etc.).
        *   Unrestricted Coupled Cluster methods (UCCD, UCCSD, etc.).
        *   Other methods that explicitly use spin-orbital based ERIs.
    *   It is typically called after a UHF calculation has converged and the MO coefficients `c` are available.

---
*This documentation was auto-generated by a script.*
