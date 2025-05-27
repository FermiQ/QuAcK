# `include/parameters.h` - Global Parameters and Constants

## Overview

The `parameters.h` file is a central header file in the QuAcK codebase. It defines a collection of global integer and double precision parameters that are accessible throughout the Fortran parts of the program. These parameters establish fixed dimensions, limits, physical constants, conversion factors, and mathematical constants used in various modules.

By centralizing these values, `parameters.h` ensures consistency and facilitates easier modification if any of these fundamental quantities need to be updated. Fortran modules typically include this file via an `include 'parameters.h'` statement.

## Key Parameters and Constants

### Integer Parameters (Dimensions and Limits)

*   **`ncart = 3`**:
    *   Description: Number of Cartesian coordinates (x, y, z).
    *   Usage: Defines the dimensionality for spatial vectors and tensor components.
*   **`nspin = 2`**:
    *   Description: Number of spin components (typically alpha and beta, or spin-up and spin-down).
    *   Usage: Used in methods that distinguish between spins (e.g., UHF, GHF) for array dimensions related to spin.
*   **`nsp = 3`**:
    *   Description: Likely related to spin polarization components in UHF (e.g., alpha, beta, total or difference). The exact usage would need to be checked in context where `nsp` is used (e.g. `EJ(nsp)` in UHF).
*   **`maxEns = 4`**:
    *   Description: Maximum number of ensembles or states for a particular type of calculation (e.g. ensemble DFT).
*   **`maxCC = 5`**:
    *   Description: Could refer to a maximum number of coupled-cluster iterations, or a limit in a specific CC routine.
*   **`maxShell = 512`**:
    *   Description: Maximum number of basis function shells.
    *   Usage: Array dimensioning for basis set information.
*   **`maxL = 7`**:
    *   Description: Maximum angular momentum quantum number (L) supported (e.g., L=0 for s, L=1 for p, ..., L=6 for i-functions, L=7 for k-functions).
    *   Usage: Basis set processing and integral calculations.
*   **`n1eInt = 3`**, **`n2eInt = 4`**, **`n3eInt = 3`**, **`n4eInt = 3`**:
    *   Description: These likely refer to the number of indices for different types of integrals or tensors (e.g., 1-electron, 2-electron). For `n2eInt = 4`, this is standard for two-electron repulsion integrals (ijkl). The others might be for specific intermediate tensors.
*   **`maxK = 20`**:
    *   Description: A maximum value, possibly for K-points in periodic calculations, or iterations in a specific solver, or number of vectors in DIIS.

### Double Precision Parameters (Tolerances, Physical and Mathematical Constants)

*   **`threshold = 1d-15`**:
    *   Description: A very small number, likely used as a threshold for numerical comparisons, convergence criteria, or for considering values to be zero.
*   **`pi = acos(-1d0)`**:
    *   Description: The mathematical constant Pi ($\pi$).
*   **`HaToeV = 27.211386245988d0`**:
    *   Description: Conversion factor from Hartrees (atomic unit of energy) to electronvolts (eV).
*   **`pmtoau = 0.0188973d0`**:
    *   Description: Conversion factor from picometers (pm) to atomic units of length (Bohr).
*   **`BoToAn = 0.529177210903d0`**:
    *   Description: Conversion factor from Bohr (atomic unit of length) to Angstroms ($\AA$).
*   **`AnToBo = 1.8897261246257702d0`**:
    *   Description: Conversion factor from Angstroms ($\AA$) to Bohr. (Reciprocal of `BoToAn`).
*   **`auToD = 2.5415802529d0`**:
    *   Description: Conversion factor from atomic units of dipole moment to Debye.

### DFT Functional Constants

*   **`CxLDA = - (3d0/4d0)*(3d0/pi)**(1d0/3d0)`**:
    *   Description: Constant for the Slater exchange functional (Local Density Approximation - LDA). This is $C_x = -
rac{3}{4} (
rac{3}{\pi})^{1/3}$.
*   **`CxLSDA = - (3d0/2d0)*(3d0/(4d0*pi))**(1d0/3d0)`**:
    *   Description: Constant for the Local Spin Density Approximation (LSDA) exchange functional. This appears to be $2 	imes (
rac{3}{4\pi})^{1/3}$, which is related to the spin-polarized form. Note: Standard $X_lpha$ form has $C_X = 
rac{3}{2}(
rac{3}{4\pi})^{1/3}$. The negative sign implies it's part of the energy expression.

## Usage

Fortran modules that require any of these parameters typically include this file at the beginning of their declaration section:

```fortran
implicit none
include 'parameters.h'

! ... rest of the code ...
```

This makes all parameters defined in `parameters.h` directly available within the module.

## Dependencies and Interactions

*   This file is self-contained and does not depend on other files.
*   It is a foundational file included by many, if not most, of the Fortran source files in the `src/` directory that require these predefined constants or dimensional parameters.
