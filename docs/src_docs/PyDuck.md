# `PyDuck.py` - Main Execution Script for QuAcK

## Overview

`PyDuck.py` is the primary Python script for initiating and managing quantum chemistry calculations within the QuAcK package. It serves as the main interface for users to define molecular systems, specify calculation parameters, and launch the core Fortran executable (`QuAcK`).

The script leverages the `pyscf` library to handle initial molecular setup, such as parsing molecular geometries (XYZ files) and basis set information. It then calculates necessary one-electron (overlap, nuclear attraction, kinetic energy, dipole) and two-electron repulsion integrals (ERIs). These integrals and other molecular information are written to specific files in a structured format that the backend Fortran program expects.

Optionally, `PyDuck.py` can integrate with `pyopencap` to incorporate Complex Absorbing Potentials (CAPs) into the calculation, which is useful for studying metastable electronic states (resonances).

After preparing all necessary input data, the script executes the compiled Fortran program (`$QUACK_ROOT/bin/QuAcK`) to perform the main quantum chemistry computations. It also handles the optional generation of Molden files for visualizing molecular orbitals if requested.

## Key Components

### 1. Environment Setup and Argument Parsing
*   **`QUACK_ROOT` Environment Variable**: Checks for the `QUACK_ROOT` environment variable, which should point to the root directory of the QuAcK installation. Exits if not set.
*   **Command-Line Arguments**: Uses the `argparse` module to define and parse a comprehensive set of command-line options. These options allow users to specify:
    *   Molecular geometry (`-x`/`--xyz`)
    *   Basis set (`-b`/`--basis`)
    *   Molecular charge (`-c`/`--charge`)
    *   Spin multiplicity (`-m`/`--multiplicity`)
    *   Units for coordinates (`--bohr`)
    *   Use of Cartesian basis functions (`--cartesian`)
    *   Options for handling two-electron integrals (`--print_2e`, `--formatted_2e`, `--mmap_2e`, `--aosym_2e`)
    *   Freezing core orbitals (`-fc`/`--frozen_core`) - Note: The script currently sets `frozen_core` from args but doesn't seem to directly use it in the `gto.M` call or file writing.
    *   Working directory (`--working_dir`)
    *   Dumping Molden files (`-dm`/`--dump_molden`)

### 2. Molecular Information Processing (with `pyscf`)
*   **Reading XYZ Geometry**: Reads molecular geometry from the specified XYZ file.
*   **PySCF Molecule Creation (`gto.M`)**:
    *   Constructs a `pyscf.gto.Mole` object using the provided geometry, basis set, charge, and spin multiplicity.
    *   Handles basis set definition, loading from a file if `use_cap` is true and `pyopencap` is available, or directly using the input basis string otherwise.
    *   Sets coordinate units (Angstrom or Bohr) and Cartesian/spherical basis functions.
*   **Electron and Orbital Information**:
    *   Extracts the number of alpha and beta electrons.
    *   Writes basic molecular information (number of atoms, electrons) to `input/molecule`.
    *   Calculates and writes the nuclear repulsion energy to `int/ENuc.dat`.
    *   Determines the number of basis functions (`norb`) and writes it to `int/nBas.dat`.

### 3. Integral Calculation and File Writing (with `pyscf`)
*   **One-Electron Integrals**:
    *   Calculates overlap (`int1e_ovlp`), nuclear attraction (`int1e_nuc`), kinetic energy (`int1e_kin`), and dipole (`int1e_r`) integrals using `mol.intor()`.
    *   Writes these matrices to `int/Ov.dat`, `int/Nuc.dat`, `int/Kin.dat`, `int/x.dat`, `int/y.dat`, and `int/z.dat` respectively, using the custom `write_matrix_to_file` function.
*   **Two-Electron Integrals (ERIs)**:
    *   Optionally calculates and writes ERIs to disk if `args.print_2e` is true.
    *   Supports different formats and symmetries:
        *   Formatted text file (`int/ERI.dat`) if `args.formatted_2e` is true (uses `write_tensor_to_file`).
        *   Binary file with 8-fold symmetry (`int/ERI_chem.bin`) if `args.aosym_2e` is true.
        *   Default binary file (`int/ERI.bin`), with an option for memory-mapped file (`args.mmap_2e`) to handle large integral sets. Integrals are transposed from `pyscf`'s "chemist" order (0,1,2,3) to "physicist" order (0,2,1,3) before writing to the binary file if not memory-mapped.
*   **Helper Functions for Writing Data**:
    *   `write_matrix_to_file(matrix, size, file, cutoff=1e-15)`: Writes a 2D matrix to a file in a sparse format (i+1, j+1, value).
    *   `write_tensor_to_file(tensor, size, file_name, cutoff=1e-15)`: Writes a 4D tensor (ERIs) to a file in a sparse format (i+1, j+1, k+1, l+1, value).

### 4. CAP Integration (with `pyopencap`)
*   **Conditional Import**: Attempts to import `pyopencap`. Sets `use_cap` flag accordingly.
*   **CAP Setup**: If `pyopencap` is available and CAP data for the molecule exists in `cap_data/`:
    *   Reads CAP parameters (onset coordinates, eta_opt).
    *   Creates a Psi4-formatted basis string using `create_psi4_basis` for `pyopencap`.
    *   Initializes `pyopencap.System` and `pyopencap.CAP` objects.
    *   Calculates AO CAP matrix and writes it to `int/CAP.dat`.
    *   `create_psi4_basis(basis_dict)`: Converts `pyscf` basis dictionary to a Psi4-formatted string, saved to `input/basis_psi4`.

### 5. Execution of Fortran Backend
*   **Memory Management**: Deletes large integral arrays and calls `gc.collect()` to free memory before running the Fortran code.
*   **Subprocess Call**: Executes the main QuAcK Fortran program located at `$QUACK_ROOT/bin/QuAcK`, passing the `working_dir` as an argument.

### 6. Molden File Generation
*   **Optional Output**: If `args.dump_molden` is true:
    *   Reads real and imaginary parts of MO coefficients from `real_MOs.dat` and `imag_MOs.dat` (presumably generated by the Fortran code).
    *   Uses `pyscf.tools.molden.from_mo` to generate `real_MOs.molden` and `imag_MOs.molden`.

## Important Variables/Constants (Command-Line Arguments)

*   `QUACK_ROOT`: Environment variable specifying the project's root directory.
*   `-b`/`--basis` (str, required): Path to the basis set file.
*   `-x`/`--xyz` (str, required): Name of the XYZ geometry file (without extension).
*   `-c`/`--charge` (int, default: 0): Total charge of the molecule.
*   `-m`/`--multiplicity` (int, default: 1): Spin multiplicity.
*   `--bohr` (flag, default: Angstrom): Specifies if XYZ coordinates are in Bohr.
*   `--cartesian` (flag, default: False): Use Cartesian basis functions instead of spherical.
*   `--print_2e` (bool, default: True): Whether to calculate and write 2-electron integrals.
*   `--formatted_2e` (flag, default: False): Write 2e integrals in formatted text; otherwise binary.
*   `--mmap_2e` (flag, default: False): Use memory-mapping for 2e integrals (reduces RAM usage).
*   `--aosym_2e` (flag, default: False): Use 8-fold symmetry for 2e integrals.
*   `--working_dir` (str, default: `$QUACK_ROOT`): Directory where calculations are run and files are stored.
*   `-dm`/`--dump_molden` (flag, default: False): Generate Molden files for MO visualization.

## Usage Examples

A typical command-line invocation might look like:

```bash
export QUACK_ROOT=/path/to/quack_directory
python PyDuck.py --xyz MyMolecule --basis cc-pVDZ --charge 0 --multiplicity 1
```

This command would run a calculation for `MyMolecule.xyz` using the `cc-pVDZ` basis set, for a neutral singlet state.

## Dependencies and Interactions

*   **External Libraries**:
    *   `pyscf`: Used for molecule definition, basis set parsing, and integral calculations.
    *   `numpy`: Used for numerical operations, especially array handling.
    *   `pyopencap` (optional): Used for Complex Absorbing Potential calculations.
*   **Internal Components**:
    *   Calls the main Fortran executable: `$QUACK_ROOT/bin/QuAcK`.
*   **File System Interactions**:
    *   Reads:
        *   XYZ geometry file (e.g., `$QUACK_ROOT/mol/MyMolecule.xyz`).
        *   Basis set file (e.g., `$QUACK_ROOT/basis/cc-pVDZ`).
        *   CAP data file (e.g., `$QUACK_ROOT/cap_data/MyMolecule.xyz`), if `use_cap` is active.
        *   `real_MOs.dat`, `imag_MOs.dat` (if `--dump_molden` is true, expects these from Fortran code).
    *   Writes:
        *   `input/molecule`: Basic molecular parameters.
        *   `input/basis_psi4`: Psi4-formatted basis for `pyopencap`.
        *   `input/eta_opt.dat`: CAP eta parameter.
        *   `int/ENuc.dat`: Nuclear repulsion energy.
        *   `int/nBas.dat`: Number of basis functions.
        *   `int/Ov.dat`, `int/Nuc.dat`, `int/Kin.dat`, `int/x.dat`, `int/y.dat`, `int/z.dat`: One-electron integrals.
        *   `int/CAP.dat` (if `use_cap` is active): CAP matrix in AO basis.
        *   `int/ERI.dat` or `int/ERI.bin` or `int/ERI_chem.bin`: Two-electron integrals.
        *   `real_MOs.molden`, `imag_MOs.molden` (if `--dump_molden` is true).
    *   Creates directories: `input/`, `int/` within the `working_dir`.
    *   Deletes files: Clears previous integral files before writing new ones. Removes `real_MOs.dat`, `imag_MOs.dat` after generating Molden files.
