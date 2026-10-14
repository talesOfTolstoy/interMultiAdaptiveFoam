# interMultiAdaptiveFoam

`interMultiAdaptiveFoam` is an OpenFOAM solver for transient, multi-region, incompressible two-phase simulations with heat transfer, phase change, and adaptive mesh refinement.

It couples one or more fluid regions with one or more solid regions, advances the fluid and solid energy equations, and includes custom interface, surface-tension, and adaptive-mesh components.

This is a research and development codebase. Validate the solver against your own benchmark cases before using it for production decisions or published results.

## Features

- Multi-region fluid/solid solution loop.
- Incompressible two-phase volume-of-fluid (VOF) interface tracking.
- Fluid and solid energy equations with phase-change source terms.
- Constant and temperature-dependent surface-tension models.
- Contact-angle boundary-condition support.
- AMR indicator based on interface cells, velocity-gradient information, and the fluid-side solid boundary.
- Custom 2D and axisymmetric dynamic-refinement implementations (adapted from Takehiro Shiraishi and Ahmad Baniabedalruhman's code).
- PISO pressure–velocity coupling and automatic time-step control.

## Requirements

- A working OpenFOAM installation with `wmake` available on `PATH`.
- A compatible C++ compiler and, for parallel runs, an MPI toolchain.
- A multi-region OpenFOAM case with the required meshes, fields, properties, numerical schemes, and solver controls.

The custom AMR sources are kept in `AMR/` and are built against the OpenFOAM API available in the selected installation. The repository does not currently declare a single supported OpenFOAM release or include an automated compatibility matrix, so check the build and a representative case before using a different release.

## Build

Source your OpenFOAM environment, move to this directory, and run the build script:

```bash
source /path/to/OpenFOAM/etc/bashrc
cd /path/to/interMultiAdaptiveFoam
./Allmake
```

`Allmake` builds the components in this order:

1. `interfaceProperties/` → `libinterfaceProperties`
2. `AMR/dynamicMesh/` → `libdynamicMesh_2D+Axi`
3. `AMR/dynamicFvMesh/` → `libdynamicFvMesh_2D+Axi`
4. the project root → `interMultiAdaptiveFoam`

The libraries and executable are installed into the OpenFOAM user locations (`$FOAM_USER_LIBBIN` and `$FOAM_USER_APPBIN`). To clean all four build targets, run:

```bash
./Allwclean
```

## Running a case

The executable expects a standard multi-region OpenFOAM case. At minimum, the case must define the fluid and solid regions in `constant/regionProperties` and provide the corresponding meshes, thermophysical properties, initial fields, numerical schemes, and solver controls.

Run from the case directory:

```bash
interMultiAdaptiveFoam
```

For parallel execution, use the normal OpenFOAM workflow with `decomposePar`, `mpirun`, and `reconstructPar`.

### Adaptive mesh refinement

The custom AMR library contains separate implementations for ordinary 2D and axisymmetric refinement:

- `dynamicRefineFvMesh2D`
- `dynamicRefineFvMeshAxi`

Each implementation reads its refinement controls from a case `dynamicMeshDict`. The axisymmetric implementation additionally requires a positive `targetCellSize`. The solver writes an `amrIndicator` field; cases should configure the refinement field and thresholds to match the intended physics and mesh.

## Repository layout

```text
.
├── interMultiAdaptiveFoam.C  # Solver entry point
├── fluid/                    # VOF, momentum, pressure, energy, and AMR indicator code
├── solid/                    # Solid-region energy equation code
├── interfaceProperties/      # Interface and surface-tension library
├── AMR/                      # Custom 2D and axisymmetric refinement libraries
├── include/                  # Shared time-step and region utilities
├── Make/                     # OpenFOAM build configuration
├── Allmake                   # Build all libraries and the solver
└── Allwclean                 # Clean all project build products
```

Generated OpenFOAM platform directories, object files, dependency files, and `lnInclude` trees are excluded by `.gitignore`.

## Status and limitations

This repository contains solver and library source code, not a complete validation case or automated regression suite. Results should be checked for conservation, interface capturing, phase-change behaviour, and mesh-adaptation quality for each application.

## License and attribution

Copyright (C) 2026 Borja Diaz-Guardamino Sobrin for the integration and
original contributions in this repository.

This project is licensed under the GNU General Public License, version 3 or
later. See [LICENSE](LICENSE). The repository incorporates OpenFOAM-derived
code and other upstream components; their existing copyright and license
notices remain applicable and must be preserved.

