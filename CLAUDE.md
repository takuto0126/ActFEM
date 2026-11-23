# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is ActFEM

ActFEM is a 3D finite element method (FEM) code written in Fortran for geophysical electromagnetic inversion. It supports:
- **Active source EM** (controlled-source, Bz amplitude/phase responses)
- **Magnetotelluric (MT)** forward and inversion (impedance tensor, tipper)
- **Joint inversion** of active source + MT + tipper simultaneously

The primary target application is volcanic structure imaging (Aso volcano, Japan).

## Build system

Each subdirectory under `src/` has its own `Makefile`. Before compiling, source the Intel oneAPI environment:

```bash
source /opt/intel/oneapi/setvars.sh --force
```

Build an executable by `cd`-ing into the relevant `src/` subdirectory and running `make`. Clean with `make clean`.

**Key executables and where to build them:**

| Executable | Source dir | Purpose |
|---|---|---|
| `ebfem_bxyz.exe` | `src/solver/` | Active source forward modeling |
| `n_ebfem_3DMT.exe` | `src/src_3DMT/` | MT 3D forward modeling |
| `ebfem_inv_ap.exe` | `src/src_inv_ap/` | Active source inversion (MPI) |
| `ebfem_inv_joint.exe` | `src/src_inv_joint/` | Joint inversion ACTIVE+MT+Tipper (MPI) |
| `meshgen1.exe`, `meshgen2.exe`, `mkline.exe`, `mkface.exe`, etc. | `src/src_mesh/` | Mesh generation tools |
| `slice.exe`, `change_model2cond.exe`, `modmodel.exe`, `n_model2cond.exe` | `src/src_post/` | Post-processing |

**Compilers:**
- Serial: `ifx` (Intel Fortran)
- MPI: `mpiifort` (Intel MPI + Fortran) for joint/AP inversion
- Cluster ISM variant uses `Makefile_ISM` with different MKL link flags

## Running executables

Executables read a control file via stdin:

```bash
# Forward active source
export OMP_NUM_THREADS=8
time ${src}/ebfem_bxyz.exe <<EOF
active.ctl
EOF

# MT forward
time ${src}/n_ebfem_3DMT.exe <<EOF
mt.ctl
EOF

# Joint inversion (MPI, ijoint: 1=ACTIVE only, 2=MT only, 3=both)
export OMP_NUM_THREADS=4
time mpiexec -n 6 ${src}/ebfem_inv_joint.exe <<EOF | tee result_inv/inv.log
3
active.ctl
mt.ctl
joint.ctl
EOF
```

See the `run.sh` scripts in `Joint_test/` subdirectories for working examples.

## Control file format (`.ctl`)

Lines starting with `##` are comments. Parameter values follow the `!` delimiter on each line. The label before `!` is human-readable documentation only; the parser reads positional order.

```
iflag_map          !2          ← value is "2"
UTM zone           !52S
```

## Mesh generation workflow

Mesh generation uses `meshgen.sh` scripts in mesh directories (e.g., `Joint_test/mesh_joint/`):

1. `meshgen1.exe` — generates 2D background mesh `.pos` file
2. `gmsh -2 nakadake2d.geo -bgm nakadake2d.pos -format msh2` — 2D Gmsh triangulation
3. `meshgen2.exe` — generates 3D background mesh `.pos` file
4. `gmsh -3 nakadake3d.geo -bgm nakadake3d.pos -format msh2` — 3D Gmsh tetrahedralization
5. `mkline.exe` — generates `lineinfo.dat` (edge information)
6. `mkface.exe` — generates `faceinfo.dat` (face information, needed for MT/joint)

The mesh control file (e.g., `aso.ctl`) specifies geometry bounds, element sizes (`sizein`, `sizebo`), topography files, and observatory/source locations.

## Source code architecture

**`src/common/`** — Shared Fortran modules included via `VPATH` in all Makefiles:
- `m_param.f90` — `param_forward`, `param_cond`, `param_source` types; central parameter container
- `m_mesh_type.f90` — `mesh` type for tetrahedral 3D mesh
- `m_matrix.f90` — Sparse matrix types (`real_crs_matrix`, `complex_crs_matrix`, `global_matrix`)
- `m_modelpart.f90` — `model`, `modelpara` types for conductivity model parameterization
- `m_constants.f90` — Physical constants
- `m_outresp.f90` — Response data types (`respdata`, `respmt`, `resptip`)

**`src/common_mpi/`** — MPI communication modules:
- `m_shareformpi_mt.f90` / `m_shareformpi_joint.f90` — frequency-domain MPI distribution

**`src/solver/`** — Core edge-based FEM solver (active source):
- Uses PARDISO (Intel MKL) for sparse linear systems
- COCG iterative solver as alternative

**`src/src_2D/`** — 2D TM mode solver used as boundary conditions for MT

**`src/src_3DMT/`** — 3D MT forward; calls 2D solver for each surface

**`src/src_inv_joint/`** — Main joint inversion (most actively developed):
- `n_inv_joint.f90` — Main program entry point
- `forward_joint.f90` (module `forward_joint_inv`) — Runs FEM forward for ACTIVE and/or MT at one frequency
- `m_jacobian_joint.f90` — Jacobian computation for both data types
- `m_param_jointinv.f90` — `param_joint` type; inversion hyperparameters and data file pointers
- `m_modelroughness_joint.f90` — Regularization: Smoothest Model (SM, type=1), Minimum Support (MS, type=2), Minimum Support Gradient (MSG, type=3)
- `m_solveCM_ap.f90` — Data-space solver (Conjugate Gradient on normal equations)
- `solvePARDISOjointinv.f90` — Model-space solver using PARDISO
- `m_freq_mpi_joint.f90` — MPI frequency distribution

**Inversion hyperparameters (`ialphaflag` in `joint.ctl`):**
- 1: L-curve method
- 2: Cooling strategy with `alpha_init` and `factor` (α ← α × 10^factor)
- 3: Minami cooling (data-space)
- 4: Modified Grayver cooling

## Test directories

- `Fwd_test/` — Standalone forward modeling tests (flat mesh, MPI variants)
- `Joint_test/` — Full workflow tests including mesh generation, forward runs, and inversion:
  - `mesh_joint/`, `mesh_light/` — Meshes for full/light test cases
  - `structure/`, `structure_light/` — Initial conductivity models
  - `fwd_active/`, `fwd_3DMT/` — Forward results used as synthetic data for inversion tests
  - `inv_act/`, `inv_mt/`, `inv_joint/` — Inversion tests (active-only, MT-only, joint)
  - `inv_*_icomb-1/` — Variants using `icombine = -1` (inner/outer block parameterization)

## Data file conventions

- Active source: amplitude (`.dat`) and phase (`.dat`) files per observatory per source
- MT impedance: complex impedance tensor files, errors separate
- Tipper: `Tx`, `Ty` components
- Conductivity model: `.msh` files with per-element conductivity values
- `condflag = 0`: homogeneous half-space; `condflag = 1`: conductivity from file

## Post-processing

```bash
# Extract a depth slice from model
slice.exe <<EOF
slice.ctl
EOF

# Convert model parameters to conductivity
change_model2cond.exe <<EOF
change.ctl
EOF
```

`view_slice.sh` in test directories shows typical GMT plotting workflows.
