# Decisions log

Scientific and engineering judgment calls made during the project, with reasoning.
Append new entries at the bottom; don't edit past entries except to correct errors.

## 2026-09-22 — CAMB installed via conda-forge, not Cobaya's built-in installer

Cobaya normally manages its own copy of CAMB (and other theory codes/likelihood data)
under a `--packages-path` directory via `cobaya-install`, building from source. Instead,
the `nuproj` environment installs CAMB as a precompiled conda-forge binary alongside
Cobaya (which is pip-only and not on conda-forge).

**Why:** conda-forge's CAMB build requires no local Fortran compiler, so the environment
comes up with zero `apt`/`sudo` steps. Cobaya can be pointed at an already-installed
CAMB (rather than downloading/building its own) via each YAML config's `theory` block.
If a config instead needs Cobaya to manage CAMB itself (e.g., a specific CAMB fork or
patch), switch that config to Cobaya's installer and note it in that config's comments.

**How to apply:** default to the conda-forge CAMB for all configs unless a specific run
needs a non-standard CAMB build.

## 2026-09-22 — No apt packages installed for initial setup

**Why:** the requested package list (cobaya, camb, getdist, numpy, scipy, matplotlib,
mpi4py, pytest, jupyter) is fully satisfiable from conda-forge (precompiled) + one pip
install (cobaya), so no `gfortran`/`build-essential`/`openmpi` were needed. mpi4py's
conda-forge build pulls in its own MPI implementation as a dependency.

**How to apply:** if a later step adds compiled likelihoods (e.g., real Planck/ACT/SPT
`clik`-based codes), those typically need a C/Fortran toolchain and possibly a
system MPI — revisit then with:
`sudo apt update && sudo apt install -y build-essential gfortran libopenmpi-dev openmpi-bin`
(the user runs this; not run automatically).
