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

## 2026-09-22 — This machine's shell env can silently redirect installs to the wrong place

While setting up `nuproj`, `conda run -n nuproj pip install cobaya` did **not** install
into the `nuproj` env — it silently installed cobaya, getdist, and this package into an
unrelated project's virtualenv (`/workspace/imbh-galactic-nuclei/.venv`), because this
session's shell already had that other project's `VIRTUAL_ENV`/`PATH` set ahead of
anything `conda run`/`conda activate` prepends. That pollution has been removed (its
`pip` was used to uninstall the packages it never should have had). Separately, even
with the *correct* interpreter, `~/.local/lib/python3.11/site-packages` (a user-site
directory shared by every Python 3.11 on this machine) was shadowing the conda-forge
`numpy` inside `nuproj` with a different build — a real risk for a package like CAMB
that ships compiled extensions built against a specific NumPy ABI.

**Why this matters:** silent env misdirection could mean a run "succeeds" against the
wrong package versions, or writes into the wrong project entirely, without any error.

**How to apply / fix in place:**
- `nuproj`'s own `activate.d`/`deactivate.d` hooks now set `PYTHONNOUSERSITE=1` while
  the env is active, so `conda activate nuproj` alone fixes the user-site shadowing.
- `conda activate nuproj` on this machine does **not** reliably win the `PATH` race
  against the ambient shell state. After activating, also run:
  `export PATH="$(conda info --base)/envs/nuproj/bin:$PATH"`
  or just call the interpreter by its absolute path,
  `/home/astro/.conda/envs/nuproj/bin/python`, which is what was used to verify this
  setup and what any launch command from this project should use to be safe.
- Before trusting any new-environment install on this machine, verify with
  `python -c "import X; print(X.__file__)"` that the reported file path is actually
  inside the intended env, not a lookalike sitting elsewhere.
