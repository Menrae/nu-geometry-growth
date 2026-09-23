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

## 2026-09-23 — Ordering diagnostics use GetDist's boundary-corrected KDE, not a custom one

`src/nuproj/diagnostics.py` needs a smooth density estimate for the cosmological posterior
on Sigma m_nu (for the overlap coefficient and the Bayes-factor's Z_H integral). Sigma m_nu
has a hard physical prior boundary at 0, and the oscillation-informed reference distribution
has its own hard boundary at the ordering's mass floor (~0.059 eV for NO, ~0.10 eV for IO) —
an uncorrected KDE would smear probability mass across both boundaries and bias the estimate
right where these diagnostics are most sensitive (near the floors).

**Why GetDist's own KDE instead of writing one:** `MCSamples.get1DDensity(name, ...)`, given
a hard `ranges` limit, already applies the Jones (1993)/Jones & Foster (1996) linear boundary
kernel correction by default (`boundary_correction_order=1`) — confirmed by reading GetDist's
own source (`getdist/mcsamples.py::get1DDensityGridData`). This is the same, already-vetted
machinery the cosmology community uses for exactly this class of problem (e.g. Planck's own
posteriors on `tau`, `r`, `Sigma m_nu`), and GetDist is already a project dependency. Writing
and validating a from-scratch boundary-corrected KDE (e.g. the reflection method) would
duplicate well-tested code for no accuracy benefit, and risks a subtly-wrong implementation
that's harder to trust than reusing GetDist's.

**Note:** `.get1DDensity(...).P` is normalized to peak=1 (GetDist's plotting convention), not
to integrate to 1 — `diagnostics.py` always renormalizes onto its own evaluation grid via
`np.trapezoid` before using it as a probability density. Also note this differs slightly from
Intertwined's own footnote 12, which describes a plain weighted *histogram* for `P_cosmo`
(no smoothing) in their implementation of the overlap coefficient; we use a smooth
boundary-corrected KDE instead for both `P_cosmo` and `P_osc^H`, which avoids histogram
bin-width sensitivity and is more standard practice — validated in `tests/test_diagnostics.py`
against independent analytic ground truths (truncated-Gaussian CDFs/quantiles via `scipy.stats`,
and an overlap/Bayes-factor check against `scipy.integrate.quad` on the exact pdfs), not just
against itself.

**How to apply:** any future density estimate on a physically-bounded cosmological parameter
in this project should set a hard `ranges` limit on the `MCSamples` object rather than
building an unbounded KDE and hoping the tails are small.

## 2026-09-23 — Oscillation-informed reference distribution built by Monte Carlo, not inversion

`oscillation_reference_samples()` draws the lightest-mass prior (`U(0, 2) eV`) and maps it
through the exact NO/IO mass relation (Eqs. 5.4-5.5) via direct simulation, then feeds the
resulting Sigma m_nu samples through the same boundary-corrected KDE as the cosmological
chain, rather than analytically inverting the mass relation to get a closed-form density.

**Why:** the forward map is monotonic but the closed-form inverse and its Jacobian, while
not hard, add a second, differently-shaped implementation to maintain and re-verify; sampling
and reusing the same KDE path as `P_cosmo` keeps the two distributions treated identically
end-to-end (same boundary correction, same interpolation grid) and is what
`tests/test_diagnostics.py::test_log_bayes_factor_matches_analytic_integral` checks *against*
(that test builds the analytic inversion independently, specifically so the Monte-Carlo path
in the library proper is verified by a genuinely separate calculation, not by itself).

**How to apply:** if `log_bayes_factor`/`overlap_coefficient` ever need to run inside a tight
loop (e.g. per-model-per-dataset in a later analysis script), consider caching
`oscillation_reference_samples()`'s output per hierarchy rather than redrawing it every call —
it doesn't depend on the cosmological chain at all.

## 2026-09-23 — Warm-up chain choice: DESI DR2 `base_mnu`, Planck2018-Plik CMB variant

Downloaded the DESI DR2 `base_mnu` chain into `chains/public/desi_dr2_base_mnu_cmb_bao/`
(gitignored; regenerate with the `curl` loop documented in `notebooks/01_warmup.ipynb`'s
first cell, or browse the source directory at
`https://data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/cobaya/base_mnu/`) to run
`notebooks/01_warmup.ipynb` against. Three CMB-likelihood variants exist for this model
(differing in which high-ℓ Planck likelihood is used: Plik, NPIPE-CamSpec, or Hillipop).
Picked the **Plik TTTEEE** variant specifically because it's the same underlying Planck 2018
Plik likelihood Intertwined itself uses (their Table 1/§2.2 dataset list), making this the
closest available public chain to Intertwined's own `CMB+DESI` baseline for a qualitative
comparison, even though it still differs from Intertwined's exact combination (their CMB
likelihood also adds ACT DR6, SPT-3G and BICEP/Keck; this chain has none of those; and it
has no supernova sample at all, whereas Intertwined always includes one).

**Why this matters for interpretation:** the resulting bound (~0.069 eV) is looser than
Intertwined's ~0.061 eV baseline, and every downstream diagnostic (tail probability,
overlap, Bayes factor) shifts in the direction that a looser bound implies. That's expected
from the missing datasets, not a bug -- see the notebook's closing markdown cell for the
full reasoning per diagnostic.
