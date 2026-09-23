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

## 2026-09-23 — Reduced-likelihood combination: Planck-lite CMB + tau prior + Planck-2018 lensing (native) + DESI DR2 BAO + Pantheon+

Surveyed what's actually installed in the `nuproj` Cobaya (v3.6.2) for the "reduced
Cobaya+CAMB likelihood" CLAUDE.md calls for, rather than assuming from memory (rule 1).
Checked the package's bundled `likelihoods/` directory directly plus `pip show` for
separately-distributed likelihood packages (HiLLiPoP, Lollipop, Planck PR4 lensing,
ACT DR6 lensing — none of those four are installed).

**Chosen combination:**
- **High-ell CMB:** `planck_2018_highl_plik.TTTEEE_lite_native` — the plik-lite
  compressed spectrum, kept as full TT+TE+EE (not TT-only). "Native" here means a pure
  Python/dataset implementation with no dependency on the compiled `clik` C library,
  which the full (non-lite) Plik likelihood needs along with a ~1.4 GB data file.
  Speed ~200 evaluations/s.
- **Low-ell:** a Gaussian prior on `tau`, `N(0.07, 0.02)`, in place of a real low-l EE
  likelihood (Cobaya ships this exact preset as `gauss_prior` in its `cosmo_input`
  module). At low multipoles the EE signal mainly pins down the reionization optical
  depth `tau`, which is otherwise degenerate with the primordial amplitude `As` (only
  `As * exp(-2*tau)` is well constrained by the rest of the spectrum) — substituting a
  literature Gaussian for the full low-l EE likelihood is a standard fast-MCMC shortcut
  that costs precision on `tau` itself but has negligible effect on `Sigma m_nu`, which
  is what this project actually needs. It also avoids needing the low-l `clik` data
  products.
- **CMB lensing:** `planck_2018_lensing.native` (no `clik` dependency, speed ~50/s).
- **BAO:** `bao.desi_dr2.desi_bao_all` — matches the proposal's DESI DR2 requirement.
- **SNe:** `sn.pantheonplus` — matches the proposal's Pantheon+ requirement.

All five components are bundled with cobaya 3.6.2 already and need no compiled `clik`
library and no extra `pip install` — so no apt/compiler step is needed for this
combination (consistent with the 2026-09-22 "no apt packages" entry above).

**Why the lensing choice matters here specifically:** this project's whole point is
separating geometry from growth, and CMB lensing reconstruction is the main growth-channel
information source in this pipeline (see `docs/notes/split_method.md`). `planck_2018_lensing`
is measurably weaker than the PR4 (Planck) + ACT DR6 lensing combination that both anchor
papers (Intertwined, Loverde & Weiner) actually use — but PR4/ACT DR6 lensing aren't
installed here and come from separate pip packages (`planckpr4lensing`, `act_dr6_lenslike`)
whose build requirements haven't been checked yet. Asked the user directly which to use
given this tradeoff; **decision: start with Planck 2018 lensing (native) now** to get the
full model ladder (LambdaCDM+mnu, w0waCDM+mnu, LambdaCDM+curvature+mnu) running end-to-end
with zero extra install risk, and revisit PR4+ACT DR6 lensing as a planned upgrade once the
pipeline is validated — not a dropped idea.

**How to apply:** all `configs/*.yaml` for this project's reduced-likelihood runs should
use exactly these five components until/unless the PR4+ACT DR6 lensing upgrade happens,
at which point update this entry and add a corresponding note in `split_method.md` about
tightened growth-channel information changing the quantitative (though hopefully not
qualitative) comparison to Table 4.

## 2026-09-23 — First two ladder configs: priors, neutrino/DE treatment, proposal covmat, MPI sizing

Wrote `configs/lcdm_mnu.yaml` and `configs/w0wa_mnu.yaml` (rungs 1 and 2 of the model
ladder) using the five likelihood components from the entry above. Every sampled
parameter's prior is copied directly from Intertwined Table 1 (`docs/notes/
intertwined_summary.md` Sec. 4, arXiv:2607.01226v2 p.7) except two intentional
deviations, both commented in-line in the YAML itself so they're visible at the point
they apply, not just here:

- **`tau`** uses the Gaussian prior `N(0.07, 0.02)` from the likelihood-choice entry
  above, rather than Intertwined's flat `U(0.01, 0.8)` against a real low-ell EE
  likelihood we don't have installed.
- **`nnu` (Neff) is fixed at 3.044**, CAMB/cobaya's standard value accounting for
  non-instantaneous neutrino decoupling, rather than Intertwined's rounded fixed value
  of 3.04. The 0.004 difference is far smaller than anything either likelihood can
  resolve and has no bearing on Sigma m_nu.

Everything else — `ombh2`, `omch2`, `theta_MC_100` (Cobaya's `theta_MC_100` prior of
`[0.5, 10]` for "100 theta_s" turns out to already be the literal convention Intertwined's
Table 1 uses, not just a coincidental match — this is a standard CosmoMC-family
parameterization both share), `logA`, `ns`, `mnu` (`num_massive_neutrinos: 3`,
`neutrino_hierarchy: degenerate`, `U(0,5) eV`), and, for `w0wa_mnu.yaml`, `w`/`wa` (CPL
dark energy, `dark_energy_model: ppf` since `w0` is allowed below -1) — match Table 1
and footnote 3 exactly, param-name-for-param-name against cobaya's own bundled
`cosmo_input` presets (verified by reading `cobaya/cosmo_input/input_database.py`
directly rather than assuming syntax from memory).

**Sampler:** `Rminus1_stop: 0.01` — stricter than Intertwined's own `R-1 <~ 0.02`
(`intertwined_summary.md` Sec. 1), per an explicit request for this project's own runs.
Noted as a deliberate deviation (stricter, not an error) so tighter error bars than
Table 4's aren't later mistaken for a discrepancy.

**Proposal covmat:** both configs seed `mcmc: covmat:` from
`chains/public/desi_dr2_base_mnu_cmb_bao/chain.covmat` — the real public DESI DR2
`base_mnu` posterior downloaded for the warm-up notebook. Its column names
(`ombh2, omch2, theta_MC_100, tau, mnu, logA, ns`) are exactly this project's own
parameter names (both configs deliberately use CAMB's short aliases `ombh2`/`omch2`
rather than cobaya's canonical `omegabh2`/`omegach2`, specifically so this covmat's
header matches verbatim with no renaming needed — confirmed by reading cobaya's covmat
-loading code in `cobaya/sampler.py`, which matches by literal string against each
parameter's name and any declared `renames`). **Why bother:** starting from a real
near-neighbor posterior's covariance means the sampler doesn't have to learn the
parameter degeneracies (e.g. the omch2-theta_MC_100-mnu geometric degeneracy) from
scratch, which is normally what costs the most burn-in time. `w0wa_mnu.yaml`'s `w`/`wa`
aren't in that covmat (the DESI chain is LambdaCDM+mnu only) — cobaya fills only the
shared parameters and falls back to the reference-distribution widths for `w`/`wa`,
logged automatically as "Missing proposal covariance for params ['w', 'wa']" when
`cobaya-run --test` is run (confirmed — see below). **Caveat:** this covmat lives under
the gitignored `chains/public/`, so a fresh clone of this repo needs to regenerate it
first (curl loop in `notebooks/01_warmup.ipynb`'s first cell) before either config's
`--test` will find it.

**Validation:** ran `cobaya-run configs/lcdm_mnu.yaml --test` and `cobaya-run
configs/w0wa_mnu.yaml --test` — both report "Test initialization successful." Cobaya's
own built-in speed-measurement (`[model] Measuring speeds...`, part of `--test`) gives
the actual per-component evaluation rates directly, rather than a hand-rolled timer:

| Component | lcdm_mnu (evals/s) | w0wa_mnu (evals/s) |
|---|---|---|
| `planck_2018_highl_plik.TTTEEE_lite_native` | 844 | 1110 |
| `planck_2018_lensing.native` | 1600 | 4330 |
| `bao.desi_dr2.desi_bao_all` | 576 | 1350 |
| `sn.pantheonplus` | 117 | 90.7 |
| `camb` (Cls from existing transfer functions — fast block: `logA`, `ns`) | 8.13 | 9.23 |
| `camb.transfers` (new transfer functions — **slow block**: `ombh2,omch2,theta,tau,mnu`[,`w`,`wa`]) | **2.64** | **2.47** |

The likelihoods themselves are never the bottleneck (>100/s each); every evaluation is
gated by CAMB's transfer-function recomputation whenever a "slow" parameter changes
(cobaya's own oversampling scheme already exploits this: nuisance parameters like
`A_planck` get 8-9 free re-tries per slow step since they don't need CAMB at all).

**Runtime estimate to `R-1 < 0.01`, given 10 cores (`nproc`) on this machine:**
Rule-of-thumb from CosmoMC/Cobaya practice (not from either anchor paper — this is
generic MCMC-convergence folklore, not a citable published number, and is flagged here
as an estimate, not a fact): reaching `R-1 ~ 0.01` for a similarly-sized posterior
(8 free parameters for `lcdm_mnu`, 10 for `w0wa_mnu`) with a *good* starting covmat
typically takes on the order of 5,000-50,000 accepted samples per chain; at a
Metropolis acceptance rate of roughly 25% (typical when well-tuned), that's
~20,000-200,000 *proposed* slow-block evaluations per chain. Dividing by the measured
slow-block rate above:

- **`lcdm_mnu`** (2.64 evals/s): ~2-21 hours per chain.
- **`w0wa_mnu`** (2.47 evals/s, plus 2 extra free parameters so likely nearer the high
  end of the range): ~2-23 hours per chain, plausibly longer.

Running N chains under MPI does **not** shorten this per-chain estimate — each chain
still needs to individually rack up that many steps for the Gelman-Rubin statistic to
mean anything; more chains just let you obtain that many *independent* chains at the
same wall-clock time (and MPI ranks are what R-1 is computed across in the first
place). **Chose 4 MPI processes** as the default in `scripts/run_chain.sh` usage below:
enough for a robust R-1 estimate, while leaving 6 of this machine's 10 cores free for
other work (notebooks, editing) while a chain runs in the background. Both estimated
runtimes are far past CLAUDE.md's 10-minute auto-launch limit, so neither chain is
started here — only the launch commands are provided, to be run by the user in `tmux`.

**How to apply:** if the real run's `R-1` trace (printed periodically to
`chains/<name>/<name>.log`) is converging much faster or slower than this estimate once
it's actually running, note the observed rate here so future ladder rungs (curvature)
can be estimated from real data instead of the folklore range above.
