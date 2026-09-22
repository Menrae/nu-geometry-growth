# nu-geometry-growth

## Project

Undergraduate research project asking: when extending the cosmological model roughly
doubles the 95% upper bound on the total neutrino mass, is it the **geometric**
constraint or the **growth** constraint that relaxes?

Two anchor papers (full details in `docs/papers/`):
- **Giarè, Lee & Di Valentino, "Intertwined Constraints in Extended Cosmologies"**
  (arXiv:2607.01226, preprint, unreviewed — treat numbers as provisional and recheck
  against the latest version before citing). Their Table 4 shows the mass bound rising
  from ~0.061 eV (LambdaCDM) to 0.113–0.184 eV as dark energy, curvature, and other
  sectors are opened, while sigma8 barely moves.
- **Loverde & Weiner, "Massive neutrinos and cosmic composition"** (arXiv:2410.00090)
  provides the method: split a neutrino mass constraint into what it gets from geometry
  vs. growth of structure.

**Research question:** apply the Loverde & Weiner split to three rungs of the
Intertwined model ladder (LambdaCDM+mass, w0waCDM+mass, LambdaCDM+curvature+mass) and
see which channel relaxes.

- **H1 (core deliverable):** going LambdaCDM → w0waCDM, the geometry-only bound does
  most of the relaxing; the growth-only bound barely moves.
- **H2 (extension):** curvature relaxes the bound through the same geometric channel.
- **H3 (extension):** the supernova sample choice matters in the split even though it
  barely affects the full bound.

Approach: reduced Cobaya+CAMB likelihood (not the full Planck/ACT/SPT-3G likelihoods
Intertwined uses) + DESI DR2 BAO + Pantheon-Plus. The week-4–6 gate is reproducing the
*direction* of the Table 4 trend with this lighter likelihood — if that fails, the
project needs to change before going further.

## Environment

Conda env `nuproj` (Python 3.11): `conda activate nuproj`. Spec lives in
`environment.yml` (conda-forge for camb/getdist/numpy/scipy/matplotlib/mpi4py/pytest/
jupyter; cobaya via pip since it isn't on conda-forge). Rebuild with:
```
mamba env create -f environment.yml
```

**Known gotcha on this machine:** shell sessions here can already have another
project's virtualenv on `PATH`/`VIRTUAL_ENV`, which `conda activate`/`conda run` don't
reliably override. After activating, confirm with `which python` — if it's not under
`.../envs/nuproj/bin`, run
`export PATH="$(conda info --base)/envs/nuproj/bin:$PATH"`, or just call the
interpreter by its absolute path (`.../envs/nuproj/bin/python`). See
`docs/notes/decisions.md` (2026-09-22 entry) for how this was found and why the env's
`activate.d` hook also sets `PYTHONNOUSERSITE=1` to stop a shared `~/.local`
site-packages directory from shadowing conda-forge packages like NumPy.

## Layout

```
configs/       Cobaya YAML run configs — see configs/README.md
src/nuproj/    Python package (pip install -e .)
scripts/       Launch + post-processing scripts — see scripts/README.md
tests/         pytest
notebooks/     exploratory notebooks
docs/notes/    decisions.md — judgment-call log
docs/papers/   source PDFs for every cited published number
figures/       output plots
chains/        gitignored — Cobaya chain output
packages/      gitignored — Cobaya's --packages-path data
```

## Standing rules

1. **Never invent or recall from memory numbers from published papers.** Every
   published value must come from a file in `docs/papers/` and be cited with
   table/section and arXiv version.
2. **All MCMC runs are defined by YAML files in `configs/`.** No ad hoc runs.
3. **Never launch a run expected to take more than 10 minutes.** Write the launch
   command and the user will run it in tmux.
4. **Ask before modifying any Boltzmann code source** (CAMB or otherwise).
5. **Log scientific judgment calls in `docs/notes/decisions.md`** with the reasoning,
   as they're made.
6. **The user is an undergraduate learning this.** Briefly explain non-obvious physics
   or statistics choices when making them.
