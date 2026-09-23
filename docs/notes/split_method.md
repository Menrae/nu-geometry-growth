# Loverde & Weiner's geometry/growth split — how it actually works

Source: `docs/papers/2410.00090v2_loverde-weiner.pdf` (arXiv:2410.00090v2, JCAP 12 (2024) 048).
Read directly: Sec. II (pp.5–17, formalism) and Sec. III (pp.17–31, results) in full, plus a
full-text search for a code release. Everything below cites the page/section/equation it
comes from. **No code was written for this note — reading only, per instructions.**

---

## 1. How the split works, and what each piece controls

There is **no single "split parameter."** Geometry and growth are isolated by *which data
and which nuisance parameter enter a given run*, not by a new physical parameter that lives
in the Boltzmann code itself. Three distinct mechanisms, used in different combinations:

### (a) Geometry-only from the primary CMB — the `Asmear` trick (§II C 3 p.16, §III A 2 p.21–22)

Sample a nuisance parameter `Asmear ~ U(0, 2)` and marginalize over it, **while excluding
any CMB lensing-reconstruction (4-point, `C^κκ_ℓ`) likelihood** from the run.

`Asmear` rescales *only* the degree to which gravitational lensing smears the acoustic peaks
in the primary TT/TE/EE power spectra. Footnote 14 (p.16) is explicit and important here:

> "Typically, the parameter referred to as `A_L` rescales the lensing spectrum consistently
> in its impact on temperature and polarization anisotropies as well as the predicted lensing
> spectrum itself... Here by `Asmear` we denote only the former effect... Marginalizing over
> `Asmear` also yields CMB constraints on late-time cosmology that effectively only derive
> from geometry (i.e., θs)."

So `Asmear` is a **deliberately restricted version of the standard `A_L`/`Alens` nuisance
parameter** already used in Planck-style analyses — restricted to touch *only* the
peak-smearing contribution, explicitly decoupled from the actual lensing-potential
prediction. Marginalizing it out removes essentially all late-time-structure information
from the primary CMB, leaving only `θs` (the geometric quantity) plus the early-time shape
parameters `{ωb, ωc, ns, As, zreion}` that calibrate `θs`'s interpretation (Sec. II A).

**What this controls in the Boltzmann code:** only the *lensed* Cl output used to compare
against TT/TE/EE data — i.e., the interpolation between unlensed and lensed spectra. It does
not touch the background expansion, the transfer function, or the raw lensing-potential
spectrum (`C_ℓ^φφ`) the code computes — that stays available, unscaled, for any lensing
likelihood used elsewhere.

### (b) Growth-only from CMB lensing — the "geometry-free prior" trick (§III B 1, p.26–27)

This is **not a code parameter at all** — it's a likelihood/prior construction done outside
the Boltzmann code:

1. Take (or run) a full Planck-posterior chain.
2. Fit an approximate multivariate Gaussian (mean + covariance) to the marginal posterior
   over `{ωb, ωc, zreion, ns, As}`.
3. **Explicitly drop the `θs` dimension** from that Gaussian (marginalize it out of the
   approximation).
4. Use the resulting 5-parameter Gaussian as an external prior, combined only with actual
   CMB lensing-reconstruction data (ACT DR6 + Planck PR4 lensing) — no primary CMB, no BAO,
   no SNe in this run.

The authors flag this themselves, in their own words (p.27):

> "While this procedure is rather ad hoc, it simply serves as a means to take Planck's
> information on the shape of the CMB anisotropy spectra as a prior for CMB lensing
> likelihoods."

Without `θs` pinning the distance to last scattering, the matter fraction is only weakly
constrained (~5–9% precision) and floats to compensate for neutrino free-streaming
suppression — producing a *much* steeper degeneracy, `Ωm ∝ (1+fν)^8` (p.27), compared to the
`Ωm ∝ (1+fν)^5` degeneracy that holds when `θs` **is** fixed (Eq. 2.20, p.11). This steeper
slope is itself evidence that this construction is doing what it's supposed to: growth (via
lensing) plus weak leftover early-time information, with the strong geometric anchor removed.

**What this controls:** nothing in the code — it's a prior substitution at the likelihood
level, applied *after* the Boltzmann code has already produced its usual outputs.

### (c) BAO and uncalibrated SNe are always pure geometry (§II B, pp.7–12)

Both are background-only observables by construction — distances as a function of redshift,
with no route to constrain perturbation growth at all. They're used with their standard
likelihoods, no special treatment needed. (SNe distances are always marginalized over their
fiducial magnitude, i.e. *not* SH0ES-calibrated, except where explicitly noted, p.18.)

### The neutrino-mass parameter itself

`Mν` (mass sum) is sampled `U(0, 1.5 eV)`, using the **degenerate hierarchy** approximation
(3 equal-mass neutrinos) throughout, "per standard practice" (p.5, footnote-adjacent text) —
except Fig. 2's small-mass illustration, which explicitly uses the normal hierarchy because
the degenerate approximation "underestimates the suppression for the minimum mass sum"
(Fig. 2 caption, p.16).

---

## 2. Boltzmann code, and whether/how it was modified

- **CLASS, not CAMB.** Explicit statement (p.18): *"We use CLASS [108, 109] and model
  nonlinear structure growth with HMCODE-2016 [110] (as used by CLASS, which has yet to
  update to the 2020 version [111])."* CAMB appears only in two footnotes (pp.12, 14) as an
  aside about a shared multipole-approximation convention (`L(ℓ) ≈ ℓ+1/2`) — they did not use
  it.

- **No modification to CLASS's internals is described.** `Asmear` is implemented as an
  ordinary nuisance parameter applied to CLASS's *lensed* Cl output (the same conceptual
  trick as the long-standing `A_L`/`Alens` parameter in CosmoMC-family codes), not a
  patched/forked Boltzmann solver.

- **Sampler: emcee, not Cobaya.** (p.18: *"We employ parameter sampling methods (using
  emcee [112–114]) and likelihood implementations as described in Ref. [58]."*) Ref. [58] is
  a companion methods paper with an overlapping author — Baryakhtar, Simon & Weiner (2024),
  *"Varying-constant cosmology from hyperlight, coupled scalars,"* arXiv:2405.10358 — likely
  where the actual pipeline/likelihood code is documented in more detail.

- **No code release found.** I searched the full PDF text for "github", "zenodo", "gitlab",
  "data availab*", "code availab*", "publicly available", and "public repository" — no hits
  anywhere in the paper. I also checked arXiv:2405.10358's abstract page (the cited methods
  reference) for a code link in its comments field — none found there either. I'm not aware
  of a public repository for this pipeline; worth asking the group directly (see closing
  list) rather than assuming one exists unfound.

**Practical implication for us:** we're reproducing a *method*, not reusing *their* code —
we're on CAMB/Cobaya, they were on CLASS/emcee. That's fine for a qualitative reproduction,
but it means every implementation choice below is ours to make and verify, not a port of an
existing script.

---

## 3. Which datasets inform which channel

| Channel | Datasets used |
| --- | --- |
| **Geometry** (always) | BAO: SDSS (eBOSS DR12/DR16, DR7 MGS), DESI DR1. Uncalibrated SNe: Pantheon, Pantheon+, DES 5YR, Union3. Planck PR3 CMB T/E, marginalized over `Asmear`, lensing excluded. |
| **Growth** (structure) | CMB lensing reconstruction: Planck PR3 lensing, ACT DR6 + Planck PR4 lensing — combined either with the θs-stripped "geometry-free" Planck prior (growth-only, §III B 1) or with the full Planck likelihood (joint geometry+growth, §III B 2). |

BAO and SNe **never** carry growth information in this framework — they are background-only
by construction, full stop. All of the growth information in this paper comes specifically
from CMB lensing reconstruction (the 4-point measurement of `C_ℓ^κκ`), not from any other
structure probe (no galaxy weak lensing, no RSD, no cluster counts are used here).

---

## 4. Key quantitative results to reproduce qualitatively

Not exact numbers — the *patterns*, since we're on a different code/sampler/reduced
likelihood set.

- **Fig. 3 (p.20), geometric tension:** BAO (SDSS, DESI DR1) prefer low/negative-leaning
  `Mν`; uncalibrated SNe (DES 5YR, Pantheon+) prefer `Mν ~ 0.2–0.4 eV`. Both are individually
  consistent with Planck's `θs`-defined line in the `(Ωm, rd√ωm)` plane, but in disjoint
  regions of it — this is the geometric root of the whole tension.

- **Fig. 5 (p.23), geometry-only (`Asmear`-marginalized, no lensing), 95th-percentile `Mν`
  and fraction of posterior below the NO/IO oscillation floors (0.0588/0.099 eV):**
  - **+ DESI DR1** → 0.198 eV (45% / 68%)
  - **+ SDSS** → 0.295 eV (21% / 38%)
  - **+ Pantheon+** → 0.555 eV (4.5% / 8.8%); text separately reports a positive detection,
    `Mν/eV ≈ 0.30 +0.16/-0.15`
  - **+ DES 5YR** → 0.680 eV (0.4% / 1.1%); detection `Mν/eV = 0.44 ± 0.15`, a ~3σ pull away
    from zero driven by DES's larger inferred `Ωm`
  (Column order NO-then-IO verified directly against the rendered figure, not just text
  extraction — smaller floor's tail probability is always the smaller number, as expected.)

- **Fig. 8 (p.27), growth-only (CMB lensing + θs-stripped prior):** steeper degeneracy
  `Ωm ∝ (1+fν)^8`, vs. `(1+fν)^5` when `θs` is held fixed (Eq. 2.20) — direct evidence the
  geometric anchor has been removed. 95th-percentile `Mν ≈ 0.4 eV` from lensing alone.

- **Fig. 9/11 (p.28, 30), joint geometry+growth:** combining both channels *narrows* the
  BAO-vs-SNe disagreement relative to either channel alone. **DES 5YR combo** →
  `Mν/eV = 0.141 +0.089/-0.081` (still nonzero); **DESI combo** → posterior nearly
  incompatible with either mass hierarchy.

- **Headline claim (abstract):** geometry is "at least as important" as structure
  suppression in current neutrino-mass bounds; DESI's tight bound in particular is
  essentially a **geometric** effect, not a growth-suppression effect.

---

## 5. Implementation options in CAMB via Cobaya

**A. Use CAMB's built-in `Alens`-style lensing rescaling + dataset selection only.**
Sample a lensing-smearing rescaling parameter (Cobaya exposes CAMB's `Alens`), exclude any
lensing-reconstruction likelihood for the geometry run; include the lensing-reconstruction
likelihood for the growth run.

- *Difficulty:* low — config-only, no new code.
- *Risk — must verify before trusting this:* whether CAMB's `Alens`, as exposed through
  Cobaya, actually leaves the raw lensing-potential Cl (`Cl['pp']`) **untouched** when a
  lensing likelihood is active in the same run. Footnote 14 requires exactly that decoupling
  (rescale the smearing, not the reported lensing spectrum). If CAMB's `Alens` couples both
  together (which is the *conventional* `A_L` behavior the paper explicitly contrasts itself
  against), this option silently breaks the geometry/growth firewall and both "channels"
  would leak into each other.

**B. Custom Cobaya theory/likelihood wrapper implementing `Asmear` exactly as defined.**
Call CAMB for both unlensed and lensed Cl, manually form
`Cl_smeared = Cl_unlensed + Asmear·(Cl_lensed − Cl_unlensed)` for use against TT/TE/EE data,
while separately exposing the true, unscaled lensing-potential Cl to any lensing likelihood.
Pure-Python wrapper around CAMB's outputs — no CAMB source modification.

- *Difficulty:* medium — requires learning Cobaya's theory-component API and CAMB's
  lensed/unlensed Cl internals.

- *Risk:* moderate but bounded, since this mimics the well-established `Alens`
  interpolation trick from CosmoMC-family codes rather than inventing new physics.

**C. Post-processing / Gaussian-prior-from-chain — a direct reproduction of the paper's own
growth-only method.**
Run one full chain, fit a Gaussian to the posterior over `{ωb, ωc, τ, ns, As}` (dropping the
geometric direction), then run a lensing-only chain using that fit as an external prior
(Cobaya supports arbitrary external priors natively — no custom theory code needed).

- *Difficulty:* low–medium — two ordinary Cobaya runs plus a small prior-fitting script.
- *Risk:* inherits the same Gaussianization/θs-dropping approximation the authors themselves
  called ad hoc; requires an existing full-data chain to bootstrap the prior from.

My inclination: **A for geometry** (cheap, standard, likely fine once verified) **+ C for
growth** (matches the paper's own method, avoids writing custom theory code), with **B** as
the fallback if A's `Alens` behavior doesn't cleanly decouple on inspection.

---

## Things I'm unsure about — for your supervisor

1. Whether CAMB's `Alens` (as exposed through Cobaya) actually decouples peak-smearing from
   the raw lensing-potential Cl the way the paper's `Asmear` does. This needs to be checked
   against CAMB's docs/source before Option A above can be trusted — I have not verified this
   myself and don't want to guess at CAMB internals.
2. They sample `zreion` (reionization redshift) rather than `τ` (optical depth) directly,
   unlike Intertwined's `τ ∈ [0.01, 0.8]`. Does this parametrization choice matter enough for
   reproducing their qualitative results that we should match it, or is `τ` fine?
3. Should we literally reproduce their "ad hoc" Gaussian-prior growth-isolation (Option C),
   or is designing a cleaner method acceptable, given we're after qualitative agreement with
   their trends rather than their exact numbers?
4. They use CLASS + HMCODE-2016 for nonlinear structure; our pipeline is CAMB (default
   nonlinear model is more likely HMcode2020). Is that mismatch acceptable here, or does it
   risk washing out the effect we're trying to isolate?
5. No code release was found for this paper or its cited methods companion
   (arXiv:2405.10358). Is it worth asking the group directly whether Loverde/Weiner's actual
   pipeline exists privately, rather than re-deriving the geometry-free-prior trick (Option C)
   from scratch and risking a subtle implementation bug?
6. Should our configs match Intertwined's explicit `neutrino_hierarchy: degenerate` CAMB
   setting for consistency across both source papers — given this paper's own footnote that
   the degenerate approximation matters "only for mass sums close to the minimum" (p.5)?
