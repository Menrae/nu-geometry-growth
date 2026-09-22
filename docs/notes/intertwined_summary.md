# Intertwined (Giarè, Lee & Di Valentino) — extracted results and methodology

Source: `docs/papers/2607.01226v2_intertwined.pdf` (arXiv:2607.01226v2, submitted
10 Aug 2026 — v1 was 1 Jul 2026). Also referenced: `docs/papers/2410.00090v2_loverde-weiner.pdf`
(Loverde & Weiner, arXiv:2410.00090v2, JCAP 12 (2024) 048, 16 Dec 2024 — v1 was 30 Sep 2024).

Everything below is transcribed directly from the PDFs (page/section/table cited inline).
Nothing here is from memory of these papers.

---

## 1. Table 4 — neutrino-sector results (p.21)

All models use CAMB + Cobaya, Gelman-Rubin R−1 ≲ 0.02 convergence. **Free-Σmν runs use
`num_massive_neutrinos: 3`, `neutrino_hierarchy: degenerate`** (three equal-mass neutrinos) —
see §4 below; this is *not* the same construction used for the oscillation-informed reference
distributions in the ordering diagnostics (§2).

Columns: `Σmν` = 95% CL upper bound (or fixed value); `ln B_NO/IO` = Eq. 5.3 log Bayes factor;
`P(Σmν<osc)` = Eq. 5.7 tail probability below the oscillation floor; `OVL` = Eq. 5.8 overlap
coefficient. NO floor ≈ 0.059 eV, IO floor ≈ 0.10 eV (Eq. 5.1).

| Model | Dataset | Neff | ΔNeff/σ | Σmν [eV] | ln B_NO/IO | P(Σmν<osc) NO | P(Σmν<osc) IO | OVL NO | OVL IO |
|---|---|---|---|---|---|---|---|---|---|
| ΛCDM + Σmν | CMB+DESI+PP | 3.04 (fixed) | – | < 0.0608 | 3.329 ± 0.154 | 0.942 | 0.999 | 0.013 | 0.001 |
| ΛCDM + Σmν | CMB+DESI+DD | 3.04 (fixed) | – | < 0.0613 | 3.169 ± 0.131 | 0.941 | 0.998 | 0.014 | 0.002 |
| ΛCDM + Neff | CMB+DESI+PP | 3.00 ± 0.10 | −0.393 | 0.06 (fixed) | – | – | – | – | – |
| ΛCDM + Neff | CMB+DESI+DD | 3.00 ± 0.11 | −0.456 | 0.06 (fixed) | – | – | – | – | – |
| ΛCDM + Σmν + Neff | CMB+DESI+PP | 2.97 ± 0.10 | −0.734 | < 0.0587 | 3.068 ± 0.127 | 0.950 | 0.998 | 0.013 | 0.001 |
| ΛCDM + Σmν + Neff | CMB+DESI+DD | 2.97 ± 0.11 | −0.728 | < 0.0592 | 3.139 ± 0.135 | 0.949 | 0.998 | 0.014 | 0.002 |
| ΛCDM + r+αs+βs+Ωk+Σmν+Neff | CMB+DESI+PP | 2.93 ± 0.18 | −0.655 | < 0.122 | 0.936 ± 0.012 | 0.679 | 0.900 | 0.036 | 0.031 |
| ΛCDM + r+αs+βs+Ωk+Σmν+Neff | CMB+DESI+DD | 2.92 ± 0.18 | −0.725 | < 0.122 | 0.936 ± 0.013 | 0.676 | 0.900 | 0.035 | 0.030 |
| w0waCDM + Σmν | CMB+DESI+PP | 3.04 (fixed) | – | < 0.113 | 1.048 ± 0.012 | 0.709 | 0.921 | 0.031 | 0.026 |
| w0waCDM + Σmν | CMB+DESI+DD | 3.04 (fixed) | – | < 0.123 | 0.915 ± 0.009 | 0.644 | 0.889 | 0.035 | 0.029 |
| w0waCDM + Neff | CMB+DESI+PP | 2.90 ± 0.11 | −1.286 | 0.06 (fixed) | – | – | – | – | – |
| w0waCDM + Neff | CMB+DESI+DD | 2.88 ± 0.11 | −1.442 | 0.06 (fixed) | – | – | – | – | – |
| w0waCDM + Σmν + Neff | CMB+DESI+PP | 2.89 ± 0.11 | −1.337 | < 0.0992 | 1.360 ± 0.020 | 0.769 | 0.953 | 0.029 | 0.022 |
| w0waCDM + Σmν + Neff | CMB+DESI+DD | 2.88 ± 0.11 | −1.447 | < 0.113 | 1.054 ± 0.014 | 0.703 | 0.920 | 0.032 | 0.026 |
| w0waCDM + r+αs+βs+Ωk+Σmν+Neff | CMB+DESI+PP | 2.90 +0.18/−0.20 | −0.757 | < 0.171 | 0.477 ± 0.007 | 0.502 | 0.752 | 0.050 | 0.046 |
| w0waCDM + r+αs+βs+Ωk+Σmν+Neff | CMB+DESI+DD | 2.88 +0.18/−0.20 | −0.882 | < 0.184 | 0.387 ± 0.006 | 0.443 | 0.696 | 0.052 | 0.049 |

**There is no row with curvature opened alone (or curvature + Σmν alone).** Ωk only appears
bundled with r, αs, βs, and Neff in the "kitchen-sink" row. See §5 (proposal check) — this
matters directly for H2.

---

## 2. Ordering diagnostics — exact definitions (§5.2–5.3, pp.22–25)

### Log Bayes factor, Eq. 5.3–5.6 (§5.2, p.22–24)

```
ln B_NO/IO = ln Z_NO − ln Z_IO
Z_H = ∫₀^∞ dΣ  P_cosmo(Σ) · P_osc^H(Σ),   H ∈ {NO, IO},  Σ ≡ Σmν
```

`Z_H` is the "overlap evidence" between the cosmological posterior `P_cosmo(Σ)` (from the
actual MCMC chain, degenerate-hierarchy prior) and an oscillation-informed reference
distribution `P_osc^H(Σ)` built as follows (footnote 9, p.24):

- Prior: **uniform in the *lightest* neutrino mass**, `m_lightest ∈ [0, 2] eV` — not a
  degenerate/equal-mass assumption. This is a genuine 3-eigenstate spectrum.
- For NO (lightest = m1): `Σmν = m1 + √(m1²+Δm²21) + √(m1²+Δm²31)`, with
  `Δm²21 = 7.42×10⁻⁵ eV²`, `Δm²31 = 2.517×10⁻³ eV²` (Eq. 5.4).
- For IO (lightest = m3): `Σmν = m3 + √(m3²+|Δm²31|) + √(m3²+|Δm²31|+Δm²21)`, with
  `|Δm²31| = 2.498×10⁻³ eV²` (Eq. 5.5).
- Both resulting `P_osc^NO`, `P_osc^IO` are normalized over `Σ ∈ [0, 5] eV`, matching the
  cosmological prior range (Table 1).
- `ln B_NO/IO` uncertainty is from bootstrap resampling of the weighted cosmological chain.

Interpretation scale (Jeffreys/Trotta): `|ln B| < 1` inconclusive, `1–2.5` weak,
`2.5–5` moderate-to-strong, `≥5` strong/very strong (p.24).

### Tail probability, Eq. 5.7 (§5.3, p.25)

```
P^H(Σmν<osc) = ∫₀^(Σmν)_min^H dΣ  P_cosmo(Σ)
```

Fraction of the cosmological posterior below the oscillation floor for ordering H
(`(Σmν)_min^NO ≈ 0.059 eV`, `(Σmν)_min^IO ≈ 0.10 eV`). Computed directly from the MCMC chain
as the posterior weight of samples below the floor (footnote 11).

### Overlap coefficient, Eq. 5.8 (§5.3, p.25)

```
OVL^H = ∫₀^∞ dΣ  min[ P_cosmo(Σ), P_osc^H(Σ) ]
```

Both distributions normalized over the same range, so `0 ≤ OVL^H ≤ 1` by construction.
Computed numerically as a weighted-histogram common area (footnote 12).

**Key nuance for our own configs:** the diagnostics compare a *degenerate-hierarchy*
cosmological posterior against a *non-degenerate, splitting-informed* oscillation prior.
These are deliberately different neutrino treatments serving different purposes — the
cosmological fit only needs the total mass, so degenerate is a standard simplification
(footnote 3, p.6), while the oscillation reference distributions need the real splittings to
mean anything as "what oscillation data imply."

---

## 3. Ωm and σ8 — Table 7 (p.37) and Figure 12 (p.40)

68% CL, CAMB+Cobaya. Full table (H0, Ωm, σ8, S8, rd all transcribed since they're one table):

| Model | Dataset | H0 | Ωm | σ8 | S8 | rd [Mpc] |
|---|---|---|---|---|---|---|
| ΛCDM | CMB+DESI+PP | 68.17 ± 0.25 | 0.3040 ± 0.0033 | 0.8182 ± 0.0042 | 0.8236 ± 0.0062 | 147.48 ± 0.17 |
| ΛCDM | CMB+DESI+DD | 68.13 ± 0.25 | 0.3044 ± 0.0034 | 0.8183 ± 0.0042 | 0.8243 ± 0.0061 | 147.47 ± 0.17 |
| ΛCDM + Σmν | CMB+DESI+PP | 68.35 +0.28/−0.24 | 0.3021 +0.0032/−0.0036 | 0.8241 +0.0053/−0.0046 | 0.8270 ± 0.0064 | 147.41 ± 0.18 |
| ΛCDM + Σmν | CMB+DESI+DD | 68.33 ± 0.26 | 0.3024 ± 0.0034 | 0.8241 +0.0053/−0.0046 | 0.8274 ± 0.0064 | 147.41 ± 0.18 |
| ΛCDM + Neff | CMB+DESI+PP | 67.91 ± 0.68 | 0.3048 ± 0.0037 | 0.8165 ± 0.0063 | 0.8230 ± 0.0066 | 147.9 ± 1.1 |
| ΛCDM + Neff | CMB+DESI+DD | 67.85 ± 0.68 | 0.3052 ± 0.0038 | 0.8161 ± 0.0063 | 0.8231 ± 0.0067 | 148.0 ± 1.1 |
| ΛCDM + Σmν + Neff | CMB+DESI+PP | 67.89 ± 0.68 | 0.3033 ± 0.0038 | 0.8209 ± 0.0065 | 0.8254 ± 0.0067 | 148.2 ± 1.1 |
| ΛCDM + Σmν + Neff | CMB+DESI+DD | 67.86 ± 0.68 | 0.3037 ± 0.0038 | 0.8207 ± 0.0067 | 0.8257 ± 0.0067 | 148.2 ± 1.1 |
| ΛCDM + Ωk | CMB+DESI+PP | 68.57 ± 0.32 | 0.3043 ± 0.0033 | 0.8230 ± 0.0047 | 0.8288 ± 0.0066 | 147.05 ± 0.26 |
| ΛCDM + Ωk | CMB+DESI+DD | 68.56 ± 0.31 | 0.3045 ± 0.0033 | 0.8231 ± 0.0047 | 0.8293 ± 0.0066 | 147.04 ± 0.26 |
| ΛCDM + r | CMB+DESI+PP | 68.17 ± 0.25 | 0.3039 ± 0.0033 | 0.8179 ± 0.0041 | 0.8231 ± 0.0061 | 147.49 ± 0.17 |
| ΛCDM + r | CMB+DESI+DD | 68.14 ± 0.24 | 0.3043 ± 0.0033 | 0.8178 ± 0.0042 | 0.8237 ± 0.0063 | 147.47 ± 0.18 |
| ΛCDM + r+αs | CMB+DESI+PP | 68.14 ± 0.25 | 0.3041 ± 0.0033 | 0.8176 ± 0.0041 | 0.8232 ± 0.0062 | 147.52 ± 0.18 |
| ΛCDM + r+αs | CMB+DESI+DD | 68.12 ± 0.25 | 0.3044 ± 0.0034 | 0.8178 ± 0.0040 | 0.8237 ± 0.0062 | 147.51 ± 0.18 |
| ΛCDM + r+αs+βs | CMB+DESI+PP | 68.10 ± 0.25 | 0.3046 ± 0.0034 | 0.8201 +0.0043/−0.0048 | 0.8264 ± 0.0066 | 147.53 ± 0.18 |
| ΛCDM + r+αs+βs | CMB+DESI+DD | 68.09 ± 0.25 | 0.3047 ± 0.0033 | 0.8201 ± 0.0046 | 0.8265 ± 0.0067 | 147.52 ± 0.18 |
| ΛCDM + r+αs+βs+Ωk+Σmν+Neff | CMB+DESI+PP | 68.01 ± 0.96 | 0.3059 ± 0.0047 | 0.8238 ± 0.0093 | 0.8318 ± 0.0080 | 148.2 ± 1.7 |
| ΛCDM + r+αs+βs+Ωk+Σmν+Neff | CMB+DESI+DD | 67.94 ± 0.94 | 0.3065 +0.0043/−0.0048 | 0.8236 ± 0.0092 | 0.8324 ± 0.0080 | 148.3 ± 1.7 |
| w0waCDM | CMB+DESI+PP | 67.63 ± 0.61 | 0.3109 ± 0.0058 | 0.8187 ± 0.0073 | 0.8334 ± 0.0069 | 147.26 ± 0.20 |
| w0waCDM | CMB+DESI+DD | 67.44 ± 0.54 | 0.3129 ± 0.0053 | 0.8176 ± 0.0070 | 0.8350 ± 0.0069 | 147.23 ± 0.20 |
| w0waCDM + Σmν | CMB+DESI+PP | 67.65 ± 0.60 | 0.3104 ± 0.0057 | 0.8200 ± 0.0081 | 0.8340 ± 0.0071 | 147.26 ± 0.20 |
| w0waCDM + Σmν | CMB+DESI+DD | 67.44 ± 0.54 | 0.3127 ± 0.0054 | 0.8183 ± 0.0078 | 0.8353 ± 0.0072 | 147.24 ± 0.20 |
| w0waCDM + Neff | CMB+DESI+PP | 66.94 ± 0.81 | 0.3122 ± 0.0058 | 0.8150 ± 0.0078 | 0.8313 ± 0.0070 | 148.8 ± 1.2 |
| w0waCDM + Neff | CMB+DESI+DD | 66.65 ± 0.78 | 0.3146 ± 0.0054 | 0.8133 ± 0.0075 | 0.8327 ± 0.0070 | 148.9 ± 1.2 |
| w0waCDM + Σmν + Neff | CMB+DESI+PP | 66.92 ± 0.81 | 0.3117 ± 0.0059 | 0.8166 ± 0.0083 | 0.8323 ± 0.0071 | 148.8 ± 1.2 |
| w0waCDM + Σmν + Neff | CMB+DESI+DD | 66.66 ± 0.77 | 0.3142 ± 0.0056 | 0.8147 ± 0.0083 | 0.8337 ± 0.0072 | 148.9 ± 1.2 |
| w0waCDM + Ωk | CMB+DESI+PP | 67.76 ± 0.61 | 0.3111 ± 0.0056 | 0.8194 ± 0.0074 | 0.8343 ± 0.0070 | 147.11 ± 0.26 |
| w0waCDM + Ωk | CMB+DESI+DD | 67.55 ± 0.58 | 0.3129 ± 0.0054 | 0.8180 ± 0.0071 | 0.8353 ± 0.0070 | 147.12 ± 0.26 |
| w0waCDM + r | CMB+DESI+PP | 67.65 ± 0.61 | 0.3107 ± 0.0057 | 0.8184 ± 0.0074 | 0.8328 ± 0.0068 | 147.27 ± 0.20 |
| w0waCDM + r | CMB+DESI+DD | 67.45 ± 0.54 | 0.3127 ± 0.0052 | 0.8172 ± 0.0070 | 0.8342 ± 0.0070 | 147.25 ± 0.21 |
| w0waCDM + r+αs | CMB+DESI+PP | 67.63 ± 0.59 | 0.3108 ± 0.0056 | 0.8183 ± 0.0073 | 0.8328 ± 0.0068 | 147.31 ± 0.20 |
| w0waCDM + r+αs | CMB+DESI+DD | 67.44 ± 0.55 | 0.3127 ± 0.0053 | 0.8171 ± 0.0070 | 0.8342 ± 0.0068 | 147.29 ± 0.20 |
| w0waCDM + r+αs+βs | CMB+DESI+PP | 67.59 ± 0.60 | 0.3113 ± 0.0057 | 0.8208 ± 0.0075 | 0.8360 ± 0.0072 | 147.32 ± 0.21 |
| w0waCDM + r+αs+βs | CMB+DESI+DD | 67.42 ± 0.54 | 0.3131 ± 0.0053 | 0.8201 ± 0.0070 | 0.8377 ± 0.0072 | 147.29 ± 0.21 |
| w0waCDM + r+αs+βs+Ωk+Σmν+Neff | CMB+DESI+PP | 67.1 ± 1.0 | 0.3135 ± 0.0061 | 0.817 ± 0.010 | 0.8351 ± 0.0082 | 148.5 ± 1.9 |
| w0waCDM + r+αs+βs+Ωk+Σmν+Neff | CMB+DESI+DD | 66.8 ± 1.0 | 0.3155 ± 0.0059 | 0.8148 ± 0.0099 | 0.8355 ± 0.0084 | 148.7 ± 1.9 |

**Text summary (p.35–36):** ΛCDM-branch models cluster at Ωm ≃ 0.303–0.306; w0waCDM-branch
models sit systematically higher, Ωm ≃ 0.311–0.316. σ8 stays in ≃0.816–0.824 (ΛCDM branch)
and ≃0.813–0.821 (w0waCDM branch) — i.e. the two branches barely differ in σ8 even though they
separate clearly in Ωm. Fig. 12 shows exactly this: in the (σ8, Ωm) plane the two branches
separate along Ωm with σ8 "broadly similar"; in the (S8, Ωm) plane the w0waCDM branch shifts to
higher S8 "because of the higher matter density" (p.39–40). Fig. 12 itself only plots baseline
ΛCDM vs. w0waCDM for CMB+DESI+DD (not the full model list).

---

## 4. Datasets, likelihoods, and priors (§2, pp.5–8)

### Priors — Table 1 (p.7), all flat

| Parameter | Prior | Description |
|---|---|---|
| Ωbh² ≡ ωb | [0.005, 0.1] | Physical baryon density |
| Ωch² ≡ ωc | [0.001, 0.99] | Physical cold-DM density |
| 100θs | [0.5, 10] | Angular size of sound horizon at recombination |
| τ | [0.01, 0.8] | Optical depth to reionization |
| ln(10¹⁰As) | [1.61, 3.91] | Amplitude of primordial scalar power spectrum |
| ns | [0.8, 1.2] | Scalar spectral index |
| Σmν [eV] | [0, 5] | Total neutrino mass |
| Neff | [0.05, 10] | Effective number of relativistic species |
| Ωk | [−0.3, 0.3] | Spatial curvature |
| r | [0, 0.5] | Tensor-to-scalar ratio |
| αs | [−1, 1] | Running of ns |
| βs | [−1, 1] | Running of the running |
| w0 | [−3, 1] | Present-day DE equation of state |
| wa | [−3, 2] | First-order DE time variation |

Neutrino treatment (footnote 3, p.6): **`num_massive_neutrinos: 3`,
`neutrino_hierarchy: degenerate`** for every free-Σmν run.

### Datasets/likelihoods (§2.2, pp.6–8)

- **Planck low-T**: Planck-2018 low-ℓ TT, Commander likelihood, ℓ<30.
- **Planck low-E**: Planck-2018 low-ℓ EE, SRoll2 likelihood, ℓ<30.
- **Planck (TT-TE-EE)**: Planck-2018 mid-ℓ TT/TE/EE, ℓmax = (1000, 600, 600).
- **ACT (TT-TE-EE)**: ACT DR6 high-ℓ TT/TE/EE, ℓ≥600.
- **SPT (TT-TE-EE)**: SPT-3G D1 high-ℓ TT/TE/EE, ℓ≥400, TT cut ℓmax=3000, TE/EE cut ℓmax=4000.
- **Planck-ACT-SPT Lensing**: ACT DR6 lensing + Planck lensing + SPT-3G MUSE lensing, combined.
- **BICEP-Keck**: BK18 BB spectrum, 20 ≲ ℓ ≲ 330.
- **BAO**: DESI DR2 BAO (galaxies, quasars, Lyα), 14 redshift bins, full covariance.
- **Pantheon-Plus (PP)**: 1701 light curves, 1550 objects, z up to ≃2.26.
- **DES-Dovekie (DD)**: recalibrated DES-5yr SNe (Dovekie cross-calibration), ~1600 DES SNe +
  ~200 low-z external SNe.

Two baseline combinations used throughout:
- **CMB+DESI+DD** = Planck+ACT+SPT+lensing+BK18 + DESI DR2 BAO + DES-Dovekie.
- **CMB+DESI+PP** = identical, with Pantheon-Plus in place of DES-Dovekie.

---

## 5. Comparison against the proposal

**Overall the proposal is accurate.** Every specific number and claim checked below traces
cleanly to the paper, with one significant gap and one presentation nit.

### ✅ Confirmed correct
- Oscillation floors 0.059 eV (NO) / 0.10 eV (IO) — matches Eq. 5.1 exactly.
- Table 4 headline bounds: ΛCDM+Σmν ≃0.061 eV, w0waCDM+Σmν → 0.113–0.123 eV, fully-open →
  0.171–0.184 eV — all match.
- Weights: ~94% (ΛCDM, NO floor), ~68% (fully-open ΛCDM-branch, PP), ~64% (w0waCDM, DD), ~44%
  (fully-open w0waCDM, DD) — all match Table 4's `P(Σmν<osc)` column exactly (see note below on
  which dataset each proposal row is drawn from).
- Loverde & Weiner is indeed cited as **[222]** in Intertwined, exactly where the proposal says
  (confirmed in Intertwined's bibliography, p.92, and in the Ωm-mismatch citation list on p.20).
- "Matter density shifts upward in evolving-DE models while σ8 barely changes": confirmed by
  Table 7 (Ωm ≃0.304 → ≃0.311–0.316 across branches; σ8 stays ≃0.816–0.824 vs. ≃0.813–0.821).
- Loverde & Weiner's claims as summarized in the proposal (CMB forces Ωm up steeply with
  neutrino mass; lensing observables depend on matter *fraction* not density; DESI's bound is
  mostly geometric; SNe pull the opposite direction from BAO because they prefer larger Ωm) all
  match their abstract essentially verbatim.

### ⚠️ Gap that affects H2 directly
**There is no "ΛCDM + curvature + mass" (or w0waCDM + curvature + mass) row in Table 4.**
Intertwined only ever opens Ωk for the neutrino-mass results bundled together with r, αs, βs,
and Neff (the "kitchen-sink" row), never in isolation the way w0waCDM+Σmν isolates the DE
channel. (Table 7 *does* have a standalone ΛCDM+Ωk row, but that run holds Σmν fixed — it's not
a neutrino-mass constraint at all.)

This means the proposal's three-rung ladder — "LambdaCDM + mass, w0waCDM + mass, and LambdaCDM
+ curvature + mass" — has no published number for the third rung to validate against. Running
a clean ΛCDM+Ωk+Σmν model is still a perfectly good thing to do (it's a new result, not a
reproduction), but the week-4–6 "qualitative agreement with Table 4" gate as written can only
apply to rungs 1 and 2. **H2 will need its own validation logic, or be reframed as "does the
isolated curvature channel behave like the kitchen-sink row's curvature+inflation+Neff channel"
rather than expecting a directly comparable published bound.** Worth raising this explicitly
before committing time to H2 over H3.

### 📝 Minor presentation nit
The proposal's summary table blends CMB+DESI+PP and CMB+DESI+DD numbers row-to-row without
labeling which is which (row 1 is close either way; row 2's 0.122 eV/68% is the **PP** column;
row 3's 0.123 eV/64% is the **DD** column; row 4's 0.184 eV/44% is the **DD** column). Not
wrong — every number is real — just worth knowing so you don't accidentally compare across rows
as if they were the same dataset combination.

### ℹ️ Useful nuance not in the proposal
The proposal doesn't say what neutrino-mass prior/hierarchy assumption to use. Intertwined's
own cosmological fits use **degenerate hierarchy** (footnote 3, p.6) — that's the standard
choice and what our own `configs/` should default to for the main chains. The ordering
diagnostics (Bayes factor, tail probability, overlap) additionally build *separate*
oscillation-informed reference distributions using the real mass splittings and a uniform prior
on the lightest mass in [0,2] eV (footnote 9, p.24) — but that's only for those three specific
diagnostics, not the main P mν posterior. Good to know precisely because the proposal's Step 1
("rebuild the ordering diagnostics ... on a public chain") needs exactly this distinction to get
right.

### Version check (v1 → v2, both papers)
Diffed v1 against v2 for both papers directly. **No numerical or definitional changes** in
either paper's core content:
- Intertwined: Table 4, Table 7, and the Section 5.2/5.3 diagnostic definitions are
  byte-identical text between v1 and v2. The only differences found are bibliography citation
  renumbering (e.g. Jeffreys/Trotta scale ref [327,328]→[338,339]) from new references added
  elsewhere, accounting for the +1 page in v2.
- Loverde & Weiner: near-every page has small wording/copyedit changes (consistent with
  JCAP referee revisions) and citation renumbering, but the abstract, core equations (e.g. Eq.
  A1, `100θ_CMB,⊥(ad) ~ N(1.06061, 4.2×10⁻⁴)`), and the appendix robustness table are unchanged
  in substance — the appendix table just moved pages because an Acknowledgments section was
  added before it.

Nothing here should change how the proposal cites either paper, other than pinning both to v2
as the current version.
