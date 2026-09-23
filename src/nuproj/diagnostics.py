"""Neutrino-mass ordering diagnostics.

Implements the three ordering diagnostics and the 95% upper bound defined in
``docs/notes/intertwined_summary.md`` (Sec. 2), which transcribes Sec. 5.2-5.3
and footnotes 9-12 of Giare, Lee & Di Valentino (2026), arXiv:2607.01226v2
(``docs/papers/2607.01226v2_intertwined.pdf``):

- :func:`tail_probability` -- Eq. 5.7, the posterior weight below an
  oscillation-implied mass floor.
- :func:`upper_limit` -- the 95% one-sided upper credible bound on Sigma m_nu.
- :func:`overlap_coefficient` -- Eq. 5.8, the overlap between the cosmological
  posterior and an oscillation-informed reference distribution.
- :func:`log_bayes_factor` -- Eq. 5.3-5.6, the log Bayes factor between the
  normal and inverted mass orderings.

All functions accept weighted samples as a ``getdist.MCSamples`` instance, a
``(values, weights)`` tuple, or a bare 1D array (uniform weights).

See ``docs/notes/decisions.md`` (2026-09-23 entries) for why density
estimation here uses GetDist's boundary-corrected KDE rather than a
from-scratch implementation, and why the oscillation-informed reference
distribution is built by Monte Carlo rather than inverting the hierarchy
relation analytically.
"""

from __future__ import annotations

from typing import Literal, Optional, Tuple, Union

import numpy as np

from getdist import MCSamples

HierarchyName = Literal["NO", "IO"]
SampleInput = Union[MCSamples, Tuple[np.ndarray, np.ndarray], np.ndarray]

# Neutrino mass-squared splittings, footnote 9 of docs/papers/2607.01226v2 (p.24):
# uniform prior on the lightest mass, splittings fixed to these central values.
DM21_SQ_DEFAULT = 7.42e-5  # eV^2, Delta m^2_21
DM31_SQ_NO_DEFAULT = 2.517e-3  # eV^2, Delta m^2_31 for normal ordering
DM31_SQ_IO_DEFAULT = 2.498e-3  # eV^2, |Delta m^2_31| for inverted ordering


def _validate_hierarchy(hierarchy: str) -> str:
    hierarchy = hierarchy.upper()
    if hierarchy not in ("NO", "IO"):
        raise ValueError(f"hierarchy must be 'NO' or 'IO', got {hierarchy!r}")
    return hierarchy


def _extract_samples(
    samples: SampleInput,
    weights: Optional[np.ndarray] = None,
    param: str = "mnu",
) -> Tuple[np.ndarray, np.ndarray]:
    """Normalize any accepted input into (values, weights) float arrays."""
    if isinstance(samples, MCSamples):
        values = np.asarray(samples.samples[:, samples.index[param]], dtype=float)
        w = np.asarray(samples.weights, dtype=float)
        return values, w
    if isinstance(samples, tuple):
        values, w = samples
        values = np.asarray(values, dtype=float)
        w = np.ones_like(values) if w is None else np.asarray(w, dtype=float)
        return values, w
    values = np.asarray(samples, dtype=float)
    w = np.ones_like(values) if weights is None else np.asarray(weights, dtype=float)
    return values, w


def oscillation_floor(
    hierarchy: HierarchyName,
    dm21_sq: float = DM21_SQ_DEFAULT,
    dm31_sq: Optional[float] = None,
) -> float:
    """Minimum total neutrino mass sum for ``hierarchy`` (lightest mass -> 0).

    NO: Sigma_min = sqrt(dm21_sq) + sqrt(dm31_sq)            [Eq. 5.4 at m1=0]
    IO: Sigma_min = sqrt(dm31_sq) + sqrt(dm31_sq + dm21_sq)  [Eq. 5.5 at m3=0]
    """
    hierarchy = _validate_hierarchy(hierarchy)
    if hierarchy == "NO":
        d31 = DM31_SQ_NO_DEFAULT if dm31_sq is None else dm31_sq
        return float(np.sqrt(dm21_sq) + np.sqrt(d31))
    d31 = DM31_SQ_IO_DEFAULT if dm31_sq is None else dm31_sq
    return float(np.sqrt(d31) + np.sqrt(d31 + dm21_sq))


def oscillation_reference_samples(
    hierarchy: HierarchyName,
    n_samples: int = 200_000,
    m_lightest_max: float = 2.0,
    dm21_sq: float = DM21_SQ_DEFAULT,
    dm31_sq: Optional[float] = None,
    seed: Optional[int] = None,
) -> np.ndarray:
    """Monte Carlo draws of Sigma m_nu implied by a uniform-in-lightest-mass
    prior, ``m_lightest ~ U(0, m_lightest_max)``, mapped through the exact
    NO/IO relation (Eqs. 5.4-5.5), per footnote 9.
    """
    hierarchy = _validate_hierarchy(hierarchy)
    rng = np.random.default_rng(seed)
    m1 = rng.uniform(0.0, m_lightest_max, size=n_samples)
    if hierarchy == "NO":
        d31 = DM31_SQ_NO_DEFAULT if dm31_sq is None else dm31_sq
        return m1 + np.sqrt(m1**2 + dm21_sq) + np.sqrt(m1**2 + d31)
    d31 = DM31_SQ_IO_DEFAULT if dm31_sq is None else dm31_sq
    return m1 + np.sqrt(m1**2 + d31) + np.sqrt(m1**2 + d31 + dm21_sq)


def tail_probability(
    samples: SampleInput,
    floor: float,
    *,
    weights: Optional[np.ndarray] = None,
    param: str = "mnu",
) -> float:
    """Eq. 5.7: fraction of the posterior with Sigma m_nu below ``floor``.

    Computed as the empirical weighted CDF directly from the samples (no
    density estimation), matching footnote 11's "computed directly from the
    MCMC chain as the total posterior weight of samples" prescription.
    """
    values, w = _extract_samples(samples, weights, param)
    mask = values < floor
    return float(w[mask].sum() / w.sum())


def upper_limit(
    samples: SampleInput,
    cl: float = 0.95,
    *,
    weights: Optional[np.ndarray] = None,
    param: str = "mnu",
) -> float:
    """Weighted one-sided upper credible bound (e.g. the 95% CL upper limit)."""
    values, w = _extract_samples(samples, weights, param)
    order = np.argsort(values)
    values_sorted = values[order]
    w_sorted = w[order]
    cum = np.cumsum(w_sorted) / w_sorted.sum()
    idx = min(int(np.searchsorted(cum, cl)), len(values_sorted) - 1)
    return float(values_sorted[idx])


def _boundary_kde(
    values: np.ndarray,
    weights: np.ndarray,
    lower: Optional[float],
    upper: Optional[float],
) -> Tuple[np.ndarray, np.ndarray]:
    """Boundary-corrected 1D KDE via GetDist (Jones 1993 linear boundary
    kernel correction, GetDist's default ``boundary_correction_order=1``,
    applied automatically whenever ``ranges`` gives a hard limit). See
    docs/notes/decisions.md for why this is used instead of a custom KDE.

    Returns (x, P) with P *not* normalized to integrate to 1 (GetDist
    normalizes its density to peak=1 for plotting) -- callers must
    renormalize on their own evaluation grid.
    """
    mc = MCSamples(
        samples=values,
        weights=weights,
        names=["x"],
        ranges={"x": (lower, upper)},
        settings={"smooth_scale_1D": -1.0},
    )
    density = mc.get1DDensity("x")
    return density.x, density.P


def overlap_of_samples(
    values_a: np.ndarray,
    weights_a: np.ndarray,
    values_b: np.ndarray,
    weights_b: np.ndarray,
    *,
    bounds_a: Tuple[Optional[float], Optional[float]] = (0.0, None),
    bounds_b: Tuple[Optional[float], Optional[float]] = (0.0, None),
    grid_range: Tuple[float, float] = (0.0, 5.0),
    num_points: int = 4096,
) -> float:
    """Overlap coefficient integral(min(P_a, P_b)) dx between two weighted
    1D sample sets, each boundary-corrected-KDE-estimated on its own hard
    range, evaluated on a shared grid over ``grid_range``. This is the
    general primitive behind :func:`overlap_coefficient` (Eq. 5.8).
    """
    x_a, p_a = _boundary_kde(values_a, weights_a, *bounds_a)
    x_b, p_b = _boundary_kde(values_b, weights_b, *bounds_b)
    grid = np.linspace(grid_range[0], grid_range[1], num_points)
    p_a_grid = np.interp(grid, x_a, p_a, left=0.0, right=0.0)
    p_b_grid = np.interp(grid, x_b, p_b, left=0.0, right=0.0)
    p_a_grid /= np.trapezoid(p_a_grid, grid)
    p_b_grid /= np.trapezoid(p_b_grid, grid)
    return float(np.trapezoid(np.minimum(p_a_grid, p_b_grid), grid))


def overlap_coefficient(
    samples: SampleInput,
    hierarchy: HierarchyName,
    *,
    weights: Optional[np.ndarray] = None,
    param: str = "mnu",
    prior_upper: float = 5.0,
    n_osc_samples: int = 200_000,
    seed: Optional[int] = 0,
) -> float:
    """Eq. 5.8: overlap between the cosmological posterior and the
    oscillation-informed reference distribution for ``hierarchy``.
    """
    hierarchy = _validate_hierarchy(hierarchy)
    values, w = _extract_samples(samples, weights, param)
    floor = oscillation_floor(hierarchy)
    osc_samples = oscillation_reference_samples(hierarchy, n_samples=n_osc_samples, seed=seed)
    osc_samples = osc_samples[osc_samples <= prior_upper]  # renormalize over [0, prior_upper], footnote 9
    return overlap_of_samples(
        values,
        w,
        osc_samples,
        np.ones_like(osc_samples),
        bounds_a=(0.0, None),
        bounds_b=(floor, None),
        grid_range=(0.0, prior_upper),
    )


def log_bayes_factor(
    samples: SampleInput,
    *,
    weights: Optional[np.ndarray] = None,
    param: str = "mnu",
    prior_upper: float = 5.0,
    n_osc_samples: int = 200_000,
    seed: Optional[int] = 0,
    n_bootstrap: int = 200,
) -> Tuple[float, float]:
    """Eq. 5.3-5.6: log Bayes factor ln(B_NO/IO) = ln Z_NO - ln Z_IO, with
    Z_H = <P_osc^H(Sigma)>_cosmo (footnote 10: evaluate the oscillation-
    informed weight at each cosmological sample and average over the
    posterior). Uncertainty from bootstrap resampling of the weighted
    cosmological samples (footnote 10).

    Returns (ln_B, uncertainty).
    """
    values, w = _extract_samples(samples, weights, param)
    grid = np.linspace(0.0, prior_upper, 4096)
    density_fn = {}
    for h in ("NO", "IO"):
        floor = oscillation_floor(h)
        osc_samples = oscillation_reference_samples(h, n_samples=n_osc_samples, seed=seed)
        osc_samples = osc_samples[osc_samples <= prior_upper]
        x_osc, p_osc = _boundary_kde(osc_samples, np.ones_like(osc_samples), floor, None)
        p_osc_grid = np.interp(grid, x_osc, p_osc, left=0.0, right=0.0)
        p_osc_grid /= np.trapezoid(p_osc_grid, grid)
        density_fn[h] = lambda sigma, g=grid, p=p_osc_grid: np.interp(sigma, g, p, left=0.0, right=0.0)

    z_no = np.average(density_fn["NO"](values), weights=w)
    z_io = np.average(density_fn["IO"](values), weights=w)
    ln_b = float(np.log(z_no) - np.log(z_io))

    rng = np.random.default_rng(seed)
    probs = w / w.sum()
    n = len(values)
    boot = np.empty(n_bootstrap)
    for i in range(n_bootstrap):
        idx = rng.choice(n, size=n, replace=True, p=probs)
        resampled = values[idx]
        zb_no = np.mean(density_fn["NO"](resampled))
        zb_io = np.mean(density_fn["IO"](resampled))
        boot[i] = np.log(zb_no) - np.log(zb_io)

    return ln_b, float(np.std(boot))
