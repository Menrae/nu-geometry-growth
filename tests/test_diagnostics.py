"""Tests for nuproj.diagnostics against synthetic distributions with
analytically known answers, per docs/notes/intertwined_summary.md Sec. 2.
"""

import numpy as np
import pytest
from scipy import stats

from nuproj.diagnostics import (
    DM21_SQ_DEFAULT,
    DM31_SQ_IO_DEFAULT,
    DM31_SQ_NO_DEFAULT,
    log_bayes_factor,
    oscillation_floor,
    oscillation_reference_samples,
    overlap_coefficient,
    overlap_of_samples,
    tail_probability,
    upper_limit,
)

RNG_SEED = 12345


def _truncnorm(loc, scale, lower=0.0, upper=np.inf):
    a = (lower - loc) / scale
    b = (upper - loc) / scale
    return stats.truncnorm(a, b, loc=loc, scale=scale)


# ---------------------------------------------------------------------------
# oscillation_floor: pure closed-form, exact-testable
# ---------------------------------------------------------------------------


def test_oscillation_floor_no_matches_hand_calculation():
    expected = np.sqrt(DM21_SQ_DEFAULT) + np.sqrt(DM31_SQ_NO_DEFAULT)
    assert oscillation_floor("NO") == pytest.approx(expected)
    # matches the paper's quoted ~0.059 eV (Eq. 5.1) to 3% (rounding in the paper)
    assert oscillation_floor("NO") == pytest.approx(0.059, rel=0.03)


def test_oscillation_floor_io_matches_hand_calculation():
    expected = np.sqrt(DM31_SQ_IO_DEFAULT) + np.sqrt(DM31_SQ_IO_DEFAULT + DM21_SQ_DEFAULT)
    assert oscillation_floor("IO") == pytest.approx(expected)
    assert oscillation_floor("IO") == pytest.approx(0.10, rel=0.03)


def test_oscillation_floor_io_above_no():
    assert oscillation_floor("IO") > oscillation_floor("NO")


def test_oscillation_floor_rejects_bad_hierarchy():
    with pytest.raises(ValueError):
        oscillation_floor("XX")


def test_oscillation_floor_lowercase_accepted():
    assert oscillation_floor("no") == oscillation_floor("NO")


# ---------------------------------------------------------------------------
# tail_probability: empirical weighted CDF vs. analytic truncated-Gaussian CDF
# ---------------------------------------------------------------------------


def test_tail_probability_matches_truncnorm_cdf():
    dist = _truncnorm(loc=0.08, scale=0.04)
    rng = np.random.default_rng(RNG_SEED)
    values = dist.rvs(size=200_000, random_state=rng)
    floor = 0.0588
    expected = dist.cdf(floor)
    got = tail_probability(values, floor)
    assert got == pytest.approx(expected, abs=0.005)


def test_tail_probability_zero_below_support():
    dist = _truncnorm(loc=0.2, scale=0.02)
    values = dist.rvs(size=50_000, random_state=RNG_SEED)
    assert tail_probability(values, 0.0) == pytest.approx(0.0, abs=1e-9)


def test_tail_probability_one_above_support():
    dist = _truncnorm(loc=0.02, scale=0.01)
    values = dist.rvs(size=50_000, random_state=RNG_SEED)
    assert tail_probability(values, 1.0) == pytest.approx(1.0, abs=1e-9)


def test_tail_probability_respects_weights():
    # A population that's 90% concentrated below the floor and 10% above,
    # built explicitly with weights so the *unweighted* proportions would
    # give the wrong answer if weights were ignored.
    below = np.full(10, 0.01)
    above = np.full(10, 0.20)
    values = np.concatenate([below, above])
    weights = np.concatenate([np.full(10, 9.0), np.full(10, 1.0)])
    got = tail_probability(values, 0.1, weights=weights)
    assert got == pytest.approx(0.9)


def test_tail_probability_accepts_value_weight_tuple():
    values = np.array([0.01, 0.02, 0.2])
    weights = np.array([1.0, 1.0, 1.0])
    got = tail_probability((values, weights), 0.1)
    assert got == pytest.approx(2.0 / 3.0)


# ---------------------------------------------------------------------------
# upper_limit: weighted quantile vs. analytic truncated-Gaussian ppf
# ---------------------------------------------------------------------------


def test_upper_limit_matches_truncnorm_ppf():
    dist = _truncnorm(loc=0.06, scale=0.03)
    rng = np.random.default_rng(RNG_SEED)
    values = dist.rvs(size=200_000, random_state=rng)
    expected = dist.ppf(0.95)
    got = upper_limit(values, cl=0.95)
    assert got == pytest.approx(expected, abs=0.005)


def test_upper_limit_weighted_matches_duplication():
    rng = np.random.default_rng(RNG_SEED)
    dist = _truncnorm(loc=0.05, scale=0.02)
    base = dist.rvs(size=5_000, random_state=rng)
    # weight=3 for every sample should equal literally tripling the sample
    weighted_result = upper_limit(base, cl=0.95, weights=np.full(base.shape, 3.0))
    duplicated_result = upper_limit(np.tile(base, 3), cl=0.95)
    assert weighted_result == pytest.approx(duplicated_result)


# ---------------------------------------------------------------------------
# overlap_of_samples: boundary-corrected KDE overlap vs. quad on known pdfs
# ---------------------------------------------------------------------------


def test_overlap_of_samples_identical_distributions_near_one():
    dist = _truncnorm(loc=0.08, scale=0.03)
    rng = np.random.default_rng(RNG_SEED)
    a = dist.rvs(size=100_000, random_state=rng)
    b = dist.rvs(size=100_000, random_state=np.random.default_rng(RNG_SEED + 1))
    ovl = overlap_of_samples(
        a, np.ones_like(a), b, np.ones_like(b),
        bounds_a=(0.0, None), bounds_b=(0.0, None), grid_range=(0.0, 0.5),
    )
    assert ovl > 0.9


def test_overlap_of_samples_disjoint_near_zero():
    dist_a = _truncnorm(loc=0.02, scale=0.005)
    dist_b = _truncnorm(loc=0.3, scale=0.01, lower=0.2)
    rng = np.random.default_rng(RNG_SEED)
    a = dist_a.rvs(size=100_000, random_state=rng)
    b = dist_b.rvs(size=100_000, random_state=np.random.default_rng(RNG_SEED + 1))
    ovl = overlap_of_samples(
        a, np.ones_like(a), b, np.ones_like(b),
        bounds_a=(0.0, None), bounds_b=(0.2, None), grid_range=(0.0, 0.5),
    )
    assert ovl < 0.02


def test_overlap_of_samples_matches_analytic_integral():
    # Ground truth via scipy.integrate.quad on the exact truncated-normal pdfs,
    # independent of this module's KDE machinery.
    dist_a = _truncnorm(loc=0.08, scale=0.03)
    dist_b = _truncnorm(loc=0.12, scale=0.04)
    from scipy.integrate import quad

    exact, _ = quad(lambda x: min(dist_a.pdf(x), dist_b.pdf(x)), 0.0, 0.5, limit=200)

    rng = np.random.default_rng(RNG_SEED)
    a = dist_a.rvs(size=200_000, random_state=rng)
    b = dist_b.rvs(size=200_000, random_state=np.random.default_rng(RNG_SEED + 1))
    got = overlap_of_samples(
        a, np.ones_like(a), b, np.ones_like(b),
        bounds_a=(0.0, None), bounds_b=(0.0, None), grid_range=(0.0, 0.5),
    )
    assert got == pytest.approx(exact, abs=0.03)


# ---------------------------------------------------------------------------
# oscillation_reference_samples: sanity checks on the physical construction
# ---------------------------------------------------------------------------


def test_oscillation_reference_samples_respect_floor():
    for hierarchy in ("NO", "IO"):
        samples = oscillation_reference_samples(hierarchy, n_samples=50_000, seed=RNG_SEED)
        assert samples.min() >= oscillation_floor(hierarchy) - 1e-12


def test_oscillation_reference_samples_deterministic_with_seed():
    a = oscillation_reference_samples("NO", n_samples=1000, seed=42)
    b = oscillation_reference_samples("NO", n_samples=1000, seed=42)
    np.testing.assert_array_equal(a, b)


# ---------------------------------------------------------------------------
# overlap_coefficient / log_bayes_factor: physics-facing wrappers
# ---------------------------------------------------------------------------


def test_overlap_coefficient_bounded_between_zero_and_one():
    dist = _truncnorm(loc=0.07, scale=0.03)
    values = dist.rvs(size=100_000, random_state=RNG_SEED)
    for hierarchy in ("NO", "IO"):
        ovl = overlap_coefficient(values, hierarchy, n_osc_samples=100_000, seed=0)
        assert 0.0 <= ovl <= 1.0


def test_overlap_coefficient_higher_for_no_when_posterior_is_tight_and_low():
    # A posterior compressed near the NO floor (like Intertwined's minimal
    # LambdaCDM+mnu row) should overlap the NO reference distribution more
    # than the IO one, which requires masses beyond 0.10 eV.
    dist = _truncnorm(loc=0.065, scale=0.01)
    values = dist.rvs(size=200_000, random_state=RNG_SEED)
    ovl_no = overlap_coefficient(values, "NO", n_osc_samples=200_000, seed=0)
    ovl_io = overlap_coefficient(values, "IO", n_osc_samples=200_000, seed=0)
    assert ovl_no > ovl_io


def test_log_bayes_factor_prefers_no_when_posterior_between_floors():
    # Posterior concentrated strictly between the NO and IO floors: P_osc^NO
    # has support there but P_osc^IO does not, so NO should be strongly
    # preferred (positive ln B), matching Table 4's minimal LambdaCDM+mnu row.
    no_floor = oscillation_floor("NO")
    io_floor = oscillation_floor("IO")
    mid = 0.5 * (no_floor + io_floor)
    dist = _truncnorm(loc=mid, scale=0.005, lower=no_floor)
    values = dist.rvs(size=200_000, random_state=RNG_SEED)
    ln_b, sigma = log_bayes_factor(values, n_osc_samples=200_000, seed=0, n_bootstrap=50)
    assert ln_b > 1.0
    assert sigma >= 0.0


def test_log_bayes_factor_matches_analytic_integral():
    # Independent ground truth: invert the exact NO/IO mass relation on a
    # dense grid (not Monte Carlo) to get an analytic oscillation density,
    # then integrate against a known cosmological pdf via quad.
    from scipy.integrate import quad

    def analytic_osc_density(hierarchy, upper=5.0, m1_max=2.0, n_grid=20_000):
        m1_grid = np.linspace(0.0, m1_max, n_grid)
        if hierarchy == "NO":
            d21, d31 = DM21_SQ_DEFAULT, DM31_SQ_NO_DEFAULT
            sigma_grid = m1_grid + np.sqrt(m1_grid**2 + d21) + np.sqrt(m1_grid**2 + d31)
        else:
            d21, d31 = DM21_SQ_DEFAULT, DM31_SQ_IO_DEFAULT
            sigma_grid = m1_grid + np.sqrt(m1_grid**2 + d31) + np.sqrt(m1_grid**2 + d31 + d21)
        floor = sigma_grid[0]

        def raw_density(sigma):
            sigma = np.atleast_1d(sigma).astype(float)
            m1 = np.interp(sigma, sigma_grid, m1_grid, left=np.nan, right=np.nan)
            if hierarchy == "NO":
                dsigma_dm1 = 1 + m1 / np.sqrt(m1**2 + d21) + m1 / np.sqrt(m1**2 + d31)
            else:
                dsigma_dm1 = 1 + m1 / np.sqrt(m1**2 + d31) + m1 / np.sqrt(m1**2 + d31 + d21)
            dens = (1.0 / m1_max) / dsigma_dm1
            dens = np.where((sigma >= floor) & (sigma <= sigma_grid[-1]), dens, 0.0)
            return np.nan_to_num(dens)

        norm, _ = quad(lambda s: raw_density(s)[0], floor, upper, limit=200)
        return floor, lambda sigma: raw_density(sigma) / norm

    dist = _truncnorm(loc=0.08, scale=0.02)

    z = {}
    for hierarchy in ("NO", "IO"):
        floor, p_osc = analytic_osc_density(hierarchy)
        integrand = lambda s, p_osc=p_osc: dist.pdf(s) * p_osc(s)[0]
        z[hierarchy], _ = quad(integrand, floor, 5.0, limit=200)
    ln_b_exact = np.log(z["NO"]) - np.log(z["IO"])

    values = dist.rvs(size=300_000, random_state=RNG_SEED)
    ln_b_got, _ = log_bayes_factor(values, n_osc_samples=300_000, seed=0, n_bootstrap=20)

    assert ln_b_got == pytest.approx(ln_b_exact, abs=0.15)
