"""Seeded statistical checks against independent small hafnian probabilities.

These are diagnostics, not a proof or a performance benchmark. The analytical
theorems do not depend on Monte Carlo agreement.
"""
from __future__ import annotations
import math
from collections import Counter
import numpy as np
from gaussian_loops import GaussianLoopSampler, exact_probabilities


def run(samples: int = 12000) -> dict:
    C = np.array([[0.10, 0.21, 0.07], [0.21, 0.04, 0.14], [0.07, 0.14, 0.09]])
    D = np.exp(1j * np.array([0.4, 1.2, -0.7]))
    cases = {
        "one_mode": np.array([[0.55]]),
        "two_mode_pair": np.array([[0, 0.5], [0.5, 0]]),
        "positive_dense": C,
        "gauge_phases": C * D[:, None] * D[None, :],
        "signed_real": np.array([[0.4, 0.2], [0.2, -0.4]]),
        "complex_phases": C * np.exp(1j * np.array(
            [[0.2, 0.5, 1.4], [0.5, -0.7, 0.9], [1.4, 0.9, 0.3]])),
    }
    results = {}
    for case_number, (name, B) in enumerate(cases.items()):
        sampler = GaussianLoopSampler(B, seed=1700 + case_number)
        exact = exact_probabilities(B, 14)
        weighted = Counter()
        total = squared = photons = 0.0
        for _ in range(samples):
            n, weight = sampler.draw_cosine()
            assert int(n.sum()) % 2 == 0
            weighted[tuple(map(int, n))] += weight
            total += weight
            squared += weight*weight
            photons += int(n.sum())
        probability_error = max(abs(weighted[n] / samples - p) for n, p in exact.items())
        second = math.exp(sampler.gamma)
        # Conservative fixed tolerance expressed in natural sampling scales.
        assert probability_error < 7 * math.sqrt(second / samples)
        assert abs(total / samples - 1) < 7 * math.sqrt(second / samples)
        assert abs(squared / samples - second) < 7 * math.exp(2*sampler.delta) / math.sqrt(samples)
        singular = np.linalg.svd(sampler.C, compute_uv=False)
        photon_variance = float(2*np.sum(singular**2/(1-singular**2)**2))
        assert abs(photons / samples - sampler.mean_C) < 7*math.sqrt(photon_variance/samples) + 1e-12
        record = {
            "samples": samples,
            "enumerated_target_mass_through_14": sum(exact.values()),
            "delta": sampler.delta,
            "gamma": sampler.gamma,
            "mean_weight": total / samples,
            "empirical_second_moment": squared / samples,
            "predicted_second_moment": second,
            "maximum_probability_absolute_error": probability_error,
            "empirical_reference_photons": photons / samples,
            "predicted_reference_photons": sampler.mean_C,
        }
        if name == "signed_real":
            # Exact r is available because all holonomies are signs.
            r = sampler.delta/2
            signed_sums = Counter()
            positive_sums = Counter()
            signs_total = 0.0
            for _ in range(samples):
                n, w = sampler.draw_signed(negative_mass=r)
                signed_sums[tuple(map(int, n))] += w
                signs_total += w
                positive_sums[tuple(map(int, sampler.sample_positive_part()))] += 1
            signed_error = max(abs(signed_sums[n]/samples-p) for n,p in exact.items())
            assert signed_error < 7 * math.exp(2*r) / math.sqrt(samples)
            assert abs(signs_total/samples-1) < 7 * math.exp(2*r) / math.sqrt(samples)
            # An observed finite-event discrepancy is bounded by global TV.
            bound = math.expm1(2*r)/2
            positive_event_error = max(abs(positive_sums[n]/samples-p) for n,p in exact.items())
            assert positive_event_error < bound + 7/math.sqrt(samples)
            record.update({"negative_mass": r,
                           "signed_estimator_maximum_probability_error": signed_error,
                           "positive_part_maximum_single_event_error": positive_event_error,
                           "positive_part_analytical_TV_bound": bound})
        results[name] = record
    return {"status": "passed", "cases": results}
