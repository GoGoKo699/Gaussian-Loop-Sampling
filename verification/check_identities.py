"""Deterministic independent checks of normalization, phases, and loop moments."""
from __future__ import annotations
import itertools
import math
import numpy as np
from gaussian_loops import GaussianLoopSampler, exact_probabilities


def run() -> dict:
    C = np.array([[0.18, 0.23], [0.23, 0.13]])
    B = C * np.exp(1j * np.array([[0.3, 1.4], [1.4, -0.6]]))
    sampler = GaussianLoopSampler(B, seed=17)
    T = sampler.T
    H = np.block([[np.zeros_like(B), B], [B.conj(), np.zeros_like(B)]])
    z = np.array([0.37, 0.81])
    lifted_z = np.tile(z, 2)
    max_length = 8
    loop_total = loop_delta = loop_gamma = phase_log_pgf = positive_log_pgf = 0.0
    enumerated = 0
    for length in range(2, max_length + 1, 2):
        # Enumerate all ordered rooted bipartite walks, including periodic ones.
        trace_C = trace_H = 0j
        for layer in (0, 1):
            for modes in itertools.product(range(2), repeat=length):
                walk = [mode + 2 * ((layer + j) % 2) for j, mode in enumerate(modes)]
                edges = list(zip(walk, walk[1:] + walk[:1]))
                reference = math.prod(T[u, v] for u, v in edges)
                complex_weight = math.prod(H[u, v] for u, v in edges)
                cosine = float((complex_weight / reference).real)
                mass = reference / (2 * length)
                factor = math.prod(lifted_z[v] for v in walk) - 1
                loop_total += mass
                loop_delta += mass * (1 - cosine)
                loop_gamma += mass * (1 - cosine)**2
                phase_log_pgf += mass * cosine * factor
                positive_log_pgf += mass * factor
                trace_C += reference
                trace_H += complex_weight
                enumerated += 1
        assert abs(trace_C - np.trace(np.linalg.matrix_power(T, length))) < 1e-12
        assert abs(trace_H - np.trace(np.linalg.matrix_power(H, length))) < 1e-12

    tail = (sampler.m * sampler.rho ** (max_length + 2) /
            ((max_length + 2) * (1 - sampler.rho**2)))
    assert abs(loop_total - sampler.lambda_C) <= tail + 1e-12
    assert abs(loop_delta - sampler.delta) <= 2 * tail + 1e-12
    assert abs(loop_gamma - sampler.gamma) <= 4 * tail + 1e-12
    Z = np.diag(z)
    log_pgf_B = 0.5 * (np.linalg.slogdet(np.eye(2) - B.conj().T @ B)[1]
                      - np.linalg.slogdet(np.eye(2) - B.conj().T @ Z @ B @ Z)[1])
    log_pgf_C = 0.5 * (np.linalg.slogdet(np.eye(2) - C @ C)[1]
                      - np.linalg.slogdet(np.eye(2) - C @ Z @ C @ Z)[1])
    assert abs(phase_log_pgf - log_pgf_B) <= tail + 1e-12
    assert abs(positive_log_pgf - log_pgf_C) <= tail + 1e-12
    bubble_mass = -0.5 * np.log1p(-sampler.q).sum()
    assert abs(bubble_mass - sampler.lambda_C) < 2e-14

    # Independent closed-form single-mode and pair-squeezed probabilities.
    t = 0.53
    one = exact_probabilities([[t]], 14)
    for k in range(8):
        target = math.sqrt(1 - t*t) * math.comb(2*k, k) * (t*t/4)**k
        assert abs(one[(2*k,)] - target) < 1e-14
    pair = exact_probabilities([[0, t], [t, 0]], 14)
    for n, p in pair.items():
        target = (1 - t*t) * t**(2*n[0]) if n[0] == n[1] else 0
        assert abs(p - target) < 1e-14

    D = np.exp(1j * np.array([0.7, -0.2]))
    gauge = GaussianLoopSampler(C * D[:, None] * D[None, :], seed=191)
    gauge_error = max(gauge.delta, gauge.gamma)
    assert gauge_error < 1e-13
    for _ in range(200):
        _, weight = gauge.draw_cosine()
        assert abs(weight - 1) < 2e-14

    # For real signed B, every cosine is +/-1 and r=Delta/2 exactly.
    signed = GaussianLoopSampler([[0.4, 0.2], [0.2, -0.4]], seed=193)
    estimate = signed.negative_mass_estimate(64)
    exact_r = signed.delta / 2
    assert estimate.lower <= exact_r <= estimate.upper
    assert abs(signed.gamma - 4*exact_r) < 1e-13
    assert 4*exact_r <= signed.gamma + 1e-13

    # Validate model boundaries: a physical B can have unstable abs(B).
    physical_but_unstable_reference = np.array([[0.6, 0.6], [0.6, -0.6]])
    assert np.linalg.norm(physical_but_unstable_reference, 2) < 1
    try:
        GaussianLoopSampler(physical_but_unstable_reference)
    except ValueError:
        pass
    else:
        raise AssertionError("An unstable reference was accepted")
    vacuum = GaussianLoopSampler(np.zeros((3, 3)), seed=1)
    assert np.array_equal(vacuum.sample_reference(), np.zeros(3))
    assert vacuum.draw_cosine()[1] == 1
    return {
        "ordered_rooted_walks_checked": enumerated,
        "maximum_walk_length": max_length,
        "loop_intensity_tail_bound": tail,
        "delta_truncation_error": abs(loop_delta - sampler.delta),
        "gamma_truncation_error": abs(loop_gamma - sampler.gamma),
        "phase_log_pgf_truncation_error": abs(phase_log_pgf - log_pgf_B),
        "bubble_logdet_error": abs(bubble_mass - sampler.lambda_C),
        "gauge_parameter_error": gauge_error,
        "signed_real_negative_mass": exact_r,
        "fourier_truncation_error_bound": estimate.truncation_error,
        "status": "passed",
    }
