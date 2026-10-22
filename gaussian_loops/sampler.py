"""Reference implementation of the loop constructions in science/.

All samples use genuine Poisson loops: a Poisson number of logarithmic-series
clusters, each containing first-return excursions. Cosines are applied once
per complete cluster. Applying them to individual excursions changes the
second-moment formula.

This is numerical reference code, not interval-certified or exact random-bit
software. In particular, analytical tail certificates below exclude floating-
point rounding error. No hafnian is used by this sampler.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterator

import numpy as np


class WorkLimitExceeded(RuntimeError):
    """An optional per-draw work limit was exceeded; no sample is returned."""


@dataclass(frozen=True)
class Loop:
    """Projected vertex visits and holonomy of one complete Poisson loop."""
    visits: tuple[int, ...]
    phase: float

    @property
    def cosine(self) -> float:
        return math.cos(self.phase)


@dataclass(frozen=True)
class NegativeMassEstimate:
    """Fourier estimate of r; truncation_error does not include roundoff."""
    value: float
    truncation_error: float
    lower: float
    upper: float
    harmonics: int


def _log_partition(matrix: np.ndarray) -> float:
    singular = np.linalg.svd(matrix, compute_uv=False)
    if np.any(singular >= 1):
        raise ArithmeticError("The matrix is numerically outside the stable domain")
    return float(-0.5 * np.log1p(-singular * singular).sum())


def _nonnegative_roundoff(value: float, scale: float, name: str) -> float:
    if value < -1e-10 * max(1.0, scale):
        raise ArithmeticError(f"Numerically negative {name}: {value}")
    return max(0.0, value)


class GaussianLoopSampler:
    """Prepare a loop soup for symmetric B with stable C = abs(B).

    sample_reference() samples P_C. For complex B, draw_cosine() estimates
    bounded observables of P_B, while sample_positive_part() samples the
    explicitly approximate distribution Q_+. Neither API mislabels Q_+ or
    P_C as an exact draw from P_B.

    Matrix symmetry is checked to absolute tolerance 1e-12; accepted tiny
    discrepancies are symmetrized. stability_tol is a numerical guard, not
    an analytical condition. max_steps, when set, raises on excessive work;
    repeatedly discarding failed draws would condition and bias the result.
    """

    def __init__(self, B, seed: int | None = None, *,
                 stability_tol: float = 1e-12, max_steps: int | None = None):
        matrix = np.asarray(B, dtype=np.complex128)
        if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or not len(matrix):
            raise ValueError("B must be a nonempty square matrix")
        if not np.all(np.isfinite(matrix)):
            raise ValueError("B must have finite entries")
        if not np.allclose(matrix, matrix.T, rtol=0, atol=1e-12):
            raise ValueError("B must be complex symmetric, B = B.T")
        if not 0 < stability_tol < 1:
            raise ValueError("stability_tol must lie strictly between zero and one")
        if max_steps is not None and (not isinstance(max_steps, int) or max_steps < 1):
            raise ValueError("max_steps must be a positive integer or None")
        self.B = (matrix + matrix.T) / 2
        self.C = np.abs(self.B)
        self.m = len(matrix)
        self.rho = float(np.linalg.eigvalsh(self.C)[-1])
        if self.rho >= 1 - stability_tol:
            raise ValueError("Entrywise magnitude C must be stably below norm one")
        self.max_steps = max_steps
        self.rng = np.random.default_rng(seed)

        zero = np.zeros_like(self.C)
        self.T = np.block([[zero, self.C], [self.C, zero]])
        lifted_B = np.block([[zero, self.B], [self.B.conj(), zero]])
        self.edge_angles = np.angle(lifted_B)
        self.dimension = 2 * self.m
        self.neighbors = [np.flatnonzero(row > 0) for row in self.T]
        self.h = np.zeros_like(self.T)
        self.q = np.zeros(self.dimension)

        # Reverse Schur-complement construction. G is the inverse on the
        # current trailing principal submatrix. Dense arithmetic O(m^3).
        G = np.empty((0, 0))
        for i in reversed(range(self.dimension)):
            v = self.T[i + 1:, i]
            g = G @ v
            q = float(self.T[i, i] + v @ g)
            if q < -1e-12 or q >= 1:
                raise ArithmeticError("Unstable first-return probability")
            self.q[i] = max(0.0, q)
            self.h[i, i] = 1.0
            self.h[i, i + 1:] = g
            s = 1.0 - q
            G = np.block([
                [np.array([[1 / s]]), (g / s)[None, :]],
                [(g / s)[:, None], G + np.outer(g, g) / s],
            ])

        self.lambda_C = _log_partition(self.C)
        self.lambda_B = _log_partition(self.B)
        self.delta = _nonnegative_roundoff(
            self.lambda_C - self.lambda_B, self.lambda_C, "Delta")
        doubled_phase = self.C * np.exp(2j * np.angle(self.B))
        self.lambda_doubled_phase = _log_partition(doubled_phase)
        self.gamma = _nonnegative_roundoff(
            2 * self.delta - 0.5 * (self.lambda_C - self.lambda_doubled_phase),
            self.lambda_C, "Gamma")
        singular_C = np.linalg.svd(self.C, compute_uv=False)
        singular_B = np.linalg.svd(self.B, compute_uv=False)
        self.mean_C = float(np.sum(singular_C**2 / (1 - singular_C**2)))
        self.mean_B = float(np.sum(singular_B**2 / (1 - singular_B**2)))

    def iter_loops(self) -> Iterator[Loop]:
        """Generate a fresh independent reference soup, one complete loop at a time.

        Visits are sparse lists. Neighbor scans cost O(d) per visit when each
        lifted vertex has at most d neighbors. Memory for a realized loop is
        proportional to its length. Exhaust this iterator to obtain one soup.
        """
        steps = 0
        for i in range(self.dimension):
            q = self.q[i]
            if q == 0:
                continue
            clusters = int(self.rng.poisson(-0.5 * math.log1p(-q)))
            h = self.h[i]
            for _ in range(clusters):
                repeats = int(self.rng.logseries(q))
                visits: list[int] = []
                angle = 0.0
                for _ in range(repeats):
                    current = i
                    while True:
                        steps += 1
                        if self.max_steps is not None and steps > self.max_steps:
                            raise WorkLimitExceeded("The draw exceeded max_steps")
                        visits.append(current % self.m)
                        adjacent = self.neighbors[current]
                        eligible = adjacent[np.searchsorted(adjacent, i):]
                        weights = self.T[current, eligible] * h[eligible]
                        if np.any(weights < -1e-12):
                            raise ArithmeticError("Negative transition weight")
                        cumulative = np.cumsum(np.maximum(weights, 0.0))
                        if not len(cumulative) or cumulative[-1] <= 0:
                            raise ArithmeticError("A bridge has no positive continuation")
                        threshold = self.rng.random() * cumulative[-1]
                        index = min(int(np.searchsorted(cumulative, threshold, side="right")),
                                    len(eligible) - 1)
                        nxt = int(eligible[index])
                        angle = math.remainder(angle + self.edge_angles[current, nxt],
                                               2 * math.pi)
                        if nxt == i:
                            break
                        current = nxt
                yield Loop(tuple(visits), angle)

    @staticmethod
    def _accumulate(counts: np.ndarray, loop: Loop) -> None:
        np.add.at(counts, np.asarray(loop.visits, dtype=np.intp), 1)

    def sample_reference(self) -> np.ndarray:
        """Draw P_C (the target law when B is nonnegative or gauge equivalent)."""
        counts = np.zeros(self.m, dtype=np.int64)
        for loop in self.iter_loops():
            self._accumulate(counts, loop)
        return counts

    def draw_cosine(self) -> tuple[np.ndarray, float]:
        """Return (reference counts, signed weight) for E_B f = E[weight*f].

        The ideal second moment is exp(Gamma). A logarithmic product prevents
        premature overflow in exp(Delta); the final float may still overflow
        when an individual weight is outside the floating-point range.
        """
        counts = np.zeros(self.m, dtype=np.int64)
        sign = 1
        log_weight = self.delta
        for loop in self.iter_loops():
            self._accumulate(counts, loop)
            cosine = loop.cosine
            if cosine == 0:
                sign = 0
                log_weight = -math.inf
            elif sign:
                sign *= 1 if cosine > 0 else -1
                log_weight += math.log(abs(cosine))
        return counts, float(sign * math.exp(log_weight))

    def sample_positive_part(self, *, maximum_length: int | None = None) -> np.ndarray:
        """Draw Q_+ by retaining loops with probability max(0, cos(phi)).

        If maximum_length is supplied, also discard longer loops. The full
        Q_+ has TV error at most min(1, expm1(2*r)/2). The truncated law needs
        the additional tail certificate in science/estimators.md. This code
        still generates the discarded loops and so is not a runtime cutoff.
        """
        if maximum_length is not None and (
                not isinstance(maximum_length, int) or maximum_length < 0):
            raise ValueError("maximum_length must be a nonnegative integer or None")
        counts = np.zeros(self.m, dtype=np.int64)
        for loop in self.iter_loops():
            if maximum_length is not None and len(loop.visits) > maximum_length:
                continue
            if self.rng.random() < max(0.0, loop.cosine):
                self._accumulate(counts, loop)
        return counts

    def draw_signed(self, *, negative_mass: float) -> tuple[np.ndarray, float]:
        """Thin by abs(cos(phi)) and return a sign with magnitude exp(2*r).

        The caller supplies r. Replacing r by r_hat multiplies the expected
        answer by exp(2*(r_hat-r)); it is not exactly unbiased. Use the Fourier
        estimate and its truncation bound to budget this error. In exact
        arithmetic the second moment with the true r is exp(4*r).
        """
        if not math.isfinite(negative_mass) or negative_mass < 0:
            raise ValueError("negative_mass must be a finite nonnegative number")
        counts = np.zeros(self.m, dtype=np.int64)
        sign = 1
        for loop in self.iter_loops():
            cosine = loop.cosine
            if self.rng.random() < abs(cosine):
                self._accumulate(counts, loop)
                if cosine < 0:
                    sign = -sign
        return counts, float(sign * math.exp(2 * negative_mass))

    def negative_mass_estimate(self, harmonics: int = 32) -> NegativeMassEstimate:
        """Compute r from K determinant harmonics, with a Fourier-tail bound.

        O(K*m^3) dense arithmetic. The bound excludes numerical determinant
        and phase errors; use certified arithmetic for a rigorous certificate.
        """
        if not isinstance(harmonics, int) or harmonics < 1:
            raise ValueError("harmonics must be a positive integer")
        total = 2 * self.lambda_C / math.pi
        angle = np.angle(self.B)
        for k in range(1, harmonics + 1):
            harmonic = self.C * np.exp((2j * k) * angle)
            total += (4 / math.pi) * ((-1) ** (k - 1)) * _log_partition(harmonic) / (4*k*k - 1)
        value = 0.5 * (total - self.lambda_B)
        error = self.lambda_C / (math.pi * (2 * harmonics + 1))
        return NegativeMassEstimate(value, error, max(0.0, value - error),
                                    min(self.lambda_C, max(0.0, value + error)), harmonics)
