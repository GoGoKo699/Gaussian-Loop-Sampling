"""Gaussian photon counts from classical Markov loops.

The implementation uses floating-point arithmetic. See science/finite-precision.md
for a separate constructive bit-complexity theorem.
"""

from .sampler import GaussianLoopSampler, Loop, NegativeMassEstimate, WorkLimitExceeded
from .exact import exact_probabilities

__all__ = [
    "GaussianLoopSampler", "Loop", "NegativeMassEstimate", "WorkLimitExceeded",
    "exact_probabilities",
]
