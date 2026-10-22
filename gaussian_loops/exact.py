"""Small-system hafnian probabilities, independent of the loop sampler.

Exponential reference calculation only; never used by the sampling algorithm.
"""
from functools import lru_cache
import math
import numpy as np


def exact_probabilities(B, max_photons: int) -> dict[tuple[int, ...], float]:
    """Enumerate even sectors through max_photons, without renormalizing the tail."""
    matrix = np.asarray(B, dtype=complex)
    if (matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or not len(matrix)
            or not np.all(np.isfinite(matrix))
            or not np.allclose(matrix, matrix.T, atol=1e-12, rtol=0)):
        raise ValueError("B must be a finite, nonempty complex symmetric matrix")
    if not isinstance(max_photons, int) or max_photons < 0:
        raise ValueError("max_photons must be a nonnegative integer")
    singular = np.linalg.svd(matrix, compute_uv=False)
    if np.any(singular >= 1):
        raise ValueError("B must have norm less than one")
    normalization = float(np.exp(0.5 * np.log1p(-singular**2).sum()))
    m = len(matrix)

    @lru_cache(None)
    def haf(n: tuple[int, ...]) -> complex:
        if not any(n):
            return 1 + 0j
        i = next(i for i, value in enumerate(n) if value)
        remaining = list(n)
        remaining[i] -= 1
        value = 0j
        for j, multiplicity in enumerate(remaining):
            if multiplicity:
                remaining[j] -= 1
                value += multiplicity * matrix[i, j] * haf(tuple(remaining))
                remaining[j] += 1
        return value

    def compositions(total: int, length: int):
        if length == 1:
            yield (total,)
        else:
            for first in range(total + 1):
                for rest in compositions(total - first, length - 1):
                    yield (first,) + rest

    result = {}
    for total in range(0, max_photons + 1, 2):
        for n in compositions(total, m):
            result[n] = float(normalization * abs(haf(n))**2 /
                              math.prod(math.factorial(x) for x in n))
    return result
