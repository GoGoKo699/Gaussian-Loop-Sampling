"""Independent graph, determinant, and trace checks for science/topology.md.

Only small graph matrices are dense. The 2,184-vertex LPS construction uses
an adjacency list; it performs no large determinant or eigendecomposition.
"""

from __future__ import annotations

from collections import deque
import math
import json
import numpy as np


def adjacency(n, edges):
    a = np.zeros((n, n))
    for i, j in edges:
        a[i, j] = a[j, i] = 1
    return a


def symmetric_phases(a, rng):
    b = np.zeros(a.shape, dtype=complex)
    for i, j in zip(*np.nonzero(np.triu(a, 1))):
        b[i, j] = b[j, i] = a[i, j] * np.exp(1j * rng.uniform(-math.pi, math.pi))
    return b


def lift(b):
    z = np.zeros_like(b)
    return np.block([[z, b], [b.conj().T, z]])


def nonbacktracking(h):
    edges = list(zip(*np.nonzero(h)))
    n = np.zeros((len(edges), len(edges)), dtype=complex)
    for ei, (a, b) in enumerate(edges):
        for fi, (c, d) in enumerate(edges):
            if b == c and a != d:
                n[ei, fi] = h[c, d]
    return n


def real_logdet(a):
    sign, logabs = np.linalg.slogdet(a)
    assert abs(sign - 1) < 2e-12, sign
    return logabs


def regular_checks(a, rng):
    m = len(a)
    d = int(a.sum(axis=1)[0])
    assert np.all(a.sum(axis=1) == d)
    h0, hp = lift(a), lift(symmetric_phases(a, rng))
    n0, np_ = nonbacktracking(h0), nonbacktracking(hp)
    vertices, undirected_edges = 2 * m, m * d
    max_bass = max_delta = max_trace = 0.0
    for t in [0.2, 0.5, 0.8]:
        s = t / d
        u = 2 * s / (1 + math.sqrt(1 - 4 * (d - 1) * s * s))
        p = d - 1
        assert p * u < 1
        for h, n in [(h0, n0), (hp, np_)]:
            lhs = np.linalg.det(np.eye(len(n)) - u * n)
            rhs = (1 - u * u) ** (undirected_edges - vertices) * np.linalg.det(
                np.eye(vertices) - u * h + p * u * u * np.eye(vertices)
            )
            max_bass = max(max_bass, abs(lhs - rhs) / max(1, abs(lhs), abs(rhs)))
        delta_direct = 0.5 * (
            real_logdet(np.eye(vertices) - s * hp)
            - real_logdet(np.eye(vertices) - s * h0)
        )
        delta_bass = 0.5 * (
            real_logdet(np.eye(len(np_)) - u * np_)
            - real_logdet(np.eye(len(n0)) - u * n0)
        )
        max_delta = max(max_delta, abs(delta_bass - delta_direct))
        r0, rp = np.eye(len(n0)), np.eye(len(np_))
        partial = 0j
        for length in range(1, 181):
            r0, rp = r0 @ (u * n0), rp @ (u * np_)
            term = 0.5 * (np.trace(r0) - np.trace(rp)) / length
            assert term.real > -2e-12
            assert abs(term.imag) < 2e-12
            partial += term
        max_delta = max(max_delta, abs(delta_direct - partial))

    # Recover all edge-operator power traces from adjacency eigenvalues.
    eigenvalues = np.linalg.eigvalsh(h0)
    roots = np.array([
        root
        for lam in eigenvalues
        for root in np.roots([1, -lam, d - 1])
    ]) / (d - 1)
    current = np.eye(len(n0))
    for length in range(1, 15):
        current = current @ (n0 / (d - 1))
        predicted = (np.sum(roots ** length)
                     + (undirected_edges - vertices) * (1 + (-1) ** length) / (d - 1) ** length)
        observed = np.trace(current)
        max_trace = max(max_trace, abs(predicted - observed) / max(1, abs(observed)))
    assert max_bass < 2e-12, max_bass
    assert max_delta < 2e-12, max_delta
    assert max_trace < 2e-12, max_trace
    return {
        "relative_bass_determinant_error": float(max_bass),
        "absolute_delta_identity_error": float(max_delta),
        "normalized_power_trace_error": float(max_trace),
    }


def unequal_weight_check(rng):
    a = adjacency(5, [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 1)])
    b = symmetric_phases(a, rng)
    for i, j in zip(*np.nonzero(np.triu(a, 1))):
        b[i, j] *= rng.uniform(0.1, 0.9)
        b[j, i] = b[i, j]
    h = lift(b)
    n = nonbacktracking(h)
    max_error = 0
    for u in [0.1, 0.3, 0.55]:
        vertex = np.eye(len(h), dtype=complex)
        edge_product = 1.0
        for i, j in zip(*np.nonzero(np.triu(abs(h), 1))):
            mag2 = abs(h[i, j]) ** 2
            denominator = 1 - u * u * mag2
            edge_product *= denominator
            vertex[i, j] = -u * h[i, j] / denominator
            vertex[j, i] = -u * h[j, i] / denominator
            vertex[i, i] += u * u * mag2 / denominator
            vertex[j, j] += u * u * mag2 / denominator
        lhs = np.linalg.det(np.eye(len(n)) - u * n)
        rhs = edge_product * np.linalg.det(vertex)
        max_error = max(max_error, abs(lhs - rhs) / max(1, abs(lhs), abs(rhs)))
    assert max_error < 2e-12, max_error
    return float(max_error)


def lps_5_13_check():
    ell = 13
    iota = next(x for x in range(ell) if x * x % ell == ell - 1)
    generators = []
    for sign in [-1, 1]:
        generators.extend([
            (1 + sign * 2 * iota, 0, 0, 1 - sign * 2 * iota),
            (1, sign * 2, -sign * 2, 1),
            (1, sign * 2 * iota, sign * 2 * iota, 1),
        ])

    def normalize(g):
        g = tuple(x % ell for x in g)
        inverse = pow(next(x for x in g if x), -1, ell)
        return tuple(x * inverse % ell for x in g)

    def multiply(g, h):
        a, b, c, d = g
        e, f, k, l = h
        return normalize((a * e + b * k, a * f + b * l, c * e + d * k, c * f + d * l))

    generators = [normalize(g) for g in generators]
    identity = normalize((1, 0, 0, 1))
    index = {identity: 0}
    vertices = [identity]
    adjacency_list = []
    for g in vertices:
        neighbors = []
        for generator in generators:
            h = multiply(g, generator)
            if h not in index:
                index[h] = len(vertices)
                vertices.append(h)
            neighbors.append(index[h])
        adjacency_list.append(neighbors)
    m = ell * (ell * ell - 1)
    assert len(vertices) == m
    assert all(len(set(row)) == 6 and v not in row for v, row in enumerate(adjacency_list))
    assert all(v in adjacency_list[w] for v, row in enumerate(adjacency_list) for w in row)
    colors = [-1] * m
    colors[0] = 0
    queue = deque([0])
    while queue:
        v = queue.popleft()
        for w in adjacency_list[v]:
            if colors[w] == -1:
                colors[w] = 1 - colors[v]
                queue.append(w)
            assert colors[w] == 1 - colors[v]
    assert -1 not in colors
    # Cayley vertex-transitivity allows a single-root shortest-cycle calculation.
    distance = [-1] * m
    parent = [-1] * m
    distance[0] = 0
    queue = deque([0])
    girth = math.inf
    while queue:
        v = queue.popleft()
        for w in adjacency_list[v]:
            if distance[w] == -1:
                parent[w] = v
                distance[w] = distance[v] + 1
                queue.append(w)
            elif parent[v] != w and parent[w] != v:
                girth = min(girth, distance[v] + distance[w] + 1)
    assert girth == 8, girth
    return {
        "prime": ell,
        "vertices": m,
        "degree": 6,
        "connected": True,
        "bipartite": True,
        "girth": girth,
    }


def finite_size_certificates() -> dict:
    """Evaluate the proof's finite-instance bounds, not empirical distances."""
    d, p, t = 6, 5, 0.5
    s = t / d
    u = 2 * s / (1 + math.sqrt(1 - 4 * p * s * s))
    q = p * u
    exponent = -1 / 3 - 4 * math.log(q) / (3 * math.log(p))
    q_threshold = p ** (-0.25)
    t_threshold = d * q_threshold / (p + q_threshold * q_threshold)
    rounded_tv_bounds = {13: 0.081927, 53: 0.004609, 113: 0.001377}
    table = []
    for ell, rounded_bound in rounded_tv_bounds.items():
        m = ell * (ell * ell - 1)
        h = 2 * math.ceil((4 * math.log(ell, p) - math.log(4, p)) / 2)
        delta_bound = (
            4 * q ** h / (1 - q * q)
            + 4 * m * (q / math.sqrt(p)) ** h / (1 - q * q / p)
            + (2 * m * (d - 2) + 4) * u ** h / (1 - u * u)
        ) / h
        tv_bound = -math.expm1(-delta_bound)
        assert 0 < tv_bound <= rounded_bound
        assert rounded_bound - tv_bound < 1e-6
        table.append({
            "prime": ell,
            "modes": m,
            "guaranteed_girth_lower_bound": h,
            "delta_upper_bound": delta_bound,
            "tv_upper_bound": tv_bound,
            "rounded_up_tv_upper_bound": rounded_bound,
        })
    assert abs(t_threshold - 0.7366044601543223) < 1e-14
    assert abs(exponent - 0.3615547196555067) < 1e-14
    return {
        "degree": d,
        "coupling_t": t,
        "nonbacktracking_q": q,
        "inverse_power_exponent": exponent,
        "strict_inverse_power_coupling_threshold": t_threshold,
        "table": table,
    }


def run() -> dict:
    rng = np.random.default_rng(922723)
    graphs = {
        "C6": adjacency(6, [(i, (i + 1) % 6) for i in range(6)]),
        "K3,3": adjacency(6, [(i, j) for i in range(3) for j in range(3, 6)]),
        "Petersen": adjacency(10, [(i, (i + 1) % 5) for i in range(5)]
                                + [(i, i + 5) for i in range(5)]
                                + [(i + 5, ((i + 2) % 5) + 5) for i in range(5)]),
    }
    results = {name: regular_checks(a, rng) for name, a in graphs.items()}
    return {
        "small_graph_checks": results,
        "unequal_weight_bass_relative_error": unequal_weight_check(rng),
        "explicit_lps_graph": lps_5_13_check(),
        "finite_size_certificates": finite_size_certificates(),
        "status": "passed",
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))

