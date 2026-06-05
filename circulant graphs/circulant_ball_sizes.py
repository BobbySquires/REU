"""
Compute ball sizes beta(v; r) = |B(v; r)| for the circulant graph C_{p^3}(1, p),
where p is prime and p < 100.

Graph: C_{p^3}(1, p) on n = p^3 vertices (0..n-1).
Edges: i ~ i±1 (mod n) and i ~ i±p (mod n).
Vertex-transitive, so use vertex 0 as center.
"""

from collections import deque


def sieve(limit):
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for j in range(i * i, limit + 1, i):
                is_prime[j] = False
    return [i for i in range(2, limit + 1) if is_prime[i]]


def bfs_from_zero(n, p):
    dist = [-1] * n
    dist[0] = 0
    q = deque([0])
    while q:
        v = q.popleft()
        d = dist[v]
        for nb in ((v + 1) % n, (v - 1) % n, (v + p) % n, (v - p) % n):
            if dist[nb] == -1:
                dist[nb] = d + 1
                q.append(nb)
    return dist


def compute_ball_sizes(p):
    n = p ** 3
    dist = bfs_from_zero(n, p)
    diam = max(dist)
    sigma = [0] * (diam + 1)
    for d in dist:
        sigma[d] += 1
    beta = []
    running = 0
    for r in range(diam + 1):
        running += sigma[r]
        beta.append(running)
    return diam, sigma, beta


def print_table(p, diam, sigma, beta):
    n = p ** 3
    print(f"p = {p}, n = p^3 = {n}, diameter = {diam}")
    print(f"{'r':<6}{'sigma(r)':<12}{'beta(r)':<12}")
    for r in range(diam + 1):
        print(f"{r:<6}{sigma[r]:<12}{beta[r]:<12}")
    assert beta[diam] == n, f"beta(diam) = {beta[diam]} != n = {n}"
    print()


def main():
    primes = sieve(99)

    all_results = []

    for p in primes:
        diam, sigma, beta = compute_ball_sizes(p)
        print_table(p, diam, sigma, beta)
        all_results.append((p, p ** 3, diam, sigma, beta))

    # Summary table: show beta(r) for small r up to the max diameter across all p,
    # but cap display at r where all graphs have beta(r) = n (or max diam).
    max_diam = max(r[2] for r in all_results)

    # Header
    header_cols = ["p", "n=p^3", "diam"] + [f"beta({r})" for r in range(max_diam + 1)]
    col_widths = [6, 10, 6] + [max(8, len(f"beta({r})") + 1) for r in range(max_diam + 1)]

    def fmt_row(vals):
        return "".join(str(v).ljust(w) for v, w in zip(vals, col_widths))

    print("=" * 80)
    print("SUMMARY TABLE")
    print("=" * 80)
    print(fmt_row(header_cols))
    print("-" * sum(col_widths[:3 + max_diam + 1]))

    for p, n, diam, sigma, beta in all_results:
        # Pad beta with n for r > diam
        beta_padded = beta + [n] * (max_diam - diam)
        row = [p, n, diam] + beta_padded
        print(fmt_row(row))


if __name__ == "__main__":
    main()
