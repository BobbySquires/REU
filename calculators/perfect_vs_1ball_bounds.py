"""
Volume lower bound on the burning number of K(2k+1, k).

For a burning sequence of length b, the b sources have radii b-1, b-2, ..., 0
at the final turn. Their total ball capacity must cover all |V| vertices, so:

    β(0) + β(1) + ... + β(b-1) >= C(2k+1, k)

Defining F(t) = β(0) + ... + β(t) - C(2k+1, k), this condition is F(b-1) >= 0.

Since F is strictly increasing (F(t) - F(t-1) = β(t) > 0), F(0) < 0, and
F(k) > 0, there is a unique t0 in {0, ..., k-1} where F first becomes
non-negative. Then b-1 >= t0+1, so b >= t0+2.

The trivial upper bound is b <= k+1: a single source lit at turn 1 has k
turns to spread by turn k+1, and since the diameter of K(2k+1, k) is k it
reaches every vertex. So one source suffices in k+1 turns.
"""

from math import comb


def sigma(k: int, r: int) -> int:
    """Shell size σ(r) = C(k, floor(r/2)) * C(k+1, ceil(r/2))."""
    if r < 0 or r > k:
        return 0
    return comb(k, r // 2) * comb(k + 1, (r + 1) // 2)


def beta(k: int, r: int) -> int:
    """Ball size β(r) = sum of σ(0)..σ(r)."""
    if r < 0:
        return 0
    if r > k:
        return comb(2 * k + 1, k)
    return sum(sigma(k, l) for l in range(r + 1))


def F(k: int, t: int) -> int:
    """F(t) = β(0) + β(1) + ... + β(t) - |V|."""
    return sum(beta(k, r) for r in range(t + 1)) - comb(2 * k + 1, k)


def find_t0(k: int) -> int:
    """Find t0, the last integer in {0,...,k-1} with F(t0) < 0."""
    for t in range(k - 1, -1, -1):
        if F(k, t) < 0:
            return t
    return -1


def print_table(k: int, t_max: int) -> None:
    n_verts = comb(2 * k + 1, k)
    t0 = find_t0(k)
    lb = t0 + 2
    ub = k + 1

    print(f"\nK({2*k+1}, {k})   |V| = C({2*k+1},{k}) = {n_verts}   diameter = {k}")
    print(f"Upper bound: b <= {ub}   |   Volume lower bound: b >= {lb}")
    print()
    print("  Columns:")
    print("    t      : largest ball has radius t, which corresponds with a perfect burning sequence for t+1 rounds")
    print("    β(t)   : size of a ball of radius t")
    print("    Σβ     : β(0) + β(1) + ... + β(t)  (total ball capacity at time t+1)")
    print("    F(t)   : Σβ - |V|  (negative = graph not yet covered, zero/positive = covered)")
    print()
    print(f"  {'t':>4}  {'β(t)':>12}  {'Σβ':>14}  {'F(t)':>14}  {'sign':>6}")
    print(f"  {'-'*4}  {'-'*12}  {'-'*14}  {'-'*14}  {'-'*6}")

    cumsum = 0
    for t in range(min(t_max, k) + 1):
        b = beta(k, t)
        cumsum += b
        f_val = cumsum - n_verts
        sign = ">= 0" if f_val >= 0 else "< 0"
        marker = "  <-- t0" if t == t0 else ""
        print(f"  {t:>4}  {b:>12}  {cumsum:>14}  {f_val:>14}  {sign:>6}{marker}")

    print()
    print(f"  t0 = {t0}  (last t with F(t) < 0)")
    print(f"  Lower bound: b >= t0 + 2 = {lb}")
    print(f"  Upper bound: b <= k + 1  = {ub}")
    if lb == ub:
        print(f"  Bounds match: b(K({2*k+1},{k})) = {lb}")


if __name__ == "__main__":
    k = int(input("Enter k (graph parameter for K(2k+1, k)): "))
    t_max = int(input(f"Enter t_max to evaluate F up to (max meaningful value is k={k}): "))
    print_table(k, t_max)
