"""
Two-source upper bound for burning K(2k+1, k).

The trivial upper bound is b(K(2k+1,k)) <= k+1 (one source, diameter k).
This script computes a tighter upper bound by considering two sources A, B
burned at turns 1 and 2. After b turns, their balls have radii b-1 and b-2,
covering:

    |B(A; t) ∪ B(B; t-1)| = β(t) + β(t-1) - |B(A;t) ∩ B(B;t-1)|

where t = b-1. To maximize coverage (and thus minimize the number of turns),
we choose A, B to minimize the intersection. The intersection size depends
on s = |A ∩ B|, and we define:

    I(s, t) = |B(A; t) ∩ B(B; t-1)|   for |A ∩ B| = s.

The two-source upper bound is then b <= t* + 1, where

    t* = min{ t : β(t) + β(t-1) - min_s I(s, t) >= C(2k+1, k) }.

References:
    See burning-kneser.pdf, "Tightening the upper bound" section.
"""

from math import comb


# ---------------------------------------------------------------------------
# Ball and shell sizes (from odd_kneser_ball_formula.py)
# ---------------------------------------------------------------------------

def sigma(k: int, r: int) -> int:
    """Shell size |∂B(v; r)| in K(2k+1, k)."""
    if r < 0 or r > k:
        return 0
    return comb(k, r // 2) * comb(k + 1, (r + 1) // 2)


def beta(k: int, r: int) -> int:
    """Ball size |B(v; r)| in K(2k+1, k)."""
    if r < 0:
        return 0
    if r > k:
        return comb(2 * k + 1, k)
    return sum(sigma(k, l) for l in range(r + 1))


# ---------------------------------------------------------------------------
# Distance in K(2k+1, k)
# ---------------------------------------------------------------------------

def dist_K(k: int, overlap: int) -> int:
    """
    Distance between two vertices in K(2k+1, k) with given overlap.

    By Lemma 1: d(A, B) = min{2(k - s), 2s + 1} where s = |A ∩ B|.
    """
    return min(2 * (k - overlap), 2 * overlap + 1)


# ---------------------------------------------------------------------------
# Intersection count I(s, t)
# ---------------------------------------------------------------------------

def intersection_count(k: int, s: int, t: int) -> int:
    """
    Count |B(A; t) ∩ B(B; t-1)| in K(2k+1, k) where |A ∩ B| = s.

    A vertex C in the intersection has k elements partitioned as:
        x   elements from A ∩ B         (pool size: s)
        sA  elements from A \\ B         (pool size: k - s)
        sB  elements from B \\ A         (pool size: k - s)
        rem elements from (A ∪ B)^c     (pool size: s + 1)

    where rem = k - x - sA - sB, and we require:
        |A ∩ C| = x + sA  =>  d(A, C) <= t
        |B ∩ C| = x + sB  =>  d(B, C) <= t - 1
        0 <= rem <= s + 1
    """
    total = 0
    ks = k - s  # size of A \ B and B \ A

    for x in range(s + 1):
        cx = comb(s, x)

        for sA in range(ks + 1):
            # Check d(A, C) <= t
            if dist_K(k, x + sA) > t:
                continue
            csA = comb(ks, sA)

            for sB in range(ks + 1):
                # Check d(B, C) <= t - 1
                if dist_K(k, x + sB) > t - 1:
                    continue

                rem = k - x - sA - sB
                if rem < 0 or rem > s + 1:
                    continue

                total += cx * csA * comb(ks, sB) * comb(s + 1, rem)

    return total


# ---------------------------------------------------------------------------
# Find s* minimizing I(s, t)
# ---------------------------------------------------------------------------

def find_s_star(k: int, t: int) -> tuple[int, int]:
    """
    Find s* = argmin_s I(s, t) over s in {0, 1, ..., k}.

    Returns (s_star, I_min).
    """
    best_s = 0
    best_val = intersection_count(k, 0, t)

    for s in range(1, k + 1):
        val = intersection_count(k, s, t)
        if val < best_val:
            best_val = val
            best_s = s

    return best_s, best_val


# ---------------------------------------------------------------------------
# Volume lower bound (Theorem 2)
# ---------------------------------------------------------------------------

def volume_lower_bound(k: int) -> int:
    """
    Compute the volume lower bound on b(K(2k+1, k)).

    Find t0 = largest non-negative integer with F(t0) < 0, where
    F(t) = sum_{r=0}^{t} β(r) - C(2k+1, k).

    Then b >= t0 + 2.
    """
    target = comb(2 * k + 1, k)
    t0 = -1
    cumsum = 0
    for t in range(k + 1):
        cumsum += beta(k, t)
        if cumsum - target < 0:
            t0 = t
    return t0 + 2


# ---------------------------------------------------------------------------
# Two-source upper bound
# ---------------------------------------------------------------------------

def two_source_upper_bound(k: int) -> tuple[int, int, int]:
    """
    Compute the two-source upper bound on b(K(2k+1, k)).

    Find t* = min{ t : β(t) + β(t-1) - I(s*, t) >= C(2k+1, k) }.
    The burning number satisfies b <= t* + 1.

    Returns (upper_bound, s_star, t_star).
    """
    target = comb(2 * k + 1, k)

    for t in range(1, k + 1):
        s_star, I_min = find_s_star(k, t)
        coverage = beta(k, t) + beta(k, t - 1) - I_min
        if coverage >= target:
            return t + 1, s_star, t

    # Fallback: trivial bound
    return k + 1, None, k


# ---------------------------------------------------------------------------
# Analysis at t = k-1 (the critical threshold)
# ---------------------------------------------------------------------------

def analyze_critical_threshold(k: int) -> dict:
    """
    Analyze the two-source strategy at t = k-1.

    This is the critical case: if two balls of radii k-1 and k-2 can cover
    the graph, then b <= k (improving on the trivial k+1).

    Returns a dict with all relevant quantities.
    """
    target = comb(2 * k + 1, k)
    t = k - 1

    bk1 = beta(k, t)       # β(k-1)
    bk2 = beta(k, t - 1)   # β(k-2)
    raw = bk1 + bk2        # max possible coverage (zero overlap)
    slack = raw - target    # how much overlap we can tolerate

    result = {
        'k': k, 'target': target, 'beta_k1': bk1, 'beta_k2': bk2,
        'raw': raw, 'slack': slack,
        's_star': None, 'I_min': None, 'ratio': None, 'improves': False,
    }

    if slack <= 0:
        # Two balls can't even theoretically cover the graph
        return result

    # Find s* minimizing I(s, k-1)
    s_star, I_min = find_s_star(k, t)
    ratio = I_min / slack

    result.update({
        's_star': s_star, 'I_min': I_min,
        'ratio': ratio, 'improves': ratio < 1,
    })
    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    k_max = int(input("Enter max k (e.g. 20): "))
    print()

    # --- Table 1: bounds comparison ---
    print("=" * 62)
    print("TABLE 1: Burning number bounds for K(2k+1, k)")
    print("=" * 62)
    header = f"{'k':>4}  {'|V|':>12}  {'vol LB':>6}  {'triv UB':>7}  {'gap':>4}"
    print(header)
    print("-" * len(header))

    for k in range(1, k_max + 1):
        n_verts = comb(2 * k + 1, k)
        lb = volume_lower_bound(k)
        ub = k + 1
        gap = ub - lb

        print(f"{k:>4}  {n_verts:>12}  {lb:>6}  {ub:>7}  {gap:>4}")

    # --- Table 2: two-source analysis at t = k-1 ---
    print()
    print("=" * 72)
    print("TABLE 2: Two-source analysis at t = k-1")
    print("  Can two balls of radii k-1 and k-2 cover K(2k+1, k)?")
    print("  Need: I(s*, k-1) < slack = β(k-1) + β(k-2) - C(2k+1, k)")
    print("=" * 72)
    header2 = f"{'k':>4}  {'slack':>12}  {'s*':>4}  {'I_min':>12}  {'ratio':>8}  {'verdict':>16}"
    print(header2)
    print("-" * len(header2))

    for k in range(1, k_max + 1):
        info = analyze_critical_threshold(k)

        if info['slack'] <= 0:
            print(f"{k:>4}  {'n/a':>12}  {'--':>4}  {'--':>12}  {'--':>8}  "
                  f"{'impossible':>16}")
        else:
            verdict = "IMPROVES!" if info['improves'] else "no improvement"
            print(
                f"{k:>4}  {info['slack']:>12}  {info['s_star']:>4}  "
                f"{info['I_min']:>12}  {info['ratio']:>8.3f}  {verdict:>16}"
            )

    # --- Table 3: s* pattern ---
    print()
    print("=" * 50)
    print("TABLE 3: s* pattern (for k where slack > 0)")
    print("=" * 50)
    print(f"{'k':>4}  {'s*':>4}  {'k//2':>5}  {'match?':>7}")
    print("-" * 25)
    for k in range(9, k_max + 1):
        info = analyze_critical_threshold(k)
        if info['s_star'] is not None:
            match = "YES" if info['s_star'] == k // 2 else "no"
            print(f"{k:>4}  {info['s_star']:>4}  {k//2:>5}  {match:>7}")

    # --- Extended ratio check ---
    print()
    print("Extended ratio check (using s = k//2 to skip full search):")
    for k_ext in [50, 100, 200, 300]:
        if k_ext > k_max:
            t = k_ext - 1
            target = comb(2 * k_ext + 1, k_ext)
            raw = beta(k_ext, t) + beta(k_ext, t - 1)
            slack = raw - target
            if slack > 0:
                I_val = intersection_count(k_ext, k_ext // 2, t)
                print(f"  k={k_ext:>4}  ratio={I_val/slack:.5f}")