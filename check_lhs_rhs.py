import math
from math import comb


def lhs(k):
    """Compute C(k+1, ceil(k/2)) * C(k, floor(k/2))"""
    return math.comb(k + 1, math.ceil(k / 2)) * math.comb(k, math.floor(k / 2))


def rhs(k):
    """Compute sum_{i=0}^{k-2} (k-1-i) * C(k+1, ceil(i/2)) * C(k, floor(i/2))"""
    total = 0
    for i in range(k - 1):  # i = 0 to k-2 inclusive
        coeff = k - 1 - i
        term = coeff * math.comb(k + 1, math.ceil(i / 2)) * math.comb(k, math.floor(i / 2))
        total += term
    return total


def closed_form_intersection(k: int) -> int:
    """
    Compute:
        sum_{x=0}^{floor(k/2)}  sum_{s_a in S_a, s_b in S_b}
            C(floor(k/2), x)
            * C(k - floor(k/2), s_a)
            * C(k - floor(k/2), s_b)
            * C(floor(k/2) + 1, k - (s_a + s_b + x))

    where
        S_a = [0, ceil((k-1)/2) - x]  ∪  [k - floor((k-1)/2) - x, k - x]
        S_b = [0, ceil((k-2)/2) - x]  ∪  [k - floor((k-2)/2) - x, k - x]

    Identities used: ceil((k-1)/2) = floor(k/2),  ceil((k-2)/2) = floor((k-1)/2).
    """
    half_k = k // 2           # floor(k/2) = s = ceil((k-1)/2)
    ks = k - half_k           # ceil(k/2) — pool size for s_a and s_b
    floor_a = (k - 1) // 2   # floor((k-1)/2) = ceil((k-2)/2)
    floor_b = (k - 2) // 2   # floor((k-2)/2)

    total = 0
    for x in range(half_k + 1):
        # S_a = [0, half_k - x] ∪ [k - floor_a - x, k - x]
        S_a = (
            set(range(0, half_k - x + 1)) |
            set(range(max(0, k - floor_a - x), k - x + 1))
        )
        # S_b = [0, floor_a - x] ∪ [k - floor_b - x, k - x]  (ceil_b = floor_a)
        S_b = (
            set(range(0, floor_a - x + 1)) |
            set(range(max(0, k - floor_b - x), k - x + 1))
        )

        cx = comb(half_k, x)
        for sa in S_a:
            if sa > ks:
                continue
            csa = comb(ks, sa)
            for sb in S_b:
                if sb > ks:
                    continue
                bottom = k - (sa + sb + x)
                if bottom < 0 or bottom > half_k + 1:
                    continue
                total += cx * csa * comb(ks, sb) * comb(half_k + 1, bottom)

    return total


def check_range(k_min=1, k_max=30):
    col = 22
    header = (
        f"{'k':>4} | {'LHS':>{col}} | {'RHS':>{col}} | "
        f"{'C(k-1,2)':>10} | {'formula(k)':>{col}} | "
        f"{'C(k-1,2)*formula':>{col}} | {'RHS - coeff*formula':>{col}} | "
        f"{'LHS - (RHS-coeff*formula)':>{col}} | {'LHS >= RHS-coeff*formula':>24}"
    )
    print(header)
    print("-" * len(header))

    holds_for = []
    fails_for = []

    for k in range(k_min, k_max + 1):
        left = lhs(k)
        right = rhs(k)
        coeff = comb(k - 1, 2)
        formula = closed_form_intersection(k)
        adjusted_rhs = right - coeff * formula
        diff = left - adjusted_rhs
        holds = left >= adjusted_rhs
        status = "TRUE" if holds else "FALSE"

        print(
            f"{k:>4} | {left:>{col}} | {right:>{col}} | "
            f"{coeff:>10} | {formula:>{col}} | "
            f"{coeff * formula:>{col}} | {adjusted_rhs:>{col}} | "
            f"{diff:>{col}} | {status:>24}"
        )

        if holds:
            holds_for.append(k)
        else:
            fails_for.append(k)

    print()
    print(f"Inequality holds for k in: {holds_for}")
    print(f"Inequality fails for k in: {fails_for}")


if __name__ == "__main__":
    check_range(k_min=1, k_max=1000)
