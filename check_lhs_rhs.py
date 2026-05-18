import math
from fractions import Fraction

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

def check_range(k_min=1, k_max=30):
    print(f"{'k':>4} | {'LHS':>20} | {'RHS':>20} | {'LHS - RHS':>20} | {'LHS >= RHS':>10}")
    print("-" * 85)

    holds_for = []
    fails_for = []

    for k in range(k_min, k_max + 1):
        left = lhs(k)
        right = rhs(k)
        diff = left - right
        holds = left >= right
        status = "TRUE" if holds else "FALSE"

        print(f"{k:>4} | {left:>20} | {right:>20} | {diff:>20} | {status:>10}")

        if holds:
            holds_for.append(k)
        else:
            fails_for.append(k)

    print()
    print(f"Inequality holds for k in: {holds_for}")
    print(f"Inequality fails for k in: {fails_for}")

if __name__ == "__main__":
    check_range(k_min=1, k_max=30)
