"""
Computes the sum of ball sizes from r=1 to r=t in K(2k+1, k).

That is, for a given k and t:

    S(k, t) = sum_{r=1}^{t} |B(r)|

where each |B(r)| is computed using the explicit formula from
odd_kneser_ball_explicit.py.
"""

from odd_kneser_ball_explicit import ball_size


def sum_of_balls(k: int, t: int) -> int | None:
    """
    Compute S(k, t) = sum_{r=1}^{t} |B(r)| in K(2k+1, k).

    Parameters
    ----------
    k : int  -- graph parameter
    t : int  -- upper limit of summation

    Returns
    -------
    int -- S(k, t), or None if t is invalid
    """
    if t < 1:
        print("Error: t must be at least 1.")
        return None
    t_eff = min(t, k)
    if t > k:
        print(f"Note: the diameter of K({2*k+1}, {k}) is {k}, "
              f"so balls of radius > {k} are the entire graph. Capping t at {k}.")
    return sum(ball_size(k, r) for r in range(1, t_eff + 1))


if __name__ == "__main__":
    k = int(input("Enter k (graph parameter for K(2k+1, k)): "))
    t = int(input("Enter t (upper limit of summation): "))

    t_eff = min(t, k)
    if t > k:
        print(f"Note: the diameter of K({2*k+1}, {k}) is {k}, capping t at {k}.")

    print(f"\nS(k={k}, t) for t = 1..{t_eff}:")
    for t_prime in range(1, t_eff + 1):
        print(sum_of_balls(k, t_prime))
