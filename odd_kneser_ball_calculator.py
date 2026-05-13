"""
Computes |B(v, r)| in any odd Kneser graph K(2k+1, k).

Vertices are k-subsets of [2k+1]. Two vertices are adjacent iff disjoint.
Because K(2k+1, k) is vertex-transitive, all r-balls have the same size
regardless of the choice of center v.

Algorithm (three steps):

  Step 1 — Find the overlap s.
    Apply Lemma 4: for vertices v, x with |v ∩ x| = s,
        d(v, x) = min(2s+1, 2(k-s)).
    Since 2s+1 is odd and 2(k-s) is even, each distance r determines s
    uniquely:
        r odd  =>  r = 2s+1  =>  s = (r-1)/2
        r even =>  r = 2(k-s) =>  s = k - r/2

  Step 2 — Count k-sets of [2k+1] with overlap exactly s with v.
    - Choose s elements from the k elements of v to be the overlap: C(k, s)
    - Choose the remaining k-s elements of x from outside v,
      i.e. from the k+1 elements of [2k+1] \\ v:          C(k+1, k-s)
    Multiplying gives the boundary layer size:
        |∂B(v, r)| = C(k, s) * C(k+1, k-s)

  Step 3 — Iterate.
    Repeat Steps 1-2 for r' = r-1, r-2, ..., 0. Then, since
        B(v, r) = disjoint union of ∂B(v, l) for l = 0..r,
    sum the layers:
        |B(v, r)| = sum_{l=0}^{r} |∂B(v, l)|

The diameter of odd Kneser graphs K(2k+1, k) is always k, so r is valid for 0 <= r <= k.
"""

from math import comb


def boundary_layer_size(k: int, r: int) -> int:
    """
    Compute |∂B(v, r)|, the number of vertices at distance exactly r
    from a fixed vertex v in K(2k+1, k).

    Implements Steps 1 and 2 of the algorithm.

    Parameters
    ----------
    k : int  -- graph parameter (vertices are k-subsets of [2k+1])
    r : int  -- distance (0 <= r <= k)

    Returns
    -------
    int -- |∂B(v, r)|; returns 0 if r is outside [0, k]
    """
    if r == 0:
        return 1                    # ∂B(v, 0) = {v}
    if r < 0 or r > k:
        return 0

    # Step 1: recover overlap s from distance r (Lemma 4)
    if r % 2 == 1:                  # r = 2s+1  =>  s = (r-1)/2
        s = (r - 1) // 2
    else:                           # r = 2(k-s)  =>  s = k - r/2
        s = k - r // 2

    # Step 2: C(k, s) * C(k+1, k-s)
    #   C(k, s)     -- ways to choose the s overlapping elements from v
    #   C(k+1, k-s) -- ways to choose the k-s non-overlapping elements
    #                  from [2k+1] \ v  (which has k+1 elements)
    return comb(k, s) * comb(k + 1, k - s)


def ball_size(k: int, r: int) -> int | None:
    """
    Compute |B(v, r)|, the number of vertices within distance r from a
    fixed vertex v in K(2k+1, k).

    Implements Step 3 of the algorithm:
        |B(v, r)| = sum_{l=0}^{r} |∂B(v, l)|

    Parameters
    ----------
    k : int  -- graph parameter
    r : int  -- ball radius

    Returns
    -------
    int -- |B(v, r)|, or the full graph size if r > k, or None if r < 0
    """
    if r < 0:
        print(f"Error: negative radii are not defined.")
        return None
    if r > k:
        total = comb(2 * k + 1, k)
        print(f"Note: the diameter of K({2*k+1}, {k}) is {k}, "
              f"so B(v, {r}) is the entire graph ({total} vertices).")
        return total
    return sum(boundary_layer_size(k, l) for l in range(r + 1))


def print_distance_distribution(k: int) -> None:
    """Print the full distance distribution and cumulative ball sizes for K(2k+1, k)."""
    n = 2 * k + 1
    total_vertices = comb(n, k)
    diameter = k

    print(f"K({n}, {k})  |V| = C({n},{k}) = {total_vertices}  diameter = {diameter}")
    print(f"  {'r':>4}  {'|dB(v,r)|':>12}  {'|B(v,r)|':>12}  {'s':>4}")
    print(f"  {'-'*4}  {'-'*12}  {'-'*12}  {'-'*4}")

    running = 0
    for r in range(diameter + 1):
        layer = boundary_layer_size(k, r)
        running += layer
        if r == 0:
            s_label = "-"
        elif r % 2 == 1:
            s_label = str((r - 1) // 2)
        else:
            s_label = str(k - r // 2)
        print(f"  {r:>4}  {layer:>12}  {running:>12}  {s_label:>4}")

    assert running == total_vertices, (
        f"Layer counts sum to {running}, expected C({n},{k}) = {total_vertices}"
    )
    print()


if __name__ == "__main__":
    import sys

    if len(sys.argv) == 3:
        # Usage: python kneser_ball_sizes.py <k> <r>
        k_arg = int(sys.argv[1])
        r_arg = int(sys.argv[2])
        result = ball_size(k_arg, r_arg)
        print(f"|B(v, {r_arg})| in K({2*k_arg+1}, {k_arg}) = {result}")
    else:
        for k in range(2, 8):
            print_distance_distribution(k)
