"""
Ball sizes in Kneser graphs K(2k+1, k).

Vertices are k-subsets of [2k+1]. Two vertices are adjacent iff disjoint.

Distance formula (from the 2005 diameter paper):
    d(A, B) = min(2s+1, 2(k-s))  where s = |A ∩ B|

Since 2s+1 is odd and 2(k-s) is even, each distance x corresponds to
exactly one overlap value s:
    x odd  -> s = (x-1)/2
    x even -> s = k - x/2

Number of k-subsets of [2k+1] with overlap exactly s with a fixed k-set v:
    C(k, s) * C(k+1, k-s)        (hypergeometric coefficient, with n-k = k+1)

The diameter of K(2k+1, k) is k.
"""

from math import comb


def vertices_at_distance(k: int, x: int) -> int:
    """
    Number of vertices at distance x from a fixed vertex in K(2k+1, k).

    Parameters
    ----------
    k : int  -- the graph parameter (vertices are k-subsets of [2k+1])
    x : int  -- distance (0 <= x <= k)

    Returns
    -------
    int -- number of vertices at that distance (0 if x is out of range)
    """
    if x == 0:
        return 1
    if x < 0 or x > k:
        return 0

    # Recover the unique overlap s from distance x
    if x % 2 == 1:          # odd distance: 2s+1 = x
        s = (x - 1) // 2
    else:                   # even distance: 2(k-s) = x
        s = k - x // 2

    # Count k-subsets of [2k+1] with |intersection with v| = s
    # = C(k, s) * C(k+1, k-s)
    return comb(k, s) * comb(k + 1, k - s)


def ball_size(k: int, r: int) -> int:
    """
    Number of vertices within distance r from a fixed vertex in K(2k+1, k).

    Parameters
    ----------
    k : int  -- the graph parameter
    r : int  -- ball radius

    Returns
    -------
    int -- |B(v, r)|
    """
    return sum(vertices_at_distance(k, x) for x in range(r + 1))


def print_distance_distribution(k: int) -> None:
    """Print the full distance distribution and cumulative ball sizes for K(2k+1, k)."""
    n = 2 * k + 1
    total_vertices = comb(n, k)
    diameter = k

    print(f"K({n}, {k})  |V| = C({n},{k}) = {total_vertices}  diameter = {diameter}")
    print(f"  {'dist':>4}  {'|layer|':>10}  {'|ball|':>10}  {'s (overlap)':>12}")
    print(f"  {'-'*4}  {'-'*10}  {'-'*10}  {'-'*12}")

    running = 0
    for x in range(diameter + 1):
        layer = vertices_at_distance(k, x)
        running += layer
        if x == 0:
            s_label = "—"
        elif x % 2 == 1:
            s_label = str((x - 1) // 2)
        else:
            s_label = str(k - x // 2)
        print(f"  {x:>4}  {layer:>10}  {running:>10}  {s_label:>12}")

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
        # Print full distance distributions for small cases
        for k in range(2, 8):
            print_distance_distribution(k)
