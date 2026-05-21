"""
Computes a table of values for |B(r)| in any odd Kneser graph K(2k+1, k) up to maximal ball size k,
which is the diameter of the graph, using the explicit formula for each layer.

Vertices are k-subsets of [2k+1]. Two vertices are adjacent iff disjoint.
Because K(2k+1, k) is vertex-transitive, all r-balls have the same size
regardless of the choice of center v.

Explicit boundary formula (derived in kneser_ball_algorithm.tex):

    |∂B(l)| = C(k, floor(l/2)) * C(k+1, ceil(l/2))

This unifies the even and odd cases:
    l even:  floor(l/2) = ceil(l/2) = l/2
             => C(k, l/2) * C(k+1, l/2)
    l odd:   floor(l/2) = (l-1)/2,  ceil(l/2) = (l+1)/2
             => C(k, (l-1)/2) * C(k+1, (l+1)/2)

The full ball size is then:

    |B(r)| = sum_{l=0}^{r} C(k, floor(l/2)) * C(k+1, ceil(l/2))

The diameter of K(2k+1, k) is k, so r is valid for 0 <= r <= k.
"""

from math import comb


def boundary_layer_size(k: int, l: int) -> int:
    """
    Compute |∂B(l)|, the number of vertices at distance exactly l
    from a fixed vertex in K(2k+1, k).

    Uses the explicit formula:
        |∂B(l)| = C(k, floor(l/2)) * C(k+1, ceil(l/2))

    Parameters
    ----------
    k : int  -- graph parameter (vertices are k-subsets of [2k+1])
    l : int  -- distance (0 <= l <= k)

    Returns
    -------
    int -- |∂B(l)|; returns 0 if l is outside [0, k]
    """
    if l < 0 or l > k:
        return 0
    return comb(k, l // 2) * comb(k + 1, (l + 1) // 2)


def ball_size(k: int, r: int) -> int | None:
    """
    Compute |B(r)|, the number of vertices within distance r from a
    fixed vertex in K(2k+1, k).

    Uses the explicit formula:
        |B(r)| = sum_{l=0}^{r} C(k, floor(l/2)) * C(k+1, ceil(l/2))

    Parameters
    ----------
    k : int  -- graph parameter
    r : int  -- ball radius

    Returns
    -------
    int -- |B(r)|, or the full graph size if r > k, or None if r < 0
    """
    if r < 0:
        print("Error: negative radii are not defined.")
        return None
    if r > k:
        total = comb(2 * k + 1, k)
        print(f"Note: the diameter of K({2*k+1}, {k}) is {k}, "
              f"so B(r={r}) is the entire graph ({total} vertices).")
        return total
    return sum(boundary_layer_size(k, l) for l in range(r + 1))


def print_distance_distribution(k: int) -> None:
    """Print the full distance distribution for K(2k+1, k)."""
    n = 2 * k + 1
    total_vertices = comb(n, k)

    print(f"\nK({n}, {k})  |V| = C({n},{k}) = {total_vertices}  diameter = {k}")
    print()
    print("  Columns:")
    print("    r       : distance (radius) from the center vertex")
    print("    s       : overlap |v ∩ x| for vertices x at distance r from v")
    print("    σ(r)    : number of vertices at distance exactly r  [= C(k, floor(r/2)) * C(k+1, ceil(r/2))]")
    print("    β(r)    : cumulative number of vertices within distance r  [= sum of σ(0)..σ(r)]")
    print()
    print(f"  {'r':>4}  {'s':>4}  {'σ(r)':>12}  {'β(r)':>12}")
    print(f"  {'-'*4}  {'-'*4}  {'-'*12}  {'-'*12}")

    running = 0
    for r in range(k + 1):
        layer = boundary_layer_size(k, r)
        running += layer
        s = (r - 1) // 2 if r % 2 == 1 else k - r // 2
        print(f"  {r:>4}  {s:>4}  {layer:>12}  {running:>12}")

    print(f"\n  β({k}) = {ball_size(k, k)}")



if __name__ == "__main__":
    k = int(input("Enter k (graph parameter for K(2k+1, k)): "))
    print_distance_distribution(k)
