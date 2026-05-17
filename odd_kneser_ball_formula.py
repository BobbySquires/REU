"""
Computes |B(r)| in any odd Kneser graph K(2k+1, k) using the explicit formula.

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


def print_distance_distribution(k: int, r: int) -> None:
    """Print the distance distribution up to radius r for K(2k+1, k)."""
    n = 2 * k + 1
    total_vertices = comb(n, k)

    print(f"\nK({n}, {k})  |V| = C({n},{k}) = {total_vertices}  diameter = {k}")
    print(f"  {'l':>4}  {'floor(l/2)':>10}  {'ceil(l/2)':>9}  {'|dB(l)|':>10}  {'|B(l)|':>10}")
    print(f"  {'-'*4}  {'-'*10}  {'-'*9}  {'-'*10}  {'-'*10}")

    running = 0
    for l in range(min(r, k) + 1):
        layer = boundary_layer_size(k, l)
        running += layer
        print(f"  {l:>4}  {l//2:>10}  {(l+1)//2:>9}  {layer:>10}  {running:>10}")

    print(f"\n  |B(r={r})| = {ball_size(k, r)}")


def print_ball_sequence(k: int, r: int) -> None:
    """Print the sequence of ball sizes |B(0)|, |B(1)|, ..., |B(r)| as a list."""
    r_max = min(r, k)
    sequence = [ball_size(k, l) for l in range(r_max + 1)]
    print(f"\nBall size sequence for K({2*k+1}, {k}), r = 0..{r_max}:")
    print(sequence)


def flat_ball_sequence(k_max: int = 10) -> list[int]:
    """
    For each k from 1 to k_max, compute |B(0)|, |B(1)|, ..., |B(k)| in
    K(2k+1, k), then concatenate all sequences into a single flat list.

    Returns
    -------
    list[int] -- [|B(0)|, |B(1)| (k=1), |B(0)|, ..., |B(2)| (k=2), ...]
    """
    return [sum(boundary_layer_size(k, l) for l in range(r + 1))
            for k in range(1, k_max + 1) for r in range(k + 1)]


if __name__ == "__main__":
    k = int(input("Enter k (graph parameter for K(2k+1, k)): "))
    r = int(input("Enter r (ball radius): "))
    print("Output options:")
    print("  1 - Full table")
    print("  2 - Ball size sequence only")
    choice = input("Choose (1 or 2): ").strip()

    if choice == "2":
        print_ball_sequence(k, r)
    else:
        print_distance_distribution(k, r)
