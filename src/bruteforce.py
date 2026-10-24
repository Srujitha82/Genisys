"""Brute-force baseline and correctness checker for the VP-tree.

Interface the VP-tree (Members 2 and 3) must provide so the checker can test it:

    tree.knn(query, k)            -> list of (distance, index), nearest first
    tree.range_search(query, r)   -> list of (distance, index), any order

`query` is a packed sequence (int from encode_dna). `index` is the position of the
sequence in the `sequences` list the tree was built from. Distances are ints.
"""
import heapq
import random
from typing import List, Sequence, Tuple

from .sequence import hamming_distance

Result = Tuple[int, int]


def knn_bruteforce(sequences: Sequence[int], query: int, k: int) -> List[Result]:
    """k nearest sequences by scanning everything. O(N log k)."""
    if k <= 0:
        return []
    scored = ((hamming_distance(query, s), i) for i, s in enumerate(sequences))
    return heapq.nsmallest(k, scored)  # sorted by (distance, index)


def range_bruteforce(sequences: Sequence[int], query: int, r: int) -> List[Result]:
    """All sequences with distance <= r, sorted by (distance, index). O(N)."""
    out = [(d, i) for i, s in enumerate(sequences) if (d := hamming_distance(query, s)) <= r]
    out.sort()
    return out


def count_distance_computations(sequences: Sequence[int]) -> int:
    """Brute force always computes N distances per query (baseline for benchmarks)."""
    return len(sequences)


def check_correctness(tree, sequences: Sequence[int], queries: Sequence[int],
                      ks: Sequence[int] = (1, 5, 10), radii: Sequence[int] = (0, 1, 2, 3),
                      verbose: bool = True) -> bool:
    """Compare the tree against brute force for every query.

    kNN: compare the sorted DISTANCE lists (ties may legitimately return different
    sequences at equal distance, so ids are not compared).
    Range: compare the exact SET of indices.
    Returns True only if everything matches.
    """
    ok = True
    for qi, q in enumerate(queries):
        for k in ks:
            expected = [d for d, _ in knn_bruteforce(sequences, q, k)]
            got = [d for d, _ in tree.knn(q, k)]
            if expected != got:
                ok = False
                if verbose:
                    print(f"[MISMATCH] kNN query#{qi} k={k}\n  expected={expected}\n  got     ={got}")
        for r in radii:
            expected_ids = {i for _, i in range_bruteforce(sequences, q, r)}
            got_ids = {i for _, i in tree.range_search(q, r)}
            if expected_ids != got_ids:
                ok = False
                if verbose:
                    print(f"[MISMATCH] range query#{qi} r={r}\n"
                          f"  missing={sorted(expected_ids - got_ids)} extra={sorted(got_ids - expected_ids)}")
    if verbose:
        print("All results match brute force." if ok else "Correctness check FAILED.")
    return ok


def sample_queries(sequences: Sequence[int], n: int = 20, seed: int = 0) -> List[int]:
    """Pick n random library sequences as queries (reproducible with seed)."""
    rng = random.Random(seed)
    return [rng.choice(sequences) for _ in range(n)]
