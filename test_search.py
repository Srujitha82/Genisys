"""Member 3: tests for kNN and range search, checked against brute force."""
import os
import random

import pytest

from src.bruteforce import check_correctness, knn_bruteforce, range_bruteforce, sample_queries
from src.sequence import encode_dna, load_library
from src.simulate import generate_library, simulate_read
from src.vptree import VPTree

SAMPLE = os.path.join(os.path.dirname(__file__), "..", "data", "sample_reference.csv")


def lib(n, length=32, seed=0):
    return [encode_dna(s) for s in generate_library(n, length, seed=seed)]


def noisy_queries(seqs_str, n, errors, seed=0):
    rng = random.Random(seed)
    return [encode_dna(simulate_read(rng.choice(seqs_str), errors, rng)) for _ in range(n)]


@pytest.mark.parametrize("leaf", [1, 4, 8, 16])
def test_matches_bruteforce_on_sample_library(leaf):
    seqs, _ = load_library(SAMPLE)
    tree = VPTree(seqs, leaf_size=leaf)
    assert check_correctness(tree, seqs, sample_queries(seqs, 40), verbose=False)


def test_matches_bruteforce_on_larger_library_with_noisy_reads():
    strs = generate_library(2000, 32, seed=1)
    seqs = [encode_dna(s) for s in strs]
    tree = VPTree(seqs, leaf_size=8)
    for errors in (0, 1, 2, 3, 6):
        qs = noisy_queries(strs, 15, errors, seed=errors)
        assert check_correctness(tree, seqs, qs, ks=(1, 3, 10), radii=(0, 1, 2, 4, 8), verbose=False)


def test_queries_that_are_not_in_the_library():
    strs = generate_library(500, 24, seed=2)
    seqs = [encode_dna(s) for s in strs]
    tree = VPTree(seqs)
    qs = [encode_dna("".join(random.Random(i).choice("ACGT") for _ in range(24))) for i in range(20)]
    assert check_correctness(tree, seqs, qs, verbose=False)


def test_knn_results_sorted_and_exact_match_first():
    seqs = lib(300, seed=3)
    tree = VPTree(seqs)
    res = tree.knn(seqs[17], 5)
    assert len(res) == 5 and res[0] == (0, 17)
    assert [d for d, _ in res] == sorted(d for d, _ in res)


def test_knn_edge_cases():
    seqs = lib(20, 16, seed=4)
    tree = VPTree(seqs, leaf_size=2)
    assert tree.knn(seqs[0], 0) == []
    assert tree.knn(seqs[0], -3) == []
    assert len(tree.knn(seqs[0], 1000)) == 20  # k > N returns everything
    assert VPTree([]).knn(encode_dna("ACGT"), 3) == []


def test_range_edge_cases():
    seqs = lib(50, 16, seed=5)
    tree = VPTree(seqs, leaf_size=2)
    q = seqs[3]
    assert {i for _, i in tree.range_search(q, 16)} == set(range(50))  # r = L -> all
    assert tree.range_search(q, 0) == range_bruteforce(seqs, q, 0)
    assert tree.range_search(q, -1) == []
    assert VPTree([]).range_search(q, 3) == []


def test_single_sequence_library():
    s = encode_dna("ACGTACGT")
    tree = VPTree([s])
    assert tree.knn(encode_dna("ACGTACGA"), 3) == [(1, 0)]
    assert tree.range_search(encode_dna("ACGTACGA"), 0) == []
    assert tree.range_search(encode_dna("ACGTACGA"), 1) == [(1, 0)]


def test_duplicates_and_ties():
    a, b = encode_dna("ACGTACGT"), encode_dna("TTTTTTTT")
    seqs = [a, a, a, b, b]
    tree = VPTree(seqs, leaf_size=1)
    assert len(tree.range_search(a, 0)) == 3
    assert [d for d, _ in tree.knn(a, 4)] == [d for d, _ in knn_bruteforce(seqs, a, 4)]


def test_vptree_prunes_work_on_clustered_data():
    """For the typical use (best match, or a small radius) the tree must compute far
    fewer distances than brute force. (Asking for k=5 on this synthetic data prunes
    much less: the 5th neighbour is far away, so tau stays large. That is expected.)"""
    strs = generate_library(5000, 32, seed=6)
    seqs = [encode_dna(s) for s in strs]
    tree = VPTree(seqs, leaf_size=8)
    qs = noisy_queries(strs, 30, 1, seed=1)
    for run in (lambda q: tree.knn(q, 1), lambda q: tree.range_search(q, 2)):
        total = 0
        for q in qs:
            tree.reset_counter()
            run(q)
            total += tree.n_dist
        assert total / len(qs) < 0.3 * len(seqs)  # well below brute force's N


def test_counter_counts_only_query_distances():
    seqs = lib(200, seed=7)
    tree = VPTree(seqs)
    tree.reset_counter()
    tree.knn(seqs[0], 3)
    first = tree.n_dist
    assert 0 < first <= len(seqs)
    tree.knn(seqs[0], 3)
    assert tree.n_dist == 2 * first
