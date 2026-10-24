import os

from src.bruteforce import (check_correctness, knn_bruteforce, range_bruteforce,
                            sample_queries)
from src.sequence import encode_dna, hamming_distance, load_library
from src.simulate import generate_library

SAMPLE = os.path.join(os.path.dirname(__file__), "..", "data", "sample_reference.csv")

# the worked example from the README
EXAMPLE = [encode_dna(s) for s in
           ["ACGTACGT", "ACGTACGA", "ACGAACGA", "ACCAACGA", "TTTTGGGG", "GGGGTTTT"]]


class BruteForceTree:
    """Stand-in 'tree' that is correct by construction. Replace with the real VP-tree."""
    def __init__(self, seqs):
        self.seqs = seqs

    def knn(self, q, k):
        return knn_bruteforce(self.seqs, q, k)

    def range_search(self, q, r):
        return range_bruteforce(self.seqs, q, r)


class BrokenKnnTree(BruteForceTree):
    def knn(self, q, k):
        return knn_bruteforce(self.seqs, q, k)[1:]  # drops the nearest neighbour


class BrokenRangeTree(BruteForceTree):
    def range_search(self, q, r):
        return range_bruteforce(self.seqs, q, r)[:-1]  # misses one result


def test_worked_example_knn():
    q = encode_dna("ACGTACGC")  # read with one sequencing error
    res = knn_bruteforce(EXAMPLE, q, 1)
    assert res == [(1, 0)]  # ties at distance 1 (index 0 and 1); lowest index first
    assert [d for d, _ in knn_bruteforce(EXAMPLE, q, 4)] == [1, 1, 2, 3]


def test_worked_example_range():
    q = encode_dna("ACGTACGC")
    assert {i for _, i in range_bruteforce(EXAMPLE, q, 1)} == {0, 1}
    assert {i for _, i in range_bruteforce(EXAMPLE, q, 3)} == {0, 1, 2, 3}


def test_knn_sorted_and_self_first():
    seqs = [encode_dna(s) for s in generate_library(100, 32, seed=3)]
    res = knn_bruteforce(seqs, seqs[10], 5)
    assert len(res) == 5 and res[0] == (0, 10)
    assert [d for d, _ in res] == sorted(d for d, _ in res)


def test_knn_edge_cases():
    seqs = [encode_dna(s) for s in generate_library(10, 16, seed=4)]
    assert knn_bruteforce(seqs, seqs[0], 0) == []
    assert len(knn_bruteforce(seqs, seqs[0], 50)) == 10  # k > N returns all
    assert knn_bruteforce([], encode_dna("ACGT"), 3) == []  # empty library


def test_range_edge_cases():
    seqs = [encode_dna(s) for s in generate_library(50, 16, seed=5)]
    q = seqs[0]
    assert {i for _, i in range_bruteforce(seqs, q, 16)} == set(range(50))  # r = L -> everything
    zero = range_bruteforce(seqs, q, 0)
    assert zero == [(0, 0)]
    assert range_bruteforce([], q, 2) == []
    assert range_bruteforce(seqs, q, -1) == []


def test_range_results_all_within_radius():
    seqs, _ = load_library(SAMPLE)
    q = seqs[0]
    for d, i in range_bruteforce(seqs, q, 5):
        assert d <= 5 and d == hamming_distance(q, seqs[i])


def test_duplicates_are_all_returned():
    s = encode_dna("ACGTACGT")
    seqs = [s, s, s, encode_dna("TTTTTTTT")]
    assert len(range_bruteforce(seqs, s, 0)) == 3


def test_checker_passes_on_correct_tree():
    seqs, _ = load_library(SAMPLE)
    assert check_correctness(BruteForceTree(seqs), seqs, sample_queries(seqs, 10), verbose=False)


def test_checker_detects_broken_knn():
    seqs, _ = load_library(SAMPLE)
    assert not check_correctness(BrokenKnnTree(seqs), seqs, sample_queries(seqs, 5), verbose=False)


def test_checker_detects_broken_range():
    seqs, _ = load_library(SAMPLE)
    assert not check_correctness(BrokenRangeTree(seqs), seqs, sample_queries(seqs, 5), verbose=False)


def test_sample_queries_reproducible():
    seqs = [encode_dna(s) for s in generate_library(50, 16, seed=6)]
    assert sample_queries(seqs, 10, seed=1) == sample_queries(seqs, 10, seed=1)
