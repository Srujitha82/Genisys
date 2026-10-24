import itertools
import os
import random

import pytest

from src.sequence import (decode_dna, encode_dna, hamming_distance, load_fasta, load_library,
                          popcount)
from src.simulate import generate_library, simulate_read

SAMPLE = os.path.join(os.path.dirname(__file__), "..", "data", "sample_reference.csv")


def naive_hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


# ---------- encoding ----------
def test_encode_known_value():
    assert encode_dna("ACGT") == 0b00011011  # A=00 C=01 G=10 T=11
    assert encode_dna("TTTT") == 0b11111111
    assert encode_dna("AAAA") == 0


def test_encode_decode_roundtrip():
    rng = random.Random(1)
    for _ in range(100):
        s = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 60)))
        assert decode_dna(encode_dna(s), len(s)) == s


def test_lowercase_and_whitespace_accepted():
    assert encode_dna("  acgt ") == encode_dna("ACGT")


@pytest.mark.parametrize("bad", ["", "   ", "ACGN", "ACGU", "AC GT", "ACG1", None, 123])
def test_invalid_input_raises(bad):
    with pytest.raises(ValueError):
        encode_dna(bad)


def test_expected_length_enforced():
    assert encode_dna("ACGT", expected_len=4)
    with pytest.raises(ValueError):
        encode_dna("ACGT", expected_len=5)


def test_too_long_raises():
    with pytest.raises(ValueError):
        encode_dna("A" * 300)


def test_popcount():
    assert popcount(0) == 0
    assert popcount(0b1011) == 3


# ---------- Hamming distance ----------
def test_hand_computed_examples():
    a = encode_dna("ACGTACGT")
    assert hamming_distance(a, encode_dna("ACGTACGA")) == 1
    assert hamming_distance(a, encode_dna("ACGAACGA")) == 2
    assert hamming_distance(a, encode_dna("ACCAACGA")) == 3
    assert hamming_distance(a, encode_dna("TTTTGGGG")) == 6
    assert hamming_distance(a, encode_dna("GGGGTTTT")) == 6


def test_identical_is_zero():
    s = encode_dna("ACGTACGTAC")
    assert hamming_distance(s, s) == 0


def test_all_positions_differ_gives_length():
    assert hamming_distance(encode_dna("AAAAAAAA"), encode_dna("CCCCCCCC")) == 8
    assert hamming_distance(encode_dna("ACGTACGT"), encode_dna("TGCATGCA")) == 8


def test_matches_naive_implementation_on_random_strings():
    rng = random.Random(2)
    for _ in range(500):
        n = rng.randint(1, 100)
        a = "".join(rng.choice("ACGT") for _ in range(n))
        b = "".join(rng.choice("ACGT") for _ in range(n))
        assert hamming_distance(encode_dna(a), encode_dna(b)) == naive_hamming(a, b)


def test_symmetry_and_bounds():
    rng = random.Random(3)
    for _ in range(200):
        a = encode_dna("".join(rng.choice("ACGT") for _ in range(32)))
        b = encode_dna("".join(rng.choice("ACGT") for _ in range(32)))
        d = hamming_distance(a, b)
        assert d == hamming_distance(b, a)
        assert 0 <= d <= 32


def test_triangle_inequality_random():
    rng = random.Random(4)
    for _ in range(2000):
        a, b, c = (encode_dna("".join(rng.choice("ACGT") for _ in range(16))) for _ in range(3))
        assert hamming_distance(a, c) <= hamming_distance(a, b) + hamming_distance(b, c)


def test_triangle_inequality_library():
    seqs, _ = load_library(SAMPLE)
    for a, b, c in itertools.combinations(seqs[:15], 3):
        assert hamming_distance(a, c) <= hamming_distance(a, b) + hamming_distance(b, c)


# ---------- loaders ----------
def test_load_sample_library():
    seqs, records = load_library(SAMPLE)
    assert len(seqs) == len(records) == 100
    assert records[0]["name"] == "ref_001"
    assert decode_dna(seqs[0], 32) == records[0]["sequence"]


def test_load_library_skips_bad_rows(tmp_path, capsys):
    p = tmp_path / "lib.csv"
    p.write_text("id,name,sequence\n1,a,ACGTACGT\n2,badbase,ACGTNCGT\n3,short,ACGT\n4,b,TTTTGGGG\n")
    seqs, records = load_library(str(p))
    assert [r["name"] for r in records] == ["a", "b"]
    assert "skipped 2" in capsys.readouterr().out


def test_load_fasta(tmp_path):
    p = tmp_path / "ref.fasta"
    p.write_text(">seq1 first\nACGTAC\nGTAC\n>seq2\nTTTTGGGGCC\n>seq3\nACGT\n")
    seqs, records = load_fasta(str(p))
    assert len(seqs) == 2  # seq3 has the wrong length and is skipped
    assert records[0]["sequence"] == "ACGTACGTAC"  # multi-line record joined
    assert records[1]["name"] == "seq2"


# ---------- simulation ----------
def test_simulate_read_has_exact_number_of_errors():
    rng = random.Random(5)
    ref = "ACGTACGTACGTACGTACGT"
    for e in range(0, 8):
        read = simulate_read(ref, e, rng)
        assert len(read) == len(ref) and naive_hamming(read, ref) == e


def test_generate_library_properties():
    lib = generate_library(60, length=24, seed=1)
    assert len(lib) == len(set(lib)) == 60
    assert all(len(s) == 24 and set(s) <= set("ACGT") for s in lib)
    assert lib == generate_library(60, length=24, seed=1)  # reproducible
