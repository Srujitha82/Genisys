"""DNA encoding and Hamming distance.

A DNA sequence is packed into ONE Python int, 2 bits per base:
    A=00  C=01  G=10  T=11
(first base in the most significant position). Packing lets us compare two
sequences with a single XOR plus a popcount instead of looping over bases.
The VP-tree treats these ints as opaque values and only calls `hamming_distance`.

All sequences in one library must have the SAME length (Hamming distance is only
defined for equal-length strings). Length is checked when sequences are loaded.
"""
import csv
from typing import List, Optional, Tuple

BASE_TO_BITS = {"A": 0, "C": 1, "G": 2, "T": 3}
BITS_TO_BASE = "ACGT"

MAX_LEN = 256  # longest supported sequence (bases)
# 0b0101...01 covering 2*MAX_LEN bits: one flag bit per base position
_FLAG_MASK = int("01" * MAX_LEN, 2)


def popcount(x: int) -> int:
    """Number of set bits in x."""
    try:
        return x.bit_count()  # Python 3.10+
    except AttributeError:
        return bin(x).count("1")


def encode_dna(seq: str, expected_len: Optional[int] = None) -> int:
    """Pack a DNA string (A/C/G/T, any case) into an int.

    Raises ValueError for empty input, unknown bases (including N), sequences
    longer than MAX_LEN, or a length different from `expected_len`.
    """
    if not isinstance(seq, str):
        raise ValueError("sequence must be a string")
    seq = seq.strip().upper()
    if not seq:
        raise ValueError("sequence is empty")
    if len(seq) > MAX_LEN:
        raise ValueError(f"sequence longer than {MAX_LEN} bases")
    if expected_len is not None and len(seq) != expected_len:
        raise ValueError(f"expected length {expected_len}, got {len(seq)}")
    value = 0
    for base in seq:
        if base not in BASE_TO_BITS:
            raise ValueError(f"invalid base {base!r} in {seq!r}")
        value = (value << 2) | BASE_TO_BITS[base]
    return value


def decode_dna(value: int, length: int) -> str:
    """Inverse of encode_dna (the length must be supplied, leading A's are 00)."""
    return "".join(BITS_TO_BASE[(value >> (2 * (length - 1 - i))) & 3] for i in range(length))


def hamming_distance(a: int, b: int) -> int:
    """Number of positions where two packed sequences differ.

    x = a XOR b is non-zero in the 2 bits of every mismatching base. OR-ing x with
    itself shifted right by 1 puts a 1 in the low bit of each mismatching base, the
    mask keeps only those low bits, and popcount counts them. O(L/word-size).

    Hamming distance is a true metric on equal-length strings (non-negative, zero
    iff identical, symmetric, triangle inequality), so it is valid for a VP-tree.
    """
    x = a ^ b
    return popcount((x | (x >> 1)) & _FLAG_MASK)


def load_library(path: str, expected_len: Optional[int] = None) -> Tuple[List[int], List[dict]]:
    """Load a CSV with columns id,name,sequence.

    The library length is `expected_len`, or the length of the first valid row.
    Rows with invalid bases or a different length are skipped and counted.
    Returns (sequences, records); sequences[i] belongs to records[i].
    """
    seqs: List[int] = []
    records: List[dict] = []
    skipped = 0
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            try:
                if expected_len is None:
                    expected_len = len(row["sequence"].strip())
                seqs.append(encode_dna(row["sequence"], expected_len))
                records.append(row)
            except ValueError:
                skipped += 1
    if skipped:
        print(f"[load_library] skipped {skipped} row(s) (invalid base or wrong length)")
    return seqs, records


def load_fasta(path: str, expected_len: Optional[int] = None) -> Tuple[List[int], List[dict]]:
    """Load a FASTA file (">name" header lines followed by sequence lines).

    Same rules as load_library. Records have keys id, name, sequence.
    """
    seqs: List[int] = []
    records: List[dict] = []
    skipped = 0
    entries, name, parts = [], None, []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    entries.append((name, "".join(parts)))
                name, parts = line[1:].strip(), []
            elif line:
                parts.append(line)
    if name is not None:
        entries.append((name, "".join(parts)))
    for i, (name, seq) in enumerate(entries, start=1):
        try:
            if expected_len is None:
                expected_len = len(seq)
            seqs.append(encode_dna(seq, expected_len))
            records.append({"id": str(i), "name": name, "sequence": seq.upper()})
        except ValueError:
            skipped += 1
    if skipped:
        print(f"[load_fasta] skipped {skipped} record(s) (invalid base or wrong length)")
    return seqs, records
