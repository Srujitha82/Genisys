"""Synthetic DNA data: reference libraries and reads with sequencing errors."""
import random
from typing import List


def random_sequence(length: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def simulate_read(seq: str, n_errors: int, rng: random.Random) -> str:
    """Return `seq` with exactly n_errors substitution errors (each to a different base)."""
    n_errors = min(n_errors, len(seq))
    read = list(seq)
    for pos in rng.sample(range(len(seq)), n_errors):
        read[pos] = rng.choice([b for b in "ACGT" if b != read[pos]])
    return "".join(read)


def generate_library(n: int, length: int = 32, family_size: int = 10, max_family_errors: int = 4,
                     seed: int = 0) -> List[str]:
    """n reference sequences grouped into 'families' (like related genes or barcodes).

    Each family has a random centre; members differ from it by 0..max_family_errors
    substitutions. Clustered data is realistic, and (unlike purely random DNA, where
    almost all distances are ~0.75*length) gives the VP-tree something to prune.
    """
    rng = random.Random(seed)
    centres = [random_sequence(length, rng) for _ in range(max(1, n // family_size))]
    out, seen = [], set()
    while len(out) < n:
        s = simulate_read(rng.choice(centres), rng.randint(0, max_family_errors), rng)
        if s not in seen:  # keep references unique
            seen.add(s)
            out.append(s)
    return out
