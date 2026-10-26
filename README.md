# VP-Tree DNA Read Matching

**[▶ Live interactive demo](https://srujitha82.github.io/Genisys/)** – type a read with errors and watch the VP-tree find its match.

Sequencing machines read DNA as millions of short reads, and every read contains errors (roughly 0.1 to 1% of bases are misread), so exact string matching fails to link a read to its true source. This project solves that with a **Vantage-Point Tree (VP-tree)**, a data structure for searching a metric space. Each DNA sequence is packed into a single integer (2 bits per base), and the distance between two sequences is their **Hamming distance**, the number of positions where they differ. Because Hamming distance obeys the triangle inequality, the tree can discard whole groups of reference sequences without ever comparing the read to them. Given a read, the system returns the k closest reference sequences or every sequence within a chosen number of mismatches, and for the best match or a small mismatch radius it needs many times fewer distance computations than a linear scan. A brute-force baseline with a correctness checker verifies that the tree returns exactly the same results.

## Contents

1. [Problem](#problem)
2. [Approach](#approach)
3. [Why a VP-tree?](#why-a-vp-tree)
4. [How the VP-tree works](#how-the-vp-tree-works)
5. [Complexity](#complexity)
6. [Project structure](#project-structure)
7. [Setup](#setup)
8. [Usage](#usage)
9. [Visual demo](#visual-demo)
10. [Testing](#testing)
11. [Results](#results)
12. [Edge cases handled](#edge-cases-handled)
13. [Team](#team)
14. [License](#license)

## Problem

Given a reference library of N DNA sequences of length L (barcodes, marker genes, reference k-mers) and a query read that may contain sequencing errors, return:

- the **k closest** reference sequences (kNN search), and
- all references within **Hamming distance r** (range search),

without comparing the read against every reference.

Real-world uses: barcode demultiplexing (assigning reads to samples despite barcode errors), species identification from marker genes, and CRISPR guide off-target search.

## Approach

1. **Encoding.** Each sequence is packed into one Python integer, 2 bits per base (`A=00, C=01, G=10, T=11`).
2. **Distance.** Hamming distance is the number of positions where two equal-length sequences differ. With packed integers it needs no per-base loop:

   ```
   x = a XOR b
   distance = popcount((x | x >> 1) & 0x5555...)
   ```

3. **Metric.** Hamming distance is non-negative, zero only for identical sequences, symmetric, and satisfies the triangle inequality, so it is a true metric and a VP-tree can use it.
4. **Index.** The VP-tree splits the library by distance from a vantage point, and the triangle inequality lets the search skip entire subtrees.

All sequences in a library must have the same length, because Hamming distance is only defined for equal-length strings. Reads with unknown bases (`N`) or the wrong length are rejected with a clear error.

## Why a VP-tree?

| Alternative | Why it is weaker here |
|---|---|
| Hash table / exact match | Any sequencing error causes a miss |
| Linear scan | O(N) distance computations for every read |
| Trie | An error early in the read sends the search down the wrong branch |
| k-d tree / R-tree | DNA has no coordinate axes to split on |
| Suffix array / FM-index | Used by real aligners (BWA, Bowtie) for exact seed matching inside one long genome, not for near-neighbour search over many sequences |
| BK-tree | The closest rival (it also uses integer distances), but its shape depends on the distance values and can become unbalanced |
| **VP-tree** | Balanced by the median split, works with any metric (including weighted ones), and returns exact results |

## How the VP-tree works

**Build.** Pick a vantage point, compute its distance to every other sequence, and split at the **median distance μ**. Sequences with distance ≤ μ go to the *inside* subtree, the rest to the *outside* subtree. Repeat on each side. The median split keeps the tree balanced, with depth about log₂ N.

**Search (kNN).** Keep τ, the distance of the k-th best match found so far (infinity until k matches exist). At a node with vantage point v:

1. Compute d = distance(read, v), and update the best matches and τ.
2. Search the **inside** subtree if `d − τ ≤ μ`.
3. Search the **outside** subtree if `d + τ ≥ μ`.
4. Visit the side that contains the read first, so τ shrinks quickly.

Skipping a subtree is safe because of the triangle inequality: every sequence x in the outside subtree has distance(v, x) > μ, so distance(read, x) ≥ distance(v, x) − d > μ − d. If μ − d ≥ τ, nothing there can beat the current best match. The inside case is symmetric.

**Range search** works the same way with τ fixed to the radius r.

**Worked example (length-8 sequences).** With vantage point `ACGTACGT`, the references are at distances `ACGTACGA`=1, `ACGAACGA`=2, `ACCAACGA`=3, `TTTTGGGG`=6, `GGGGTTTT`=6, so μ = 3. For the read `ACGTACGC` (distance 1 from the vantage point), τ = 1 after the first comparison. After searching the inside subtree, d + τ = 2 < μ = 3, so both outside sequences are skipped. The tree computes 4 distances instead of 6.

## Complexity

| Operation | Time |
|---|---|
| Distance (L bases, 64-bit words) | O(L / 32), effectively constant for short reads |
| Build | O(N log N) distance computations |
| Query (average, small radius) | O(log N) distance computations |
| Query (worst case) | O(N) |
| Space | O(N) |
| Brute-force kNN / range (baseline) | O(N log k) / O(N) |

Random DNA sequences differ in about 75% of positions, so pruning works well only for small radii. That matches the real use case, because most reads contain few errors. When the radius is large relative to μ, both subtrees must be searched and the cost approaches O(N).

## Project structure

```
├── src/
│   ├── sequence.py      # DNA <-> int encoding, Hamming distance, CSV/FASTA loaders
│   ├── simulate.py      # synthetic reference libraries and reads with errors
│   ├── bruteforce.py    # brute-force baseline + correctness checker
│   └── vptree.py        # VP-tree build, kNN and range search
├── tests/
│   ├── test_sequence.py
│   ├── test_bruteforce.py
│   └── test_vptree.py
├── docs/
│   └── index.html       # interactive visual demo (GitHub Pages)
├── data/
│   └── sample_reference.csv   # 100 sequences, length 32
├── make_sample_data.py
├── brute_demo.py        # CLI: brute-force nearest sequences
├── vptree_demo.py       # CLI: VP-tree search + distance counter
├── benchmark.py         # VP-tree vs brute force
├── requirements.txt
└── README.md
```

## Setup

Requires Python 3.8 or newer. The only dependency is pytest.

```bash
git clone https://github.com/Srujitha82/Genisys.git
cd Genisys
pip install -r requirements.txt
```

## Usage

### Command line (brute force)

```bash
python brute_demo.py TCCGACAGTTATGCTTTCGCGTTTGCCTAAAC --k 5
```

This read is `ref_001` with 2 sequencing errors:

```
distance 2   ref_001    TCCGACAGTGATGCTTTCGCCTTTGCCTAAAC
distance 5   ref_066    TGCGACAGTGACGCTTTCGCCGTTGCCTAAAC
distance 7   ref_045    TGCGACAGTGAGGCTTTGGCCGTTGCCTAATC
distance 7   ref_092    TGCGGCAGTCACGCTTTCGCCGTTGCCTAATC
distance 17  ref_034    AAGCACAATAAACCTCTCTGACTGGCCGAATA
```

### Command line (VP-tree)

```bash
python vptree_demo.py TCCGACAGTTATGCTTTCGCGTTTGCCTAAAC --k 3
python vptree_demo.py TCCGACAGTTATGCTTTCGCGTTTGCCTAAAC --radius 4
```

The output ends with the number of distances the tree computed compared with brute force.

### Benchmark

```bash
python benchmark.py --out results.md
```

### From Python

```python
from src.sequence import encode_dna, load_library
from src.bruteforce import knn_bruteforce, range_bruteforce
from src.vptree import VPTree          # adjust to the actual class name

seqs, records = load_library("data/sample_reference.csv")
read = encode_dna("TCCGACAGTTATGCTTTCGCGTTTGCCTAAAC")

tree = VPTree(seqs)
tree.knn(read, k=5)             # [(distance, index), ...] nearest first
tree.range_search(read, r=3)    # everything within 3 mismatches

knn_bruteforce(seqs, read, 5)   # same answer, computed by scanning everything
```

### Generate your own data

```python
import random
from src.simulate import generate_library, simulate_read

refs = generate_library(100000, length=32, seed=1)   # clustered synthetic library
read = simulate_read(refs[0], n_errors=2, rng=random.Random(1))
```

Real libraries can be loaded with `load_library` (CSV with columns `id,name,sequence`) or `load_fasta`.

## Visual demo

The page in `docs/index.html` is served by GitHub Pages and needs no installation. Open the live link at the top of this README.

- Blue nodes were compared with the read, amber dashed nodes were **pruned** without ever being compared, and green nodes are the current best matches.
- The base-pair ladder shows the read against the compared sequence, with broken red bonds at mismatches.
- Counters show the distance computations of the tree against brute force, and the current radius τ.

The page uses a small 15-sequence library so the tree stays readable. The Python code is the real implementation.

## Testing

```bash
python -m pytest
```

The tests cover DNA encoding (lowercase input, invalid bases like `N`, wrong lengths), Hamming distance against a naive implementation, the metric properties (symmetry, bounds, triangle inequality), edge cases, the loaders, and the VP-tree itself.

**Correctness checker.** `check_correctness` runs many queries and compares the tree with brute force. For kNN it compares the sorted distances (ties may legitimately return different sequences at equal distance), and for range queries it compares the exact set of results.

```python
from src.bruteforce import check_correctness, sample_queries

check_correctness(tree, seqs, sample_queries(seqs, 50))
```

The tree must provide `tree.knn(query, k)` and `tree.range_search(query, r)`, each returning `(distance, index)` pairs, where `index` is the position in the list the tree was built from.

## Results

Benchmarks (`python benchmark.py`) compare the VP-tree with brute force on synthetic clustered DNA libraries (length 32). Reads have 0 to 3 simulated sequencing errors, and the figures are averaged over 100 queries. Every answer checked was identical to brute force.

| Library size N | Search | Brute-force distances | VP-tree distances (avg) | Reduction |
|---|---|---|---|---|
| 1,000 | best match (k=1) | 1,000 | 209 | 4.8x |
| 1,000 | range r=2 | 1,000 | 228 | 4.4x |
| 1,000 | k=5 | 1,000 | 687 | 1.5x |
| 10,000 | best match (k=1) | 10,000 | 1,750 | 5.7x |
| 10,000 | range r=2 | 10,000 | 1,239 | 8.1x |
| 10,000 | k=5 | 10,000 | 7,192 | 1.4x |
| 100,000 | best match (k=1) | 100,000 | 10,244 | 9.8x |
| 100,000 | range r=2 | 100,000 | 6,604 | 15.1x |
| 100,000 | k=5 | 100,000 | 64,383 | 1.6x |

The saving grows with library size for the best match and for small-radius searches, which is the main use case. Asking for k=5 prunes much less on this data: the fifth-nearest sequence is far from the read, so the search radius τ stays large and few subtrees can be skipped. In pure Python that case can even be slower than brute force, because the tree has extra bookkeeping. This matches the theory: pruning works only when τ is small compared with the median distance μ.

## Edge cases handled

- Exact match (distance 0), and reads with 1, 2 or more errors
- A read that matches nothing within radius r
- `k` larger than the library, empty library, single sequence
- Duplicate references and ties at equal distance
- Unknown bases such as `N`, lowercase input, and wrong-length reads (rejected with `ValueError`)

## Team

| Member | GitHub | Contribution |
|---|---|---|
| Member 1 | [@Srujitha82](https://github.com/Srujitha82) | DNA encoding and Hamming distance, CSV and FASTA loaders, read simulator, brute-force baseline, correctness checker, unit tests, interactive visual demo |
| Member 2 | [@yasaswinimattaparthi26](https://github.com/yasaswinimattaparthi26) | VP-tree construction, vantage point selection, leaf buckets |
| Member 3 | [@jaladhi1926](https://github.com/jaladhi1926) | kNN and range search, benchmarks, CLI demo |

Each member committed from their own GitHub account; see the commit history and **Insights → Contributors**.

## License

MIT
