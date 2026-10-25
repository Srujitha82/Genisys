"""Benchmark: VP-tree vs brute force (Member 3).

Usage:
    python benchmark.py                       # N = 1000, 10000, 100000
    python benchmark.py --sizes 1000 5000     # custom sizes
    python benchmark.py --out results.md      # also save the markdown table

For each library size it builds a tree, simulates reads with 0..3 sequencing errors,
and reports the average number of distance computations per query (brute force always
needs N) plus the average query time. It also verifies the tree's answers against
brute force, so a speed-up is never reported for a wrong answer.
"""
import argparse
import random
import time

from src.bruteforce import knn_bruteforce, range_bruteforce
from src.sequence import encode_dna
from src.simulate import generate_library, simulate_read
from src.vptree import VPTree

parser = argparse.ArgumentParser()
parser.add_argument("--sizes", type=int, nargs="+", default=[1000, 10000, 100000])
parser.add_argument("--length", type=int, default=32)
parser.add_argument("--queries", type=int, default=100)
parser.add_argument("--verify", type=int, default=10, help="queries also run by brute force")
parser.add_argument("--leaf", type=int, default=8)
parser.add_argument("--out", default=None)
args = parser.parse_args()

SEARCHES = [("kNN k=1", lambda t, q: t.knn(q, 1), lambda s, q: knn_bruteforce(s, q, 1)),
            ("kNN k=5", lambda t, q: t.knn(q, 5), lambda s, q: knn_bruteforce(s, q, 5)),
            ("range r=2", lambda t, q: t.range_search(q, 2), lambda s, q: range_bruteforce(s, q, 2))]

rows = ["| N | search | brute-force distances | VP-tree distances (avg) | reduction | VP-tree ms/query | brute-force ms/query |",
        "|---|---|---|---|---|---|---|"]
for n in args.sizes:
    strs = generate_library(n, args.length, seed=1)
    seqs = [encode_dna(s) for s in strs]
    t0 = time.perf_counter()
    tree = VPTree(seqs, leaf_size=args.leaf)
    build = time.perf_counter() - t0
    rng = random.Random(2)
    queries = [encode_dna(simulate_read(rng.choice(strs), rng.randint(0, 3), rng)) for _ in range(args.queries)]
    print(f"N={n}: built in {build:.2f}s, height {tree.stats()['height']}")
    for name, tree_fn, brute_fn in SEARCHES:
        total, t0 = 0, time.perf_counter()
        for q in queries:
            tree.reset_counter()
            tree_fn(tree, q)
            total += tree.n_dist
        tree_ms = (time.perf_counter() - t0) / len(queries) * 1000
        checked = queries[:args.verify]
        t0 = time.perf_counter()
        expected = [brute_fn(seqs, q) for q in checked]
        brute_ms = (time.perf_counter() - t0) / max(1, len(checked)) * 1000
        for q, exp in zip(checked, expected):  # correctness check (not timed)
            got = tree_fn(tree, q)
            assert [d for d, _ in exp] == [d for d, _ in got], f"MISMATCH for {name}, N={n}"
        avg = total / len(queries)
        rows.append(f"| {n:,} | {name} | {n:,} | {avg:,.0f} | {n / avg:.1f}x | {tree_ms:.2f} | {brute_ms:.2f} |")

table = "\n".join(rows)
print("\n" + table)
print(f"\nAll verified answers matched brute force ({args.verify} queries per search).")
if args.out:
    with open(args.out, "w", encoding="utf-8") as f:
        f.write(table + "\n")
