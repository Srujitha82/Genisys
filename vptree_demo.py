"""Demo CLI (Member 3): find the closest reference sequences for a read with the VP-tree.

Usage:
    python vptree_demo.py TCCGACAGTTATGCTTTCGCGTTTGCCTAAAC --k 3
    python vptree_demo.py TCCGACAGTTATGCTTTCGCGTTTGCCTAAAC --radius 4
"""
import argparse

from src.sequence import encode_dna, load_library
from src.vptree import VPTree

p = argparse.ArgumentParser()
p.add_argument("read", help="DNA read (A/C/G/T)")
p.add_argument("--k", type=int, default=3, help="number of nearest sequences")
p.add_argument("--radius", type=int, default=None, help="range search: all within this many mismatches")
p.add_argument("--data", default="data/sample_reference.csv")
p.add_argument("--leaf", type=int, default=8)
args = p.parse_args()

seqs, records = load_library(args.data)
length = len(records[0]["sequence"])
try:
    q = encode_dna(args.read, expected_len=length)
except ValueError as e:
    raise SystemExit(f"Error: {e} (library sequences have length {length})")

tree = VPTree(seqs, leaf_size=args.leaf)
print(f"Library: {len(seqs)} sequences, tree height {tree.stats()['height']}")
tree.reset_counter()
results = tree.range_search(q, args.radius) if args.radius is not None else tree.knn(q, args.k)
for d, i in results:
    print(f"distance {d:<3} {records[i]['name']:<10} {records[i]['sequence']}")
if not results:
    print("No sequence found.")
print(f"\nVP-tree computed {tree.n_dist} distances; brute force would compute {len(seqs)}.")
