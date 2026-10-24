"""Member 1 demo: brute-force nearest reference sequences for a DNA read.
Usage: python brute_demo.py ACGT... --k 5
"""
import argparse

from src.bruteforce import knn_bruteforce
from src.sequence import decode_dna, encode_dna, load_library

p = argparse.ArgumentParser()
p.add_argument("read", help="DNA read (A/C/G/T)")
p.add_argument("--k", type=int, default=5)
p.add_argument("--data", default="data/sample_reference.csv")
args = p.parse_args()

seqs, records = load_library(args.data)
length = len(records[0]["sequence"])
try:
    q = encode_dna(args.read, expected_len=length)
except ValueError as e:
    raise SystemExit(f"Error: {e} (library sequences have length {length})")
for d, i in knn_bruteforce(seqs, q, args.k):
    print(f"distance {d:<3} {records[i]['name']:<10} {records[i]['sequence']}")
