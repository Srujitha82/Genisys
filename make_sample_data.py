"""Regenerate data/sample_reference.csv (run once; the CSV is committed).
Usage: python make_sample_data.py [n] [length]
"""
import csv
import sys

from src.simulate import generate_library

n = int(sys.argv[1]) if len(sys.argv) > 1 else 100
length = int(sys.argv[2]) if len(sys.argv) > 2 else 32
seqs = generate_library(n, length, seed=42)
with open("data/sample_reference.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["id", "name", "sequence"])
    for i, s in enumerate(seqs, start=1):
        w.writerow([i, f"ref_{i:03d}", s])
print(f"wrote {n} sequences of length {length} to data/sample_reference.csv")
