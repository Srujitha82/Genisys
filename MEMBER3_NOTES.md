# Member 3 – search, benchmarks, demo

Files you own: `src/search.py`, `tests/test_search.py`, `benchmark.py`, `vptree_demo.py`.
Needs Member 2's `src/vptree.py` (fields: `tree.root`, `tree.dist()`, `tree.reset_counter()`, node `.vp/.mu/.child/.bucket`).

What the code does
- `knn_search(tree, query, k)`: max-heap of the best k, tau = k-th best distance; visit the query's side first; skip the outside if d + tau < mu, the inside if d - tau > mu.
- `range_search(tree, query, r)`: same pruning with tau = r.
- `benchmark.py`: N = 1k / 10k / 100k, reads with 0-3 errors, prints a markdown table and verifies answers against brute force. Run `python benchmark.py --out results.md` and paste the table into the README Results section.
- `vptree_demo.py`: CLI, prints matches and "VP-tree computed X distances; brute force would compute N".

You must be able to explain: the triangle-inequality pruning proof, why visit the near side first, why k=1 and small radius prune well but k=5 barely prunes on this data (the 5th neighbour is far, so tau stays large).
