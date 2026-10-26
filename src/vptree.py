"""VP-tree construction (Member 2).

Structure
---------
Every sequence stays in the list `tree.seqs`; nodes only store INDICES into that
list, so the data is never copied and search results can report `(distance, index)`.

Internal node : vp = index of the vantage point, mu = median distance,
                child[0] = inside subtree  (distance to vp <= mu)
                child[1] = outside subtree (distance to vp >  mu)
Leaf node     : vp = None, bucket = list of indices (brute-forced during search)

Search (Member 3, src/search.py) uses exactly these fields plus `tree.dist()`.
"""
import random
from typing import Callable, List, Optional, Sequence

from .sequence import hamming_distance


class Node:
    __slots__ = ("vp", "mu", "child", "bucket")

    def __init__(self, vp: Optional[int] = None, mu: int = 0, bucket: Optional[List[int]] = None):
        self.vp = vp
        self.mu = mu
        self.child = [None, None]  # [inside, outside]
        self.bucket = bucket


class VPTree:
    def __init__(self, seqs: Sequence[int], distance: Callable[[int, int], float] = hamming_distance,
                 leaf_size: int = 8, seed: int = 0):
        if leaf_size < 1:
            raise ValueError("leaf_size must be >= 1")
        self.seqs = list(seqs)
        self.distance = distance
        self.leaf_size = leaf_size
        self.rng = random.Random(seed)
        self.n_dist = 0  # distance computations made by QUERIES (reset with reset_counter)
        holder = [None]
        self._build(list(range(len(self.seqs))), holder)
        self.root = holder[0]

    # ----- query-time distance with counter (used by Member 3's search code) -----
    def dist(self, query: int, idx: int):
        self.n_dist += 1
        return self.distance(query, self.seqs[idx])

    def reset_counter(self) -> None:
        self.n_dist = 0

    # ----- search entry points (implemented in src/search.py by Member 3) -----
    def knn(self, query, k):
        from .search import knn_search
        return knn_search(self, query, k)

    def range_search(self, query, r):
        from .search import range_search
        return range_search(self, query, r)

    # ----- build -----
    def _pick_vantage(self, items: List[int]) -> int:
        """Sample a few candidates and keep the one whose distances to a probe
        sample have the highest variance: a well-spread vantage point splits the
        data more evenly, which gives better pruning than a random pick."""
        if len(items) <= 2:
            return items[0]
        cands = self.rng.sample(items, min(5, len(items)))
        probe = self.rng.sample(items, min(20, len(items)))
        best, best_var = cands[0], -1.0
        for c in cands:
            ds = [self.distance(self.seqs[c], self.seqs[p]) for p in probe if p != c]
            if len(ds) < 2:
                continue
            m = sum(ds) / len(ds)
            var = sum((x - m) ** 2 for x in ds) / len(ds)
            if var > best_var:
                best, best_var = c, var
        return best

    def _build(self, items: List[int], holder: list) -> None:
        """Iterative build (explicit stack) so deep or degenerate data cannot hit
        Python's recursion limit."""
        stack = [(items, holder, 0)]
        while stack:
            items, hold, slot = stack.pop()
            if not items:
                continue
            if len(items) <= self.leaf_size:
                hold[slot] = Node(bucket=items)
                continue
            vp = self._pick_vantage(items)
            vp_seq = self.seqs[vp]
            ds = sorted((self.distance(vp_seq, self.seqs[i]), i) for i in items if i != vp)
            if ds[-1][0] == 0:  # everything identical to the vantage point: one bucket
                hold[slot] = Node(bucket=items)
                continue
            mu = ds[(len(ds) - 1) // 2][0]  # median distance (inside is never empty)
            inside = [i for d, i in ds if d <= mu]
            outside = [i for d, i in ds if d > mu]
            node = Node(vp=vp, mu=mu)
            hold[slot] = node
            stack.append((inside, node.child, 0))
            stack.append((outside, node.child, 1))

    # ----- statistics (handy for tests, README and the interview) -----
    def stats(self) -> dict:
        nodes = leaves = height = 0
        stack = [(self.root, 1)] if self.root else []
        while stack:
            n, depth = stack.pop()
            nodes += 1
            height = max(height, depth)
            if n.vp is None:
                leaves += 1
            else:
                for c in n.child:
                    if c is not None:
                        stack.append((c, depth + 1))
        return {"nodes": nodes, "leaves": leaves, "height": height, "size": len(self.seqs)}
