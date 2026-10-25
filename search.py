"""kNN and range search on a VPTree (Member 3).

Both searches rely on the triangle inequality. For a node with vantage point v,
radius mu, and a query q at distance d = dist(q, v), with current search radius tau:

  every x in the OUTSIDE subtree has dist(v, x) > mu  =>  dist(q, x) >= dist(v, x) - d > mu - d
  every x in the INSIDE  subtree has dist(v, x) <= mu =>  dist(q, x) >= d - dist(v, x) >= d - mu

so the outside can be skipped when d + tau < mu, and the inside when d - tau > mu.
We always visit the side containing the query first so that tau shrinks early.
"""
import heapq
from typing import List, Tuple

Result = Tuple[int, int]  # (distance, index)


def knn_search(tree, query, k: int) -> List[Result]:
    """The k nearest sequences, nearest first, as (distance, index)."""
    if k <= 0 or tree.root is None:
        return []
    heap: list = []  # max-heap of the best k so far, stored as (-distance, -index)

    def tau():
        return -heap[0][0] if len(heap) == k else float("inf")

    def consider(idx, d):
        if len(heap) < k:
            heapq.heappush(heap, (-d, -idx))
        elif d < tau():
            heapq.heapreplace(heap, (-d, -idx))

    def visit(node):
        if node is None:
            return
        if node.vp is None:  # leaf bucket: brute force
            for i in node.bucket:
                consider(i, tree.dist(query, i))
            return
        d = tree.dist(query, node.vp)
        consider(node.vp, d)
        inside, outside = node.child
        if d <= node.mu:
            visit(inside)
            if d + tau() >= node.mu:
                visit(outside)
        else:
            visit(outside)
            if d - tau() <= node.mu:
                visit(inside)

    visit(tree.root)
    return sorted((-nd, -ni) for nd, ni in heap)


def range_search(tree, query, r) -> List[Result]:
    """All sequences within distance r of the query, sorted by (distance, index)."""
    if tree.root is None or r < 0:
        return []
    out: List[Result] = []

    def visit(node):
        if node is None:
            return
        if node.vp is None:
            for i in node.bucket:
                d = tree.dist(query, i)
                if d <= r:
                    out.append((d, i))
            return
        d = tree.dist(query, node.vp)
        if d <= r:
            out.append((d, node.vp))
        inside, outside = node.child
        if d - r <= node.mu:
            visit(inside)
        if d + r >= node.mu:
            visit(outside)

    visit(tree.root)
    out.sort()
    return out
