from collections import defaultdict
from src.utils import TransportTable
import numpy as np


def calc_u_v(c: np.ndarray, x: np.ndarray) -> (np.array, np.array):
    """
    v[j]-u[i]=c[i,j] if x[i,j]>0
    v[0]=0
    To solve these equations, we construct a graph with m+n nodes associated with v, u:
        nodes[i] = u[i], 0 <= i < m
        nodes[m+j] = v[j], 0 <= j < n
    where edges represent equations, and then perform a breadth first search.

    :param c: m*n transition matrix
    :param x: m*n matrix; empty cells marked as `None`
    :return: potentials vectors u[m], v[n]
    """
    assert c.shape == x.shape
    (m, n) = c.shape

    # constructing a graph (adjacency list)
    graph = defaultdict(list)
    for i in range(m):
        for j in range(n):
            if x[i, j] is not None:
                graph[i].append(m+j)
                graph[m+j].append(i)

    # performing BFS run
    u, v = np.zeros(m), np.zeros(n)
    visited = [False for _ in range(m+n)]
    stack = [m + 0]
    visited[m + 0] = True  # v[0] := 0 <=> node no. m+0 is visited
    while stack:
        node = stack.pop()
        for adj in graph[node]:  # for all adjacent nodes:
            if not visited[adj]:
                visited[adj] = True
                stack.append(adj)
                if node < m:  # `node` in `u`, `adj` in `v`
                    i = node
                    j = adj - m
                    v[j] = c[i, j] + u[i]
                else:  # `node` in `v`, `adj` in `u`
                    i = adj
                    j = node - m
                    u[i] = v[j] - c[i, j]

    # performing a check: all components of u, v found?
    if False in visited:
        raise Exception("unable to define potentials, perhaps the problem is ill-posed")

    return u, v


def calc_alpha(x: np.ndarray, u: np.array, v: np.array) -> np.ndarray:
    m, n = len(u), len(v)
    assert x.shape == (m, n)
    """
    10 lines of code below are equivalent to:
    return np.array([np.array([v[j] - u[i] if x[i, j] is None else None for j in range(m)]) for i in range(n)])
    """
    alpha = []
    for i in range(m):
        row = []
        for j in range(n):
            if x[i, j] is None:
                row.append(v[j] - u[i])
            else:
                row.append(None)
        alpha.append(row)
    return np.array(alpha)


def solve_transportation_potentials(a: np.array, b: np.array, c: np.ndarray, x0: np.ndarray):
    u, v = calc_u_v(c, x0)
    alpha = calc_alpha(x0, u, v)
    t = TransportTable(u=u, v=v, c=c, x0=x0, alpha=alpha)
