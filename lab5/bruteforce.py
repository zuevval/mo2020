from typing import Iterator, Callable
import numpy as np  # type: ignore
from convert_utils import CanonicalForm
from utils import binomial_grid, subset_by_index


def find_min_sv(cf: CanonicalForm,
                target_function: Callable[[np.array], float],
                abs_tolerance: float = 1e-3) -> np.array:
    binom_table = binomial_grid(cf.n, cf.m)
    N = np.array(list(range(cf.n)))
    min_sv = None
    for idx in range(binom_table[-1, -1]):
        Nk = subset_by_index(N, binom_table, idx)
        AMNk = np.array([cf.A[:, i] for i in Nk]).T
        if np.abs(np.linalg.det(AMNk)) > abs_tolerance:
            xNk = np.matmul(np.linalg.inv(AMNk), cf.b)
            if np.min(xNk) < 0:
                continue
            xN = np.zeros(cf.n)
            for i in range(len(Nk)):
                xN[Nk[i]] = xNk[i]
            if min_sv is None or target_function(xN) < target_function(min_sv):
                min_sv = xN
    return min_sv


def bruteforce(cf: CanonicalForm) -> np.array:
    """
    A brute force algorithm for solving canonical linear programming problem.
    :param cf: parameters of problem in canonical form
    :return: support vector which minimizes the target function c^T*x
    """
    def target_f(x: np.array) -> float:
        return float(np.dot(x, cf.c))
    return find_min_sv(cf, target_f)
