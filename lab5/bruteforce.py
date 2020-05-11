from typing import Iterator, Callable
import numpy as np  # type: ignore
from convert_utils import CanonicalForm
from utils import binomial_grid, subset_by_index


class BruteForceResult:
    def __init__(self, x: np.array, inv_AMNk: np.array, Nk: np.array):
        self.x, self.inv_AMNk, self.Nk = x, inv_AMNk, Nk


def find_min_sv(cf: CanonicalForm,
                target_function: Callable[[np.array], float],
                abs_tolerance: float = 1e-3) -> BruteForceResult:
    binom_table = binomial_grid(cf.n, cf.m)
    N = np.array(list(range(cf.n)))
    min_sv, inv_AMNk_min, Nk_min = None, None, None
    for idx in range(binom_table[-1, -1]):
        Nk = subset_by_index(N, binom_table, idx)
        AMNk = np.array([cf.A[:, i] for i in Nk]).T
        if np.abs(np.linalg.det(AMNk)) > abs_tolerance:
            inv_AMNk = np.linalg.inv(AMNk)
            xNk = np.matmul(inv_AMNk, cf.b)
            if np.min(xNk) < 0:
                continue
            xN = np.zeros(cf.n)
            for i in range(len(Nk)):
                xN[Nk[i]] = xNk[i]
            if min_sv is None or target_function(xN) < target_function(min_sv):
                min_sv, inv_AMNk_min, Nk_min = xN, inv_AMNk, Nk
    return BruteForceResult(min_sv, inv_AMNk_min, Nk_min)


def bruteforce(cf: CanonicalForm) -> BruteForceResult:
    """
    A brute force algorithm for solving canonical linear programming problem.
    :param cf: parameters of problem in canonical form
    :return: support vector which minimizes the target function c^T*x
    and indices of its basis vectors
    """
    def target_f(x: np.array) -> float:
        return float(np.dot(x, cf.c))
    return find_min_sv(cf, target_f)
